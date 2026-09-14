"""VariMAT loader (spec §3.2, §3.3).

One VariMAT file = one sample. A cohort is many samples, so the public entry
point takes a *list* of files and ingests them into the analytics store in one
transaction.

Three things here are load-bearing and every one of them is a documented trap
(spec E04):

1. DEDUP.  Rows per distinct variant is 1.75x because the same variant is
   annotated once per transcript, under multiple gene symbols. We dedup on
   CHROM:START:REF:ALT and pick the canonical row by MANE. Skip this and every
   count is ~75% too high.

2. REVIEWABLE SUBSET.  autoACMGPrediction runs genome-wide and marks ~309k
   variants VUS in a *single* sample. A raw VUS count is meaningless. Only the
   reviewable subset (PASS + coding consequence + rare) reaches the VUS
   inventory.

3. FINGERPRINT.  Coverage has no identifier in VariMAT (CRDB is empty). The
   set of genes with any ONTARGET call is the only evidence of what the assay
   looked at, and it is what the denominator service clusters on (§3.4).

Raw intronic/intergenic rows are ~98% of the file and are dropped after
fingerprinting (§3.2 storage).
"""
from __future__ import annotations

import csv
import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

from ..config import REVIEWABLE_MAX_GNOMAD_AF
from ..reference.genes import GENE_DISEASE, resolve_symbol

# --------------------------------------------------------------- contract ----
# Spec §3.3 field mapping. Names are matched case-insensitively; the loader
# tolerates a superset of columns (a real file has 230).
COLUMNS = {
    "chrom": ["CHROM", "CHR", "CHROMOSOME"],
    "start": ["START", "POS", "POSITION"],
    "ref": ["REF", "REF_ALLELE"],
    "alt": ["ALT", "ALT_ALLELE"],
    "gene_id": ["GENE_ID"],
    "gene_name": ["GENE_NAME", "GENE", "SYMBOL"],
    "varclass": ["VARCLASS"],
    "vartype": ["VARTYPE"],
    "aa_chg": ["AA_CHG"],
    "aa_pos": ["AA_POS"],
    "prot_len": ["PROT_LEN"],
    "cdna_chg": ["CDNA_CHG"],
    "zygosity": ["ZYGOSITY"],
    "vaf": ["ALT_ALLELE_PERCENTAGE"],
    "depth": ["OVERALL_READ_DEPTH"],
    "alt_depth": ["ALT_DEPTH"],
    "filter_status": ["VARIANT_FILTER_STATUS"],
    "var_qual": ["VAR_QUAL"],
    "location": ["VARIANT_LOCATION"],
    "acmg": ["autoACMGPrediction", "AUTOACMGPREDICTION"],
    "acmg_rules": ["autoACMGRules", "AUTOACMGRULES"],
    "acmg_rules_info": ["autoACMGRulesInfo", "AUTOACMGRULESINFO"],
    "clinvar_sig": ["ClinVar_Significance", "CLINVAR_SIGNIFICANCE"],
    "clinvar_disease": ["ClinVar_Disease", "CLINVAR_DISEASE"],
    "clinvar_id": ["ClinVar_ID", "CLINVAR_ID"],
    "gnomad_af": ["gnomAD_AF", "GNOMAD_AF"],
    "gnomad_sas_af": ["gnomAD_SAS_AF", "GNOMAD_SAS_AF"],
    "ga100k_sas_af": ["GA100K_SAS_af", "GA100K_SAS_AF"],
    "gnomad_v2_sas": ["gnomADv2_AF_sas", "GNOMADV2_AF_SAS"],
    "medvardb_af": ["MedVarDb_ALT_AF", "MEDVARDB_ALT_AF"],
    "medvardb_het": ["MedVarDb_PassHetSamples"],
    "medvardb_hom": ["MedVarDb_PassHomSamples"],
    "omim_id": ["OMIM_ID"],
    "omim_disease": ["OMIM_DISEASE"],
    "mane": ["MANE"],
    "refseq_id": ["REFSEQ_ID"],
    "canonical": ["CANNONICAL_TRAS", "CANONICAL_TRANS"],
}

NA_VALUES = {"", "NA", "N/A", ".", "-", "NULL", "None", "nan", "NaN"}

# ------------------------------------------------- consequence normalisation --
# VARCLASS vocabulary -> the consequence vocabulary the UI and colour
# conventions use (spec §2.8).
_CONSEQ_EXACT = {
    "MISSENSE": "Missense",
    "NONSENSE": "Nonsense",
    "STOPGAIN": "Nonsense",
    "STOP-GAIN": "Nonsense",
    "SILENT": "Synonymous",
    "SYNONYMOUS": "Synonymous",
    "5UTR": "5' UTR",
    "3UTR": "3' UTR",
    "EXONIC-NC": "Non-coding exonic",
    "STARTLOSS": "Start loss",
    "STOPLOSS": "Stop loss",
}

# Coding / protein-altering consequences that enter the reviewable subset.
REVIEWABLE_CONSEQUENCES = {
    "Missense", "Nonsense", "Frameshift", "In-frame indel", "Splice site",
    "Synonymous", "5' UTR", "3' UTR", "Non-coding exonic", "Start loss", "Stop loss",
}

# Consequences that alter the protein — the tightest tier of the §3.2 funnel.
PROTEIN_ALTERING = {
    "Missense", "Nonsense", "Frameshift", "In-frame indel", "Splice site",
    "Start loss", "Stop loss",
}

TRUNCATING = {"Nonsense", "Frameshift", "Splice site", "Exon deletion"}


def normalise_consequence(varclass: Optional[str]) -> Optional[str]:
    if not varclass:
        return None
    v = varclass.strip().upper()
    if v in NA_VALUES:
        return None
    if v in _CONSEQ_EXACT:
        return _CONSEQ_EXACT[v]
    if v.startswith("FRAMESHIFT"):
        return "Frameshift"
    if v.startswith("INFRAME") or v.startswith("IN-FRAME"):
        return "In-frame indel"
    # "*-SS-*" splice-site forms, e.g. ACCEPTOR-SS-VARIANT, DONOR-SS-LOSS
    if "-SS-" in v or v.endswith("-SS") or "SPLICE" in v:
        return "Splice site"
    if "UTR" in v:
        return "5' UTR" if "5" in v else "3' UTR"
    if "INTRON" in v or "INTERGENIC" in v or "UPSTREAM" in v or "DOWNSTREAM" in v:
        return None  # not reviewable; dropped after fingerprinting
    return varclass.strip().title()


_VARTYPE_MAP = {
    "SNV": "SNV", "SNP": "SNV", "MNV": "SNV", "SUB": "SNV", "SUBSTITUTION": "SNV",
    "INS": "Indel", "DEL": "Indel", "INDEL": "Indel", "DELINS": "Indel",
    "CNV": "CNV", "DUP": "CNV", "AMP": "CNV",
    "SV": "SV", "BND": "SV", "INV": "SV", "TRA": "SV", "FUSION": "SV",
}


def normalise_var_class(vartype: Optional[str], consequence: Optional[str]) -> str:
    if vartype:
        v = vartype.strip().upper()
        if v in _VARTYPE_MAP:
            return _VARTYPE_MAP[v]
    if consequence == "Splice site":
        return "Splice"
    if consequence in ("Frameshift", "In-frame indel"):
        return "Indel"
    return "SNV"


_ACMG_MAP = {
    "PATHOGENIC": "Pathogenic",
    "LIKELY PATHOGENIC": "Likely pathogenic",
    "LIKELY_PATHOGENIC": "Likely pathogenic",
    "UNCERTAIN SIGNIFICANCE": "Uncertain significance",
    "VUS": "Uncertain significance",
    "UNCERTAIN_SIGNIFICANCE": "Uncertain significance",
    "LIKELY BENIGN": "Likely benign",
    "LIKELY_BENIGN": "Likely benign",
    "BENIGN": "Benign",
}


def normalise_acmg(value: Optional[str]) -> str:
    if not value:
        return "Uncertain significance"
    return _ACMG_MAP.get(value.strip().upper(), "Uncertain significance")


# ------------------------------------------------------------------ parsing --
def _num(value: Optional[str]) -> Optional[float]:
    if value is None:
        return None
    v = value.strip()
    if v in NA_VALUES:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _int(value: Optional[str]) -> Optional[int]:
    f = _num(value)
    return None if f is None else int(f)


def _truthy(value: Optional[str]) -> bool:
    if value is None:
        return False
    return value.strip().upper() in {"Y", "YES", "TRUE", "1", "MANE", "MANE_SELECT", "MANE SELECT"}


def _open_text(path: Path) -> io.TextIOBase:
    if path.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _sniff_delimiter(header_line: str) -> str:
    return "," if header_line.count(",") > header_line.count("\t") else "\t"


class ColumnMap:
    """Resolves the logical field names in COLUMNS against a file's header."""

    def __init__(self, header: Sequence[str]):
        lookup = {h.strip().upper(): h for h in header}
        self._map: Dict[str, Optional[str]] = {}
        for logical, candidates in COLUMNS.items():
            actual = None
            for cand in candidates:
                actual = lookup.get(cand.upper())
                if actual:
                    break
            self._map[logical] = actual
        self.header = list(header)

    def get(self, row: Dict[str, str], logical: str) -> Optional[str]:
        col = self._map.get(logical)
        if not col:
            return None
        val = row.get(col)
        if val is None:
            return None
        val = val.strip()
        return None if val in NA_VALUES else val

    @property
    def missing_required(self) -> List[str]:
        required = ["chrom", "start", "ref", "alt", "varclass", "filter_status"]
        return [r for r in required if not self._map.get(r)]


# -------------------------------------------------------------------- funnel --
@dataclass
class LoadFunnel:
    """The §3.2 observed funnel, reproduced per file so dedup can be verified.

    Observed on the real exome: 195,798 distinct -> 125,740 PASS -> 7,028
    coding -> 4,968 rare -> 1,124 protein-altering.
    """
    source_file: str = ""
    raw_rows: int = 0
    distinct_variants: int = 0
    pass_filter: int = 0
    coding: int = 0
    rare: int = 0
    protein_altering: int = 0
    on_target_genes: int = 0
    unresolved_symbols: Set[str] = field(default_factory=set)

    @property
    def rows_per_variant(self) -> float:
        return self.raw_rows / self.distinct_variants if self.distinct_variants else 0.0

    def as_dict(self) -> Dict[str, Any]:
        return {
            "source_file": self.source_file,
            "raw_rows": self.raw_rows,
            "distinct_variants": self.distinct_variants,
            "rows_per_variant": round(self.rows_per_variant, 3),
            "pass_filter": self.pass_filter,
            "coding": self.coding,
            "rare": self.rare,
            "protein_altering": self.protein_altering,
            "on_target_genes": self.on_target_genes,
            "unresolved_symbols": sorted(self.unresolved_symbols)[:50],
        }


@dataclass
class ParsedVariant:
    variant_key: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gene_symbol: Optional[str]
    gene_id: Optional[str]
    consequence: Optional[str]
    var_class: str
    hgvs_c: Optional[str]
    hgvs_p: Optional[str]
    aa_pos: Optional[int]
    protein_len: Optional[int]
    zygosity_raw: Optional[str]
    vaf: Optional[float]
    depth: Optional[int]
    alt_depth: Optional[int]
    filter_status: Optional[str]
    gnomad_af: Optional[float]
    gnomad_sas_af: Optional[float]
    ga100k_sas_af: Optional[float]
    clinvar_sig: Optional[str]
    clinvar_id: Optional[str]
    acmg: str
    acmg_codes: Optional[str]
    mane: bool
    on_target: bool
    reviewable: bool


TRANSITION_PAIRS = {"AG", "GA", "CT", "TC"}


def compute_qc(variants: List["ParsedVariant"]) -> Dict[str, Any]:
    """Sample QC over every PASS call, before the reviewable filter.

    This has to happen here. The reviewable subset that reaches the store is a
    handful of rows per sample — enough to report on, nowhere near enough to
    compute a Ti/Tv or a het/hom ratio from. Once the intronic rows are archived
    (§3.2 storage) the information is gone.
    """
    n_ti = n_tv = n_het = n_hom = n_x = n_x_het = 0
    vafs: List[float] = []
    n = 0
    for v in variants:
        if not v.filter_status or v.filter_status.upper() != "PASS":
            continue
        n += 1
        if len(v.ref) == 1 and len(v.alt) == 1 and v.ref != v.alt:
            if (v.ref + v.alt).upper() in TRANSITION_PAIRS:
                n_ti += 1
            else:
                n_tv += 1
        z = (v.zygosity_raw or "").lower()
        if z.startswith("het"):
            n_het += 1
            if v.vaf is not None:
                vafs.append(v.vaf)
        elif z.startswith("hom"):
            n_hom += 1
        if v.chrom in ("X", "chrX"):
            n_x += 1
            if z.startswith("het"):
                n_x_het += 1
    return {
        "qc_n_called": n,
        "qc_ti_tv": round(n_ti / n_tv, 4) if n_tv else None,
        "qc_het_hom": round(n_het / n_hom, 4) if n_hom else None,
        "qc_mean_het_vaf": round(sum(vafs) / len(vafs), 4) if vafs else None,
        "qc_x_het_rate": round(n_x_het / n_x, 4) if n_x >= 5 else None,
    }


@dataclass
class ParsedSample:
    sample_id: str
    source_file: str
    variants: List[ParsedVariant]
    fingerprint: Set[str]
    funnel: LoadFunnel
    qc: Dict[str, Any] = field(default_factory=dict)
    reference_build: Optional[str] = None
    pipeline_version: Optional[str] = None
    caller: Optional[str] = None


# Filename convention (spec §3.4): sample ID, pipeline version, caller,
# reference build. No panel or test code — that is exactly why coverage must be
# inferred.
_FILENAME_RE = re.compile(
    r"^(?P<sample>[A-Za-z0-9_.-]+?)"
    r"(?:[._](?P<pipeline>(?:LINC|PIPE)[-_]?v?[\d.]+))?"
    r"(?:[._](?P<caller>DRAGEN|GATK|STRELKA|DEEPVARIANT)[-_]?(?P<callerver>[\d.]+)?)?"
    r"(?:[._](?P<build>GRCh3[78]|hg19|hg38))?"
    r"(?:\.varimat)?(?:\.tsv|\.txt|\.csv)?(?:\.gz)?$",
    re.IGNORECASE,
)


def parse_filename(path: Path) -> Dict[str, Optional[str]]:
    m = _FILENAME_RE.match(path.name)
    if not m:
        return {"sample_id": path.stem, "pipeline_version": None,
                "caller": None, "reference_build": None}
    g = m.groupdict()
    caller = g.get("caller")
    if caller and g.get("callerver"):
        caller = "{}-{}".format(caller.upper(), g["callerver"])
    elif caller:
        caller = caller.upper()
    build = g.get("build")
    if build:
        build = {"hg19": "GRCh37", "hg38": "GRCh38"}.get(build.lower(), build)
    return {
        "sample_id": g.get("sample") or path.stem,
        "pipeline_version": g.get("pipeline"),
        "caller": caller,
        "reference_build": build,
    }


# -------------------------------------------------------------------- loader --
def parse_file(path: Path, sample_id: Optional[str] = None) -> ParsedSample:
    """Parse one VariMAT file into a deduped, classified sample record."""
    path = Path(path)
    meta = parse_filename(path)
    sid = sample_id or meta["sample_id"]
    funnel = LoadFunnel(source_file=path.name)

    # variant_key -> best row seen so far. "Best" = MANE, then longest
    # protein consequence evidence, so the canonical transcript wins (§3.2).
    best: Dict[str, ParsedVariant] = {}
    fingerprint: Set[str] = set()

    with _open_text(path) as fh:
        first = fh.readline()
        if not first:
            raise ValueError("{}: file is empty".format(path.name))
        delim = _sniff_delimiter(first)
        header = next(csv.reader([first], delimiter=delim))
        cmap = ColumnMap(header)
        missing = cmap.missing_required
        if missing:
            raise ValueError(
                "{}: not a VariMAT file — missing required columns {}".format(
                    path.name, ", ".join(missing)))

        reader = csv.DictReader(fh, fieldnames=header, delimiter=delim)
        for raw in reader:
            funnel.raw_rows += 1
            parsed = _parse_row(raw, cmap, funnel)
            if parsed is None:
                continue

            # Fingerprint BEFORE any consequence filtering: an ONTARGET call of
            # any kind is evidence the assay looked at that gene (§3.4 step 1).
            if parsed.on_target and parsed.gene_symbol:
                fingerprint.add(parsed.gene_symbol)

            prev = best.get(parsed.variant_key)
            if prev is None or _prefer(parsed, prev):
                best[parsed.variant_key] = parsed

    funnel.distinct_variants = len(best)
    funnel.on_target_genes = len(fingerprint)

    kept: List[ParsedVariant] = []
    for v in best.values():
        if v.filter_status and v.filter_status.upper() == "PASS":
            funnel.pass_filter += 1
        else:
            continue
        if v.consequence in REVIEWABLE_CONSEQUENCES:
            funnel.coding += 1
        else:
            continue
        if v.gnomad_af is None or v.gnomad_af < REVIEWABLE_MAX_GNOMAD_AF:
            funnel.rare += 1
        else:
            continue
        if v.consequence in PROTEIN_ALTERING:
            funnel.protein_altering += 1
        kept.append(v)

    return ParsedSample(
        sample_id=sid,
        source_file=str(path),
        variants=kept,
        fingerprint=fingerprint,
        funnel=funnel,
        qc=compute_qc(list(best.values())),
        reference_build=meta["reference_build"],
        pipeline_version=meta["pipeline_version"],
        caller=meta["caller"],
    )


def _prefer(new: ParsedVariant, old: ParsedVariant) -> bool:
    """Canonical-row selection for a duplicated locus (§3.2).

    MANE wins outright. Otherwise prefer the row that actually resolved to a
    known gene, then the one with a protein consequence, then the one carrying
    an amino-acid position — i.e. the most informative annotation.
    """
    if new.mane != old.mane:
        return new.mane
    if bool(new.gene_symbol) != bool(old.gene_symbol):
        return bool(new.gene_symbol)
    new_pa = new.consequence in PROTEIN_ALTERING
    old_pa = old.consequence in PROTEIN_ALTERING
    if new_pa != old_pa:
        return new_pa
    if (new.aa_pos is not None) != (old.aa_pos is not None):
        return new.aa_pos is not None
    return False


def _parse_row(raw: Dict[str, str], cmap: ColumnMap, funnel: LoadFunnel) -> Optional[ParsedVariant]:
    chrom = cmap.get(raw, "chrom")
    start = cmap.get(raw, "start")
    ref = cmap.get(raw, "ref")
    alt = cmap.get(raw, "alt")
    if not chrom or not start or ref is None or alt is None:
        return None
    try:
        pos = int(float(start))
    except ValueError:
        return None

    chrom = chrom.replace("chr", "").replace("CHR", "")
    variant_key = "{}:{}:{}:{}".format(chrom, pos, ref, alt)

    gene_id = cmap.get(raw, "gene_id")
    gene_name = cmap.get(raw, "gene_name")
    symbol = resolve_symbol(gene_id, gene_name)
    if symbol is None and gene_name:
        funnel.unresolved_symbols.add(gene_name)

    consequence = normalise_consequence(cmap.get(raw, "varclass"))
    var_class = normalise_var_class(cmap.get(raw, "vartype"), consequence)

    vaf_pct = _num(cmap.get(raw, "vaf"))          # VariMAT stores percent, e.g. 47.73
    vaf = None if vaf_pct is None else round(vaf_pct / 100.0, 4)

    location = cmap.get(raw, "location")
    on_target = (location or "").upper() == "ONTARGET"

    protein_len = _int(cmap.get(raw, "prot_len"))
    if protein_len is None and symbol:
        rec = GENE_DISEASE.get(symbol)
        protein_len = rec.protein_len if rec else None

    return ParsedVariant(
        variant_key=variant_key,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        gene_symbol=symbol,
        gene_id=gene_id,
        consequence=consequence,
        var_class=var_class,
        hgvs_c=cmap.get(raw, "cdna_chg"),
        hgvs_p=cmap.get(raw, "aa_chg"),
        aa_pos=_int(cmap.get(raw, "aa_pos")),
        protein_len=protein_len,
        zygosity_raw=cmap.get(raw, "zygosity"),
        vaf=vaf,
        depth=_int(cmap.get(raw, "depth")),
        alt_depth=_int(cmap.get(raw, "alt_depth")),
        filter_status=cmap.get(raw, "filter_status"),
        gnomad_af=_num(cmap.get(raw, "gnomad_af")),
        gnomad_sas_af=_num(cmap.get(raw, "gnomad_sas_af")),
        ga100k_sas_af=_num(cmap.get(raw, "ga100k_sas_af")),
        clinvar_sig=cmap.get(raw, "clinvar_sig"),
        clinvar_id=cmap.get(raw, "clinvar_id"),
        acmg=normalise_acmg(cmap.get(raw, "acmg")),
        acmg_codes=cmap.get(raw, "acmg_rules"),
        mane=_truthy(cmap.get(raw, "mane")),
        on_target=on_target,
        reviewable=False,   # set by the funnel pass in parse_file
    )


def discover(directory: Path) -> List[Path]:
    """Every VariMAT file under a directory. One file = one sample."""
    directory = Path(directory)
    if not directory.exists():
        return []
    out: List[Path] = []
    for pattern in ("*.tsv", "*.txt", "*.csv", "*.tsv.gz", "*.txt.gz", "*.csv.gz"):
        out.extend(sorted(directory.rglob(pattern)))
    return sorted(set(out))


def parse_many(paths: Iterable[Path]) -> Iterator[ParsedSample]:
    """Parse a set of per-sample files — the multi-sample cohort entry point."""
    for p in paths:
        yield parse_file(Path(p))
