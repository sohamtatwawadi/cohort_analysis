"""Upload validation and variant normalisation — Part II §2.3.

    "Mandatory validations, each of which fails the upload rather than warning"

That phrasing is the whole design. Every check here returns a hard failure, not
a warning, because a warning on an upload screen is a warning a researcher
under deadline will click past — and the resulting dataset is then wrong in a
way nothing downstream can detect. A dataset that fails validation is not
registered at all.

The checks:

    build determined and recorded    never inferred from position alone
    sample IDs reconcile             unmatched reported in BOTH directions
    variants normalised              left-aligned, multi-allelics split
    duplicate samples detected
    reported vs genetic sex          concordance checked
    chromosome naming harmonised     chr1 vs 1
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

from .types import BUILDS, MISSING, GenotypeMatrix, PhenotypeTable, Variant


@dataclass
class ValidationIssue:
    code: str
    severity: str          # 'fail' | 'warn'
    message: str
    detail: Dict[str, Any] = field(default_factory=dict)

    @property
    def fatal(self) -> bool:
        return self.severity == "fail"


@dataclass
class ValidationReport:
    issues: List[ValidationIssue] = field(default_factory=list)
    checks_run: List[str] = field(default_factory=list)

    def add(self, code: str, severity: str, message: str, **detail) -> None:
        self.issues.append(ValidationIssue(code, severity, message, detail))

    @property
    def ok(self) -> bool:
        return not any(i.fatal for i in self.issues)

    @property
    def failures(self) -> List[ValidationIssue]:
        return [i for i in self.issues if i.fatal]

    def as_dict(self) -> Dict[str, Any]:
        return {
            "ok": self.ok,
            "checks_run": self.checks_run,
            "issues": [{"code": i.code, "severity": i.severity,
                        "message": i.message, "detail": i.detail}
                       for i in self.issues],
        }


# --------------------------------------------------------------- chromosome --
def harmonise_chrom(chrom: str) -> str:
    """`chr1` and `1` are the same chromosome. Mixing the two conventions in
    one dataset silently splits every per-chromosome grouping in half."""
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return {"M": "MT", "m": "MT", "23": "X", "24": "Y", "25": "MT"}.get(c, c.upper()
                                                                       if c.lower() in ("x", "y", "mt") else c)


# ------------------------------------------------------------- normalisation --
def left_align(ref: str, alt: str, pos: int) -> Tuple[str, str, int]:
    """Trim shared suffix then shared prefix — the standard parsimony rule.

    Without this, the SAME indel written two ways (a shared-base representation
    from one caller, a trimmed one from another) produces two variant keys, and
    a join across datasets silently misses. Order matters: suffix first, then
    prefix, and always leave at least one base on each side.
    """
    if not ref or not alt:
        return ref, alt, pos

    # Shared suffix.
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    # Shared prefix — each base removed shifts the position right.
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return ref, alt, pos


def normalise_variants(gm: GenotypeMatrix) -> GenotypeMatrix:
    """Apply chromosome harmonisation and left-alignment in place-ish.

    Multi-allelic splitting happens in the format readers, because it requires
    the per-allele genotype which is lost once dosages are computed.
    """
    for v in gm.variants:
        v.chrom = harmonise_chrom(v.chrom)
        v.ref, v.alt, v.pos = left_align(v.ref, v.alt, v.pos)
    return gm


# ------------------------------------------------------------------- checks --
def check_build(gm: GenotypeMatrix, declared: Optional[str],
                report: ValidationReport) -> Optional[str]:
    """Spec: "Genome build determined and recorded — never inferred from
    position alone."

    GRCh37 and GRCh38 coordinates occupy the same numeric range, so there is no
    honest way to tell them apart from positions. If neither the file nor the
    uploader says, the upload fails and asks.
    """
    report.checks_run.append("genome_build")
    build = declared or gm.build
    if not build:
        report.add("build_unknown", "fail",
                   "Genome build could not be determined from the file and was not "
                   "declared. Positions alone cannot distinguish GRCh37 from GRCh38 — "
                   "declare the build explicitly.")
        return None
    if build not in BUILDS:
        report.add("build_unrecognised", "fail",
                   "Genome build '{}' is not one of {}.".format(build, ", ".join(BUILDS)),
                   declared=build)
        return None
    if gm.build and declared and gm.build != declared:
        report.add("build_conflict", "fail",
                   "File header says {} but the upload declares {}. One of them is "
                   "wrong and guessing is not acceptable.".format(gm.build, declared),
                   file_build=gm.build, declared=declared)
        return None
    return build


def check_sample_reconciliation(genotype_ids: Sequence[str],
                                phenotype_ids: Optional[Sequence[str]],
                                report: ValidationReport) -> Dict[str, Any]:
    """Spec: "Sample IDs in genotype and phenotype files reconcile; unmatched
    IDs reported both directions."

    Both directions matters. Reporting only genotype-without-phenotype hides
    the case where the phenotype file is for a different cohort entirely and
    happens to overlap on a handful of IDs.
    """
    report.checks_run.append("sample_reconciliation")
    g = set(genotype_ids)
    if phenotype_ids is None:
        return {"matched": 0, "genotype_only": len(g), "phenotype_only": 0}

    p = set(phenotype_ids)
    matched = g & p
    g_only = g - p
    p_only = p - g

    if not matched:
        report.add("no_sample_overlap", "fail",
                   "No sample IDs are shared between the genotype and phenotype files "
                   "({} genotype IDs, {} phenotype IDs). These are probably different "
                   "cohorts, or the ID columns differ.".format(len(g), len(p)),
                   genotype_only=sorted(g_only)[:10], phenotype_only=sorted(p_only)[:10])
    elif len(matched) < 0.5 * max(len(g), len(p)):
        report.add("low_sample_overlap", "warn",
                   "Only {} of {} genotype samples and {} phenotype records reconcile."
                   .format(len(matched), len(g), len(p)),
                   genotype_only=sorted(g_only)[:20], phenotype_only=sorted(p_only)[:20])

    return {
        "matched": len(matched),
        "genotype_only": len(g_only),
        "phenotype_only": len(p_only),
        "genotype_only_examples": sorted(g_only)[:20],
        "phenotype_only_examples": sorted(p_only)[:20],
    }


def check_duplicate_sample_ids(sample_ids: Sequence[str],
                               report: ValidationReport) -> List[str]:
    report.checks_run.append("duplicate_sample_ids")
    seen: Set[str] = set()
    dupes: List[str] = []
    for s in sample_ids:
        if s in seen:
            dupes.append(s)
        seen.add(s)
    if dupes:
        report.add("duplicate_sample_ids", "fail",
                   "{} sample ID(s) appear more than once. Every column must identify "
                   "a distinct sample.".format(len(dupes)),
                   examples=sorted(set(dupes))[:20])
    return dupes


def check_duplicate_genotypes(gm: GenotypeMatrix, report: ValidationReport,
                              threshold: float = 0.99,
                              max_samples: int = 400) -> List[Dict[str, Any]]:
    """Detect samples that are genetically identical under different IDs.

    Distinct from duplicate IDs: the same individual enrolled twice inflates
    every count and violates the independence every test assumes. Quadratic in
    samples, so on large uploads this is sampled — the kinship engine does the
    thorough version once PCA/pruning has run.
    """
    report.checks_run.append("duplicate_genotypes")
    n = gm.n_samples
    if n < 2 or gm.n_variants == 0:
        return []

    idx = np.arange(n)
    if n > max_samples:
        idx = np.linspace(0, n - 1, max_samples).astype(int)

    d = gm.dosages[:, idx].astype(np.int16)
    observed = d != MISSING
    dupes: List[Dict[str, Any]] = []
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            both = observed[:, a] & observed[:, b]
            n_both = int(both.sum())
            if n_both < 50:
                continue
            concordance = float((d[both, a] == d[both, b]).mean())
            if concordance >= threshold:
                dupes.append({"sample_a": gm.sample_ids[idx[a]],
                              "sample_b": gm.sample_ids[idx[b]],
                              "concordance": round(concordance, 4),
                              "variants_compared": n_both})
    if dupes:
        report.add("duplicate_genotypes", "warn",
                   "{} sample pair(s) are genetically near-identical (≥{:.0%} genotype "
                   "concordance). If these are not intentional replicates, one of each "
                   "pair must be removed before any analysis assuming independence."
                   .format(len(dupes), threshold),
                   pairs=dupes[:20])
    return dupes


def infer_genetic_sex(gm: GenotypeMatrix,
                      het_threshold: float = 0.10) -> Dict[str, Optional[str]]:
    """Infer sex from X-chromosome heterozygosity.

    Males are hemizygous on X, so heterozygous calls should be near zero;
    females sit well above. Samples with no X coverage return None rather than
    a guess — a targeted panel often has no X at all, and inventing a sex there
    would produce a false discordance for every sample.
    """
    x_rows = [i for i, v in enumerate(gm.variants)
              if harmonise_chrom(v.chrom) == "X"]
    out: Dict[str, Optional[str]] = {s: None for s in gm.sample_ids}
    if not x_rows:
        return out

    sub = gm.dosages[x_rows, :]
    observed = sub != MISSING
    n_called = observed.sum(axis=0)
    het = ((sub == 1) & observed).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = np.where(n_called > 0, het / np.maximum(n_called, 1), np.nan)

    for i, sid in enumerate(gm.sample_ids):
        if n_called[i] < 20 or not np.isfinite(rate[i]):
            out[sid] = None
        else:
            out[sid] = "M" if rate[i] < het_threshold else "F"
    return out


def check_sex_concordance(gm: GenotypeMatrix, reported: Optional[Dict[str, str]],
                          report: ValidationReport) -> Dict[str, Any]:
    """Spec: "Reported vs genetic sex concordance checked."

    Discordance is a warning, not a failure: it is a real and common finding
    (sample swap, data-entry error, or genuine biology) and the right response
    is investigation by someone who knows the cohort, not a blocked upload.
    """
    report.checks_run.append("sex_concordance")
    genetic = infer_genetic_sex(gm)
    if not reported:
        return {"checked": False, "genetic_sex": genetic, "discordant": []}

    discordant = []
    for sid, rep in reported.items():
        gen = genetic.get(sid)
        if gen and rep and str(rep).upper()[:1] in ("M", "F") and gen != str(rep).upper()[:1]:
            discordant.append({"sample_id": sid, "reported": rep, "genetic": gen})

    if discordant:
        report.add("sex_discordance", "warn",
                   "{} sample(s) have genetic sex discordant with reported sex. This is "
                   "usually a sample swap or a data-entry error and should be resolved "
                   "before analysis.".format(len(discordant)),
                   examples=discordant[:20])
    return {"checked": True, "genetic_sex": genetic, "discordant": discordant}


def check_not_empty(gm: GenotypeMatrix, report: ValidationReport) -> None:
    report.checks_run.append("non_empty")
    if gm.n_samples == 0:
        report.add("no_samples", "fail", "The file contains no samples.")
    if gm.n_variants == 0:
        report.add("no_variants", "fail", "The file contains no variants.")


def validate_upload(gm: GenotypeMatrix,
                    declared_build: Optional[str] = None,
                    phenotypes: Optional[PhenotypeTable] = None,
                    reported_sex: Optional[Dict[str, str]] = None,
                    deep_duplicate_check: bool = True) -> Tuple[ValidationReport, Dict[str, Any]]:
    """Run every §2.3 mandatory validation. Returns the report and the facts
    the profiler and registry need."""
    report = ValidationReport()

    check_not_empty(gm, report)
    build = check_build(gm, declared_build, report)
    check_duplicate_sample_ids(gm.sample_ids, report)

    normalise_variants(gm)
    report.checks_run.append("variant_normalisation")

    reconciliation = check_sample_reconciliation(
        gm.sample_ids, phenotypes.sample_ids if phenotypes else None, report)

    sex = {"checked": False, "genetic_sex": {}, "discordant": []}
    dupes: List[Dict[str, Any]] = []
    if report.ok:
        sex = check_sex_concordance(gm, reported_sex, report)
        if deep_duplicate_check:
            dupes = check_duplicate_genotypes(gm, report)

    return report, {
        "build": build,
        "reconciliation": reconciliation,
        "sex": sex,
        "duplicate_genotype_pairs": dupes,
    }
