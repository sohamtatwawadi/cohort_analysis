"""Emit VariMAT-shaped fixture files plus a matching clinical sidecar.

    python -m backend.tools.make_varimat_fixtures data/varimat --samples 12

One file per sample, because that is how VariMAT arrives and a cohort is many
samples. The files deliberately reproduce the properties that make the real
loader hard (spec §3.2):

  * MULTI-TRANSCRIPT DUPLICATION — each variant is emitted 1-3 times under
    different transcripts and sometimes different gene symbols, with exactly one
    row flagged MANE. Target ~1.75 rows per distinct variant, matching the real
    exome. A loader that does not dedup will read ~75% too many variants, and
    this fixture is what catches that.
  * GENOME-WIDE VUS NOISE — a large block of intronic/intergenic rows carrying
    autoACMGPrediction = "Uncertain significance", which is what makes a raw VUS
    count meaningless.
  * NON-PASS ROWS — so the PASS stage of the funnel has something to drop.
  * GENE SYMBOL ALIASES — e.g. GBA1 for GBA, so the alias map is exercised.

These are fixtures, not clinical data. They exist so the ingest path can be run
and asserted before real VariMAT files are available.
"""
from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

from ..app.reference.genes import GENE_DISEASE, GENE_LIST
from ..app.reference.tests import TEST_CODES
from ..app.ingest.clinical import SIDECAR_COLUMNS

# A realistic subset of the 230 VariMAT columns — every one the loader reads,
# plus enough filler that the column-matching logic is genuinely tested.
HEADER = [
    "CHROM", "START", "END", "REF", "ALT", "GENE_ID", "GENE_NAME", "VARCLASS",
    "VARTYPE", "AA_CHG", "AA_POS", "PROT_LEN", "CDNA_CHG", "ZYGOSITY",
    "ALT_ALLELE_PERCENTAGE", "OVERALL_READ_DEPTH", "ALT_DEPTH",
    "VARIANT_FILTER_STATUS", "VAR_QUAL", "VARIANT_LOCATION",
    "autoACMGPrediction", "autoACMGRules", "autoACMGRulesInfo",
    "ClinVar_Significance", "ClinVar_Disease", "ClinVar_ID",
    "gnomAD_AF", "gnomAD_SAS_AF", "GA100K_SAS_af", "gnomADv2_AF_sas",
    "MedVarDb_ALT_AF", "MedVarDb_PassHetSamples", "MedVarDb_PassHomSamples",
    "OMIM_ID", "OMIM_DISEASE", "ClinVar_gene", "MANE", "REFSEQ_ID",
    "CANNONICAL_TRAS", "CRDB",
]

ALIASES = {"GBA": "GBA1", "STK11": "LKB1", "PTEN": "MMAC1", "SLC26A4": "PDS"}
AAS = "ACDEFGHIKLMNPQRSTVWY"
ANCESTRIES = ["South Asian", "South Asian", "Middle Eastern", "European", "East Asian"]
CONSENTS = [
    "Full research + secondary findings",
    "Full research + secondary findings",
    "Clinical only, secondary findings declined",
    "Clinical only, research declined",
    "Withdrawn",
]


def _row(**kw):
    r = {c: "NA" for c in HEADER}
    r["CRDB"] = ""            # empty in real files — why coverage must be inferred
    r["VARIANT_LOCATION"] = "ONTARGET"
    r.update(kw)
    return r


def _multi_transcript(base, rng, n_tx=None):
    """Expand one logical variant into its per-transcript annotation rows.

    Multi-transcript annotation applies to EVERY call, not only the coding ones
    — intronic and intergenic rows are annotated against overlapping transcripts
    just the same. Duplicating only the interesting rows would leave the
    file-wide rows-per-variant ratio near 1.0 and the fixture would not
    reproduce the 1.75x inflation the dedup step exists to remove.
    """
    n_tx = n_tx or rng.choice([1, 2, 2, 2, 3])
    mane_idx = rng.randrange(n_tx)
    out = []
    for t in range(n_tx):
        r = dict(base)
        r["MANE"] = "MANE_Select" if t == mane_idx else "NA"
        r["CANNONICAL_TRAS"] = "Y" if t == mane_idx else "N"
        r["REFSEQ_ID"] = "NM_{:06d}.{}".format(rng.randrange(999999), t + 1)
        if t != mane_idx:
            r["GENE_ID"] = "NA"
        out.append(r)
    return out


def emit_sample(path: Path, sample_id: str, test_code: str, rng: random.Random) -> int:
    test = TEST_CODES[test_code]
    rows = []

    # --- reportable coding variants, each duplicated across transcripts -------
    n_coding = rng.randint(14, 26)
    for _ in range(n_coding):
        gene = rng.choice(test.genes)
        rec = GENE_DISEASE[gene]
        chrom = rec.chrom
        pos = rng.randint(1_000_000, 200_000_000)
        ref, alt = rng.choice("ACGT"), rng.choice("ACGT")
        if ref == alt:
            alt = "ACGT".replace(ref, "")[0]

        varclass = rng.choice(
            ["MISSENSE"] * 6 + ["NONSENSE", "FRAMESHIFT-DEL", "INFRAME-INS",
                                "ACCEPTOR-SS-VARIANT", "SILENT"])
        aa_pos = rng.randint(1, min(rec.protein_len or 800, 4000))
        acmg = rng.choice(["Uncertain significance"] * 6 +
                          ["Pathogenic", "Likely pathogenic", "Likely benign", "Benign"])
        is_plp = acmg in ("Pathogenic", "Likely pathogenic")
        gnomad = round(10 ** -(3.2 + rng.random() * 2.6), 8) if is_plp \
            else round(10 ** -(1.6 + rng.random() * 2.4), 8)
        zyg = rng.choice(["Heterozygous"] * 8 + ["Homozygous"] * 2)
        depth = rng.randint(40, 300)
        vaf = round(rng.uniform(95, 100), 2) if zyg == "Homozygous" \
            else round(rng.uniform(38, 62), 2)

        # 1-3 transcript rows for this ONE variant. Exactly one is MANE.
        n_tx = rng.choice([1, 2, 2, 2, 3])
        mane_idx = rng.randrange(n_tx)
        for t in range(n_tx):
            # Off-MANE rows sometimes carry an alias or a neighbouring gene
            # symbol — the real cause of the 1.75x row inflation.
            symbol = gene
            if t != mane_idx and gene in ALIASES and rng.random() < 0.5:
                symbol = ALIASES[gene]
            rows.append(_row(
                CHROM=chrom, START=pos, END=pos + 1, REF=ref, ALT=alt,
                GENE_ID=rec.gene_id if t == mane_idx else "NA",
                GENE_NAME=symbol,
                VARCLASS=varclass if t == mane_idx else rng.choice([varclass, "INTRONIC"]),
                VARTYPE="SNV" if "FRAMESHIFT" not in varclass else "DEL",
                AA_CHG="p.{}{}{}".format(rng.choice(AAS), aa_pos, rng.choice(AAS)),
                AA_POS=aa_pos, PROT_LEN=rec.protein_len or "NA",
                CDNA_CHG="c.{}{}>{}".format(aa_pos * 3, ref, alt),
                ZYGOSITY=zyg,
                ALT_ALLELE_PERCENTAGE=vaf,
                OVERALL_READ_DEPTH=depth,
                ALT_DEPTH=int(depth * vaf / 100),
                VARIANT_FILTER_STATUS="PASS",
                VAR_QUAL=rng.randint(200, 900),
                autoACMGPrediction=acmg,
                autoACMGRules="PVS1,PM2" if is_plp else "PM2,BP4",
                ClinVar_Significance=acmg if rng.random() < 0.6 else "NA",
                ClinVar_Disease=rec.condition,
                ClinVar_ID="VCV{:09d}".format(rng.randrange(999999999)),
                gnomAD_AF=gnomad,
                gnomAD_SAS_AF=round(gnomad * rng.uniform(0.6, 2.2), 8),
                GA100K_SAS_af=round(gnomad * rng.uniform(0.5, 2.8), 8),
                OMIM_DISEASE=rec.condition,
                ClinVar_gene=gene,
                MANE="MANE_Select" if t == mane_idx else "NA",
                REFSEQ_ID="NM_{:06d}.{}".format(rng.randrange(999999), t + 1),
                CANNONICAL_TRAS="Y" if t == mane_idx else "N",
            ))

    # --- rows the funnel must drop -------------------------------------------
    # Non-PASS calls.
    for _ in range(rng.randint(10, 20)):
        gene = rng.choice(test.genes)
        rec = GENE_DISEASE[gene]
        rows.extend(_multi_transcript(_row(
            CHROM=rec.chrom, START=rng.randint(1_000_000, 200_000_000), END=0,
            REF="A", ALT="G", GENE_ID=rec.gene_id, GENE_NAME=gene,
            VARCLASS="MISSENSE", VARTYPE="SNV",
            VARIANT_FILTER_STATUS=rng.choice(["LowQual", "LowDepth", "StrandBias"]),
            autoACMGPrediction="Uncertain significance",
            OVERALL_READ_DEPTH=rng.randint(5, 18)), rng))

    # Genome-wide intronic/intergenic noise, all called VUS by autoACMG. This is
    # the block that makes an unfiltered VUS count meaningless (spec §3.2).
    for _ in range(rng.randint(180, 320)):
        gene = rng.choice(GENE_LIST)
        rec = GENE_DISEASE[gene]
        rows.extend(_multi_transcript(_row(
            CHROM=rec.chrom, START=rng.randint(1_000_000, 200_000_000), END=0,
            REF="C", ALT="T", GENE_ID=rec.gene_id, GENE_NAME=gene,
            VARCLASS=rng.choice(["INTRONIC", "INTERGENIC", "UPSTREAM", "DOWNSTREAM"]),
            VARTYPE="SNV", VARIANT_FILTER_STATUS="PASS",
            autoACMGPrediction="Uncertain significance",
            gnomAD_AF=round(rng.uniform(0.02, 0.4), 6),
            OVERALL_READ_DEPTH=rng.randint(30, 200)), rng))

    # Common coding variants — PASS and coding, dropped by the AF threshold.
    for _ in range(rng.randint(20, 40)):
        gene = rng.choice(test.genes)
        rec = GENE_DISEASE[gene]
        rows.extend(_multi_transcript(_row(
            CHROM=rec.chrom, START=rng.randint(1_000_000, 200_000_000), END=0,
            REF="G", ALT="A", GENE_ID=rec.gene_id, GENE_NAME=gene,
            VARCLASS="MISSENSE", VARTYPE="SNV", VARIANT_FILTER_STATUS="PASS",
            autoACMGPrediction="Likely benign",
            gnomAD_AF=round(rng.uniform(0.02, 0.45), 6),
            AA_POS=rng.randint(1, 300), PROT_LEN=rec.protein_len or "NA",
            OVERALL_READ_DEPTH=rng.randint(60, 300),
            ZYGOSITY="Heterozygous",
            ALT_ALLELE_PERCENTAGE=round(rng.uniform(40, 60), 2)), rng))

    rng.shuffle(rows)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return len(rows)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("directory", nargs="?", default="data/varimat")
    p.add_argument("--samples", type=int, default=12)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--clinical", default="data/clinical")
    args = p.parse_args(argv)

    rng = random.Random(args.seed)
    out = Path(args.directory)
    out.mkdir(parents=True, exist_ok=True)
    # Filenames encode a randomly chosen pipeline/caller/build, so regenerating
    # with a different seed leaves the previous generation behind and the
    # directory accumulates two files per sample.
    for stale in out.glob("*.varimat.tsv"):
        stale.unlink()
    clinical_dir = Path(args.clinical)
    clinical_dir.mkdir(parents=True, exist_ok=True)

    codes = list(TEST_CODES)
    sidecar = []
    total = 0

    # A few multi-member families, so probands-only has something to suppress
    # even in the fixture set.
    family_of = {}
    fam_n = 0

    for i in range(1, args.samples + 1):
        sample_id = "MG{:05d}".format(1000 + i)
        code = codes[i % len(codes)]
        test = TEST_CODES[code]
        pipeline = rng.choice(["LINC-v2.8.3", "LINC-v2.9.0", "LINC-v2.9.1"])
        caller = rng.choice(["DRAGEN-4.2", "GATK-4.5"])
        build = "GRCh38" if rng.random() < 0.85 else "GRCh37"

        # Filename carries sample, pipeline, caller, build — and NO test code.
        fname = "{}.{}.{}.{}.varimat.tsv".format(sample_id, pipeline, caller, build)
        n = emit_sample(out / fname, sample_id, code, rng)
        total += n

        if i % 4 == 1:
            fam_n += 1
            family_of[i] = ("FAM-F{:03d}".format(fam_n), True)
        else:
            fam = "FAM-F{:03d}".format(max(1, fam_n))
            family_of[i] = (fam, False) if i % 4 == 2 else (
                "FAM-F{:03d}".format(fam_n), False)
        family_id, is_proband = family_of[i]

        sidecar.append({
            "sample_id": sample_id,
            "subject_id": "SJ-" + sample_id,
            "family_id": family_id,
            "mrn": "MRN{}".format(600000 + i),
            "sex": rng.choice(["F", "M"]),
            "age": rng.randint(1, 72),
            "ancestry": rng.choice(ANCESTRIES),
            "consent_class": rng.choice(CONSENTS),
            "relation": "Proband" if is_proband else rng.choice(["Mother", "Father", "Sibling"]),
            "is_proband": "true" if is_proband else "false",
            "indication": test.indication,
            "affected_status": rng.choice(["Affected", "Unaffected", "At risk / unaffected"]),
            "family_history": rng.choice(["Positive", "Negative / unknown"]),
            "referral_source": rng.choice(["Medical genetics", "Oncology", "Paediatrics"]),
            "phenotype_hpo": "",
            "test_code": code,
            "sample_type": rng.choice(["Blood (EDTA)", "Saliva", "Buccal swab"]),
            "collection_date": "2026-{:02d}-{:02d}".format(rng.randint(1, 8), rng.randint(1, 28)),
            "pipeline_version": pipeline,
            "reference_build": build,
            "mean_depth": rng.randint(60, 240),
            "pct_bases_20x": round(rng.uniform(94, 99.8), 1),
            "qc_status": "Pass" if rng.random() < 0.95 else "Low coverage",
        })

    sc_path = clinical_dir / "subjects.csv"
    with sc_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=SIDECAR_COLUMNS)
        w.writeheader()
        for rec in sidecar:
            w.writerow({k: rec.get(k, "") for k in SIDECAR_COLUMNS})

    print("wrote {} VariMAT files ({:,} rows) to {}".format(args.samples, total, out))
    print("wrote clinical sidecar: {}".format(sc_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
