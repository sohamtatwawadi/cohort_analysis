"""Zygosity & inheritance — spec §G04.

    "Separate carrier from affected. A homozygous P/LP in a recessive gene and
     a heterozygous P/LP in the same gene are different clinical entities."

    biallelic     = zygosity in (Homozygous, Compound heterozygous)
                    AND inheritance in (AR, AR/AD)
    AR carrier    = zygosity = Heterozygous AND inheritance = AR   -- NOT affected
    hemizygous    = COMPUTED: chrX/chrY variant AND subject sex = M
    compound het  = COMPUTED: >=2 distinct P/LP variants, same gene, same subject

Done-when: "Het, hom, compound-het and hemizygous are never summed into one
number." Nothing in this module returns a total across zygosity states, and the
heatmap is a cross-tab precisely so the states stay separated.

The two computed states are produced by ingest.store.derive_zygosity over the
whole store; this module reads them. They are never read from ZYGOSITY, which
carries Het/Hom only (spec §3.5).
"""
from __future__ import annotations

from typing import Any, Dict, List

from .. import db
from .carrier import PLP_SQL, _obs_where
from .cohort import OBS_JOIN
from . import denominator as den

ZYGOSITY_ORDER = ["Heterozygous", "Homozygous", "Compound heterozygous", "Hemizygous"]
INHERITANCE_ORDER = ["AD", "AR", "AR/AD", "XLR", "XLD"]


def overview(criteria: Dict[str, Any]) -> Dict[str, Any]:
    obs = _obs_where(criteria)
    where = """WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
               AND {obs} AND {plp}""".format(obs=obs, plp=PLP_SQL)

    kpis = db.row("""
        SELECT
          COUNT(DISTINCT CASE WHEN f.zygosity IN ('Homozygous','Compound heterozygous')
                               AND gd.inheritance IN ('AR','AR/AD')
                              THEN f.subject_id END) AS biallelic_subjects,
          COUNT(DISTINCT CASE WHEN f.zygosity = 'Heterozygous'
                               AND gd.inheritance = 'AR'
                              THEN f.subject_id END) AS ar_carriers,
          COUNT(DISTINCT CASE WHEN f.zygosity = 'Hemizygous'
                              THEN f.subject_id END) AS hemizygous_subjects,
          COUNT(DISTINCT CASE WHEN f.zygosity = 'Compound heterozygous'
                              THEN f.subject_id END) AS compound_het_subjects,
          COUNT(DISTINCT CASE WHEN f.zygosity = 'Homozygous'
                              THEN f.subject_id END) AS homozygous_subjects,
          COUNT(DISTINCT CASE WHEN gd.inheritance = 'AD' AND f.zygosity = 'Heterozygous'
                              THEN f.subject_id END) AS ad_het_subjects
        {join} {where}
    """.format(join=OBS_JOIN, where=where)) or {}

    grid = db.rows("""
        SELECT gd.inheritance, f.zygosity, COUNT(*) AS n,
               COUNT(DISTINCT f.subject_id) AS subjects
        {join} {where}
        GROUP BY 1, 2
    """.format(join=OBS_JOIN, where=where))

    matrix = {i: {z: 0 for z in ZYGOSITY_ORDER} for i in INHERITANCE_ORDER}
    for r in grid:
        if r["inheritance"] in matrix and r["zygosity"] in matrix[r["inheritance"]]:
            matrix[r["inheritance"]][r["zygosity"]] = int(r["subjects"])

    return {
        "kpis": {k: int(v or 0) for k, v in kpis.items()},
        "heatmap": {
            "rows": INHERITANCE_ORDER,
            "cols": ZYGOSITY_ORDER,
            "matrix": matrix,
            "max": max([max(r.values()) for r in matrix.values()] + [1]),
        },
        "by_gene": by_gene(criteria),
        "compound_het": compound_het_pairs(criteria),
        "hemizygous": hemizygous_rows(criteria),
        "compound_het_footer": (
            "Flagged for trio or long-read phase confirmation before reporting "
            "as biallelic."),
        "hemizygous_footer": (
            "Hemizygous is derived from chromosome and subject sex. ZYGOSITY in "
            "the source file carries Het/Hom only; reading it directly "
            "misclassifies X-linked males as carriers."),
    }


def by_gene(criteria: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Per-gene zygosity split with the §G04 interpretation flags."""
    obs = _obs_where(criteria)
    rows = db.rows("""
        SELECT f.gene_symbol AS gene, gd.inheritance, gd.condition, gd.validity,
               COUNT(DISTINCT CASE WHEN f.zygosity = 'Heterozygous'
                                   THEN f.subject_id END) AS het,
               COUNT(DISTINCT CASE WHEN f.zygosity = 'Homozygous'
                                   THEN f.subject_id END) AS hom,
               COUNT(DISTINCT CASE WHEN f.zygosity = 'Compound heterozygous'
                                   THEN f.subject_id END) AS chet,
               COUNT(DISTINCT CASE WHEN f.zygosity = 'Hemizygous'
                                   THEN f.subject_id END) AS hemi,
               COUNT(DISTINCT f.subject_id) AS subjects
        {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {plp}
        GROUP BY 1, 2, 3, 4
        ORDER BY subjects DESC
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL))

    out = []
    for r in rows:
        het, hom = int(r["het"]), int(r["hom"])
        chet, hemi = int(r["chet"]), int(r["hemi"])
        inh = r["inheritance"]
        biallelic = hom + chet

        # Spec §G04 interpretation flags.
        if inh in ("AR", "AR/AD") and biallelic:
            flag, tone = "biallelic — consistent with affected", "red"
        elif inh in ("AR", "AR/AD"):
            flag, tone = "heterozygous carriers only", "grey"
        elif inh.startswith("XL") and hemi:
            flag, tone = "hemizygous male", "amber"
        elif inh == "AD":
            flag, tone = "dominant — het sufficient", "blue"
        else:
            flag, tone = "", "grey"

        out.append({
            "gene": r["gene"], "condition": r["condition"], "inheritance": inh,
            "validity": r["validity"],
            # Deliberately four separate columns. Never summed (§G04 done-when).
            "het": het, "hom": hom, "compound_het": chet, "hemizygous": hemi,
            "biallelic": biallelic, "subjects": int(r["subjects"]),
            "flag": flag, "flag_tone": tone,
        })
    return out


def compound_het_pairs(criteria: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Subjects with >=2 distinct P/LP variants in one gene.

    Phase is unconfirmed — two variants in one gene may be on the same allele
    (in cis), which is not biallelic disease. The table exists to route these
    to trio or long-read confirmation, not to report them as solved.
    """
    obs = _obs_where(criteria)
    rows = db.rows("""
        SELECT f.subject_id, f.gene_symbol AS gene, gd.condition, gd.inheritance,
               s.sex, s.affected_status, s.family_id,
               STRING_AGG(DISTINCT COALESCE(f.hgvs_p, f.hgvs_c), ' + ') AS variants,
               COUNT(DISTINCT f.variant_key) AS n_variants,
               STRING_AGG(DISTINCT i.segregation, '; ') AS segregation
        {join}
        JOIN subject s ON s.subject_id = f.subject_id
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {plp}
          AND f.zygosity = 'Compound heterozygous'
        GROUP BY 1, 2, 3, 4, 5, 6, 7
        HAVING COUNT(DISTINCT f.variant_key) >= 2
        ORDER BY f.gene_symbol, f.subject_id
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL))
    for r in rows:
        r["phase"] = "Unconfirmed"
    return rows


def hemizygous_rows(criteria: Dict[str, Any]) -> List[Dict[str, Any]]:
    """X/Y-linked calls in male subjects, showing what the file said vs what
    was derived — the §E04 trap made visible rather than merely avoided."""
    obs = _obs_where(criteria)
    return db.rows("""
        SELECT f.subject_id, f.gene_symbol AS gene, f.chrom, s.sex,
               gd.inheritance, gd.condition, f.hgvs_p, f.hgvs_c,
               f.zygosity_raw AS file_says, f.zygosity AS derived,
               i.classification, s.affected_status
        {join}
        JOIN subject s ON s.subject_id = f.subject_id
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {plp}
          AND f.zygosity = 'Hemizygous'
        ORDER BY f.gene_symbol, f.subject_id
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL))
