"""Diagnostic yield and carrier rates — spec §G02, §G03.

Diagnostic yield is "the number clients will scrutinise hardest", and the whole
difficulty is in the fourth condition:

    solved(subject) = EXISTS observation WHERE
                        classification in (Pathogenic, Likely pathogenic)
                    AND gene_disease_validity in (Definitive, Strong)
                    AND reportable = true
                    AND gene is phenotype-relevant for the subject's indication

Drop the relevance gate and a broad exome counts any P/LP as diagnostic — which
pushed yield to 81% in the spec authors' testing. A P/LP in an unrelated organ
system is a secondary finding, not a diagnosis. Calibrated expectation: ~41.5%
overall, 40-60% targeted panels, 23-36% exome/genome.

Carrier rate uses the gene-specific denominator from denominator.py, never the
cohort size.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from .. import config, db
from ..reference.genes import GENE_DISEASE
from . import denominator as den
from .cohort import OBS_JOIN, observation_filter_sql

PLP_SQL = "i.classification IN ('Pathogenic','Likely pathogenic')"
VUS_SQL = "i.classification = 'Uncertain significance'"

# What makes a finding SECONDARY (spec §G09). Shared by §G02's KPI and §G09 so
# the two can never disagree.
#
# Three conditions, and each one is load-bearing:
#   1. the gene is on the ACMG SF v3.3 list
#   2. the gene is NOT relevant to the subject's own indication — a P/LP in
#      MYBPC3 found on a cardiomyopathy panel is the PRIMARY result, not an
#      incidental. "Secondary" is defined relative to why the test was ordered.
#   3. the subject's test can actually report secondary findings
# Consent is applied on top of this, at subject level.
SECONDARY_FINDING_SQL = """
    gd.gene_sets LIKE '%SF%'
    AND NOT EXISTS (
        SELECT 1 FROM indication_gene ig
        WHERE ig.indication = s.indication AND ig.gene_symbol = f.gene_symbol)
    AND EXISTS (
        SELECT 1 FROM run r2 JOIN test_code tc2 ON tc2.code = r2.test_code
        WHERE r2.subject_id = f.subject_id AND tc2.sf_capable)
"""
SF_CONSENT_SQL = "s.consent_class = 'Full research + secondary findings'"


def _obs_where(criteria: Dict[str, Any]) -> str:
    """Observation-level scoping shared by every germline analytic.

    Only `reportable_only` is carried through. The other genomic criteria
    already selected cohort MEMBERSHIP in §G01 step 8; re-applying them to the
    observations would hide a subject's other findings and make per-subject
    result classification disagree with the subject list (spec §G12 done-when).
    """
    params: List[Any] = []
    sql = observation_filter_sql(
        {"reportable_only": bool(criteria.get("reportable_only"))}, params)
    assert not params, "reportable_only must not parameterise"
    return sql


# ------------------------------------------------------------------- yield ---
def diagnostic_yield(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """§G02 diagnostic yield, with the phenotype-relevance gate applied."""
    obs = _obs_where(criteria)
    total = int(db.scalar("SELECT COUNT(*) FROM cohort_subject", default=0))

    solved = int(db.scalar("""
        SELECT COUNT(DISTINCT f.subject_id) {join}
        JOIN subject s ON s.subject_id = f.subject_id
        JOIN indication_gene ig ON ig.indication = s.indication
                               AND ig.gene_symbol = f.gene_symbol
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs}
          AND {plp}
          AND gd.validity IN ('Definitive','Strong')
          AND i.reportable
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL), default=0))

    # VUS-only: not solved, but carrying at least one uncertain result.
    vus_only = int(db.scalar("""
        SELECT COUNT(DISTINCT cs.subject_id) FROM cohort_subject cs
        WHERE cs.subject_id NOT IN ({solved_ids})
          AND cs.subject_id IN (
            SELECT f.subject_id {join}
            WHERE {obs} AND {vus})
    """.format(solved_ids=_solved_subject_sql(obs), join=OBS_JOIN, obs=obs, vus=VUS_SQL),
        default=0))

    return {
        "solved": solved,
        "vus_only": vus_only,
        "negative": max(0, total - solved - vus_only),
        "total": total,
        "rate": den.rate(solved, total),
        "vus_only_rate": den.rate(vus_only, total),
        "caveat": ("Yield is computed against the subject's stated indication. A "
                   "P/LP finding in an unrelated organ system is counted as a "
                   "secondary finding, not a diagnosis."),
    }


def _solved_subject_sql(obs: str) -> str:
    return """
        SELECT f.subject_id {join}
        JOIN subject s ON s.subject_id = f.subject_id
        JOIN indication_gene ig ON ig.indication = s.indication
                               AND ig.gene_symbol = f.gene_symbol
        WHERE {obs} AND {plp}
          AND gd.validity IN ('Definitive','Strong') AND i.reportable
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL)


def subject_results(criteria: Dict[str, Any]) -> Dict[str, str]:
    """subject_id -> 'P/LP' | 'VUS-only' | 'Negative'.

    §G12 done-when: "The result classification matches the yield computation in
    G02 exactly." It does, because both read this one definition.
    """
    obs = _obs_where(criteria)
    solved = {r["subject_id"] for r in db.rows("""
        SELECT DISTINCT subject_id FROM ({}) t
        WHERE subject_id IN (SELECT subject_id FROM cohort_subject)
    """.format(_solved_subject_sql(obs)))}
    vus = {r["subject_id"] for r in db.rows("""
        SELECT DISTINCT f.subject_id {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {vus}
    """.format(join=OBS_JOIN, obs=obs, vus=VUS_SQL))}

    out: Dict[str, str] = {}
    for r in db.rows("SELECT subject_id FROM cohort_subject"):
        sid = r["subject_id"]
        out[sid] = "P/LP" if sid in solved else "VUS-only" if sid in vus else "Negative"
    return out


def yield_by(dimension: str, criteria: Dict[str, Any]) -> List[Dict[str, Any]]:
    """§G03 yield stratified by test code or self-reported ancestry."""
    col = {"test_code": "r.test_code", "ancestry": "s.ancestry",
           "indication": "s.indication"}.get(dimension)
    if not col:
        raise ValueError("unsupported yield dimension: {}".format(dimension))
    obs = _obs_where(criteria)

    rows = db.rows("""
        WITH solved AS ({solved}),
        base AS (
            SELECT cs.subject_id, {col} AS bucket
            FROM cohort_subject cs
            JOIN subject s ON s.subject_id = cs.subject_id
            JOIN run r ON r.subject_id = cs.subject_id)
        SELECT COALESCE(CAST(bucket AS VARCHAR), '(unspecified)') AS bucket,
               COUNT(DISTINCT base.subject_id) AS subjects,
               COUNT(DISTINCT CASE WHEN base.subject_id IN (SELECT subject_id FROM solved)
                                   THEN base.subject_id END) AS solved
        FROM base GROUP BY 1 ORDER BY subjects DESC
    """.format(solved=_solved_subject_sql(obs), col=col))

    return [{
        "bucket": r["bucket"],
        "subjects": int(r["subjects"]),
        "solved": int(r["solved"]),
        "rate": den.rate(int(r["solved"]), int(r["subjects"])),
    } for r in rows]


# ----------------------------------------------------------- carrier rates ---
def carrier_table(criteria: Dict[str, Any], plp_only: bool = True) -> List[Dict[str, Any]]:
    """§G03 carrier-rate table. One row per gene with any finding.

        carrier_rate = P/LP subjects / subjects whose run assays and reports
                       the gene
    """
    obs = _obs_where(criteria)
    dens = den.gene_denominators()

    rows = db.rows("""
        SELECT f.gene_symbol AS gene,
               COUNT(DISTINCT CASE WHEN {plp} THEN f.subject_id END) AS plp_subjects,
               COUNT(DISTINCT CASE WHEN {plp} THEN f.family_id END) AS plp_families,
               COUNT(DISTINCT f.subject_id) AS any_subjects,
               COUNT(DISTINCT CASE WHEN {plp} THEN f.variant_key END) AS unique_plp_variants,
               COUNT(DISTINCT f.variant_key) AS unique_variants,
               COUNT(*) FILTER (WHERE {plp}) AS plp_observations,
               COUNT(*) FILTER (WHERE {vus}) AS vus_observations,
               COUNT(*) FILTER (WHERE f.var_class = 'CNV') AS cnv_observations
        {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject) AND {obs}
        GROUP BY f.gene_symbol
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL, vus=VUS_SQL))

    out: List[Dict[str, Any]] = []
    for r in rows:
        g = r["gene"]
        d = dens.get(g, {})
        rec = GENE_DISEASE.get(g)
        plp_subjects = int(r["plp_subjects"] or 0)
        if plp_only and not plp_subjects:
            continue
        assayed = int(d.get("assayed", 0))
        out.append({
            "gene": g,
            "condition": rec.condition if rec else "",
            "mondo": rec.mondo if rec else "",
            "inheritance": rec.inheritance if rec else "",
            "validity": rec.validity if rec else "",
            "penetrance": rec.penetrance if rec else "",
            "plp_subjects": plp_subjects,
            "plp_families": int(r["plp_families"] or 0),
            "any_subjects": int(r["any_subjects"] or 0),
            "unique_variants": int(r["unique_variants"] or 0),
            "unique_plp_variants": int(r["unique_plp_variants"] or 0),
            "plp_observations": int(r["plp_observations"] or 0),
            "vus_observations": int(r["vus_observations"] or 0),
            # CNV coverage cannot be inferred — absence is uninformative, so
            # these are reported as counts with no rate (spec §3.4, E03.4).
            "cnv_observations": int(r["cnv_observations"] or 0),
            "assayed": assayed,
            "provisional": int(d.get("provisional", 0)),
            "denominator_label": "{}{}".format(
                assayed, " +{}?".format(d["provisional"]) if d.get("provisional") else ""),
            "carrier_rate": den.rate(plp_subjects, assayed),
            "any_rate": den.rate(int(r["any_subjects"] or 0), assayed),
            "warn": bool(d.get("warn")),
            "warn_reason": d.get("warn_reason", ""),
            # Families < subjects means the same allele is being counted more
            # than once from one pedigree — independence is violated (§G07).
            "family_clustered": int(r["plp_families"] or 0) < plp_subjects,
        })
    out.sort(key=lambda r: (-r["carrier_rate"]["pct"], -r["plp_subjects"], r["gene"]))
    return out


WARN_LEGEND = ("⚠ = denominator below {}, or >{:.0%} of the denominator has "
               "provisional scope. Rate shown but not fit for external quotation."
               ).format(config.DENOM_MIN, config.PROVISIONAL_MAX)

ANCESTRY_CAVEAT = ("Differences reflect referral pattern and panel choice as much "
                   "as population genetics. Not a population frequency estimate.")


# --------------------------------------------------------------- dashboard ---
def dashboard(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """§G02 · L1 subject-level summary."""
    obs = _obs_where(criteria)
    total = int(db.scalar("SELECT COUNT(*) FROM cohort_subject", default=0))
    y = diagnostic_yield(criteria)

    fam = db.row("""
        WITH f AS (SELECT family_id, COUNT(*) n FROM cohort_subject GROUP BY family_id)
        SELECT COUNT(*) AS families,
               COALESCE(SUM(CASE WHEN n > 1 THEN 1 ELSE 0 END), 0) AS multi
        FROM f
    """) or {}

    counts = db.row("""
        SELECT COUNT(*) FILTER (WHERE {plp}) AS plp_findings,
               COUNT(DISTINCT CASE WHEN {plp} THEN f.variant_key END) AS plp_unique,
               COUNT(*) FILTER (WHERE {vus}) AS vus_findings,
               COUNT(*) AS all_findings
        {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject) AND {obs}
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL, vus=VUS_SQL)) or {}

    # Secondary findings: consented subjects only, and the consent test is on
    # the cohort query, not on display (spec §G09).
    sf = db.row("""
        SELECT COUNT(DISTINCT f.subject_id) AS subjects,
               COUNT(*) AS findings
        {join}
        JOIN subject s ON s.subject_id = f.subject_id
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {plp}
          AND {secondary}
          AND {consent}
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL,
               secondary=SECONDARY_FINDING_SQL, consent=SF_CONSENT_SQL)) or {}

    return {
        "kpis": {
            "subjects": total,
            "families": int(fam.get("families") or 0),
            "multi_member_families": int(fam.get("multi") or 0),
            "yield": y["rate"],
            "vus_only": y["vus_only_rate"],
            "plp_findings": int(counts.get("plp_findings") or 0),
            "plp_unique_variants": int(counts.get("plp_unique") or 0),
            "vus_findings": int(counts.get("vus_findings") or 0),
            "all_findings": int(counts.get("all_findings") or 0),
            "sf_subjects": int(sf.get("subjects") or 0),
            "sf_findings": int(sf.get("findings") or 0),
        },
        "yield": y,
        "outcome_stack": [
            {"label": "P/LP — diagnostic", "value": y["solved"], "tone": "plp"},
            {"label": "VUS only", "value": y["vus_only"], "tone": "vus"},
            {"label": "Negative", "value": y["negative"], "tone": "neg"},
        ],
        "zygosity_stack": db.rows("""
            SELECT f.zygosity AS label, COUNT(*) AS value
            {join}
            WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
              AND {obs} AND {plp}
            GROUP BY 1 ORDER BY value DESC
        """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL)),
        "consent_stack": db.rows("""
            SELECT s.consent_class AS label, COUNT(*) AS value
            FROM cohort_subject cs JOIN subject s USING (subject_id)
            GROUP BY 1 ORDER BY value DESC
        """),
        "age_stack": db.rows("""
            SELECT s.age_bucket AS label, COUNT(*) AS value
            FROM cohort_subject cs JOIN subject s USING (subject_id)
            GROUP BY 1 ORDER BY label
        """),
        "ancestry_stack": db.rows("""
            SELECT s.ancestry AS label, COUNT(*) AS value
            FROM cohort_subject cs JOIN subject s USING (subject_id)
            GROUP BY 1 ORDER BY value DESC
        """),
        "yield_by_indication": yield_by("indication", criteria),
        "accrual": db.rows("""
            SELECT strftime(sm.collection_date, '%Y-%m') AS month, COUNT(*) AS value
            FROM cohort_subject cs
            JOIN run r ON r.subject_id = cs.subject_id
            JOIN sample sm ON sm.sample_id = r.sample_id
            GROUP BY 1 ORDER BY 1 DESC LIMIT 24
        """)[::-1],
        "carrier_top": carrier_table(criteria)[:12],
    }
