"""Gene-specific denominator service — spec §3.7, §3.4, C01.

    denominator(gene) = COUNT(runs in cohort WHERE run_scope covers gene)
    frequency(gene)   = altered / denominator(gene)     -- NEVER / cohort_size

This is the single most important module in the tool. From §E04:

    Trap: using cohort size as denominator
    Consequence: frequencies understated, plausible, quoted
    Guard: gene-specific denominator service

"Plausible" is the dangerous word. A rate computed on cohort size is not
obviously wrong — it is just quietly, consistently too low, and it gets quoted.

Germline unit of analysis is the SUBJECT, not the sample (spec Part B), so the
denominator counts distinct subjects whose run assays and reports the gene.

Three buckets, and the distinction between the second and third is the point:

    assayed      a reportable scope row whose confidence is declared or inferred
    provisional  a reportable scope row whose confidence is 'unknown' — we know
                 the gene is in scope but not how the scope was established
                 (typically a test code whose registry entry is incomplete).
                 Counted in numerators, flagged in denominators.
    not assayed  no reportable scope row at all. NOT wild-type, simply
                 unobserved — and non-reportable rows land here too, because a
                 gene the test does not report is a gene the cohort has no
                 result for.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from .. import config, db
from ..reference.genes import GENE_DISEASE, GENE_LIST

KNOWN = ("declared", "inferred")


def gene_denominators() -> Dict[str, Dict[str, Any]]:
    """Per-gene denominator over the resolved cohort. Keyed by gene symbol."""
    cohort_n = int(db.scalar("SELECT COUNT(*) FROM cohort_subject", default=0))

    rows = db.rows("""
        SELECT rs.gene_symbol AS gene,
               COUNT(DISTINCT CASE WHEN rs.confidence IN ('declared','inferred')
                                   THEN cr.subject_id END) AS assayed,
               COUNT(DISTINCT CASE WHEN rs.confidence = 'unknown'
                                   THEN cr.subject_id END) AS provisional,
               COUNT(DISTINCT CASE WHEN rs.confidence = 'declared'
                                   THEN cr.subject_id END) AS declared,
               COUNT(DISTINCT CASE WHEN rs.confidence = 'inferred'
                                   THEN cr.subject_id END) AS inferred
        FROM cohort_run cr
        JOIN run_scope rs ON rs.run_id = cr.run_id
        WHERE rs.reportable
        GROUP BY rs.gene_symbol
    """)

    out: Dict[str, Dict[str, Any]] = {}
    for g in GENE_LIST:
        out[g] = {"gene": g, "assayed": 0, "provisional": 0, "declared": 0,
                  "inferred": 0, "not_assayed": cohort_n, "cohort_n": cohort_n}
    for r in rows:
        g = r["gene"]
        if g not in out:
            continue
        assayed = int(r["assayed"] or 0)
        provisional = int(r["provisional"] or 0)
        out[g].update({
            "assayed": assayed,
            "provisional": provisional,
            "declared": int(r["declared"] or 0),
            "inferred": int(r["inferred"] or 0),
            "not_assayed": max(0, cohort_n - assayed - provisional),
        })

    for g, d in out.items():
        d.update(_annotate(d, cohort_n))
    return out


def _annotate(d: Dict[str, Any], cohort_n: int) -> Dict[str, Any]:
    assayed = d["assayed"]
    provisional = d["provisional"]
    total = assayed + provisional
    # Spec §3.7:
    #   denominator_warn = denominator < 30
    #                      OR unknown_scope / (denominator + unknown_scope) > 0.25
    low_n = assayed < config.DENOM_MIN
    high_prov = (provisional / total) > config.PROVISIONAL_MAX if total else False
    reasons: List[str] = []
    if low_n:
        reasons.append("denominator below {}".format(config.DENOM_MIN))
    if high_prov:
        reasons.append("more than {:.0%} of the denominator has provisional scope"
                       .format(config.PROVISIONAL_MAX))
    return {
        "warn": low_n or high_prov,
        "warn_low_n": low_n,
        "warn_provisional": high_prov,
        "warn_reason": " and ".join(reasons),
        "coverage_pct": round(100.0 * total / cohort_n, 1) if cohort_n else 0.0,
        # 'mixed' rather than collapsing to 'inferred' the moment one subject
        # arrives by inference: 99 declared and 1 inferred is a different claim
        # from 100 inferred, and the inspector is where that distinction is
        # supposed to be visible.
        "confidence": ("declared" if d["declared"] and not d["inferred"]
                       else "inferred" if d["inferred"] and not d["declared"]
                       else "mixed" if d["declared"] and d["inferred"]
                       else "unknown"),
    }


def rate(numerator: int, denominator: int) -> Dict[str, Any]:
    """Every percentage in the UI renders with its counts (spec §3.7).

    Returning the counts alongside the percentage is what makes that structural
    rather than a convention a screen can forget. `label` is the string the UI
    shows: "32 / 250 (12.8%)", never a bare 13%.
    """
    pct = (100.0 * numerator / denominator) if denominator else 0.0
    return {
        "n": numerator, "d": denominator, "pct": round(pct, 1),
        "label": "{} / {}".format(numerator, denominator),
        "display": "{} / {} ({}%)".format(numerator, denominator, round(pct, 1))
        if denominator else "{} / 0 (—)".format(numerator),
    }


def coverage_inspector() -> Dict[str, Any]:
    """C01 · Denominator / coverage inspector.

    "Prove the numbers are honest. Answers 'where does 32 out of 250 come
     from?' in one click."
    """
    cohort_n = int(db.scalar("SELECT COUNT(*) FROM cohort_subject", default=0))
    dens = gene_denominators()

    observed = {
        r["gene"]: int(r["n"]) for r in db.rows("""
            SELECT f.gene_symbol AS gene, COUNT(DISTINCT f.subject_id) AS n
            FROM finding f
            WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
            GROUP BY f.gene_symbol
        """)
    }

    genes = []
    for g in GENE_LIST:
        d = dens[g]
        rec = GENE_DISEASE[g]
        genes.append({
            "gene": g, "condition": rec.condition, "inheritance": rec.inheritance,
            "assayed": d["assayed"], "provisional": d["provisional"],
            "not_assayed": d["not_assayed"], "declared": d["declared"],
            "inferred": d["inferred"], "coverage_pct": d["coverage_pct"],
            "confidence": d["confidence"], "warn": d["warn"],
            "warn_reason": d["warn_reason"],
            "subjects_with_finding": observed.get(g, 0),
        })
    genes.sort(key=lambda r: (-r["coverage_pct"], r["gene"]))

    thresholds = db.meta_get("scope_thresholds")
    return {
        "cohort_n": cohort_n,
        "genes": genes,
        "explain": ("Cohort has {} subjects. Each gene's denominator is the subset "
                    "whose test code actually assays and reports it. This is why a "
                    "single cohort has many denominators.".format(cohort_n)),
        "thresholds_json": thresholds,
        "panels": db.rows("""
            SELECT COALESCE(r.implied_panel_id, '(unresolved)') AS implied_panel,
                   COALESCE(r.test_code, '(no test code)') AS test_code,
                   COUNT(*) AS runs
            FROM run r WHERE r.run_id IN (SELECT run_id FROM cohort_run)
            GROUP BY 1, 2 ORDER BY runs DESC
        """),
    }


def small_cell(n: int) -> str:
    """Spec §3.7: any count < 5 renders as '<5' in exports and client-visible
    views. A single rare-condition family in a named ancestry and city is
    re-identifying (spec §G09)."""
    if n <= 0:
        return "0"
    return "<{}".format(config.SMALL_CELL_THRESHOLD) if n < config.SMALL_CELL_THRESHOLD else str(n)
