"""Germline analysis modules G05-G09.

    G05  Gene-disease association   — yield inflation from weak-validity genes
    G06  Phenotype (HPO) analysis   — hypothesis generation, never reclassification
    G07  Population frequency       — founder candidates and QA signals
    G08  VUS inventory              — the curation work queue
    G09  Secondary findings         — ACMG SF v3.3, consent-enforced at query level
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from .. import config, db
from ..reference.genes import GENE_DISEASE, GENE_LIST, SF_GENES
from .carrier import (PLP_SQL, SECONDARY_FINDING_SQL, SF_CONSENT_SQL,
                      VUS_SQL, _obs_where)
from .cohort import OBS_JOIN
from . import denominator as den

# Anchor date for staleness. The store's most recent collection, not wall-clock
# today, so a saved cohort's priority scores are reproducible (spec E03.6).
ANCHOR = "(SELECT MAX(collection_date) FROM sample)"


# ============================================================ G05 gene-disease
def gene_disease(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """§G05. "Prevents yield inflation from weak-validity genes."

        any_plp     = DISTINCT subjects with >=1 P/LP (all validity levels)
        established = DISTINCT subjects with >=1 P/LP in Definitive|Strong genes
        inflation % = (any_plp - established) / established

    "That inflation figure is the point of the screen — it shows what
     unfiltered yield would have claimed."
    """
    obs = _obs_where(criteria)
    base = """WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
              AND {obs} AND {plp}""".format(obs=obs, plp=PLP_SQL)

    counts = db.row("""
        SELECT COUNT(DISTINCT f.subject_id) AS any_plp,
               COUNT(DISTINCT CASE WHEN gd.validity IN ('Definitive','Strong')
                                   THEN f.subject_id END) AS established
        {join} {base}
    """.format(join=OBS_JOIN, base=base)) or {}

    any_plp = int(counts.get("any_plp") or 0)
    established = int(counts.get("established") or 0)
    inflation = round(100.0 * (any_plp - established) / established, 1) if established else 0.0

    conditions = db.rows("""
        SELECT gd.condition, gd.mondo_id, gd.inheritance, gd.validity, gd.penetrance,
               STRING_AGG(DISTINCT f.gene_symbol, ', ') AS genes,
               COUNT(DISTINCT f.subject_id) AS subjects,
               COUNT(DISTINCT f.family_id) AS families,
               COUNT(DISTINCT f.variant_key) AS variants
        {join} {base}
        GROUP BY 1, 2, 3, 4, 5
        ORDER BY subjects DESC
    """.format(join=OBS_JOIN, base=base))

    validity_dist = db.rows("""
        SELECT gd.validity AS label, COUNT(DISTINCT f.subject_id) AS value
        {join} {base} GROUP BY 1
    """.format(join=OBS_JOIN, base=base))
    penetrance_dist = db.rows("""
        SELECT gd.penetrance AS label, COUNT(DISTINCT f.subject_id) AS value
        {join} {base} GROUP BY 1
    """.format(join=OBS_JOIN, base=base))

    # Full curated reference table, annotated with this cohort's observations.
    observed = {r["gene"]: r for r in db.rows("""
        SELECT f.gene_symbol AS gene,
               COUNT(DISTINCT f.subject_id) AS plp_subjects,
               COUNT(DISTINCT f.variant_key) AS plp_variants
        {join} {base} GROUP BY 1
    """.format(join=OBS_JOIN, base=base))}

    reference = []
    for g in GENE_LIST:
        rec = GENE_DISEASE[g]
        hit = observed.get(g, {})
        reference.append({
            "gene": g, "condition": rec.condition, "mondo": rec.mondo,
            "inheritance": rec.inheritance, "validity": rec.validity,
            "penetrance": rec.penetrance, "gene_sets": list(rec.sets),
            "plp_subjects": int(hit.get("plp_subjects") or 0),
            "plp_variants": int(hit.get("plp_variants") or 0),
            "counts_toward_yield": rec.validity in ("Definitive", "Strong"),
        })
    reference.sort(key=lambda r: (-r["plp_subjects"], r["gene"]))

    return {
        "kpis": {
            "any_plp": any_plp,
            "established": established,
            "weak_only": max(0, any_plp - established),
            "inflation_pct": inflation,
            "conditions": len(conditions),
        },
        "inflation_label": "{} → {} subjects".format(any_plp, established),
        "conditions": conditions,
        "validity_dist": validity_dist,
        "penetrance_dist": penetrance_dist,
        "reference": reference,
        "guidance": ("Moderate-validity, low-penetrance findings belong in a separate "
                     "reporting category from definitive high-penetrance genes. Pooling "
                     "them is what makes a yield number unusable for counselling."),
    }


# =============================================================== G06 phenotype
def phenotype(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """§G06. Gene x HPO co-occurrence and a phenotype-informed review queue.

    Both footers are mandatory and non-negotiable: this is hypothesis
    generation, and the queue never mutates a classification.
    """
    obs = _obs_where(criteria)

    terms = db.rows("""
        SELECT s.phenotype_hpo AS term, COUNT(*) AS subjects
        FROM cohort_subject cs JOIN subject s USING (subject_id)
        WHERE s.phenotype_hpo IS NOT NULL
        GROUP BY 1 ORDER BY subjects DESC LIMIT 10
    """)
    term_list = [t["term"] for t in terms]

    genes = [r["gene"] for r in db.rows("""
        SELECT f.gene_symbol AS gene, COUNT(DISTINCT f.subject_id) AS n
        {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {plp}
        GROUP BY 1 ORDER BY n DESC LIMIT 14
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL))]

    cells = db.rows("""
        SELECT f.gene_symbol AS gene, s.phenotype_hpo AS term,
               COUNT(DISTINCT f.subject_id) AS n
        {join}
        JOIN subject s ON s.subject_id = f.subject_id
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {plp}
        GROUP BY 1, 2
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL))

    matrix = {g: {t: 0 for t in term_list} for g in genes}
    for c in cells:
        if c["gene"] in matrix and c["term"] in matrix[c["gene"]]:
            matrix[c["gene"]][c["term"]] = int(c["n"])

    # Review queue: VUS observations where the gene's condition text overlaps
    # the subject's HPO term — candidate PP4 evidence, surfaced for a curator.
    queue = db.rows("""
        SELECT f.finding_id, f.subject_id, f.gene_symbol AS gene,
               COALESCE(f.hgvs_p, f.hgvs_c) AS variant, gd.condition,
               s.phenotype_hpo AS term, s.indication, f.gnomad_af,
               i.acmg_codes, i.interpreted_at
        {join}
        JOIN subject s ON s.subject_id = f.subject_id
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {vus}
          AND EXISTS (
            SELECT 1 FROM indication_gene ig
            WHERE ig.indication = s.indication AND ig.gene_symbol = f.gene_symbol)
        ORDER BY f.gene_symbol, f.subject_id
        LIMIT ?
    """.format(join=OBS_JOIN, obs=obs, vus=VUS_SQL), [config.TABLE_RENDER_CAP])

    return {
        "terms": terms,
        "genes": genes,
        "matrix": matrix,
        "max": max([max(r.values()) for r in matrix.values()] + [1]) if matrix else 1,
        "queue": queue,
        "heatmap_footer": ("Sparse cells are expected — this is hypothesis generation, "
                           "not an association test."),
        "queue_footer": ("Phenotype match raises PP4 consideration. Surfaced for curator "
                         "review — no automatic reclassification."),
    }


# ========================================================= G07 population freq
def population_frequency(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """§G07. Internal allele frequency against gnomAD.

        internal_AF = carrier_subjects / (2 x assayed_subjects)  -- assumes het
        ratio       = internal_AF / gnomAD_AF
        rows        = variants seen in >=2 subjects, sorted by ratio DESC

    Finds founder candidates — and implausibly high internal frequencies, which
    are a QA signal, not biology.
    """
    obs = _obs_where(criteria)
    dens = den.gene_denominators()

    rows = db.rows("""
        SELECT f.variant_key, f.gene_symbol AS gene,
               MAX(COALESCE(f.hgvs_p, f.hgvs_c)) AS variant,
               MAX(f.consequence) AS consequence,
               MAX(i.classification) AS classification,
               MAX(f.gnomad_af) AS gnomad_af,
               MAX(f.gnomad_sas_af) AS gnomad_sas_af,
               MAX(f.ga100k_sas_af) AS ga100k_sas_af,
               MAX(f.clinvar_sig) AS clinvar_sig,
               MAX(gd.condition) AS condition,
               COUNT(DISTINCT f.subject_id) AS carriers,
               COUNT(DISTINCT f.family_id) AS families,
               STRING_AGG(DISTINCT s.ancestry, ', ') AS ancestries
        {join}
        JOIN subject s ON s.subject_id = f.subject_id
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject) AND {obs}
        GROUP BY f.variant_key, f.gene_symbol
        HAVING COUNT(DISTINCT f.subject_id) >= 2
    """.format(join=OBS_JOIN, obs=obs))

    dominant = db.row("""
        SELECT s.ancestry, COUNT(*) n FROM cohort_subject cs
        JOIN subject s USING (subject_id)
        GROUP BY 1 ORDER BY n DESC LIMIT 1
    """) or {"ancestry": "unspecified"}

    out: List[Dict[str, Any]] = []
    for r in rows:
        assayed = dens.get(r["gene"], {}).get("assayed", 0)
        carriers = int(r["carriers"])
        families = int(r["families"])
        # Assumes heterozygosity — stated in the footer, not hidden in code.
        internal_af = carriers / (2.0 * assayed) if assayed else None
        gnomad = r["gnomad_af"]
        ratio = (internal_af / gnomad) if (internal_af and gnomad) else None
        is_plp = r["classification"] in ("Pathogenic", "Likely pathogenic")

        flags: List[Dict[str, str]] = []
        if is_plp and gnomad and gnomad > 0.01:
            flags.append({"label": "AF too high for P/LP", "tone": "red",
                          "why": "classification review indicated"})
        if ratio and ratio > 100:
            flags.append({"label": "check for artefact", "tone": "amber",
                          "why": "internal AF >100x gnomAD"})
        if families < carriers:
            flags.append({"label": "family clustered", "tone": "grey",
                          "why": "violates independence — {} carriers in {} families"
                                 .format(carriers, families)})
        if is_plp and ratio and ratio > 5 and families == carriers:
            flags.append({"label": "founder · {}".format(r["ancestries"].split(",")[0].strip()),
                          "tone": "teal", "why": "recurrent, unrelated, enriched over gnomAD"})

        out.append({
            **r,
            "carriers": carriers, "families": families,
            "assayed": assayed,
            "internal_af": round(internal_af, 6) if internal_af else None,
            "ratio": round(ratio, 1) if ratio else None,
            "flags": flags,
            "is_plp": is_plp,
        })
    out.sort(key=lambda r: -(r["ratio"] or 0))

    return {
        "rows": out[:config.TABLE_RENDER_CAP],
        "total_rows": len(out),
        "kpis": {
            "recurrent_variants": len(out),
            "founder_candidates": sum(1 for r in out
                                      if any(f["label"].startswith("founder") for f in r["flags"])),
            "af_too_high": sum(1 for r in out
                               if any(f["label"] == "AF too high for P/LP" for f in r["flags"])),
            "family_clustered": sum(1 for r in out
                                    if any(f["label"] == "family clustered" for f in r["flags"])),
        },
        "banner": ("This cohort is {}-predominant and referral-selected. Internal AF is "
                   "not a population frequency estimate and must not be quoted as one."
                   .format(dominant["ancestry"])),
        "footer": ("internal_AF = carrier subjects / (2 × assayed subjects). The factor "
                   "of 2 assumes every carrier is heterozygous; homozygous and hemizygous "
                   "carriers make this an underestimate."),
        "sources": ("gnomAD_SAS_AF, GA100K_SAS_af (Indian 100K genomes) and "
                    "gnomADv2_AF_sas. For a South-Asian-predominant cohort GA100K is "
                    "the better comparator."),
    }


# ============================================================= G08 VUS inventory
def vus_inventory(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """§G08. Every uncertain result ranked by review priority.

        priority = evidence_delta x 0.6
                 + subjects_affected x 7
                 + MIN(age_days / 12, 30)

    CRITICAL SCOPING (spec §G08): restricted to the reviewable subset per §3.2.
    autoACMGPrediction returns 309,331 VUS for a single exome; unfiltered, this
    screen is unusable. `finding.reviewable` is set by the loader funnel.
    """
    obs = _obs_where(criteria)

    rows = db.rows("""
        SELECT f.variant_key, f.gene_symbol AS gene,
               MAX(COALESCE(f.hgvs_p, f.hgvs_c)) AS variant,
               MAX(f.consequence) AS consequence,
               MAX(gd.condition) AS condition,
               MAX(gd.validity) AS validity,
               COUNT(DISTINCT f.subject_id) AS subjects,
               COUNT(DISTINCT f.family_id) AS families,
               MAX(f.gnomad_af) AS gnomad_af,
               MAX(f.clinvar_sig) AS clinvar_sig,
               MAX(i.acmg_codes) AS acmg_codes,
               MAX(i.evidence_delta) AS evidence_delta,
               MIN(i.interpreted_at) AS oldest_interpretation,
               SUM(CASE WHEN i.curated_flag THEN 1 ELSE 0 END) AS curated,
               DATE_DIFF('day', MIN(i.interpreted_at), {anchor}) AS age_days
        {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {vus}
          AND f.reviewable
        GROUP BY f.variant_key, f.gene_symbol
    """.format(join=OBS_JOIN, obs=obs, vus=VUS_SQL, anchor=ANCHOR))

    out = []
    for r in rows:
        delta = float(r["evidence_delta"] or 0)
        subjects = int(r["subjects"])
        age_days = int(r["age_days"] or 0)
        # The three components stay visible per row so the score can be argued
        # with (spec §G08 done-when: "score components are visible per row").
        c_evidence = round(delta * 0.6, 1)
        c_subjects = subjects * 7
        c_age = round(min(age_days / 12.0, 30), 1)
        out.append({
            **r,
            "subjects": subjects,
            "families": int(r["families"]),
            "age_days": age_days,
            "curated": int(r["curated"] or 0),
            "evidence_delta": delta,
            "score": round(c_evidence + c_subjects + c_age, 1),
            "score_parts": {"evidence": c_evidence, "subjects": c_subjects, "age": c_age},
            "new_evidence": delta >= 40,
            "stale": age_days > 180,
        })
    out.sort(key=lambda r: -r["score"])

    by_gene = db.rows("""
        SELECT f.gene_symbol AS label, COUNT(DISTINCT f.variant_key) AS value
        {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {vus} AND f.reviewable
        GROUP BY 1 ORDER BY value DESC LIMIT 16
    """.format(join=OBS_JOIN, obs=obs, vus=VUS_SQL))

    return {
        "rows": out[:config.TABLE_RENDER_CAP],
        "total_rows": len(out),
        "kpis": {
            "unique_vus": len(out),
            "subjects_affected": int(db.scalar("""
                SELECT COUNT(DISTINCT f.subject_id) {join}
                WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
                  AND {obs} AND {vus} AND f.reviewable
            """.format(join=OBS_JOIN, obs=obs, vus=VUS_SQL), default=0)),
            "stale": sum(1 for r in out if r["stale"]),
            "new_evidence": sum(1 for r in out if r["new_evidence"]),
            "curated": sum(1 for r in out if r["curated"]),
            "high_priority": sum(1 for r in out if r["score"] >= 40),
        },
        "by_gene": by_gene,
        "age_histogram": _histogram([r["age_days"] for r in out], bins=12),
        "footer": ("Priority is a triage heuristic, not a classification. Nothing here "
                   "is reclassified without curator sign-off."),
        "scoping_note": ("Scoped to the reviewable subset: PASS filter, coding "
                         "consequence, gnomAD AF < {}. A raw autoACMG VUS count is "
                         "meaningless — it runs genome-wide."
                         .format(config.REVIEWABLE_MAX_GNOMAD_AF)),
    }


def _histogram(values: List[float], bins: int = 12) -> List[Dict[str, Any]]:
    if not values:
        return []
    lo, hi = min(values), max(values)
    width = (hi - lo) / bins or 1
    out = [{"lo": round(lo + i * width), "hi": round(lo + (i + 1) * width), "value": 0}
           for i in range(bins)]
    for v in values:
        idx = min(bins - 1, int((v - lo) / width))
        out[idx]["value"] += 1
    return out


# ======================================================== G09 secondary findings
def secondary_findings(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """§G09 · ACMG SF v3.3.

        eligible    = subjects with explicit secondary-findings consent
        sf_capable  = eligible subjects whose test can report SF
        sf_findings = P/LP in an ACMG SF v3.3 gene, eligible subjects only
        sf_rate     = DISTINCT subjects with SF / sf_capable

    CONSENT ENFORCEMENT IS THE DEFINING REQUIREMENT. Subjects who declined are
    excluded from the cohort QUERY, not filtered from display — see §G01 step 7
    in cohort.py, which is where the gate actually lives. This module computes
    over whatever the cohort contains and additionally re-applies the consent
    predicate, so a caller who reaches it without setting sf_only still cannot
    obtain a rate that includes a declining subject.
    """
    obs = _obs_where(criteria)
    consent_ok = SF_CONSENT_SQL

    eligible = int(db.scalar("""
        SELECT COUNT(*) FROM cohort_subject cs JOIN subject s USING (subject_id)
        WHERE {}
    """.format(consent_ok), default=0))

    sf_capable = int(db.scalar("""
        SELECT COUNT(DISTINCT cs.subject_id)
        FROM cohort_subject cs
        JOIN subject s USING (subject_id)
        JOIN run r ON r.subject_id = cs.subject_id
        LEFT JOIN test_code tc ON tc.code = r.test_code
        WHERE {} AND COALESCE(tc.sf_capable, FALSE)
    """.format(consent_ok), default=0))

    declined = int(db.scalar("""
        SELECT COUNT(*) FROM subject s
        WHERE s.consent_class IN ('Clinical only, secondary findings declined',
                                  'Clinical only, research declined')
    """, default=0))

    # Numerator population must be a SUBSET of the denominator population
    # (sf_capable), or the rate can exceed 100%. SECONDARY_FINDING_SQL carries
    # the sf_capable test on the subject's run, so the two agree by construction.
    sf_where = """
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject)
          AND {obs} AND {plp}
          AND {secondary}
          AND {consent}
    """.format(obs=obs, plp=PLP_SQL, secondary=SECONDARY_FINDING_SQL, consent=consent_ok)

    sf_subjects = int(db.scalar("""
        SELECT COUNT(DISTINCT f.subject_id) {join}
        JOIN subject s ON s.subject_id = f.subject_id {where}
    """.format(join=OBS_JOIN, where=sf_where), default=0))

    by_gene = db.rows("""
        SELECT f.gene_symbol AS gene, MAX(gd.condition) AS condition,
               MAX(gd.inheritance) AS inheritance, MAX(gd.penetrance) AS penetrance,
               COUNT(DISTINCT f.subject_id) AS subjects,
               COUNT(DISTINCT f.family_id) AS families,
               COUNT(DISTINCT f.variant_key) AS variants
        {join}
        JOIN subject s ON s.subject_id = f.subject_id {where}
        GROUP BY 1 ORDER BY subjects DESC
    """.format(join=OBS_JOIN, where=sf_where))
    for r in by_gene:
        # Small-cell suppression before the value leaves the service, not in
        # the template (spec §G09, §3.7).
        r["subjects_display"] = den.small_cell(int(r["subjects"]))
        r["families_display"] = den.small_cell(int(r["families"]))

    queue = db.rows("""
        SELECT f.subject_id, f.gene_symbol AS gene,
               COALESCE(f.hgvs_p, f.hgvs_c) AS variant, gd.condition,
               i.classification, f.zygosity, s.age, s.sex, s.indication,
               r.test_code, i.interpreted_at
        {join}
        JOIN subject s ON s.subject_id = f.subject_id
        JOIN run r ON r.subject_id = f.subject_id
        {where}
        ORDER BY f.gene_symbol, f.subject_id
        LIMIT ?
    """.format(join=OBS_JOIN, where=sf_where), [config.TABLE_RENDER_CAP])

    return {
        "kpis": {
            "eligible": eligible,
            "sf_capable": sf_capable,
            "sf_subjects": sf_subjects,
            "sf_genes": len(by_gene),
            "declined": declined,
            "rate": den.rate(sf_subjects, sf_capable),
        },
        "by_gene": by_gene,
        "queue": queue,
        "consent_posture": db.rows("""
            SELECT s.consent_class AS label, COUNT(*) AS value
            FROM cohort_subject cs JOIN subject s USING (subject_id)
            GROUP BY 1 ORDER BY value DESC
        """),
        "banner": ("{} subjects declined secondary findings and are excluded from this "
                   "analysis entirely — not merely suppressed in display. The cohort "
                   "cannot be constructed to include them.".format(declined)),
        "suppression_note": ("Counts below {} render as <{}. A single rare-condition "
                             "family in a named ancestry and city is re-identifying."
                             .format(config.SMALL_CELL_THRESHOLD,
                                     config.SMALL_CELL_THRESHOLD)),
        "sf_gene_count": len(SF_GENES),
    }
