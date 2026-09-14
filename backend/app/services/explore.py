"""Drill-path modules G10-G13 — the L4/L5/L2/L3 levels.

    G10  Gene analysis     (L4)  per-gene carrier rate, zygosity, ancestry split
    G11  Variant analysis  (L5)  every observation, with hand-off to interpretation
    G12  Subject list      (L2)  one row per subject, with the result classification
    G13  Assay list / QC   (L3)  provenance and QC per run

G11 contains NO interpretation logic. It reads and links out (spec §S05/§G11).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from .. import config, db
from ..reference.genes import GENE_DISEASE
from .carrier import PLP_SQL, VUS_SQL, _obs_where, subject_results
from .cohort import OBS_JOIN
from . import denominator as den


# ================================================================ G10 · gene ==
def gene_detail(gene: str, criteria: Dict[str, Any]) -> Dict[str, Any]:
    """§G10 · L4. Done-when: "the side-panel rate matches G03 exactly" — it does,
    because both call denominator.gene_denominators()."""
    obs = _obs_where(criteria)
    rec = GENE_DISEASE.get(gene)
    if not rec:
        return {"gene": gene, "error": "unknown gene"}

    dens = den.gene_denominators().get(gene, {})
    assayed = int(dens.get("assayed", 0))
    scope = """WHERE f.gene_symbol = ?
               AND f.subject_id IN (SELECT subject_id FROM cohort_subject)
               AND {obs}""".format(obs=obs)

    summary = db.row("""
        SELECT COUNT(DISTINCT CASE WHEN {plp} THEN f.subject_id END) AS plp_subjects,
               COUNT(DISTINCT CASE WHEN {plp} THEN f.family_id END) AS plp_families,
               COUNT(DISTINCT CASE WHEN {vus} THEN f.subject_id END) AS vus_subjects,
               COUNT(DISTINCT f.subject_id) AS any_subjects,
               COUNT(DISTINCT f.variant_key) AS unique_variants,
               COUNT(*) AS observations
        {join} {scope}
    """.format(join=OBS_JOIN, scope=scope, plp=PLP_SQL, vus=VUS_SQL), [gene]) or {}

    findings = db.rows("""
        SELECT f.variant_key, MAX(COALESCE(f.hgvs_p, f.hgvs_c)) AS variant,
               MAX(f.hgvs_c) AS hgvs_c, MAX(f.hgvs_p) AS hgvs_p,
               MAX(f.aa_pos) AS aa_pos, MAX(f.consequence) AS consequence,
               MAX(f.var_class) AS var_class,
               MAX(i.classification) AS classification,
               MAX(i.acmg_codes) AS acmg_codes,
               MAX(f.clinvar_sig) AS clinvar_sig, MAX(f.gnomad_af) AS gnomad_af,
               COUNT(DISTINCT f.subject_id) AS subjects,
               COUNT(DISTINCT f.family_id) AS families,
               STRING_AGG(DISTINCT f.zygosity, ', ') AS zygosities
        {join} {scope}
        GROUP BY f.variant_key
        ORDER BY subjects DESC, variant
    """.format(join=OBS_JOIN, scope=scope), [gene])

    classification_stack = db.rows("""
        SELECT i.classification AS label, COUNT(*) AS value
        {join} {scope} GROUP BY 1 ORDER BY value DESC
    """.format(join=OBS_JOIN, scope=scope), [gene])

    zygosity_stack = db.rows("""
        SELECT f.zygosity AS label, COUNT(*) AS value
        {join} {scope} AND {plp} GROUP BY 1 ORDER BY value DESC
    """.format(join=OBS_JOIN, scope=scope, plp=PLP_SQL), [gene])

    # Per-ancestry, per-gene denominators (spec §G10 done-when). The
    # denominator for each ancestry is that ancestry's subjects assayed FOR
    # THIS GENE, not that ancestry's subjects in the cohort.
    ancestry_rows = db.rows("""
        WITH assayed AS (
            SELECT s.ancestry, COUNT(DISTINCT cr.subject_id) AS n
            FROM cohort_run cr
            JOIN subject s ON s.subject_id = cr.subject_id
            JOIN run_scope rs ON rs.run_id = cr.run_id
            WHERE rs.gene_symbol = ? AND rs.reportable
              AND rs.confidence IN ('declared','inferred')
            GROUP BY 1),
        carriers AS (
            SELECT s.ancestry, COUNT(DISTINCT f.subject_id) AS n
            {join}
            JOIN subject s ON s.subject_id = f.subject_id
            {scope} AND {plp}
            GROUP BY 1)
        SELECT a.ancestry, a.n AS assayed, COALESCE(c.n, 0) AS carriers
        FROM assayed a LEFT JOIN carriers c USING (ancestry)
        ORDER BY a.n DESC
    """.format(join=OBS_JOIN, scope=scope, plp=PLP_SQL), [gene, gene])

    plp_subjects = int(summary.get("plp_subjects") or 0)
    return {
        "gene": gene,
        "gene_disease": {
            "condition": rec.condition, "mondo": rec.mondo,
            "inheritance": rec.inheritance, "validity": rec.validity,
            "penetrance": rec.penetrance, "gene_sets": list(rec.sets),
            "protein_len": rec.protein_len, "chrom": rec.chrom,
            "counts_toward_yield": rec.validity in ("Definitive", "Strong"),
        },
        "kpis": {
            "plp_subjects": plp_subjects,
            "plp_families": int(summary.get("plp_families") or 0),
            "vus_subjects": int(summary.get("vus_subjects") or 0),
            "any_subjects": int(summary.get("any_subjects") or 0),
            "unique_variants": int(summary.get("unique_variants") or 0),
            "observations": int(summary.get("observations") or 0),
            "carrier_rate": den.rate(plp_subjects, assayed),
        },
        "denominator": dens,
        "findings": findings,
        "classification_stack": classification_stack,
        "zygosity_stack": zygosity_stack,
        "ancestry": [{
            "ancestry": r["ancestry"],
            "rate": den.rate(int(r["carriers"]), int(r["assayed"])),
            "warn": int(r["assayed"]) < config.DENOM_MIN,
        } for r in ancestry_rows],
        "lollipop": _lollipop(gene, rec.protein_len, findings, rec.domains),
        "ancestry_caveat": ("Differences reflect referral pattern and panel choice as "
                            "much as population genetics. Not a population frequency "
                            "estimate."),
    }


def _lollipop(gene: str, protein_len: Optional[int], findings: List[Dict[str, Any]],
              domains) -> Optional[Dict[str, Any]]:
    """Grouped by (aa_pos, consequence). None when the gene has no protein
    model — §G10 edge case: hide the plot, keep the panel."""
    if not protein_len:
        return None
    groups: Dict[str, Dict[str, Any]] = {}
    for f in findings:
        if f["aa_pos"] is None or f["var_class"] in ("CNV", "SV"):
            continue
        key = "{}|{}".format(f["aa_pos"], f["consequence"])
        g = groups.setdefault(key, {
            "aa_pos": int(f["aa_pos"]), "consequence": f["consequence"],
            "count": 0, "plp": 0, "variant": f["variant"]})
        g["count"] += int(f["subjects"])
        if f["classification"] in ("Pathogenic", "Likely pathogenic"):
            g["plp"] += int(f["subjects"])
    pts = sorted(groups.values(), key=lambda g: g["aa_pos"])
    if not pts:
        return None
    mx = max(p["count"] for p in pts)
    for p in pts:
        p["x"] = p["aa_pos"] / protein_len
        p["label_it"] = p["count"] >= 0.55 * mx
    return {"protein_len": protein_len, "points": pts, "max": mx,
            "domains": [{"start": d[0], "end": d[1], "name": d[2]} for d in (domains or ())]}


def gene_pills(criteria: Dict[str, Any], limit: int = 16) -> List[Dict[str, Any]]:
    obs = _obs_where(criteria)
    return db.rows("""
        SELECT f.gene_symbol AS gene, COUNT(DISTINCT f.subject_id) AS subjects
        {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject) AND {obs}
        GROUP BY 1 ORDER BY subjects DESC, gene LIMIT ?
    """.format(join=OBS_JOIN, obs=obs), [limit])


# ============================================================= G11 · variants ==
def variant_list(criteria: Dict[str, Any], gene: Optional[str] = None,
                 limit: Optional[int] = None) -> Dict[str, Any]:
    """§G11 · L5. Every observation. No interpretation logic — read and link out."""
    obs = _obs_where(criteria)
    limit = limit or config.TABLE_RENDER_CAP
    params: List[Any] = []
    gene_pred = ""
    if gene:
        gene_pred = "AND f.gene_symbol = ?"
        params.append(gene)

    total = int(db.scalar("""
        SELECT COUNT(*) {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject) AND {obs} {gp}
    """.format(join=OBS_JOIN, obs=obs, gp=gene_pred), params, default=0))

    rows = db.rows("""
        SELECT f.finding_id, f.subject_id, f.family_id, f.gene_symbol AS gene,
               COALESCE(f.hgvs_p, f.hgvs_c) AS variant, f.hgvs_c, f.hgvs_p,
               f.variant_key, f.chrom, f.pos, f.consequence, f.var_class,
               f.zygosity, f.zygosity_raw, f.vaf, f.depth, f.gnomad_af,
               f.gnomad_sas_af, f.ga100k_sas_af, f.clinvar_sig, f.clinvar_id,
               i.classification, i.acmg_codes, i.curated_flag, i.reportable,
               i.interpreted_at, i.inherited_from, i.segregation,
               gd.condition, gd.mondo_id, gd.inheritance, gd.validity, gd.penetrance,
               gd.gene_sets LIKE '%SF%' AS is_sf_gene,
               s.relation, s.sex, s.age, s.ancestry, s.indication,
               s.affected_status, s.consent_class,
               r.test_code, r.assay_version, r.pipeline_version, r.reference_build,
               r.mean_depth, sm.sample_type,
               (SELECT COUNT(DISTINCT f2.subject_id) FROM finding f2
                WHERE f2.variant_key = f.variant_key
                  AND f2.subject_id IN (SELECT subject_id FROM cohort_subject)
               ) AS recurrence
        {join}
        JOIN subject s ON s.subject_id = f.subject_id
        JOIN run r ON r.run_id = f.run_id
        JOIN sample sm ON sm.sample_id = r.sample_id
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject) AND {obs} {gp}
        ORDER BY
          CASE i.classification
            WHEN 'Pathogenic' THEN 0 WHEN 'Likely pathogenic' THEN 1
            WHEN 'Uncertain significance' THEN 2 WHEN 'Likely benign' THEN 3
            ELSE 4 END,
          f.gene_symbol, f.pos
        LIMIT ?
    """.format(join=OBS_JOIN, obs=obs, gp=gene_pred), params + [limit])

    return {
        "rows": rows,
        "total": total,
        "shown": len(rows),
        "truncated": total > len(rows),
        "gene": gene,
        "footer": ("Showing first {} of {} observations.".format(len(rows), total)
                   if total > len(rows) else "{} observations.".format(total)),
        "handoff_note": ("This module reads and links out. Classification happens in "
                         "the interpretation view; nothing here changes a call."),
    }


def variant_detail(finding_id: str, criteria: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """§G11 detail modal. Includes reportability and SF consent status per
    observation — done-when requires both to be visible."""
    obs = _obs_where(criteria)
    row = db.row("""
        SELECT f.*, i.classification, i.acmg_codes, i.curated_flag, i.reportable,
               i.interpreted_at, i.inherited_from, i.segregation, i.kb_snapshot_id,
               gd.condition, gd.mondo_id, gd.inheritance, gd.validity, gd.penetrance,
               gd.gene_sets,
               s.relation, s.sex, s.age, s.ancestry, s.indication, s.family_id AS fam,
               s.affected_status, s.family_history, s.consent_class, s.phenotype_hpo,
               r.run_id, r.test_code, r.assay_version, r.pipeline_version,
               r.reference_build, r.caller, r.mean_depth, r.pct_bases_20x,
               r.qc_status, r.implied_panel_id, sm.sample_type, sm.collection_date,
               tc.name AS test_name, tc.sf_capable, tc.scope_complete
        {join}
        JOIN subject s ON s.subject_id = f.subject_id
        JOIN run r ON r.run_id = f.run_id
        JOIN sample sm ON sm.sample_id = r.sample_id
        LEFT JOIN test_code tc ON tc.code = r.test_code
        WHERE f.finding_id = ? AND {obs}
    """.format(join=OBS_JOIN, obs=obs), [finding_id])
    if not row:
        return None

    row["family_size"] = int(db.scalar(
        "SELECT COUNT(*) FROM subject WHERE family_id = ?", [row["fam"]], default=0))
    row["recurrence"] = int(db.scalar("""
        SELECT COUNT(DISTINCT subject_id) FROM finding
        WHERE variant_key = ?
          AND subject_id IN (SELECT subject_id FROM cohort_subject)
    """, [row["variant_key"]], default=0))
    row["is_sf_gene"] = "SF" in (row.get("gene_sets") or "")
    row["sf_consent"] = row["consent_class"] == "Full research + secondary findings"
    row["sf_returnable"] = bool(row["is_sf_gene"] and row["sf_consent"]
                                and row.get("sf_capable"))
    row["scope_status"] = "Complete" if row.get("scope_complete") else "Provisional"
    return row


# ============================================================== G12 · subjects ==
def subject_list(criteria: Dict[str, Any], limit: Optional[int] = None) -> Dict[str, Any]:
    """§G12 · L2. Done-when: "The result classification matches the yield
    computation in G02 exactly." It reads carrier.subject_results()."""
    obs = _obs_where(criteria)
    limit = limit or config.TABLE_RENDER_CAP
    results = subject_results(criteria)

    rows = db.rows("""
        SELECT s.subject_id, s.family_id, s.relation, s.is_proband, s.indication,
               s.age, s.age_bucket, s.sex, s.ancestry, s.affected_status,
               s.family_history, s.consent_class, s.referral_source, s.phenotype_hpo,
               r.test_code, r.reference_build, r.pipeline_version, r.qc_status,
               sm.collection_date,
               (SELECT COUNT(*) FROM subject s2
                WHERE s2.family_id = s.family_id) AS family_size
        FROM cohort_subject cs
        JOIN subject s USING (subject_id)
        JOIN run r ON r.subject_id = s.subject_id
        JOIN sample sm ON sm.sample_id = r.sample_id
        ORDER BY s.subject_id
        LIMIT ?
    """, [limit])

    ids = [r["subject_id"] for r in rows]
    findings = _findings_by_subject(ids, obs) if ids else {}

    for r in rows:
        f = findings.get(r["subject_id"], {})
        r["findings"] = f.get("n", 0)
        r["plp_findings"] = f.get("plp", 0)
        r["vus_findings"] = f.get("vus", 0)
        r["key_genes"] = f.get("genes", "")
        r["result"] = results.get(r["subject_id"], "Negative")

    total = int(db.scalar("SELECT COUNT(*) FROM cohort_subject", default=0))
    return {"rows": rows, "total": total, "shown": len(rows),
            "truncated": total > len(rows),
            "result_counts": _tally(results.values())}


def _findings_by_subject(ids: List[str], obs: str) -> Dict[str, Dict[str, Any]]:
    rows = db.rows("""
        SELECT f.subject_id, COUNT(*) AS n,
               COUNT(*) FILTER (WHERE {plp}) AS plp,
               COUNT(*) FILTER (WHERE {vus}) AS vus,
               STRING_AGG(DISTINCT CASE WHEN {plp} THEN f.gene_symbol END, ', ') AS genes
        {join}
        WHERE f.subject_id IN ({ph}) AND {obs}
        GROUP BY 1
    """.format(join=OBS_JOIN, obs=obs, plp=PLP_SQL, vus=VUS_SQL,
               ph=", ".join(["?"] * len(ids))), ids)
    return {r["subject_id"]: r for r in rows}


def _tally(values) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for v in values:
        out[v] = out.get(v, 0) + 1
    return out


def subject_detail(subject_id: str, criteria: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    obs = _obs_where(criteria)
    row = db.row("""
        SELECT s.*, r.run_id, r.test_code, r.assay_version, r.pipeline_version,
               r.reference_build, r.caller, r.mean_depth, r.pct_bases_20x,
               r.qc_status, r.implied_panel_id, sm.sample_type, sm.collection_date,
               tc.name AS test_name, tc.sf_capable, tc.scope_complete
        FROM subject s
        JOIN run r ON r.subject_id = s.subject_id
        JOIN sample sm ON sm.sample_id = r.sample_id
        LEFT JOIN test_code tc ON tc.code = r.test_code
        WHERE s.subject_id = ?
    """, [subject_id])
    if not row:
        return None

    row["family"] = db.rows("""
        SELECT subject_id, relation, is_proband, sex, age, affected_status,
               consent_class
        FROM subject WHERE family_id = ? ORDER BY is_proband DESC, subject_id
    """, [row["family_id"]])
    row["findings"] = db.rows("""
        SELECT f.finding_id, f.gene_symbol AS gene,
               COALESCE(f.hgvs_p, f.hgvs_c) AS variant, f.consequence, f.zygosity,
               f.vaf, f.depth, f.gnomad_af, f.clinvar_sig,
               i.classification, i.acmg_codes, i.reportable, i.curated_flag,
               gd.condition, gd.inheritance, gd.validity, gd.penetrance
        {join}
        WHERE f.subject_id = ? AND {obs}
        ORDER BY CASE i.classification
                   WHEN 'Pathogenic' THEN 0 WHEN 'Likely pathogenic' THEN 1
                   WHEN 'Uncertain significance' THEN 2 ELSE 3 END, f.gene_symbol
    """.format(join=OBS_JOIN, obs=obs), [subject_id])
    row["result"] = subject_results(criteria).get(subject_id, "Negative")
    row["scope_status"] = "Complete" if row.get("scope_complete") else "Provisional"
    return row


# ================================================================== G13 · runs ==
def run_list(criteria: Dict[str, Any], limit: Optional[int] = None) -> Dict[str, Any]:
    """§G13 · L3. Provenance and QC per run, with per-run scope status."""
    limit = limit or config.TABLE_RENDER_CAP
    obs = _obs_where(criteria)

    rows = db.rows("""
        SELECT r.run_id, r.subject_id, r.test_code, r.assay_version,
               r.pipeline_version, r.reference_build, r.caller, r.mean_depth,
               r.pct_bases_20x, r.qc_status, r.implied_panel_id, r.source_file,
               sm.sample_type, sm.collection_date,
               tc.name AS test_name, tc.scope_complete,
               (SELECT COUNT(*) FROM run_scope rs
                WHERE rs.run_id = r.run_id) AS genes_in_scope,
               (SELECT COUNT(*) FROM run_scope rs
                WHERE rs.run_id = r.run_id AND rs.confidence = 'unknown'
               ) AS genes_provisional,
               (SELECT COUNT(*) FROM run_scope rs
                WHERE rs.run_id = r.run_id AND rs.confidence = 'inferred'
               ) AS genes_inferred,
               (SELECT COUNT(*) FROM finding f WHERE f.run_id = r.run_id) AS findings
        FROM cohort_run cr
        JOIN run r USING (run_id)
        JOIN sample sm ON sm.sample_id = r.sample_id
        LEFT JOIN test_code tc ON tc.code = r.test_code
        ORDER BY r.run_id
        LIMIT ?
    """, [limit])

    for r in rows:
        prov = int(r["genes_provisional"] or 0)
        r["scope_status"] = "Provisional" if prov else "Complete"
        r["scope_source"] = ("declared" if r["test_code"]
                             else "inferred · {}".format(r["implied_panel_id"] or "?"))

    total = int(db.scalar("SELECT COUNT(*) FROM cohort_run", default=0))
    prov = db.row("""
        SELECT COUNT(DISTINCT reference_build) AS builds,
               COUNT(DISTINCT pipeline_version) AS pipelines,
               COUNT(DISTINCT caller) AS callers
        FROM run WHERE run_id IN (SELECT run_id FROM cohort_run)
    """) or {}

    warning = None
    if int(prov.get("builds") or 0) > 1 or int(prov.get("pipelines") or 0) > 1:
        warning = ("Cohort mixes {} reference builds and {} pipeline versions. "
                   "Acceptable for descriptive counts; not acceptable for frequency "
                   "comparison without stratification."
                   .format(prov.get("builds"), prov.get("pipelines")))

    return {
        "rows": rows, "total": total, "shown": len(rows),
        "truncated": total > len(rows),
        "provenance": prov, "provenance_warning": warning,
        "qc_summary": db.rows("""
            SELECT r.qc_status AS label, COUNT(*) AS value
            FROM cohort_run cr JOIN run r USING (run_id)
            GROUP BY 1 ORDER BY value DESC
        """),
    }
