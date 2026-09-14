"""Governance layer — spec C02, C03, C04, and §2.5/§2.7.

    "These four are modals, reachable from the rail on every screen. They are
     the trust layer and the strongest part of the client demo."

Output classification (§2.5) is stored as module metadata so a screen cannot
ship without one, and it is logged on every run.
"""
from __future__ import annotations

import csv
import io
import json
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from .. import config, db
from . import denominator as den

# ------------------------------------------------------- output classes (§2.5)
# CLINICAL is never produced by this module.
MODULE_CLASS: Dict[str, str] = {
    "g_builder": "OPERATIONAL",
    "g_dashboard": "RESEARCH",
    "g_carrier": "RESEARCH",
    "g_zygosity": "RESEARCH",
    "g_genedisease": "RESEARCH",
    "g_phenotype": "RESEARCH",
    "g_popfreq": "RESEARCH",
    "g_vus": "OPERATIONAL",
    "g_sf": "RESEARCH",
    "g_gene": "RESEARCH",
    "g_variants": "RESEARCH",
    "g_subjects": "RESEARCH",
    "g_runs": "OPERATIONAL",
    "c_denominator": "OPERATIONAL",
    "c_sampleqc": "OPERATIONAL",
    "c_manifest": "OPERATIONAL",
    "c_audit": "OPERATIONAL",
    "d_compile": "OPERATIONAL",
}

MODULE_TITLES: Dict[str, str] = {
    "g_builder": "Cohort builder", "g_dashboard": "Cohort dashboard",
    "g_carrier": "Carrier & diagnostic yield", "g_zygosity": "Zygosity & inheritance",
    "g_genedisease": "Gene–disease association", "g_phenotype": "Phenotype (HPO)",
    "g_popfreq": "Population frequency", "g_vus": "VUS inventory",
    "g_sf": "Secondary findings (ACMG SF v3.3)", "g_gene": "Gene analysis",
    "g_variants": "Variant analysis", "g_subjects": "Subject list",
    "g_runs": "Assay list / QC", "c_denominator": "Denominator inspector",
    "c_sampleqc": "Sample QC",
    "c_manifest": "Cohort manifest", "c_audit": "Audit trail",
    "d_compile": "Ask a cohort question",
}


def output_class(module: str) -> str:
    """A module with no declared class is a bug, not a default."""
    if module not in MODULE_CLASS:
        raise KeyError("module '{}' has no output class (spec §2.5)".format(module))
    return MODULE_CLASS[module]


# ================================================================ C02 manifest
def manifest(cohort, module: Optional[str] = None) -> Dict[str, Any]:
    """§C02 · Reproducibility as an artifact.

    criteria_hash + kb_snapshot_id together make a figure reproducible.
    member_hash detects silent membership drift. pipeline_versions and
    reference_builds expose provenance mixing. Without KB pinning, re-running
    a cohort after a ClinVar update yields a different number with no
    explanation.
    """
    builds = [r["reference_build"] for r in db.rows("""
        SELECT DISTINCT reference_build FROM run
        WHERE run_id IN (SELECT run_id FROM cohort_run) AND reference_build IS NOT NULL
        ORDER BY 1
    """)]
    pipelines = [r["pipeline_version"] for r in db.rows("""
        SELECT DISTINCT pipeline_version FROM run
        WHERE run_id IN (SELECT run_id FROM cohort_run) AND pipeline_version IS NOT NULL
        ORDER BY 1
    """)]
    test_codes = [r["test_code"] for r in db.rows("""
        SELECT DISTINCT COALESCE(test_code, '(inferred from fingerprint)') AS test_code
        FROM run WHERE run_id IN (SELECT run_id FROM cohort_run) ORDER BY 1
    """)]

    dens = den.gene_denominators()
    provisional_genes = sum(1 for d in dens.values() if d["provisional"])
    inferred_genes = sum(1 for d in dens.values() if d["confidence"] == "inferred")

    return {
        "cohort_name": cohort.name,
        "profile": config.PROFILE,
        "tenant_id": config.TENANT_ID,
        "resolved_at": cohort.resolved_at,
        # The date a relative collection window was measured back from. Without
        # it a "last 24 months" cohort cannot be reproduced after an ingest,
        # because the anchor moves with the newest sample.
        "collection_anchor": cohort.as_of,
        "criteria": cohort.criteria,
        "member_count": cohort.n_subjects,
        "family_count": cohort.n_families,
        "member_hash": cohort.member_hash,
        "criteria_hash": cohort.criteria_hash,
        "kb_snapshot_id": cohort.kb_snapshot_id,
        "pipeline_versions": pipelines,
        "reference_builds": builds,
        "test_codes": test_codes,
        "output_class": output_class(module) if module else "RESEARCH",
        "small_cell_suppression": True,
        "suppression_threshold": config.SMALL_CELL_THRESHOLD,
        "coverage_confidence": {
            "genes_total": len(dens),
            "genes_with_provisional_scope": provisional_genes,
            "genes_inferred_scope": inferred_genes,
        },
        "store": {
            "source": db.meta_get("source", "unknown"),
            "built_at": db.meta_get("built_at"),
        },
        "generated_by": config.DEFAULT_USER,
    }


# =================================================================== C03 audit
def log_run(module: str, cohort, user: Optional[str] = None) -> Dict[str, Any]:
    """§2.7. Every screen render logs timestamp, profile, module, output class,
    cohort size and criteria hash."""
    entry = {
        "analysis_id": uuid.uuid4().hex[:12],
        "snapshot_id": None,
        "module": module,
        "output_class": output_class(module),
        "executed_at": datetime.utcnow(),
        "user_id": user or config.DEFAULT_USER,
        "cohort_size": cohort.n_subjects,
        "criteria_hash": cohort.criteria_hash,
    }
    db.insert_rows("analysis_run", list(entry), [list(entry.values())])
    return entry


def audit_trail(limit: int = 200) -> Dict[str, Any]:
    rows = db.rows("""
        SELECT analysis_id, module, output_class, executed_at, user_id,
               cohort_size, criteria_hash
        FROM analysis_run ORDER BY executed_at DESC LIMIT ?
    """, [limit])
    for r in rows:
        r["module_title"] = MODULE_TITLES.get(r["module"], r["module"])
    return {
        "rows": rows,
        "total": int(db.scalar("SELECT COUNT(*) FROM analysis_run", default=0)),
        "footer": ("Every analysis run is recorded with its output class. Clinical "
                   "report templates reject any asset whose class is RESEARCH at "
                   "build time."),
    }


def snapshot(cohort) -> str:
    """Pin a cohort resolution so a figure can be reproduced later (§C02)."""
    sid = "snap-" + uuid.uuid4().hex[:10]
    db.insert_rows("cohort_snapshot",
                   ["snapshot_id", "cohort_id", "resolved_at", "member_hash",
                    "criteria_hash", "kb_snapshot_id", "member_count"],
                   [[sid, cohort.name, datetime.utcnow(), cohort.member_hash,
                     cohort.criteria_hash, cohort.kb_snapshot_id, cohort.n_subjects]])
    return sid


# ================================================================== C04 export
WATERMARK = ("RESEARCH USE ONLY — not for clinical reporting. Generated by "
             "ImpactOmics Germline Cohort Analytics.")


def _suppress(rows: List[Dict[str, Any]], count_columns: Sequence[str]) -> List[Dict[str, Any]]:
    """Small-cell suppression applied BEFORE write (spec C04 rules)."""
    out = []
    for r in rows:
        rr = dict(r)
        for col in count_columns:
            if col in rr and isinstance(rr[col], (int, float)) and rr[col] is not None:
                rr[col] = den.small_cell(int(rr[col]))
        out.append(rr)
    return out


COUNT_COLUMNS = ("subjects", "families", "plp_subjects", "plp_families",
                 "any_subjects", "carriers", "vus_subjects", "n", "value",
                 "unique_variants", "unique_plp_variants")


def export(cohort, module: str, rows: List[Dict[str, Any]],
           fmt: str = "csv", user: Optional[str] = None) -> Dict[str, Any]:
    """Produce an export that "can be handed to a third party and fully
    understood without asking us" (spec C04 done-when).

    Every format embeds the manifest. RESEARCH exports carry a visible in-file
    watermark. Suppression happens before the bytes are written, not in a view.
    """
    mf = manifest(cohort, module=module)
    cls = mf["output_class"]
    safe_rows = _suppress(rows, COUNT_COLUMNS)
    stamp = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    base = "germline_{}_{}_{}".format(module, cohort.criteria_hash, stamp)

    config.EXPORT_DIR.mkdir(parents=True, exist_ok=True)

    if fmt == "json":
        path = config.EXPORT_DIR / (base + ".json")
        path.write_text(json.dumps({
            "watermark": WATERMARK if cls == "RESEARCH" else None,
            "manifest": mf, "module": module,
            "module_title": MODULE_TITLES.get(module, module),
            "row_count": len(safe_rows), "rows": safe_rows,
        }, indent=2, default=str))

    elif fmt == "xlsx":
        path = config.EXPORT_DIR / (base + ".xlsx")
        _write_xlsx(path, safe_rows, mf, module, cls)

    else:
        fmt = "csv"
        path = config.EXPORT_DIR / (base + ".csv")
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            if cls == "RESEARCH":
                w.writerow(["# " + WATERMARK])
            for k, v in mf.items():
                w.writerow(["# {}".format(k),
                            json.dumps(v, default=str) if isinstance(v, (dict, list)) else v])
            w.writerow([])
            if safe_rows:
                cols = list(safe_rows[0])
                w.writerow(cols)
                for r in safe_rows:
                    w.writerow([r.get(c) for c in cols])

    log_run(module, cohort, user=user)
    return {"path": str(path), "filename": path.name, "format": fmt,
            "rows": len(safe_rows), "output_class": cls,
            "watermarked": cls == "RESEARCH", "manifest": mf}


def _write_xlsx(path: Path, rows: List[Dict[str, Any]], mf: Dict[str, Any],
                module: str, cls: str) -> None:
    from openpyxl import Workbook
    from openpyxl.styles import Font

    wb = Workbook()
    ws = wb.active
    ws.title = "Data"
    start = 1
    if cls == "RESEARCH":
        ws.cell(row=1, column=1, value=WATERMARK).font = Font(bold=True, color="B00020")
        start = 3
    if rows:
        cols = list(rows[0])
        for j, c in enumerate(cols, 1):
            ws.cell(row=start, column=j, value=c).font = Font(bold=True)
        for i, r in enumerate(rows, start + 1):
            for j, c in enumerate(cols, 1):
                v = r.get(c)
                ws.cell(row=i, column=j,
                        value=json.dumps(v, default=str) if isinstance(v, (dict, list)) else v)

    # Manifest on a second sheet (spec C04 formats).
    ms = wb.create_sheet("Manifest")
    ms.cell(row=1, column=1, value="Field").font = Font(bold=True)
    ms.cell(row=1, column=2, value="Value").font = Font(bold=True)
    for i, (k, v) in enumerate(mf.items(), 2):
        ms.cell(row=i, column=1, value=k)
        ms.cell(row=i, column=2,
                value=json.dumps(v, default=str) if isinstance(v, (dict, list)) else str(v))
    wb.save(path)
