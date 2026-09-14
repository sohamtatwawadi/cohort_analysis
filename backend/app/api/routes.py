"""HTTP API.

Every endpoint that renders a module resolves the cohort, logs an audit entry
(spec §2.7) and returns the module payload with its output class attached. The
browser receives computed numbers; it never computes one.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Body, HTTPException, Query
from pydantic import BaseModel, Field

from .. import config, db
from ..ingest import store, synthetic, varimat
from ..reference.genes import GENE_DISEASE, GENE_LIST, GENE_SETS
from ..reference.tests import (AFFECTED_STATUSES, AGE_BUCKETS, ANCESTRIES,
                               ASSAY_VERSIONS, CONSENT_CLASSES, FAMILY_HISTORY,
                               INDICATIONS, PIPELINE_VERSIONS, REFERENCE_BUILDS,
                               REFERRAL_SOURCES, RELATIONS, SAMPLE_TYPES,
                               TEST_CODE_LIST, TEST_CODES)
from ..services import (analysis, carrier, cohort as cohort_svc, compiler,
                        denominator, explore, governance, library, sampleqc,
                        suggest as suggest_svc, zygosity)

router = APIRouter()


class CohortRequest(BaseModel):
    criteria: Optional[Dict[str, Any]] = None
    name: str = "All germline subjects · probands"


class SaveRequest(BaseModel):
    name: str
    criteria: Dict[str, Any] = Field(default_factory=dict)
    view: str = "g_dashboard"
    description: str = ""


class QuestionRequest(BaseModel):
    question: str
    payload: Optional[Dict[str, Any]] = None


class ExportRequest(BaseModel):
    criteria: Optional[Dict[str, Any]] = None
    name: str = "All germline subjects · probands"
    module: str
    format: str = "csv"


def _resolve(req: CohortRequest):
    return cohort_svc.resolve(req.criteria, name=req.name)


def _render(module: str, req: CohortRequest, payload_fn):
    """Resolve → compute → log. The shape every module endpoint returns.

    Held under COHORT_LOCK end to end: the cohort lives in connection-scoped
    temp tables, so releasing between resolve and compute would let a concurrent
    request swap the cohort out from under this one.
    """
    with cohort_svc.cohort_scope(req.criteria, req.name) as co:
        payload = payload_fn(co)
        governance.log_run(module, co)
    return {
        "module": module,
        "module_title": governance.MODULE_TITLES[module],
        "output_class": governance.output_class(module),
        "cohort": co.as_dict(),
        "data": payload,
    }


# ------------------------------------------------------------------- meta ----
@router.get("/meta")
def meta() -> Dict[str, Any]:
    """Everything the shell needs to render: vocabularies, modules, store state."""
    return {
        "profile": config.PROFILE,
        "tenant_id": config.TENANT_ID,
        "kb_snapshot_id": db.meta_get("kb_snapshot_id", config.KB_SNAPSHOT_ID),
        "store": {
            "source": db.meta_get("source", "none"),
            "built_at": db.meta_get("built_at"),
            "counts": db.table_counts(),
            "empty": store.is_empty(),
        },
        "modules": [
            {"module": m, "title": governance.MODULE_TITLES[m],
             "output_class": governance.MODULE_CLASS[m]}
            for m in ["g_builder", "g_dashboard", "g_carrier", "g_zygosity",
                      "g_genedisease", "g_phenotype", "g_popfreq", "g_vus",
                      "g_sf", "g_gene", "g_variants"]
        ],
        "drill_path": [
            {"level": "L1", "module": "g_dashboard", "label": "Cohort"},
            {"level": "L2", "module": "g_subjects", "label": "Subject"},
            {"level": "L3", "module": "g_runs", "label": "Assay"},
            {"level": "L4", "module": "g_gene", "label": "Gene"},
            {"level": "L5", "module": "g_variants", "label": "Variant"},
        ],
        "vocabulary": {
            "indication": INDICATIONS, "sex": ["F", "M"], "age_bucket": AGE_BUCKETS,
            "ancestry": ANCESTRIES, "affected_status": AFFECTED_STATUSES,
            "family_history": FAMILY_HISTORY, "relation": RELATIONS,
            "referral_source": REFERRAL_SOURCES, "consent_class": CONSENT_CLASSES,
            "test_code": TEST_CODE_LIST, "assay_version": ASSAY_VERSIONS,
            "pipeline_version": PIPELINE_VERSIONS, "reference_build": REFERENCE_BUILDS,
            "sample_type": SAMPLE_TYPES, "gene": GENE_LIST,
            "var_class": compiler.VAR_CLASSES, "consequence": compiler.CONSEQUENCES,
            "classification": compiler.CLASSIFICATIONS, "zygosity": compiler.ZYGOSITIES,
            "inheritance": compiler.INHERITANCES, "validity": compiler.VALIDITIES,
            "penetrance": ["High", "Moderate", "Low"],
            "clinvar_sig": ["Pathogenic", "Likely pathogenic",
                            "Pathogenic/Likely pathogenic", "Uncertain significance",
                            "Conflicting interpretations", "Benign", "Not reported"],
        },
        "gene_sets": {k: v for k, v in GENE_SETS.items()},
        "test_codes": [{
            "code": t.code, "name": t.name, "version": t.assay_version,
            "indication": t.indication, "sf_capable": t.sf_capable,
            "scope_complete": t.scope_complete, "genes": len(t.genes),
        } for t in TEST_CODES.values()],
        "criteria_labels": cohort_svc.CRITERIA_LABELS,
        "blank_criteria": cohort_svc.blank(),
        "suggested_questions": compiler.SUGGESTED_QUESTIONS,
        "constants": {
            "denom_min": config.DENOM_MIN,
            "provisional_max": config.PROVISIONAL_MAX,
            "small_cell": config.SMALL_CELL_THRESHOLD,
            "render_cap": config.TABLE_RENDER_CAP,
        },
    }


# ----------------------------------------------------------------- cohort ----
@router.post("/cohort/resolve")
def cohort_resolve(req: CohortRequest) -> Dict[str, Any]:
    """§G01 builder. Returns the funnel that attributes every dropped subject."""
    # family_stats, provisional_stats and the composition queries all read the
    # cohort temp tables, so they belong inside the same critical section as the
    # resolve that created them.
    with cohort_svc.cohort_scope(req.criteria, req.name) as co:
        governance.log_run("g_builder", co)
        out = co.as_dict()
        out["family_stats"] = cohort_svc.family_stats()
        out["provisional"] = cohort_svc.provisional_stats()
        out["composition"] = {
            "indication": db.rows("""
            SELECT s.indication AS label, COUNT(*) AS value
            FROM cohort_subject cs JOIN subject s USING (subject_id)
            GROUP BY 1 ORDER BY value DESC"""),
            "test_code": db.rows("""
            SELECT COALESCE(r.test_code, '(inferred scope)') AS label, COUNT(*) AS value
            FROM cohort_run cr JOIN run r USING (run_id)
            GROUP BY 1 ORDER BY value DESC"""),
            "ancestry": db.rows("""
            SELECT s.ancestry AS label, COUNT(*) AS value
            FROM cohort_subject cs JOIN subject s USING (subject_id)
            GROUP BY 1 ORDER BY value DESC"""),
            "consent": db.rows("""
            SELECT s.consent_class AS label, COUNT(*) AS value
            FROM cohort_subject cs JOIN subject s USING (subject_id)
            GROUP BY 1 ORDER BY value DESC"""),
        }
    out["module"] = "g_builder"
    out["output_class"] = governance.output_class("g_builder")
    return out


@router.post("/cohort/suggestions")
def cohort_suggestions(req: CohortRequest) -> Dict[str, Any]:
    """Which questions are worth asking of THIS cohort, and why.

    Advice, not a gate — a lab cohort can run any descriptive analysis. The
    point is that a flat menu makes a pointless option look identical to a
    useful one until you click it.
    """
    with cohort_svc.cohort_scope(req.criteria, req.name) as co:
        return suggest_svc.suggest(co.criteria)


@router.post("/module/g_dashboard")
def m_dashboard(req: CohortRequest) -> Dict[str, Any]:
    return _render("g_dashboard", req, lambda co: carrier.dashboard(co.criteria))


@router.post("/module/g_carrier")
def m_carrier(req: CohortRequest) -> Dict[str, Any]:
    def build(co):
        return {
            "kpis": carrier.dashboard(co.criteria)["kpis"],
            "yield": carrier.diagnostic_yield(co.criteria),
            "table": carrier.carrier_table(co.criteria),
            "by_test_code": carrier.yield_by("test_code", co.criteria),
            "by_ancestry": carrier.yield_by("ancestry", co.criteria),
            "warn_legend": carrier.WARN_LEGEND,
            "ancestry_caveat": carrier.ANCESTRY_CAVEAT,
        }
    return _render("g_carrier", req, build)


@router.post("/module/g_zygosity")
def m_zygosity(req: CohortRequest) -> Dict[str, Any]:
    return _render("g_zygosity", req, lambda co: zygosity.overview(co.criteria))


@router.post("/module/g_genedisease")
def m_genedisease(req: CohortRequest) -> Dict[str, Any]:
    return _render("g_genedisease", req, lambda co: analysis.gene_disease(co.criteria))


@router.post("/module/g_phenotype")
def m_phenotype(req: CohortRequest) -> Dict[str, Any]:
    return _render("g_phenotype", req, lambda co: analysis.phenotype(co.criteria))


@router.post("/module/g_popfreq")
def m_popfreq(req: CohortRequest) -> Dict[str, Any]:
    return _render("g_popfreq", req,
                   lambda co: analysis.population_frequency(co.criteria))


@router.post("/module/g_vus")
def m_vus(req: CohortRequest) -> Dict[str, Any]:
    return _render("g_vus", req, lambda co: analysis.vus_inventory(co.criteria))


@router.post("/module/g_sf")
def m_sf(req: CohortRequest) -> Dict[str, Any]:
    return _render("g_sf", req, lambda co: analysis.secondary_findings(co.criteria))


@router.post("/module/g_gene")
def m_gene(req: CohortRequest, gene: Optional[str] = Query(None)) -> Dict[str, Any]:
    def build(co):
        pills = explore.gene_pills(co.criteria)
        target = gene or (pills[0]["gene"] if pills else None)
        return {"pills": pills,
                "detail": explore.gene_detail(target, co.criteria) if target else None}
    return _render("g_gene", req, build)


@router.post("/module/g_variants")
def m_variants(req: CohortRequest, gene: Optional[str] = Query(None)) -> Dict[str, Any]:
    def build(co):
        return {"pills": explore.gene_pills(co.criteria),
                **explore.variant_list(co.criteria, gene=gene)}
    return _render("g_variants", req, build)


@router.post("/module/g_subjects")
def m_subjects(req: CohortRequest) -> Dict[str, Any]:
    return _render("g_subjects", req, lambda co: explore.subject_list(co.criteria))


@router.post("/module/g_runs")
def m_runs(req: CohortRequest) -> Dict[str, Any]:
    return _render("g_runs", req, lambda co: explore.run_list(co.criteria))


@router.post("/detail/variant/{finding_id}")
def d_variant(finding_id: str, req: CohortRequest) -> Dict[str, Any]:
    with cohort_svc.cohort_scope(req.criteria, req.name) as co:
        row = explore.variant_detail(finding_id, co.criteria)
    if not row:
        raise HTTPException(404, "finding not in cohort")
    return row


@router.post("/detail/subject/{subject_id}")
def d_subject(subject_id: str, req: CohortRequest) -> Dict[str, Any]:
    with cohort_svc.cohort_scope(req.criteria, req.name) as co:
        row = explore.subject_detail(subject_id, co.criteria)
    if not row:
        raise HTTPException(404, "subject not found")
    return row


# ------------------------------------------------------------- governance ----
@router.post("/qc/samples")
def qc_samples(req: CohortRequest) -> Dict[str, Any]:
    """Sample QC for the resolved cohort. A contaminated sample or a swap
    invalidates every rate downstream, however careful the denominators."""
    return _render("c_sampleqc", req, lambda co: sampleqc.sample_qc(co.criteria))


@router.post("/governance/denominator")
def g_denominator(req: CohortRequest) -> Dict[str, Any]:
    return _render("c_denominator", req, lambda co: denominator.coverage_inspector())


@router.post("/governance/manifest")
def g_manifest(req: CohortRequest, module: Optional[str] = Query(None)) -> Dict[str, Any]:
    with cohort_svc.cohort_scope(req.criteria, req.name) as co:
        return governance.manifest(co, module=module)


@router.get("/governance/audit")
def g_audit(limit: int = 200) -> Dict[str, Any]:
    return governance.audit_trail(limit=limit)


@router.post("/governance/export")
def g_export(req: ExportRequest) -> Dict[str, Any]:
    with cohort_svc.cohort_scope(req.criteria, req.name) as co:
        rows = _rows_for_module(req.module, co)
        return governance.export(co, req.module, rows, fmt=req.format)


def _rows_for_module(module: str, co) -> List[Dict[str, Any]]:
    c = co.criteria
    if module == "g_carrier":
        return [{k: v for k, v in r.items() if not isinstance(v, dict)}
                | {"carrier_rate_pct": r["carrier_rate"]["pct"],
                   "carrier_rate": r["carrier_rate"]["label"]}
                for r in carrier.carrier_table(c)]
    if module == "g_vus":
        return [{k: v for k, v in r.items() if not isinstance(v, dict)}
                for r in analysis.vus_inventory(c)["rows"]]
    if module == "g_variants":
        return explore.variant_list(c, limit=100000)["rows"]
    if module == "g_subjects":
        return explore.subject_list(c, limit=100000)["rows"]
    if module == "g_runs":
        return explore.run_list(c, limit=100000)["rows"]
    if module == "g_popfreq":
        return [{k: v for k, v in r.items() if not isinstance(v, list)}
                for r in analysis.population_frequency(c)["rows"]]
    if module == "g_zygosity":
        return zygosity.by_gene(c)
    if module == "g_genedisease":
        return analysis.gene_disease(c)["reference"]
    if module == "c_denominator":
        return denominator.coverage_inspector()["genes"]
    raise HTTPException(400, "module '{}' has no export shape".format(module))


# ---------------------------------------------------------------- library ----
@router.get("/cohorts")
def c_list() -> Dict[str, Any]:
    return {"cohorts": library.list_cohorts()}


@router.post("/cohorts")
def c_save(req: SaveRequest) -> Dict[str, Any]:
    return library.save(req.name, req.criteria, view=req.view,
                        description=req.description)


@router.delete("/cohorts/{cohort_id}")
def c_delete(cohort_id: str) -> Dict[str, Any]:
    return {"deleted": library.delete(cohort_id)}


@router.get("/cohorts/{cohort_id}/verify")
def c_verify(cohort_id: str) -> Dict[str, Any]:
    return library.verify_reproducible(cohort_id)


# --------------------------------------------------------------- compiler ----
@router.get("/compile/vocabulary")
def q_vocab() -> Dict[str, Any]:
    return compiler.vocabulary()


@router.post("/compile")
def q_compile(req: QuestionRequest) -> Dict[str, Any]:
    """Compile only. Nothing executes here — the user must approve first."""
    return compiler.compile_and_validate(req.question, payload=req.payload)


@router.post("/compile/{compilation_id}/execute")
def q_execute(compilation_id: str, req: CohortRequest) -> Dict[str, Any]:
    compiler.mark_executed(compilation_id)
    with cohort_svc.cohort_scope(req.criteria, req.name) as co:
        return co.as_dict()


@router.get("/compile/corpus")
def q_corpus(limit: int = 200) -> Dict[str, Any]:
    return {"rows": compiler.corpus(limit)}


# ----------------------------------------------------------------- ingest ----
@router.post("/admin/rebuild/synthetic")
def a_synthetic(n_families: int = 640, seed: int = 20260813) -> Dict[str, Any]:
    result = synthetic.rebuild(n_families=n_families, seed=seed)
    library.seed_builtins()
    return result


@router.post("/admin/rebuild/varimat")
def a_varimat(directory: Optional[str] = None,
              clinical_dir: Optional[str] = None) -> Dict[str, Any]:
    from pathlib import Path
    try:
        result = store.rebuild_from_varimat(
            directory=Path(directory) if directory else None,
            clinical_dir=Path(clinical_dir) if clinical_dir else None)
    except FileNotFoundError as e:
        raise HTTPException(400, str(e))
    library.seed_builtins()
    return result


@router.get("/admin/varimat/discover")
def a_discover(directory: Optional[str] = None) -> Dict[str, Any]:
    from pathlib import Path
    files = varimat.discover(Path(directory) if directory else config.VARIMAT_DIR)
    return {"directory": str(directory or config.VARIMAT_DIR),
            "files": [{"name": f.name, "path": str(f),
                       "size_mb": round(f.stat().st_size / 1e6, 2)} for f in files]}
