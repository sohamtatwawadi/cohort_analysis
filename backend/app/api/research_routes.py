"""Research Mode HTTP API — Part II.

Mirrors the Lab Mode surface but for uploaded datasets. The capability gate is
enforced in the service layer (registry.can_run, called by jobs.submit), so an
endpoint cannot start an analysis the dataset does not support even if the
client calls it directly.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Body, HTTPException, Query, UploadFile, File, Form
from pydantic import BaseModel, Field

from .. import config, db
from ..research import capability, jobs, registry, store as rstore, variantset
from ..research import analyses  # noqa: F401  — registers the analyses
from ..research.formats import detect
from ..research.profile import profile_dataset
from ..research.stats import power as power_mod
from ..research.types import PhenotypeTable
from ..research.validate import validate_upload

router = APIRouter()

UPLOAD_DIR = Path(config.DB_PATH).parent / "uploads"


class ProjectRequest(BaseModel):
    name: str
    description: str = ""
    owner: Optional[str] = None
    retention_days: Optional[int] = None


class OverrideRequest(BaseModel):
    analysis: str
    justification: str
    requested_by: Optional[str] = None


class JobRequest(BaseModel):
    analysis: str
    spec: Dict[str, Any] = Field(default_factory=dict)
    submitted_by: Optional[str] = None


class PowerRequest(BaseModel):
    kind: str = "binary"
    n_cases: int = 0
    n_controls: int = 0
    n: int = 0
    maf: float = 0.01
    odds_ratio: float = 2.0
    beta_sd: float = 0.2
    alpha: Optional[float] = None
    alpha_scope: str = "genome"


class VariantSetRequest(BaseModel):
    name: str
    definition: Dict[str, Any]
    created_by: Optional[str] = None


# --------------------------------------------------------------------- meta --
@router.get("/meta")
def research_meta() -> Dict[str, Any]:
    return {
        "mode": "research",
        "ownership_notice": registry.DATA_OWNERSHIP_NOTICE,
        "consent_statement": registry.CONSENT_STATEMENT,
        "thresholds": capability.THRESHOLDS,
        "underpowered_stamp": capability.UNDERPOWERED_STAMP,
        "non_overridable": sorted(capability.NON_OVERRIDABLE),
        "analyses": jobs.registered(),
        "variant_set_presets": variantset.PRESETS,
        "supported_formats": [
            {"format": "VCF / VCF.gz", "priority": "P1", "status": "supported"},
            {"format": "PLINK 1 (BED/BIM/FAM)", "priority": "P1", "status": "supported"},
            {"format": "VariMAT", "priority": "P1", "status": "supported (Lab Mode loader)"},
            {"format": "PLINK 2 (PGEN)", "priority": "P1", "status": "not implemented"},
            {"format": "BGEN", "priority": "P2", "status": "not implemented"},
            {"format": "Hail MatrixTable", "priority": "P2", "status": "not implemented"},
            {"format": "BCF", "priority": "P1", "status": "not implemented"},
        ],
    }


# ----------------------------------------------------------------- projects --
@router.get("/projects")
def list_projects() -> Dict[str, Any]:
    """Always returns at least one project. A dataset cannot be registered
    without one, so an empty list is a dead end for the UI rather than a
    meaningful state."""
    registry.ensure_default_project()
    return {"projects": registry.list_projects(),
            "orphaned_payloads": registry.orphaned_payloads()}


@router.post("/projects")
def create_project(req: ProjectRequest) -> Dict[str, Any]:
    return registry.create_project(req.name, req.owner or config.DEFAULT_USER,
                                   req.description, req.retention_days)


# ----------------------------------------------------------------- datasets --
@router.get("/projects/{project_id}/datasets")
def list_datasets(project_id: str) -> Dict[str, Any]:
    if not registry.get_project(project_id):
        raise HTTPException(404, "unknown project")
    return {"datasets": registry.list_datasets(project_id)}


@router.post("/projects/{project_id}/datasets/upload")
async def upload_dataset(
    project_id: str,
    name: str = Form(...),
    genome_build: Optional[str] = Form(None),
    consent_attested: bool = Form(False),
    genotypes: List[UploadFile] = File(...),
    phenotypes: Optional[UploadFile] = File(None),
    annotation_source: Optional[str] = Form(None),
) -> Dict[str, Any]:
    """§2.3 ingestion pipeline, from the browser.

    Takes a LIST of genotype files because PLINK is a triple — .bed carries the
    genotypes, .bim the variants and .fam the samples, and any one of them alone
    is unreadable. A VCF is just a list of one.
    """
    if not registry.get_project(project_id):
        raise HTTPException(404, "unknown project")
    if not consent_attested:
        raise HTTPException(400, "Consent attestation is required. "
                                 + registry.CONSENT_STATEMENT)

    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    for up in genotypes:
        if not up.filename:
            continue
        dest = UPLOAD_DIR / Path(up.filename).name
        dest.write_bytes(await up.read())
        written.append(dest)
    if not written:
        raise HTTPException(400, "No genotype file received.")

    # Point the loader at the .bed of a PLINK set; detect_format resolves the
    # siblings from there. For anything else the single file is the entry point.
    entry = next((w for w in written if w.suffix.lower() == ".bed"), written[0])

    ppath = None
    if phenotypes is not None and phenotypes.filename:
        ppath = UPLOAD_DIR / Path(phenotypes.filename).name
        ppath.write_bytes(await phenotypes.read())

    return _ingest(project_id, name, entry, ppath, genome_build, consent_attested,
                   annotation_source)


class IngestPathRequest(BaseModel):
    name: str
    genotype_path: str
    phenotype_path: Optional[str] = None
    genome_build: Optional[str] = None
    consent_attested: bool = False
    annotation_source: Optional[str] = None


@router.post("/projects/{project_id}/datasets/ingest")
def ingest_from_path(project_id: str, req: IngestPathRequest) -> Dict[str, Any]:
    """Register a dataset already on disk — the path a large upload takes,
    and the one the CLI and tests use."""
    if not registry.get_project(project_id):
        raise HTTPException(404, "unknown project")
    if not req.consent_attested:
        raise HTTPException(400, "Consent attestation is required. "
                                 + registry.CONSENT_STATEMENT)
    return _ingest(project_id, req.name, Path(req.genotype_path),
                   Path(req.phenotype_path) if req.phenotype_path else None,
                   req.genome_build, req.consent_attested, req.annotation_source)


def _ingest(project_id: str, name: str, gpath: Path, ppath: Optional[Path],
            build: Optional[str], consent: bool,
            annotation_source: Optional[str]) -> Dict[str, Any]:
    if not gpath.exists():
        raise HTTPException(400, "genotype file not found: {}".format(gpath))

    fmt = detect.detect_format(gpath)
    if fmt == "unknown":
        raise HTTPException(400, "Unrecognised genotype format. Supported: VCF, "
                                 "VCF.gz, PLINK1 (BED/BIM/FAM), VariMAT.")
    try:
        gm = detect.load_any(gpath)
    except Exception as exc:
        raise HTTPException(400, "Could not read {}: {}".format(fmt, exc))

    ph = _read_phenotypes(ppath) if ppath else None

    report, facts = validate_upload(gm, declared_build=build, phenotypes=ph,
                                    deep_duplicate_check=gm.n_samples <= 400)
    if not report.ok:
        # §2.3: validations FAIL the upload rather than warning.
        raise HTTPException(422, detail={
            "message": "Upload rejected: mandatory validation failed.",
            "report": report.as_dict()})

    prof = profile_dataset(gm, ph, coverage_confidence="declared" if build else "unknown",
                           compute_genetics=gm.n_variants <= 200_000)

    ds = registry.register_dataset(
        project_id, name, prof, fmt, [str(gpath)], facts["build"],
        config.DEFAULT_USER, consent_attested=consent,
        annotation={"source": annotation_source},
        phenotype_kinds=ph.kinds if ph else None)

    rstore.save_genotypes(ds["dataset_id"], gm)
    if ph:
        rstore.save_phenotypes(ds["dataset_id"], ph)

    return {"dataset": ds, "validation": report.as_dict(),
            "profile": prof.as_dict(),
            "capabilities": registry.get_capabilities(ds["dataset_id"])}


def _read_phenotypes(path: Path) -> PhenotypeTable:
    """CSV/TSV phenotypes. Column kind is inferred, which drives model selection
    (§4.1) — a binary column read as quantitative silently becomes a linear
    regression on 0/1."""
    import csv

    import numpy as np

    delim = "\t" if path.suffix in (".tsv", ".txt") else ","
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter=delim))
    if not rows:
        raise HTTPException(400, "phenotype file is empty")

    id_col = next((c for c in rows[0]
                   if c.lower() in ("sample_id", "iid", "id", "sample")), None)
    if not id_col:
        raise HTTPException(400, "phenotype file needs a sample_id / IID column")

    sample_ids = [r[id_col] for r in rows]
    columns: Dict[str, Any] = {}
    kinds: Dict[str, str] = {}
    for col in rows[0]:
        if col == id_col:
            continue
        raw = [r.get(col, "") for r in rows]
        numeric = []
        for v in raw:
            try:
                numeric.append(float(v) if v not in ("", "NA", "NaN", ".") else float("nan"))
            except ValueError:
                numeric = None
                break
        if numeric is None:
            columns[col] = np.asarray(raw, dtype=object)
            kinds[col] = "categorical"
        else:
            arr = np.asarray(numeric, dtype=float)
            finite = arr[np.isfinite(arr)]
            uniq = set(finite.tolist())
            kinds[col] = "binary" if uniq <= {0.0, 1.0} and len(uniq) == 2 else "quantitative"
            columns[col] = arr
    return PhenotypeTable(sample_ids=sample_ids, columns=columns, kinds=kinds)


@router.get("/datasets/{dataset_id}")
def get_dataset(dataset_id: str, project_id: Optional[str] = Query(None)) -> Dict[str, Any]:
    ds = registry.get_dataset(dataset_id, project_id)
    if not ds:
        raise HTTPException(404, "unknown dataset")
    return {"dataset": ds, "profile": registry.get_profile(dataset_id),
            "capabilities": registry.get_capabilities(dataset_id)}


@router.get("/datasets/{dataset_id}/capabilities")
def get_capabilities(dataset_id: str) -> Dict[str, Any]:
    return {"capabilities": registry.get_capabilities(dataset_id)}


@router.post("/datasets/{dataset_id}/override")
def create_override(dataset_id: str, req: OverrideRequest) -> Dict[str, Any]:
    try:
        return registry.create_override(dataset_id, req.analysis, req.justification,
                                        req.requested_by or config.DEFAULT_USER)
    except PermissionError as exc:
        raise HTTPException(403, str(exc))
    except (ValueError, KeyError) as exc:
        raise HTTPException(400, str(exc))


@router.delete("/datasets/{dataset_id}")
def delete_dataset(dataset_id: str, project_id: Optional[str] = Query(None)) -> Dict[str, Any]:
    try:
        return registry.delete_dataset(dataset_id, project_id)
    except KeyError as exc:
        raise HTTPException(404, str(exc))


# --------------------------------------------------------------------- jobs --
@router.post("/datasets/{dataset_id}/jobs")
def submit_job(dataset_id: str, req: JobRequest,
               project_id: Optional[str] = Query(None)) -> Dict[str, Any]:
    ds = registry.get_dataset(dataset_id, project_id)
    if not ds:
        raise HTTPException(404, "unknown dataset")
    try:
        return jobs.submit(ds["project_id"], dataset_id, req.analysis, req.spec,
                           req.submitted_by)
    except PermissionError as exc:
        # §3.3 — a locked analysis is refused, not warned about.
        raise HTTPException(403, detail={"message": str(exc), "locked": True})
    except KeyError as exc:
        raise HTTPException(400, str(exc))


@router.get("/datasets/{dataset_id}/jobs")
def list_jobs(dataset_id: str) -> Dict[str, Any]:
    return {"jobs": jobs.list_jobs(dataset_id=dataset_id)}


@router.get("/jobs/{job_id}")
def get_job(job_id: str) -> Dict[str, Any]:
    job = jobs.get_job(job_id)
    if not job:
        raise HTTPException(404, "unknown job")
    return job


@router.get("/jobs/{job_id}/result")
def get_result(job_id: str) -> Dict[str, Any]:
    res = jobs.get_result(job_id)
    if not res:
        job = jobs.get_job(job_id)
        if not job:
            raise HTTPException(404, "unknown job")
        raise HTTPException(409, "job is {} — no result yet".format(job["status"]))
    return res


@router.post("/jobs/{job_id}/cancel")
def cancel_job(job_id: str) -> Dict[str, Any]:
    return jobs.cancel(job_id)


@router.get("/estimate")
def estimate(analysis: str, n_samples: int, n_variants: int = 1) -> Dict[str, Any]:
    return jobs.estimate_cost(analysis, n_samples, n_variants)


# -------------------------------------------------------------------- power --
@router.post("/power")
def compute_power(req: PowerRequest) -> Dict[str, Any]:
    """§5.1 — power shown BEFORE the analysis runs."""
    alpha = req.alpha or power_mod.suggested_alpha(1, req.alpha_scope)["alpha"]
    if req.kind == "binary":
        return power_mod.power_case_control(req.n_cases, req.n_controls, req.maf,
                                            req.odds_ratio, alpha)
    return power_mod.power_quantitative(req.n or (req.n_cases + req.n_controls),
                                        req.maf, req.beta_sd, alpha)


# ------------------------------------------------------------- variant sets --
@router.get("/projects/{project_id}/variant-sets")
def list_variant_sets(project_id: str) -> Dict[str, Any]:
    return {"presets": variantset.PRESETS, "sets": variantset.list_sets(project_id)}


@router.post("/projects/{project_id}/variant-sets")
def save_variant_set(project_id: str, req: VariantSetRequest) -> Dict[str, Any]:
    return variantset.save_set(project_id, req.name, req.definition, req.created_by)
