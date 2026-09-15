"""Job queue and compute lifecycle — Part II §6.

    "Research Mode analyses cannot run in a web request."

    UI / API -> Analysis service -> Job queue -> Compute engine -> Result store

A GWAS over a million variants is minutes of CPU; a web worker blocked that
long is a web worker that has timed out. Every analysis is submitted, queued,
and polled.

Two things this module enforces that are easy to skip and expensive to retrofit:

CAPABILITY IS CHECKED AT SUBMISSION (§3). `submit` calls `registry.can_run`,
so the gate sits in front of the queue rather than inside each analysis. An
analysis module cannot be invoked on a dataset that does not support it, even
by a caller that forgot to check.

THE UNDERPOWERED STAMP IS CARRIED, NOT COPIED (§3.4). If the job ran under an
override, `underpowered` is set on the job row at submission and every result
and export reads it from there. It cannot be edited off a result payload.
"""
from __future__ import annotations

import json
import platform
import sys
import threading
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from typing import Any, Callable, Dict, List, Optional

from .. import config, db
from . import registry
from .capability import UNDERPOWERED_STAMP

# Analyses are CPU-bound numpy work, which releases the GIL inside BLAS, so a
# small thread pool is adequate here. A production deployment swaps this for a
# real broker without changing the job contract.
_POOL = ThreadPoolExecutor(max_workers=2, thread_name_prefix="research-job")
_LOCK = threading.RLock()

STATUS_QUEUED = "queued"
STATUS_VALIDATING = "validating"
STATUS_RUNNING = "running"
STATUS_POST = "post-processing"
STATUS_COMPLETE = "complete"
STATUS_FAILED = "failed"
STATUS_CANCELLED = "cancelled"

_CANCELLED: set = set()

# Registered analysis implementations: name -> callable(spec, context) -> dict
_REGISTRY: Dict[str, Callable[..., Dict[str, Any]]] = {}


def register(analysis: str):
    def deco(fn):
        _REGISTRY[analysis] = fn
        return fn
    return deco


def registered() -> List[str]:
    return sorted(_REGISTRY)


# Statuses a job can still move out of. Anything here is owned by a live thread.
ACTIVE_STATUSES = (STATUS_QUEUED, STATUS_VALIDATING, STATUS_RUNNING, STATUS_POST)


def reconcile_interrupted() -> int:
    """Fail any job still marked active at startup.

    Job state lives in the database but the thread running it does not. If the
    process dies mid-job — a crash, a restart, a laptop lid — the row stays
    "running" forever with nothing left to advance it, and the UI polls a job
    that will never finish. Nothing else notices, because from the outside a
    hung job and a slow one look identical.

    Startup is the one moment we can be certain no job is genuinely in flight,
    so it is the only safe place to make this call.
    """
    stale = db.rows(
        "SELECT job_id, analysis FROM analysis_job WHERE status IN ({})".format(
            ",".join("?" * len(ACTIVE_STATUSES))), list(ACTIVE_STATUSES))
    for row in stale:
        db.execute(
            "UPDATE analysis_job SET status = ?, finished_at = ?, "
            "error = ? WHERE job_id = ?",
            [STATUS_FAILED, datetime.utcnow(),
             "Interrupted: the server stopped while this job was running. "
             "Nothing was written; submit it again.", row["job_id"]])
    return len(stale)


# ------------------------------------------------------------ reproducibility
def reproducibility_record(spec: Dict[str, Any], dataset: Dict[str, Any],
                           seed: Optional[int]) -> Dict[str, Any]:
    """§6: "every job records container image, software versions, reference
    genome, annotation release, parameters, random seed, cohort version,
    dataset release." Without this a result cannot be reproduced, only
    re-obtained and hoped to match."""
    import numpy
    import scipy

    return {
        "software": {
            "python": sys.version.split()[0],
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "engine": "in-process numpy/scipy",
        },
        "container_image": config_get("CONTAINER_IMAGE", "local-dev"),
        "reference_genome": dataset.get("genome_build"),
        "annotation": {
            "source": dataset.get("annotation_source"),
            "version": dataset.get("annotation_version"),
            "transcript_set": dataset.get("transcript_set"),
        },
        "dataset": {
            "dataset_id": dataset.get("dataset_id"),
            "n_samples": dataset.get("n_samples"),
            "n_variants": dataset.get("n_variants"),
            "uploaded_at": str(dataset.get("uploaded_at")),
        },
        "kb_snapshot_id": db.meta_get("kb_snapshot_id", config.KB_SNAPSHOT_ID),
        "parameters": spec,
        "random_seed": seed,
    }


def config_get(key: str, default: str) -> str:
    import os
    return os.environ.get(key, default)


# ------------------------------------------------------------------ submit ---
def submit(project_id: str, dataset_id: str, analysis: str,
           spec: Dict[str, Any], submitted_by: Optional[str] = None,
           seed: Optional[int] = 20260913) -> Dict[str, Any]:
    """Queue an analysis. Refuses outright when the capability matrix says no.

    The refusal is deliberately not a warning the caller can proceed past —
    §3.3: "Never let a locked analysis run 'anyway with a warning.'"
    """
    dataset = registry.get_dataset(dataset_id, project_id)
    if not dataset:
        raise KeyError("unknown dataset '{}' in project '{}'".format(dataset_id, project_id))
    if analysis not in _REGISTRY:
        raise KeyError("no implementation registered for analysis '{}'".format(analysis))

    gate = registry.can_run(dataset_id, analysis)
    if not gate["allowed"]:
        raise PermissionError(gate["reason"])

    job_id = "job-" + uuid.uuid4().hex[:10]
    now = datetime.utcnow()
    repro = reproducibility_record(spec, dataset, seed)

    db.insert_rows(
        "analysis_job",
        ["job_id", "project_id", "dataset_id", "analysis", "spec_json", "status",
         "submitted_at", "started_at", "finished_at", "submitted_by", "error",
         "log", "underpowered", "override_id", "reproducibility_json"],
        [[job_id, project_id, dataset_id, analysis, json.dumps(spec), STATUS_QUEUED,
          now, None, None, submitted_by or config.DEFAULT_USER, None, "",
          bool(gate["underpowered"]), gate.get("override_id"),
          json.dumps(repro, default=str)]])

    _POOL.submit(_run, job_id)
    return get_job(job_id)


def _set(job_id: str, **fields) -> None:
    if not fields:
        return
    with _LOCK:
        sets = ", ".join("{} = ?".format(k) for k in fields)
        db.execute("UPDATE analysis_job SET {} WHERE job_id = ?".format(sets),
                   list(fields.values()) + [job_id])


def _log(job_id: str, line: str) -> None:
    with _LOCK:
        current = db.scalar("SELECT log FROM analysis_job WHERE job_id = ?",
                            [job_id], default="") or ""
        stamped = "[{}] {}".format(datetime.utcnow().strftime("%H:%M:%S"), line)
        db.execute("UPDATE analysis_job SET log = ? WHERE job_id = ?",
                   [(current + "\n" + stamped).strip(), job_id])


def _run(job_id: str) -> None:
    job = get_job(job_id)
    if not job:
        return
    try:
        if job_id in _CANCELLED:
            _set(job_id, status=STATUS_CANCELLED, finished_at=datetime.utcnow())
            return

        _set(job_id, status=STATUS_VALIDATING, started_at=datetime.utcnow())
        _log(job_id, "validating inputs")

        fn = _REGISTRY[job["analysis"]]
        spec = json.loads(job["spec_json"])

        _set(job_id, status=STATUS_RUNNING)
        _log(job_id, "running {}".format(job["analysis"]))

        context = {
            "job_id": job_id,
            "dataset_id": job["dataset_id"],
            "project_id": job["project_id"],
            "underpowered": bool(job["underpowered"]),
            "log": lambda m: _log(job_id, m),
        }
        result = fn(spec, context)

        _set(job_id, status=STATUS_POST)
        _log(job_id, "post-processing")

        # §3.4: the stamp travels with the result rather than being applied by
        # whatever renders it.
        if job["underpowered"]:
            result = dict(result)
            result["underpowered_stamp"] = UNDERPOWERED_STAMP
            result["override_justification"] = db.scalar(
                "SELECT justification FROM capability_override WHERE override_id = ?",
                [job["override_id"]], default="")

        diagnostics = result.pop("diagnostics", None) if isinstance(result, dict) else None
        db.execute("DELETE FROM analysis_result WHERE job_id = ?", [job_id])
        db.insert_rows(
            "analysis_result",
            ["job_id", "analysis", "result_json", "diagnostics_json", "created_at",
             "output_class"],
            [[job_id, job["analysis"], json.dumps(result, default=_json_default),
              json.dumps(diagnostics, default=_json_default) if diagnostics else None,
              datetime.utcnow(), "RESEARCH"]])

        _set(job_id, status=STATUS_COMPLETE, finished_at=datetime.utcnow())
        _log(job_id, "complete")

    except Exception as exc:
        _set(job_id, status=STATUS_FAILED, finished_at=datetime.utcnow(),
             error="{}: {}".format(type(exc).__name__, exc))
        _log(job_id, "FAILED: {}".format(exc))
        _log(job_id, traceback.format_exc(limit=6))


def _json_default(o):
    import numpy as np
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        f = float(o)
        # JSON has no NaN/Infinity; emitting them produces a payload the
        # browser's JSON.parse rejects outright.
        return None if (f != f or f in (float("inf"), float("-inf"))) else f
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, datetime):
        return o.isoformat()
    return str(o)


# ------------------------------------------------------------------ reading --
def get_job(job_id: str) -> Optional[Dict[str, Any]]:
    return db.row("SELECT * FROM analysis_job WHERE job_id = ?", [job_id])


def list_jobs(project_id: Optional[str] = None, dataset_id: Optional[str] = None,
              limit: int = 100) -> List[Dict[str, Any]]:
    sql = "SELECT * FROM analysis_job WHERE 1=1"
    params: List[Any] = []
    if project_id:
        sql += " AND project_id = ?"
        params.append(project_id)
    if dataset_id:
        sql += " AND dataset_id = ?"
        params.append(dataset_id)
    sql += " ORDER BY submitted_at DESC LIMIT ?"
    params.append(limit)
    return db.rows(sql, params)


def get_result(job_id: str) -> Optional[Dict[str, Any]]:
    row = db.row("SELECT * FROM analysis_result WHERE job_id = ?", [job_id])
    if not row:
        return None
    out = {
        "job_id": job_id,
        "analysis": row["analysis"],
        "output_class": row["output_class"],
        "result": json.loads(row["result_json"]),
        "diagnostics": json.loads(row["diagnostics_json"]) if row["diagnostics_json"] else None,
        "created_at": row["created_at"],
    }
    job = get_job(job_id)
    if job:
        out["underpowered"] = bool(job["underpowered"])
        out["reproducibility"] = json.loads(job["reproducibility_json"] or "{}")
    return out


def cancel(job_id: str) -> Dict[str, Any]:
    _CANCELLED.add(job_id)
    job = get_job(job_id)
    if job and job["status"] in (STATUS_QUEUED, STATUS_VALIDATING):
        _set(job_id, status=STATUS_CANCELLED, finished_at=datetime.utcnow())
    return get_job(job_id)


def wait(job_id: str, timeout: float = 120.0, poll: float = 0.05) -> Dict[str, Any]:
    """Block until a job reaches a terminal state. For tests and the CLI — the
    HTTP API polls instead."""
    import time
    deadline = time.time() + timeout
    while time.time() < deadline:
        job = get_job(job_id)
        if job and job["status"] in (STATUS_COMPLETE, STATUS_FAILED, STATUS_CANCELLED):
            return job
        time.sleep(poll)
    raise TimeoutError("job {} did not finish within {}s".format(job_id, timeout))


def estimate_cost(analysis: str, n_samples: int, n_variants: int) -> Dict[str, Any]:
    """§6: "Estimated cost and runtime shown before submission."

    Deliberately coarse. A wrong-by-2x estimate still tells a researcher
    whether this is a coffee break or an overnight run, which is the decision
    they are actually making.
    """
    # Calibrated against observed runs on this engine: an association scan of
    # 57,333 variants x 2,400 samples took ~127s, i.e. ~9e-7 s per
    # sample-variant. Association was previously excluded from the
    # variant-scaling term, which estimated that same job at one second —
    # an estimate that wrong is worse than none, because it is the number a
    # researcher uses to decide whether to wait.
    unit = {"association": 9e-7, "gwas": 1.1e-6, "burden": 3e-6,
            "prs": 2e-7, "survival": 5e-6}.get(analysis, 9e-7)
    scans_variants = analysis in ("association", "gwas", "burden", "prs")
    work = max(1, n_samples) * max(1, n_variants if scans_variants else 1)
    seconds = max(1.0, work * unit)
    profile = ("small" if seconds < 60 else
               "memory-heavy" if seconds < 900 else "distributed")
    return {
        "estimated_seconds": round(seconds, 1),
        "estimated_label": _humanise(seconds),
        "compute_profile": profile,
        "note": "Estimate only; actual runtime depends on missingness and covariates.",
    }


def _humanise(seconds: float) -> str:
    if seconds < 90:
        return "~{:.0f} seconds".format(seconds)
    if seconds < 5400:
        return "~{:.0f} minutes".format(seconds / 60)
    return "~{:.1f} hours".format(seconds / 3600)
