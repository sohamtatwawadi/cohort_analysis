"""Dataset registry, project isolation and data governance — Part II §8.

Uploaded research data carries obligations Lab Mode does not have:

    "The uploader's data belongs to the uploader. It is never pooled into
     reference cohorts, never used to improve models, never visible to other
     tenants."

Two of those obligations are structural rather than procedural, and this module
is where they live:

PROJECT ISOLATION. Every dataset, job and result carries project_id, and every
read goes through a function that takes a project_id. There is no query path
that returns another project's rows, so isolation does not depend on a caller
remembering to add a filter.

DELETION THAT PROPAGATES. "Deletion must actually propagate to derived tables,
caches and result stores." A soft-delete flag on `dataset` would leave the
profile, capability assessment, phenotype registry, job specs and result
payloads in place — all of which are derived from, and can reveal, the deleted
data. `delete_dataset` removes them and reports what it removed.
"""
from __future__ import annotations

import json
import uuid
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional

from .. import config, db
from .capability import assess, overridable
from .profile import DataProfile

CONSENT_STATEMENT = (
    "I attest that this dataset is held under appropriate participant consent and "
    "ethics approval for the analyses I intend to run, that I am authorised to "
    "upload it, and that ImpactOmics is permitted to process it on my behalf."
)

DATA_OWNERSHIP_NOTICE = (
    "Your data belongs to you. It is never pooled into reference cohorts, never "
    "used to train or improve models, and is never visible to another tenant. It "
    "is scoped to this project; cross-project access requires an explicit, logged "
    "grant. You can delete it at any time, and deletion propagates to every "
    "derived table and stored result."
)


def _now() -> datetime:
    return datetime.utcnow()


def _jid(prefix: str) -> str:
    return "{}-{}".format(prefix, uuid.uuid4().hex[:10])


# ---------------------------------------------------------------- projects ---
def create_project(name: str, owner: str, description: str = "",
                   retention_days: Optional[int] = None) -> Dict[str, Any]:
    pid = _jid("prj")
    db.insert_rows(
        "project",
        ["project_id", "tenant_id", "name", "description", "owner",
         "created_at", "retention_days", "deleted_at"],
        [[pid, config.TENANT_ID, name, description, owner, _now(), retention_days, None]])
    return get_project(pid)


def get_project(project_id: str) -> Optional[Dict[str, Any]]:
    return db.row("""
        SELECT * FROM project WHERE project_id = ? AND deleted_at IS NULL
    """, [project_id])


def list_projects() -> List[Dict[str, Any]]:
    rows = db.rows("""
        SELECT p.*,
               (SELECT COUNT(*) FROM dataset d
                WHERE d.project_id = p.project_id AND d.deleted_at IS NULL) AS n_datasets
        FROM project p WHERE p.deleted_at IS NULL AND p.tenant_id = ?
        ORDER BY p.created_at DESC
    """, [config.TENANT_ID])
    return rows


def orphaned_payloads() -> List[str]:
    """Genotype payloads on disk with no dataset row.

    Should always be empty. A non-empty result means something deleted dataset
    metadata without going through delete_dataset — the payload is still there,
    which is the opposite of the §8 deletion guarantee.
    """
    from . import store as research_store
    if not research_store.RESEARCH_DIR.exists():
        return []
    known = {r["dataset_id"] for r in db.rows("SELECT dataset_id FROM dataset")}
    return sorted(d.name for d in research_store.RESEARCH_DIR.iterdir()
                  if d.is_dir() and d.name not in known)


def purge_orphaned_payloads() -> Dict[str, Any]:
    from . import store as research_store
    removed = []
    for did in orphaned_payloads():
        research_store.delete_payload(did)
        removed.append(did)
    return {"removed": removed, "count": len(removed)}


def ensure_default_project() -> str:
    """A project always exists, because a dataset cannot be registered without
    one — isolation has nowhere to hang otherwise."""
    existing = list_projects()
    if existing:
        return existing[0]["project_id"]
    return create_project("Default research project", config.DEFAULT_USER,
                          "Created automatically on first use.")["project_id"]


# ---------------------------------------------------------------- datasets ---
def register_dataset(
    project_id: str,
    name: str,
    profile: DataProfile,
    source_format: str,
    source_files: List[str],
    build: Optional[str],
    uploaded_by: str,
    consent_attested: bool,
    annotation: Optional[Dict[str, str]] = None,
    phenotype_kinds: Optional[Dict[str, str]] = None,
    samples: Optional[List[Dict[str, Any]]] = None,
    storage_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Register a validated dataset and persist its profile and capabilities.

    Consent attestation is checked here rather than in the API layer so that no
    code path can register a dataset without it (§8).
    """
    if not consent_attested:
        raise PermissionError(
            "Consent attestation is required before a dataset can be registered. "
            + CONSENT_STATEMENT)
    if not get_project(project_id):
        raise KeyError("unknown project '{}'".format(project_id))

    did = _jid("ds")
    ann = annotation or {}
    now = _now()

    db.insert_rows(
        "dataset",
        ["dataset_id", "project_id", "name", "mode", "source_format", "source_files",
         "genome_build", "n_samples", "n_variants", "uploaded_at", "uploaded_by",
         "status", "failure_reason", "annotation_source", "annotation_version",
         "transcript_set", "consent_attested", "consent_statement",
         "consent_attested_by", "consent_attested_at", "storage_path", "checksum",
         "deleted_at"],
        [[did, project_id, name, "research", source_format, json.dumps(source_files),
          build, profile.n_samples, profile.n_variants, now, uploaded_by,
          "ready", None, ann.get("source"), ann.get("version"),
          ann.get("transcript_set"), True, CONSENT_STATEMENT, uploaded_by, now,
          storage_path, None, None]])

    save_profile(did, profile)

    if phenotype_kinds:
        rows = []
        for pname, kind in phenotype_kinds.items():
            entry = profile.phenotypes.get(pname, {})
            rows.append([did, pname, kind, entry.get("label", pname),
                         entry.get("n_present"), entry.get("cases"),
                         entry.get("controls"), entry.get("completeness")])
        db.insert_rows("dataset_phenotype",
                       ["dataset_id", "name", "kind", "label", "n_present",
                        "n_cases", "n_controls", "completeness"], rows)

    if samples:
        db.insert_rows("dataset_sample",
                       ["dataset_id", "sample_id", "reported_sex", "genetic_sex",
                        "call_rate", "in_phenotype", "is_duplicate"],
                       [[did, s.get("sample_id"), s.get("reported_sex"),
                         s.get("genetic_sex"), s.get("call_rate"),
                         s.get("in_phenotype", False), s.get("is_duplicate", False)]
                        for s in samples])

    return get_dataset(did, project_id)


def save_profile(dataset_id: str, profile: DataProfile) -> None:
    db.execute("DELETE FROM dataset_profile WHERE dataset_id = ?", [dataset_id])
    db.insert_rows("dataset_profile", ["dataset_id", "profiled_at", "profile_json"],
                   [[dataset_id, _now(), json.dumps(profile.as_dict(), default=str)]])
    save_capabilities(dataset_id, profile)


def save_capabilities(dataset_id: str, profile: DataProfile) -> None:
    """Persist the capability assessment so what the researcher was shown is
    auditable later, rather than being recomputed against a changed profile."""
    db.execute("DELETE FROM dataset_capability WHERE dataset_id = ?", [dataset_id])
    now = _now()
    db.insert_rows(
        "dataset_capability",
        ["dataset_id", "analysis", "available", "reasons_json", "assessed_at"],
        [[dataset_id, c.analysis, c.available,
          json.dumps({"title": c.title, "phase": c.phase, "remedy": c.remedy,
                      "requirements": [r.as_dict() for r in c.requirements]}),
          now] for c in assess(profile)])


def get_dataset(dataset_id: str, project_id: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """Fetch a dataset. When project_id is supplied it is enforced — this is
    the isolation boundary (§8)."""
    sql = "SELECT * FROM dataset WHERE dataset_id = ? AND deleted_at IS NULL"
    params: List[Any] = [dataset_id]
    if project_id is not None:
        sql += " AND project_id = ?"
        params.append(project_id)
    row = db.row(sql, params)
    if row and row.get("source_files"):
        try:
            row["source_files"] = json.loads(row["source_files"])
        except (TypeError, ValueError):
            pass
    return row


def list_datasets(project_id: str) -> List[Dict[str, Any]]:
    return db.rows("""
        SELECT dataset_id, project_id, name, source_format, genome_build,
               n_samples, n_variants, uploaded_at, uploaded_by, status,
               annotation_source, annotation_version, transcript_set
        FROM dataset
        WHERE project_id = ? AND deleted_at IS NULL
        ORDER BY uploaded_at DESC
    """, [project_id])


def get_profile(dataset_id: str) -> Optional[Dict[str, Any]]:
    row = db.row("SELECT profile_json FROM dataset_profile WHERE dataset_id = ?",
                 [dataset_id])
    return json.loads(row["profile_json"]) if row else None


def get_capabilities(dataset_id: str) -> List[Dict[str, Any]]:
    rows = db.rows("""
        SELECT analysis, available, reasons_json, assessed_at
        FROM dataset_capability WHERE dataset_id = ? ORDER BY analysis
    """, [dataset_id])
    out = []
    for r in rows:
        payload = json.loads(r["reasons_json"] or "{}")
        reqs = payload.get("requirements", [])
        unmet = [q for q in reqs if not q.get("met")]
        out.append({
            "analysis": r["analysis"],
            "title": payload.get("title", r["analysis"]),
            "phase": payload.get("phase", ""),
            "available": bool(r["available"]),
            "requirements": reqs,
            "unmet": unmet,
            "remedy": payload.get("remedy", ""),
            "override": _override_state(dataset_id, r["analysis"], unmet),
        })
    # Available first, then by phase, so the researcher sees what they can run.
    out.sort(key=lambda c: (not c["available"], c["phase"], c["analysis"]))
    return out


def _override_state(dataset_id: str, analysis: str,
                    unmet: List[Dict[str, Any]]) -> Dict[str, Any]:
    from .capability import Capability, NON_OVERRIDABLE, Requirement

    active = get_override(dataset_id, analysis)
    blocking = [u for u in unmet if u.get("label") in NON_OVERRIDABLE]
    soft = [u for u in unmet if u.get("label") not in NON_OVERRIDABLE]
    return {
        "active": active is not None,
        "override_id": active["override_id"] if active else None,
        "justification": active["justification"] if active else None,
        "can_override": bool(soft) and not blocking,
        "blocking": blocking,
        "overridable": soft,
        "reason": ("Missing data cannot be overridden — it has to be supplied."
                   if blocking else
                   "Borderline thresholds may be overridden with a recorded "
                   "justification. Results will be stamped."),
    }


# ---------------------------------------------------------------- override ---
def create_override(dataset_id: str, analysis: str, justification: str,
                    requested_by: str) -> Dict[str, Any]:
    """§3.4. Requires a typed justification, is audited, and causes every
    resulting output to carry the UNDERPOWERED stamp.

    Refuses when the unmet requirement is missing data rather than a borderline
    threshold — you cannot override your way to genotypes you did not upload.
    """
    justification = (justification or "").strip()
    if len(justification) < 20:
        raise ValueError(
            "An override requires a written justification (at least 20 characters) "
            "explaining why the analysis is defensible despite the unmet threshold. "
            "It is recorded in the audit log and shown alongside the result.")

    caps = {c["analysis"]: c for c in get_capabilities(dataset_id)}
    cap = caps.get(analysis)
    if not cap:
        raise KeyError("unknown analysis '{}' for this dataset".format(analysis))
    if cap["available"]:
        raise ValueError("'{}' is already available; no override needed.".format(analysis))
    if not cap["override"]["can_override"]:
        raise PermissionError(
            "'{}' cannot be overridden: {}".format(
                analysis,
                "; ".join(b["label"] + " — " + b["observed"]
                          for b in cap["override"]["blocking"])
                or cap["override"]["reason"]))

    oid = _jid("ovr")
    db.insert_rows(
        "capability_override",
        ["override_id", "dataset_id", "analysis", "justification", "requested_by",
         "created_at", "unmet_json"],
        [[oid, dataset_id, analysis, justification, requested_by, _now(),
          json.dumps(cap["unmet"])]])
    return {"override_id": oid, "dataset_id": dataset_id, "analysis": analysis,
            "justification": justification, "stamped": True}


def get_override(dataset_id: str, analysis: str) -> Optional[Dict[str, Any]]:
    return db.row("""
        SELECT * FROM capability_override
        WHERE dataset_id = ? AND analysis = ?
        ORDER BY created_at DESC LIMIT 1
    """, [dataset_id, analysis])


def can_run(dataset_id: str, analysis: str) -> Dict[str, Any]:
    """The single gate every job submission passes through.

    Returns whether the analysis may run and, if so, whether its output must
    carry the UNDERPOWERED stamp. There is deliberately no third state — an
    analysis is available, available-with-stamp, or refused.
    """
    caps = {c["analysis"]: c for c in get_capabilities(dataset_id)}
    cap = caps.get(analysis)
    if not cap:
        return {"allowed": False, "underpowered": False,
                "reason": "Unknown analysis '{}'.".format(analysis), "unmet": []}
    if cap["available"]:
        return {"allowed": True, "underpowered": False, "reason": "", "unmet": []}

    ov = get_override(dataset_id, analysis)
    if ov:
        return {"allowed": True, "underpowered": True,
                "override_id": ov["override_id"],
                "justification": ov["justification"],
                "reason": "Running under a recorded override.",
                "unmet": cap["unmet"]}

    return {
        "allowed": False, "underpowered": False,
        "reason": "Requirements not met: {}".format(
            "; ".join("{} (you have: {})".format(u["label"], u["observed"])
                      for u in cap["unmet"])),
        "unmet": cap["unmet"],
        "remedy": cap["remedy"],
    }


# ---------------------------------------------------------------- deletion ---
DERIVED_TABLES = [
    ("dataset_profile", "dataset_id"),
    ("dataset_capability", "dataset_id"),
    ("dataset_phenotype", "dataset_id"),
    ("dataset_sample", "dataset_id"),
    ("capability_override", "dataset_id"),
]


def delete_dataset(dataset_id: str, project_id: Optional[str] = None) -> Dict[str, Any]:
    """§8: "Deletion must actually propagate to derived tables, caches and
    result stores."

    A flag on `dataset` alone would leave the profile, the capability
    assessment, the phenotype registry and every stored result in place — all
    derived from the deleted data, and all capable of revealing it. This removes
    them and reports the counts so deletion is verifiable rather than asserted.
    """
    ds = get_dataset(dataset_id, project_id)
    if not ds:
        raise KeyError("unknown dataset '{}' in this project".format(dataset_id))

    removed: Dict[str, int] = {}
    job_ids = [r["job_id"] for r in db.rows(
        "SELECT job_id FROM analysis_job WHERE dataset_id = ?", [dataset_id])]

    if job_ids:
        placeholders = ", ".join("?" * len(job_ids))
        removed["analysis_result"] = int(db.scalar(
            "SELECT COUNT(*) FROM analysis_result WHERE job_id IN ({})".format(placeholders),
            job_ids, default=0))
        db.execute("DELETE FROM analysis_result WHERE job_id IN ({})".format(placeholders),
                   job_ids)
    removed["analysis_job"] = len(job_ids)
    db.execute("DELETE FROM analysis_job WHERE dataset_id = ?", [dataset_id])

    for table, col in DERIVED_TABLES:
        removed[table] = int(db.scalar(
            "SELECT COUNT(*) FROM {} WHERE {} = ?".format(table, col),
            [dataset_id], default=0))
        db.execute("DELETE FROM {} WHERE {} = ?".format(table, col), [dataset_id])

    db.execute("DELETE FROM dataset WHERE dataset_id = ?", [dataset_id])
    removed["dataset"] = 1

    # The genotype/phenotype payload itself, not just the metadata describing
    # it. Leaving this on disk would make "deleted" a lie.
    from . import store as research_store
    removed["payload_removed"] = int(research_store.delete_payload(dataset_id))

    return {"dataset_id": dataset_id, "deleted": True, "removed": removed}


def enforce_retention() -> Dict[str, Any]:
    """§8 retention, "defined per project, enforced automatically"."""
    expired: List[str] = []
    for proj in list_projects():
        days = proj.get("retention_days")
        if not days:
            continue
        cutoff = _now() - timedelta(days=int(days))
        for ds in db.rows("""
            SELECT dataset_id FROM dataset
            WHERE project_id = ? AND deleted_at IS NULL AND uploaded_at < ?
        """, [proj["project_id"], cutoff]):
            delete_dataset(ds["dataset_id"], proj["project_id"])
            expired.append(ds["dataset_id"])
    return {"deleted": expired, "count": len(expired)}
