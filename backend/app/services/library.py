"""Saved cohort definitions — the rail's "Saved cohorts" section (§2.6).

A saved cohort stores CRITERIA, never a member list. Re-running it against the
same store and the same KB snapshot must reproduce the same member_hash
(spec E03.6: "Any saved cohort re-runs to an identical number"); re-running it
against a changed store must produce a DIFFERENT member_hash, visibly, rather
than silently drifting.

The built-in set exists so the tool opens onto real germline questions rather
than an empty builder. Each one is a question an analyst actually asks.
"""
from __future__ import annotations

import json
import uuid
from datetime import datetime
from typing import Any, Dict, List, Optional

from .. import config, db
from . import cohort as cohort_svc

BUILTIN: List[Dict[str, Any]] = [
    {
        "cohort_id": "builtin-all",
        "name": "Everyone (one per family)",
        "description": ("Every consented subject, one per family. The default "
                        "denominator for anything else."),
        "view": "g_dashboard",
        "criteria": {},
    },
    {
        "cohort_id": "builtin-hboc",
        "name": "HBOC panel · BRCA1/2 P/LP carriers",
        "description": ("Pathogenic and likely pathogenic carriers on the hereditary "
                        "breast and ovarian cancer gene set."),
        "view": "g_carrier",
        "criteria": {
            "gene_set": "HBOC (breast/ovarian)",
            "classification": ["Pathogenic", "Likely pathogenic"],
        },
    },
    {
        "cohort_id": "builtin-lynch",
        "name": "Lynch / MMR P/LP carriers",
        "description": "Mismatch-repair carriers, for surveillance-programme sizing.",
        "view": "g_carrier",
        "criteria": {
            "gene_set": "Lynch / MMR",
            "classification": ["Pathogenic", "Likely pathogenic"],
        },
    },
    {
        "cohort_id": "builtin-carrier-screen",
        "name": "Reproductive carrier screening · all",
        "description": ("Carrier-screening referrals. Carrier rate is per gene, not "
                        "per subject — a positive here is carrier status, not a "
                        "diagnosis."),
        "view": "g_carrier",
        "criteria": {"indication": ["Reproductive carrier screening"]},
    },
    {
        "cohort_id": "builtin-biallelic",
        "name": "Biallelic P/LP · recessive genes",
        "description": ("Homozygous and compound-heterozygous findings in AR genes — "
                        "consistent with affected status, not carrier status."),
        "view": "g_zygosity",
        "criteria": {
            "classification": ["Pathogenic", "Likely pathogenic"],
            "zygosity": ["Homozygous", "Compound heterozygous"],
            "inheritance": ["AR", "AR/AD"],
        },
    },
    {
        "cohort_id": "builtin-xlinked",
        "name": "X-linked · hemizygous males",
        "description": ("Derived hemizygous calls. These are affected individuals, "
                        "not carriers — and ZYGOSITY in the source file cannot tell "
                        "you that."),
        "view": "g_zygosity",
        "criteria": {
            "classification": ["Pathogenic", "Likely pathogenic"],
            "zygosity": ["Hemizygous"],
        },
    },
    {
        "cohort_id": "builtin-vus",
        "name": "All VUS · reclassification queue",
        "description": ("The curation work queue, scoped to the reviewable subset and "
                        "ranked by recontact leverage and staleness."),
        "view": "g_vus",
        "criteria": {"classification": ["Uncertain significance"]},
    },
    {
        "cohort_id": "builtin-founder",
        "name": "Founder-variant candidates",
        "description": ("Recurrent P/LP alleles enriched over gnomAD. Also surfaces "
                        "implausible internal frequencies, which are a QA signal."),
        "view": "g_popfreq",
        "criteria": {
            "founder_only": True,
            "classification": ["Pathogenic", "Likely pathogenic"],
        },
    },
    {
        "cohort_id": "builtin-sf",
        "name": "ACMG SF v3.3 · consented only",
        "description": ("Secondary findings among subjects with explicit consent. "
                        "Declining subjects are excluded at the query level."),
        "view": "g_sf",
        "criteria": {
            "sf_only": True,
            "classification": ["Pathogenic", "Likely pathogenic"],
        },
    },
    {
        "cohort_id": "builtin-exome",
        "name": "Exome & genome referrals",
        "description": ("Broad tests only. Yield here is expected to be lower than "
                        "targeted panels — that is the phenotype-relevance gate "
                        "working, not a failure."),
        "view": "g_dashboard",
        "criteria": {"test_code": ["MG-CES", "MG-WES", "MG-WGS"]},
    },
    {
        "cohort_id": "builtin-families",
        "name": "Multi-member families · probands OFF",
        "description": ("Deliberately includes related subjects, to show what family "
                        "clustering does to a carrier rate. Rates here must not be "
                        "quoted."),
        "view": "g_carrier",
        "criteria": {"probands_only": False},
    },
    {
        "cohort_id": "builtin-provisional",
        "name": "Incomplete-scope test codes",
        "description": ("Runs on test codes whose reportable scope the registry does "
                        "not declare. Their denominators are provisional."),
        "view": "g_runs",
        "criteria": {"test_code": ["MG-HEARING-96", "MG-NDD-58"]},
    },
]


def seed_builtins() -> None:
    db.execute("DELETE FROM cohort_def WHERE builtin")
    db.insert_rows(
        "cohort_def",
        ["cohort_id", "tenant_id", "profile", "name", "description",
         "criteria_json", "version", "created_by", "created_at", "builtin"],
        [[c["cohort_id"], config.TENANT_ID, config.PROFILE, c["name"],
          c["description"], json.dumps({"criteria": c["criteria"], "view": c["view"]}),
          1, "system", datetime.utcnow(), True] for c in BUILTIN],
    )


def list_cohorts() -> List[Dict[str, Any]]:
    rows = db.rows("""
        SELECT cohort_id, name, description, criteria_json, version, created_by,
               created_at, builtin
        FROM cohort_def WHERE profile = ?
        -- Built-ins sort by id so "builtin-all" leads: the first starting
        -- point a user sees should be the broadest one, not whichever row the
        -- seed happened to insert last.
        ORDER BY builtin DESC,
                 CASE WHEN builtin THEN cohort_id END ASC,
                 created_at DESC
    """, [config.PROFILE])
    for r in rows:
        payload = json.loads(r.pop("criteria_json"))
        r["criteria"] = payload.get("criteria", {})
        r["view"] = payload.get("view", "g_dashboard")
    return rows


def get(cohort_id: str) -> Optional[Dict[str, Any]]:
    for c in list_cohorts():
        if c["cohort_id"] == cohort_id:
            return c
    return None


def save(name: str, criteria: Dict[str, Any], view: str = "g_dashboard",
         description: str = "", user: Optional[str] = None) -> Dict[str, Any]:
    cid = "coh-" + uuid.uuid4().hex[:10]
    normalised = cohort_svc.normalise(criteria)
    db.insert_rows(
        "cohort_def",
        ["cohort_id", "tenant_id", "profile", "name", "description",
         "criteria_json", "version", "created_by", "created_at", "builtin"],
        [[cid, config.TENANT_ID, config.PROFILE, name, description,
          json.dumps({"criteria": normalised, "view": view}), 1,
          user or config.DEFAULT_USER, datetime.utcnow(), False]],
    )
    return {"cohort_id": cid, "name": name, "criteria": normalised, "view": view}


def delete(cohort_id: str) -> bool:
    if cohort_id.startswith("builtin-"):
        return False
    db.execute("DELETE FROM cohort_def WHERE cohort_id = ? AND NOT builtin", [cohort_id])
    return True


def verify_reproducible(cohort_id: str) -> Dict[str, Any]:
    """Resolve a saved cohort twice and assert the member hash is stable.

    Spec E03.6. Cheap to run, and it is the only thing standing between a saved
    cohort and silent membership drift.
    """
    rec = get(cohort_id)
    if not rec:
        return {"ok": False, "error": "unknown cohort"}
    a = cohort_svc.resolve(rec["criteria"], name=rec["name"])
    # Pin the anchor from the first resolve. A cohort with a collection window
    # legitimately changes membership when new samples arrive; re-running
    # against the same anchor is what actually tests reproducibility, rather
    # than reporting that moving data moved.
    b = cohort_svc.resolve(rec["criteria"], name=rec["name"], as_of=a.as_of)
    return {
        "ok": a.member_hash == b.member_hash and a.criteria_hash == b.criteria_hash,
        "cohort_id": cohort_id, "name": rec["name"],
        "member_hash": a.member_hash, "criteria_hash": a.criteria_hash,
        "member_count": a.n_subjects, "collection_anchor": a.as_of,
    }
