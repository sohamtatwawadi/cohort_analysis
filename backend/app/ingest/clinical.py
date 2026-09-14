"""Clinical / LIMS sidecar loader.

Spec §3.5 lists what VariMAT does not contain and therefore must be integrated:

    family relationships, consent class, indication, phenotype, affected status

These are not optional metadata. Consent gates cohort membership (§G01) and
family structure gates independence (probands-only, §G01). A cohort built on
genomic data alone is not a germline cohort — it is a pile of samples.

Sidecar format: one CSV or JSON file keyed by `sample_id`. Anything the lab's
LIMS can export will map onto these columns; unknown columns are ignored.
"""
from __future__ import annotations

import csv
import json
from datetime import date, datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..reference.tests import (AGE_BUCKETS, CONSENT_CLASSES, CONSENT_WITHDRAWN,
                               TEST_CODES)

SIDECAR_COLUMNS = [
    "sample_id", "subject_id", "family_id", "mrn", "sex", "age", "ancestry",
    "consent_class", "relation", "is_proband", "indication", "affected_status",
    "family_history", "referral_source", "phenotype_hpo", "test_code",
    "sample_type", "collection_date", "pipeline_version", "reference_build",
    "mean_depth", "pct_bases_20x", "qc_status",
]

_TRUE = {"1", "true", "yes", "y", "t"}


def age_bucket(age: Optional[int]) -> Optional[str]:
    if age is None:
        return None
    if age < 1:
        return AGE_BUCKETS[0]
    if age < 12:
        return AGE_BUCKETS[1]
    if age < 18:
        return AGE_BUCKETS[2]
    if age < 40:
        return AGE_BUCKETS[3]
    if age < 60:
        return AGE_BUCKETS[4]
    return AGE_BUCKETS[5]


def _coerce(rec: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(rec)
    for k in ("age", "mean_depth", "pct_bases_20x"):
        v = out.get(k)
        if v in (None, ""):
            out[k] = None
        else:
            try:
                out[k] = float(v) if k != "age" else int(float(v))
            except (TypeError, ValueError):
                out[k] = None
    pro = out.get("is_proband")
    if isinstance(pro, str):
        out["is_proband"] = pro.strip().lower() in _TRUE
    elif pro is None:
        out["is_proband"] = (out.get("relation") or "Proband") == "Proband"

    cd = out.get("collection_date")
    if isinstance(cd, str) and cd:
        try:
            out["collection_date"] = datetime.strptime(cd[:10], "%Y-%m-%d").date()
        except ValueError:
            out["collection_date"] = None

    out["age_bucket"] = out.get("age_bucket") or age_bucket(out.get("age"))

    # An unrecognised consent value must not silently become permissive.
    # Unknown consent is treated as withdrawn: the subject is excluded from
    # every cohort until the LIMS says otherwise (spec §G01 step 2).
    consent = out.get("consent_class")
    if consent not in CONSENT_CLASSES:
        out["consent_class"] = CONSENT_WITHDRAWN
        out["consent_unmapped"] = consent

    code = out.get("test_code")
    if code and code in TEST_CODES and not out.get("indication"):
        out["indication"] = TEST_CODES[code].indication
    return out


def stub_record(sample_id: str) -> Dict[str, Any]:
    """Fallback for a VariMAT file with no sidecar row.

    Deliberately conservative: consent is Withdrawn, so the sample loads and is
    inspectable but cannot enter any cohort until the clinical system supplies
    a consent class. Failing closed is the only safe default here.
    """
    return _coerce({
        "sample_id": sample_id,
        "subject_id": "SJ-" + sample_id,
        "family_id": "FAM-" + sample_id,
        "consent_class": None,          # -> Withdrawn
        "relation": "Proband",
        "is_proband": True,
        "qc_status": "Unknown",
    })


def load_sidecar(path: Path) -> Dict[str, Dict[str, Any]]:
    """Load every sidecar file under `path`, keyed by sample_id."""
    path = Path(path)
    out: Dict[str, Dict[str, Any]] = {}
    if not path.exists():
        return out

    files: List[Path] = [path] if path.is_file() else sorted(
        list(path.glob("*.csv")) + list(path.glob("*.json")))

    for f in files:
        if f.suffix == ".json":
            data = json.loads(f.read_text())
            records = data if isinstance(data, list) else data.get("subjects", [])
        else:
            with f.open(newline="", encoding="utf-8") as fh:
                records = list(csv.DictReader(fh))
        for rec in records:
            sid = (rec.get("sample_id") or "").strip()
            if not sid:
                continue
            out[sid] = _coerce(rec)
    return out


def write_sidecar_template(path: Path) -> Path:
    """Emit a header-only CSV so the lab knows exactly what to export."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        csv.writer(fh).writerow(SIDECAR_COLUMNS)
    return path
