"""Runtime configuration for the germline cohort analytics service."""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# Analytics store. Rebuildable from source at any time (spec E03.7).
DB_PATH = Path(os.environ.get("COHORT_DB", ROOT / "data" / "germline.duckdb"))

# Where `ingest varimat` looks for per-sample VariMAT files by default.
VARIMAT_DIR = Path(os.environ.get("VARIMAT_DIR", ROOT / "data" / "varimat"))

# Clinical/LIMS sidecar. VariMAT carries no subject, family, consent or
# indication data (spec §3.5) — it must come from the clinical system.
CLINICAL_DIR = Path(os.environ.get("CLINICAL_DIR", ROOT / "data" / "clinical"))

FRONTEND_DIR = ROOT / "frontend"
EXPORT_DIR = Path(os.environ.get("EXPORT_DIR", ROOT / "data" / "exports"))

TENANT_ID = int(os.environ.get("TENANT_ID", "11"))

# Knowledgebase snapshot pin. Without this, re-running a cohort after a ClinVar
# update yields a different number with no explanation (spec C02).
KB_SNAPSHOT_ID = os.environ.get("KB_SNAPSHOT_ID", "kb-2026-07-28")

PROFILE = "germline"

# --- governance constants (spec §3.7) ---------------------------------------
SMALL_CELL_THRESHOLD = 5     # counts below this render as "<5"
DENOM_MIN = 30               # denominator warning floor
PROVISIONAL_MAX = 0.25       # >25% provisional scope triggers the warning
REVIEWABLE_MAX_GNOMAD_AF = 0.01
TABLE_RENDER_CAP = 400

# The analyst identity stamped on audit entries and exports. In production this
# comes from the auth layer.
DEFAULT_USER = os.environ.get("COHORT_USER", "analyst@impactomics")
