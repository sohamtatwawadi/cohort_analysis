"""Shared fixtures. Builds one small synthetic store for the whole session."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.app import config, db  # noqa: E402

N_FAMILIES = 200
SEED = 4242


@pytest.fixture(scope="session", autouse=True)
def store(tmp_path_factory):
    """A deterministic store in a temp location — never the developer's own."""
    from backend.app.ingest import synthetic
    from backend.app.services import library

    config.DB_PATH = tmp_path_factory.mktemp("store") / "test.duckdb"
    config.EXPORT_DIR = tmp_path_factory.mktemp("exports")
    db.close()
    synthetic.rebuild(n_families=N_FAMILIES, seed=SEED)
    library.seed_builtins()
    yield
    db.close()


@pytest.fixture
def default_cohort():
    from backend.app.services import cohort
    return cohort.resolve()
