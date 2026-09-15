"""Every analysis the UI offers must actually exist.

The capability matrix and the job registry are written in different files by
different concerns — one decides what the data can support, the other decides
what the engine can do — and nothing connected them. So five analyses were
advertised as "ready to run" for months while clicking them did nothing: no
modal, no error, no log line. The card was real, the analysis was not.

These tests close that gap from both ends, and add the frontend as a third end,
since a registered analysis with no configuration form fails in exactly the same
silent way.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.app.research import analyses, capability, jobs  # noqa: E402,F401
from backend.app.research.profile import DataProfile  # noqa: E402

ANALYSE_JS = ROOT / "frontend" / "js" / "analyse.js"

# Not built, and documented as not built (§12): R3 PheWAS, R7 fine-mapping,
# colocalization and heritability. They must stay permanently unavailable —
# a card that can never turn green is honest; one that turns green and then
# does nothing is not.
NOT_BUILT = {"phewas", "colocalization", "fine_mapping", "heritability"}


def _all_capabilities():
    """A maximally-capable profile, so every gate that CAN open does."""
    p = DataProfile(
        n_samples=5000, n_variants=500000, density_class="genome_wide",
        genome_build="GRCh38",
    )
    p.phenotypes = {
        "status": {"kind": "binary", "cases": 2000, "controls": 3000},
        "age": {"kind": "quantitative"},
        "ancestry": {"kind": "categorical"},
        "sex": {"kind": "binary"},
    }
    p.ancestry = {"n_pcs": 10}
    p.family_structure_present = True
    p.n_trios = 40
    p.n_binary_phenotypes = 2
    p.n_quantitative_phenotypes = 1
    p.max_cases = 2000
    p.has_controls = True
    p.has_time_to_event = True
    p.n_coded_phenotypes = 4
    p.n_unrelated = 4800
    p.mean_call_rate = 0.99
    p.n_polymorphic = 480000
    p.coverage_confidence = "declared"
    p.annotations_present = True
    return capability.assess(p)


def test_every_advertisable_analysis_is_implemented():
    """If the matrix can ever mark it available, the engine must be able to run
    it. This is the exact defect: available=True with no registered function."""
    registered = set(jobs.registered())
    missing = sorted(
        c.analysis for c in _all_capabilities()
        if c.available and c.analysis not in registered)
    assert not missing, (
        "advertised as ready to run but not implemented: {}. Either register an "
        "implementation or stop the matrix marking it available.".format(missing))


def test_analyses_that_are_not_built_can_never_become_available():
    """The four documented gaps must stay locked even on a perfect dataset."""
    for c in _all_capabilities():
        if c.analysis in NOT_BUILT:
            assert not c.available, (
                "{} is not implemented but the matrix marked it available"
                .format(c.analysis))


def test_nothing_is_registered_that_the_matrix_never_offers():
    """The reverse leak: an implementation no card can reach is dead code."""
    offered = {c.analysis for c in _all_capabilities()}
    orphans = sorted(set(jobs.registered()) - offered)
    assert not orphans, "registered but unreachable from the UI: {}".format(orphans)


# ------------------------------------------------------------- the frontend --
def _js_keys(block_name):
    """Top-level keys of a `const <name> = { ... }` object literal."""
    src = ANALYSE_JS.read_text()
    start = src.index("const {} = {{".format(block_name))
    depth, i = 0, src.index("{", start)
    for j in range(i, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                body = src[i + 1:j]
                break
    # Only keys at nesting depth 0 of the literal.
    keys, d = [], 0
    for line in body.splitlines():
        stripped = line.strip()
        if d == 0:
            m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*:", stripped)
            if m:
                keys.append(m.group(1))
        d += line.count("{") + line.count("[") - line.count("}") - line.count("]")
    return set(keys)


@pytest.fixture(scope="module")
def registered():
    return set(jobs.registered())


def test_every_runnable_analysis_has_a_configuration_form(registered):
    """`configureModal` bails with a toast when FORMS has no entry, which is
    indistinguishable from a dead button."""
    missing = sorted(registered - _js_keys("FORMS"))
    assert not missing, (
        "no FORMS entry in analyse.js — clicking these does nothing: {}"
        .format(missing))


def test_every_runnable_analysis_has_a_result_renderer(registered):
    """Without a RESULTS entry the user gets a raw JSON dump, which was
    explicitly rejected as unusable."""
    missing = sorted(registered - _js_keys("RESULTS"))
    assert not missing, (
        "no RESULTS renderer in analyse.js — results render as JSON: {}"
        .format(missing))


# ------------------------------------------------- interrupted job recovery --
def test_a_job_left_running_by_a_crash_is_failed_at_startup(tmp_path_factory):
    """Job state lives in the database, the thread advancing it does not. Kill
    the process mid-job and the row stays "running" with nothing left to move
    it — the UI then polls a job that will never finish, and a hung job is
    indistinguishable from a slow one."""
    from datetime import datetime

    from backend.app import config, db
    from backend.app.ingest import synthetic
    from backend.app.research import jobs, registry

    config.DB_PATH = tmp_path_factory.mktemp("stale") / "s.duckdb"
    db.close()
    synthetic.rebuild(n_families=10, seed=3)
    try:
        pid = registry.ensure_default_project()
        now = datetime.utcnow()
        for i, status in enumerate(jobs.ACTIVE_STATUSES + (jobs.STATUS_COMPLETE,)):
            db.insert_rows(
                "analysis_job",
                ["job_id", "project_id", "dataset_id", "analysis", "spec_json",
                 "status", "submitted_at"],
                [["job-stale-{}".format(i), pid, "ds-none", "gwas", "{}",
                  status, now]])

        n = jobs.reconcile_interrupted()
        assert n == len(jobs.ACTIVE_STATUSES), (
            "reconciled {} but {} were active".format(n, len(jobs.ACTIVE_STATUSES)))

        rows = {r["job_id"]: r for r in db.rows(
            "SELECT job_id, status, error FROM analysis_job "
            "WHERE job_id LIKE 'job-stale-%'")}
        for i in range(len(jobs.ACTIVE_STATUSES)):
            r = rows["job-stale-{}".format(i)]
            assert r["status"] == jobs.STATUS_FAILED
            assert "Interrupted" in (r["error"] or ""), r["error"]

        # A finished job is not touched — re-failing a completed analysis would
        # discard a result that is sitting right there.
        done = rows["job-stale-{}".format(len(jobs.ACTIVE_STATUSES))]
        assert done["status"] == jobs.STATUS_COMPLETE

        assert jobs.reconcile_interrupted() == 0, "second run should be a no-op"
    finally:
        db.close()
