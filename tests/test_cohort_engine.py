"""Regression tests for four cohort-engine defects found in review.

Each of these is invisible on synthetic data or under sequential access, which
is exactly why they need a test that constructs the condition deliberately:

  1. the lock released between resolve and read  -> wrong cohort, silently
  2. funnel step 1 unjoined                      -> exclusions misattributed
  3. cohort_run ignored the assay filter         -> denominators from other tests
  4. collection anchor moved on ingest           -> false reproducibility drift
"""
from __future__ import annotations

import sys
import threading
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.app import config, db  # noqa: E402
from backend.app.ingest import synthetic  # noqa: E402
from backend.app.services import carrier, cohort, denominator  # noqa: E402


@pytest.fixture(scope="module", autouse=True)
def store(tmp_path_factory):
    config.DB_PATH = tmp_path_factory.mktemp("engine") / "e.duckdb"
    db.close()
    synthetic.rebuild(n_families=60, seed=99)
    yield
    db.close()


# ============================================================ 1. lock contract
def test_cohort_scope_is_the_supported_entry_point():
    """resolve() alone locks only its own DROP/CREATE. Anything that then reads
    the temp tables must hold the lock across both, which is what the context
    manager exists to make unmissable."""
    with cohort.cohort_scope({"gene_set": "Lynch / MMR"}) as co:
        # Inside the scope the tables belong to us.
        n = db.scalar("SELECT COUNT(*) FROM cohort_subject", default=0)
        assert n == co.n_subjects


def test_concurrent_different_cohorts_do_not_cross_contaminate():
    """The failure this prevents is silent: a request reads another request's
    cohort and returns entirely plausible numbers for the wrong people."""
    lynch = {"gene_set": "Lynch / MMR"}
    hboc = {"gene_set": "HBOC (breast/ovarian)"}

    with cohort.cohort_scope(lynch) as co:
        lynch_n = co.n_subjects
    with cohort.cohort_scope(hboc) as co:
        hboc_n = co.n_subjects
    assert lynch_n != hboc_n, "cohorts are identical; the test proves nothing"

    seen = {"lynch": set(), "hboc": set()}
    errors = []

    def hammer(label, criteria, expected):
        for _ in range(12):
            try:
                with cohort.cohort_scope(criteria) as co:
                    # Read the temp table AFTER resolving, the window the bug
                    # lived in.
                    observed = db.scalar("SELECT COUNT(*) FROM cohort_subject",
                                         default=0)
                    seen[label].add(observed)
                    assert observed == expected
            except Exception as exc:            # noqa: BLE001 — surfaced below
                errors.append("{}: {}".format(label, exc))

    threads = [threading.Thread(target=hammer, args=("lynch", lynch, lynch_n)),
               threading.Thread(target=hammer, args=("hboc", hboc, hboc_n))]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert not errors, errors
    assert seen["lynch"] == {lynch_n}, "lynch saw {}".format(seen["lynch"])
    assert seen["hboc"] == {hboc_n}, "hboc saw {}".format(seen["hboc"])


# =========================================== 2. funnel attributes every drop
def test_subject_without_sequencing_gets_its_own_funnel_step():
    """A subject registered in the clinical system but not yet sequenced used to
    vanish at whatever step followed, and be blamed on it — usually the consent
    gate, which had never touched them."""
    db.execute("""
        INSERT INTO subject (subject_id, tenant_id, sex, consent_class, family_id,
                             is_proband, indication)
        VALUES ('SJ-NORUN', 11, 'F', 'Full research + secondary findings',
                'FAM-0001', TRUE, 'Cardiac')
    """)
    try:
        before = cohort.resolve()
        steps = {s["step"]: s for s in before.funnel}
        assert "Has sequencing data" in steps

        # The un-sequenced subject is attributed to the sequencing step.
        seq_drop = steps["Has sequencing data"]["dropped"]
        assert seq_drop >= 1

        # And the consent gate drops exactly the withdrawn subjects that HAVE
        # sequencing data — no more. Asserting a flat zero here would pass or
        # fail on whether the seeded store happened to contain any, which says
        # nothing about attribution.
        withdrawn_with_run = int(db.scalar("""
            SELECT COUNT(DISTINCT s.subject_id) FROM subject s
            JOIN run r USING (subject_id)
            WHERE s.consent_class = 'Withdrawn'
        """, default=0))
        consent = next(s for s in before.funnel
                       if s["step"].startswith("Consent gate"))
        assert consent["dropped"] == withdrawn_with_run, (
            "consent gate dropped {} but only {} withdrawn subjects have a run"
            .format(consent["dropped"], withdrawn_with_run))
    finally:
        db.execute("DELETE FROM subject WHERE subject_id = 'SJ-NORUN'")


def test_funnel_drops_reconcile_with_the_final_count():
    """The funnel's promise: every excluded subject attributable to exactly one
    step. Baseline minus the sum of drops must equal the cohort."""
    for crit in ({}, {"probands_only": False}, {"indication": ["Cardiac"]},
                 {"gene_set": "Lynch / MMR"}):
        co = cohort.resolve(crit)
        steps = [s for s in co.funnel if not s.get("final")]
        baseline = steps[0]["count"]
        dropped = sum(s["dropped"] for s in steps)
        assert baseline - dropped == co.n_subjects, (
            "{}: {} - {} != {}".format(crit, baseline, dropped, co.n_subjects))


# ================================ 3. cohort_run honours the assay criteria
def test_denominator_ignores_runs_outside_the_assay_filter():
    """Filtering to a one-gene panel must not credit an exome the same subject
    happened to have. Otherwise a gene the chosen test never looks at reads as
    'tested and clear'."""
    row = db.row("SELECT subject_id FROM run WHERE test_code = 'MG-THAL-1' LIMIT 1")
    assert row, "fixture has no MG-THAL-1 run"
    sid = row["subject_id"]

    with cohort.cohort_scope({"test_code": ["MG-THAL-1"]}):
        before = denominator.gene_denominators()["BRCA1"]["assayed"]

    db.execute("INSERT INTO sample (sample_id, subject_id, sample_type) "
               "VALUES ('SP-EXTRA', ?, 'Blood (EDTA)')", [sid])
    db.execute("INSERT INTO run (run_id, sample_id, subject_id, test_code, "
               "reference_build) VALUES ('GR-EXTRA','SP-EXTRA',?, 'MG-WES','GRCh38')",
               [sid])
    db.execute("INSERT INTO run_scope SELECT 'GR-EXTRA', gene_symbol, "
               "'SNV,Indel,Splice', TRUE, 'declared' FROM gene_disease")
    try:
        with cohort.cohort_scope({"test_code": ["MG-THAL-1"]}):
            after = denominator.gene_denominators()["BRCA1"]["assayed"]
            runs = [r["test_code"] for r in db.rows("""
                SELECT r.test_code FROM cohort_run cr JOIN run r USING (run_id)
                WHERE cr.subject_id = ?""", [sid])]
        assert after == before, (
            "an out-of-filter exome leaked into the denominator: {} -> {}"
            .format(before, after))
        assert runs == ["MG-THAL-1"], runs
    finally:
        db.execute("DELETE FROM run_scope WHERE run_id = 'GR-EXTRA'")
        db.execute("DELETE FROM run WHERE run_id = 'GR-EXTRA'")
        db.execute("DELETE FROM sample WHERE sample_id = 'SP-EXTRA'")


# ================================== 4. the collection anchor can be pinned
def test_collection_window_is_reproducible_when_the_anchor_is_pinned():
    """A window cohort legitimately changes membership when new samples arrive.
    Reproducibility means re-running against the SAME anchor gives the same
    members — not that the window never moves."""
    sid = db.row("SELECT subject_id FROM subject LIMIT 1")["subject_id"]
    before = cohort.resolve({"months": 24})
    assert before.as_of, "resolve did not record the anchor it used"

    db.execute("INSERT INTO sample (sample_id, subject_id, collection_date) "
               "VALUES ('SP-FUTURE', ?, DATE '2027-06-01')", [sid])
    try:
        moved = cohort.resolve({"months": 24})
        pinned = cohort.resolve({"months": 24}, as_of=before.as_of)

        assert moved.as_of != before.as_of, "anchor did not move; test is inert"
        assert pinned.as_of == before.as_of
        assert pinned.member_hash == before.member_hash, (
            "pinning the anchor did not reproduce the cohort")
    finally:
        db.execute("DELETE FROM sample WHERE sample_id = 'SP-FUTURE'")


def test_manifest_carries_the_anchor():
    from backend.app.services import governance
    with cohort.cohort_scope({"months": 12}) as co:
        mf = governance.manifest(co, module="g_carrier")
    assert mf["collection_anchor"] == co.as_of


# ======================================================= denominator labels
def test_confidence_reports_mixed_rather_than_collapsing_to_inferred():
    """99 declared plus 1 inferred is a different claim from 100 inferred, and
    the coverage inspector is where that distinction is meant to be visible."""
    with cohort.cohort_scope():
        dens = denominator.gene_denominators()
    for gene, d in dens.items():
        if d["declared"] and d["inferred"]:
            assert d["confidence"] == "mixed", (gene, d["confidence"])
        elif d["declared"]:
            assert d["confidence"] == "declared"


def test_not_assayed_absorbs_non_reportable_scope():
    """A gene the test does not report is a gene the cohort has no result for —
    it belongs in not_assayed, not in provisional."""
    with cohort.cohort_scope() as co:
        for gene, d in denominator.gene_denominators().items():
            total = d["assayed"] + d["provisional"] + d["not_assayed"]
            assert total == co.n_subjects, (gene, d)
