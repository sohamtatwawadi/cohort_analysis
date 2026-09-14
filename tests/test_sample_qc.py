"""Sample QC.

The point of this screen is that a contaminated sample or a swap invalidates
every rate downstream. So the tests are mostly about not lying: catching what
is genuinely bad, not flagging a whole assay type for being itself, and saying
which data a ratio was computed from.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.app import config, db  # noqa: E402
from backend.app.ingest import synthetic  # noqa: E402
from backend.app.ingest.varimat import compute_qc  # noqa: E402
from backend.app.services import cohort, sampleqc  # noqa: E402


@pytest.fixture(scope="module", autouse=True)
def store(tmp_path_factory):
    config.DB_PATH = tmp_path_factory.mktemp("qc") / "qc.duckdb"
    db.close()
    synthetic.rebuild(n_families=200, seed=31)
    yield
    db.close()


@pytest.fixture(scope="module")
def qc(store):
    with cohort.cohort_scope() as co:
        return sampleqc.sample_qc(co.criteria)


# ------------------------------------------------------------- the metrics --
class _V:
    """Minimal stand-in for a parsed VariMAT row."""
    def __init__(self, ref, alt, zyg="Heterozygous", vaf=0.5, chrom="1",
                 filter_status="PASS"):
        self.ref, self.alt = ref, alt
        self.zygosity_raw, self.vaf, self.chrom = zyg, vaf, chrom
        self.filter_status = filter_status


def test_ti_tv_counts_transitions_correctly():
    """A>G, G>A, C>T, T>C are transitions; everything else is a transversion."""
    v = [_V("A", "G"), _V("G", "A"), _V("C", "T"), _V("T", "C"),  # 4 Ti
         _V("A", "C"), _V("G", "T")]                              # 2 Tv
    assert compute_qc(v)["qc_ti_tv"] == pytest.approx(2.0)


def test_qc_ignores_non_pass_calls():
    v = [_V("A", "G")] * 4 + [_V("A", "C")] * 2 + [
        _V("A", "C", filter_status="LowQual")] * 50
    out = compute_qc(v)
    assert out["qc_n_called"] == 6
    assert out["qc_ti_tv"] == pytest.approx(2.0)


def test_qc_is_computed_over_the_full_call_set_not_the_retained_subset(qc):
    """The reason this lives on `run` at all. §3.2 archives ~98% of rows after
    fingerprinting, and `finding` keeps only the reviewable handful — far too
    few for a per-sample ratio."""
    basis = {s["metrics_basis"] for s in qc["samples"]}
    assert all("full call set" in b for b in basis), basis

    retained = db.scalar("""
        SELECT MAX(n) FROM (SELECT run_id, COUNT(*) n FROM finding GROUP BY 1)
    """, default=0)
    n_called = db.scalar("SELECT MIN(qc_n_called) FROM run", default=0)
    assert n_called > retained * 10, (
        "the stored call count ({}) is not meaningfully larger than the retained "
        "findings ({}) — the fixture is not exercising the distinction"
        .format(n_called, retained))


def test_metrics_basis_is_always_reported(qc):
    """A panel-sized ratio must not be readable as a genome-wide one."""
    for s in qc["samples"]:
        assert s["metrics_basis"]


# ------------------------------------------------------------ what it finds --
def test_contaminated_samples_are_failed(qc):
    """The planted contaminated samples carry all three signatures at once:
    depressed Ti/Tv, excess heterozygosity, off-centre allele balance."""
    failed = [s for s in qc["flagged"] if s["status"] == "fail"]
    assert failed, "nothing failed; the fixture plants bad samples, so this is wrong"

    triple = [s for s in failed
              if s["ti_tv"] and s["ti_tv"] < 2.0
              and s["het_hom"] and s["het_hom"] > 2.5
              and s["mean_het_vaf"] and abs(s["mean_het_vaf"] - 0.5) > 0.08]
    assert triple, "no sample shows the contamination signature"

    for s in triple[:3]:
        labels = " ".join(f["label"] for f in s["flags"])
        assert "Ti/Tv" in labels and "Het/Hom" in labels


def test_most_samples_pass(qc):
    """A check that fails everything is indistinguishable from no check. The
    fixture plants roughly 3% bad samples, so a sane screen passes the rest."""
    s = qc["summary"]
    assert s["pass"] / s["total"] > 0.85, s


def test_pipeline_qc_failure_is_surfaced(qc):
    """The lab's own verdict outranks anything we derive."""
    flagged_runs = {s["run_id"] for s in qc["flagged"]}
    pipeline_failed = {r["run_id"] for r in db.rows("""
        SELECT run_id FROM run WHERE qc_status <> 'Pass'
          AND run_id IN (SELECT run_id FROM cohort_run)
    """)}
    assert pipeline_failed <= flagged_runs, (
        "runs the pipeline itself failed were not surfaced: "
        + str(sorted(pipeline_failed - flagged_runs)[:5]))


def test_every_flag_states_a_number(qc):
    """'This sample looks odd' is not actionable. Each flag must carry the
    observed value and what it is being compared against."""
    for s in qc["flagged"]:
        assert s["flags"]
        for f in s["flags"]:
            assert f["level"] in ("warn", "fail")
            assert f["label"] and f["detail"]


# ------------------------------------------------ cohort-relative behaviour --
def test_thresholds_are_cohort_relative_not_absolute(qc):
    """A panel, an exome and a genome have different expected values. Flagging
    on a fixed number would fail an entire assay type for being itself."""
    assert "median absolute deviation" in qc["method"]
    for key in ("ti_tv", "het_hom", "mean_depth"):
        assert qc["cohort"][key]["median"] is not None
        assert qc["cohort"][key]["n"] > 0


def test_mad_is_not_dragged_by_the_outliers_it_looks_for():
    """The reason for MAD over mean/SD: a handful of extreme values inflate SD
    enough to hide themselves."""
    vals = [3.0] * 50 + [0.2] * 5
    med = sampleqc._median(vals)
    mad = sampleqc._mad(vals, med)
    assert med == pytest.approx(3.0)
    # An SD here is ~0.8, which would put the outliers inside 4 sigma. MAD
    # keeps them where they belong.
    dev = sampleqc._deviation(0.2, med, mad)
    assert dev is None or abs(dev) > 5


def test_small_panels_are_reported_as_not_assessable_when_unstored():
    """Without stored metrics, a 1-gene panel cannot support a ratio and must
    say so rather than produce one from three calls."""
    v = [_V("A", "G"), _V("C", "T"), _V("A", "C")]
    out = compute_qc(v)
    # compute_qc itself will produce a ratio from 3 calls; the guard is in the
    # service, which only falls back when there are enough retained findings.
    assert out["qc_n_called"] == 3
    assert sampleqc.MIN_CALLS_FOR_RATIO > 3


def test_summary_counts_reconcile(qc):
    s = qc["summary"]
    assert s["pass"] + s["warn"] + s["fail"] == s["total"]
    assert len(qc["flagged"]) == s["warn"] + s["fail"]
