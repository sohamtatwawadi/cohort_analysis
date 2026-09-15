"""The score test must agree with the exact fit, or speed bought nothing.

A fast scan that ranks variants differently from the exact model is not an
optimisation, it is a different analysis wearing the same label. These tests
pin it against the GLM the scan replaced.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.app.research.stats import glm, scoretest  # noqa: E402
from backend.app.research.types import MISSING  # noqa: E402


def _sim(n=1200, seed=11, beta=0.0, n_cov=2):
    rng = np.random.default_rng(seed)
    g = rng.binomial(2, 0.25, n).astype(float)
    cov = rng.normal(0, 1, (n, n_cov))
    logit = -0.7 + beta * g + 0.3 * cov[:, 0]
    p = 1 / (1 + np.exp(-logit))
    y = (rng.random(n) < p).astype(float)
    return g, y, cov


# ------------------------------------------------------------ agreement -----
@pytest.mark.parametrize("beta", [0.0, 0.25, 0.5])
def test_score_p_tracks_the_exact_logistic_fit(beta):
    """Under the null and at modest effects the two should be close. They are
    asymptotically equivalent, not identical, so this checks the order of
    magnitude rather than equality."""
    g, y, cov = _sim(beta=beta, seed=int(100 * beta) + 3)
    out = scoretest.score_scan(g.reshape(1, -1), y, cov, kind="binary")
    design = np.column_stack([g, cov])
    exact = glm.fit_logistic(design, y, 0)
    p_exact = exact.pvalue
    p_score = out["p"][0]

    assert np.isfinite(p_score)
    # Both on the same side of significance, and within an order of magnitude.
    if p_exact < 1e-3 or p_score < 1e-3:
        ratio = np.log10(p_score) / np.log10(p_exact)
        assert 0.5 < ratio < 2.0, (p_score, p_exact)
    else:
        assert abs(p_score - p_exact) < 0.1, (p_score, p_exact)


def test_ranking_matches_the_exact_fit():
    """What a scan is actually for: the top of the list must be the same."""
    rng = np.random.default_rng(5)
    n, m = 900, 60
    cov = rng.normal(0, 1, (n, 2))
    G = rng.binomial(2, rng.uniform(0.1, 0.4, m)[:, None] * np.ones((m, n)))
    logit = -0.6 + 0.55 * G[3] + 0.5 * G[17] + 0.25 * cov[:, 0]
    y = (rng.random(n) < 1 / (1 + np.exp(-logit))).astype(float)

    fast = scoretest.score_scan(G, y, cov, kind="binary")["p"]
    slow = np.array([
        glm.fit_logistic(np.column_stack([G[i], cov]), y, 0).pvalue
        for i in range(m)])

    assert set(np.argsort(fast)[:5]) & {3, 17} == {3, 17}, "planted hits not on top"
    # Spearman correlation across the whole scan.
    rho = stats.spearmanr(fast, slow).correlation
    assert rho > 0.95, "score and exact rank variants differently (rho={:.3f})".format(rho)


def test_quantitative_outcome_matches_ols():
    rng = np.random.default_rng(9)
    n = 800
    g = rng.binomial(2, 0.3, n).astype(float)
    cov = rng.normal(0, 1, (n, 1))
    y = 1.0 + 0.4 * g + 0.5 * cov[:, 0] + rng.normal(0, 1, n)

    out = scoretest.score_scan(g.reshape(1, -1), y, cov, kind="quantitative")
    exact = glm.fit_linear(np.column_stack([g, cov]), y, 0)
    assert out["p"][0] < 1e-6 and exact.pvalue < 1e-6


# ------------------------------------------------------------ calibration ----
def test_null_p_values_are_uniform():
    """The property that makes a scan trustworthy: with no signal, p is flat.
    A scan that is even mildly anti-conservative produces a page of false hits."""
    rng = np.random.default_rng(23)
    n, m = 1500, 800
    cov = rng.normal(0, 1, (n, 2))
    y = (rng.random(n) < 0.35).astype(float)
    G = rng.binomial(2, 0.3, (m, n))

    p = scoretest.score_scan(G, y, cov, kind="binary")["p"]
    p = p[np.isfinite(p)]
    assert p.size > m * 0.9
    ks = stats.kstest(p, "uniform").pvalue
    assert ks > 0.01, "null p-values are not uniform (KS p={:.4g})".format(ks)

    lam = scoretest.genomic_inflation(p)
    assert 0.85 < lam < 1.15, "genomic inflation {:.3f} under a pure null".format(lam)


# --------------------------------------------------------------- robustness --
def test_monomorphic_variant_gets_no_p_value():
    """Dividing by a ~zero variance would manufacture a spectacular p-value for
    a variant nobody carries."""
    _, y, cov = _sim()
    G = np.zeros((1, y.size))
    out = scoretest.score_scan(G, y, cov, kind="binary")
    assert not out["usable"][0]
    assert np.isnan(out["p"][0])


def test_missing_dosages_are_imputed_not_dropped():
    g, y, cov = _sim(beta=0.4, seed=77)
    G = g.reshape(1, -1).astype(np.int8).copy()
    G[0, :50] = MISSING
    out = scoretest.score_scan(G, y, cov, kind="binary")
    assert np.isfinite(out["p"][0])
    assert out["n_used"] == y.size, "samples were dropped; n must stay constant"


def test_transposed_input_is_rejected_not_silently_analysed():
    """The defect this mirrors: a square-ish matrix passed the wrong way round
    is silently analysed as though it were correct."""
    g, y, cov = _sim(n=300)
    with pytest.raises(ValueError, match="transposed"):
        scoretest.score_scan(g.reshape(-1, 1), y, cov, kind="binary")


# ---------------------------------------------------------------- the point --
def test_it_is_actually_fast():
    """The entire justification. If this is not dramatically faster than the
    per-variant fit, the added code is not worth carrying."""
    rng = np.random.default_rng(3)
    n, m = 2000, 2000
    cov = rng.normal(0, 1, (n, 4))
    y = (rng.random(n) < 0.4).astype(float)
    G = rng.binomial(2, 0.25, (m, n))

    t0 = time.time()
    scoretest.score_scan(G, y, cov, kind="binary")
    fast = time.time() - t0

    t0 = time.time()
    for i in range(40):
        glm.fit_logistic(np.column_stack([G[i], cov]), y, 0)
    slow_per_variant = (time.time() - t0) / 40

    assert fast < slow_per_variant * m * 0.2, (
        "score scan of {} variants took {:.2f}s; {} exact fits would take about "
        "{:.0f}s".format(m, fast, m, slow_per_variant * m))
