"""Known-answer tests for the Research Mode regression kernels.

Every assertion here is against something computed independently of
`glm.py` — a closed-form result, a scipy routine, or a numerical optimiser run
on the objective function directly. A test that re-derives the answer using the
same code path it is testing proves only that the code is self-consistent, and
these numbers end up in publications.

The independent references used:

* linear      -> `scipy.stats.linregress` and an exact algebraic case
* logistic    -> direct maximisation of the Bernoulli log-likelihood with
                 `scipy.optimize.minimize`
* Firth       -> the closed-form "+½ to every cell" result, which is exactly
                 what Firth's penalty reduces to for a saturated 2x2 model
* HC3         -> the sandwich recomputed elementwise in the test
* BH          -> the hand-worked 10-p-value example that matches R's
                 `p.adjust(method="BH")`
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy import optimize, stats

from backend.app.research.stats import glm
from backend.app.research.stats.glm import (
    FitResult,
    fit_firth,
    fit_linear,
    fit_logistic,
    fit_robust_linear,
    multiple_testing,
    rank_inverse_normal,
    recommend_model,
)


# ===========================================================================
# 1. Linear regression
# ===========================================================================

def test_linear_exact_deterministic_line():
    """y = 2x + 1 with no noise. beta must be exactly 2, R^2 exactly 1."""
    x = np.arange(1.0, 21.0)
    y = 2.0 * x + 1.0
    res = fit_linear(x.reshape(-1, 1), y, 0, ["x"])

    assert res.beta == pytest.approx(2.0, abs=1e-12)
    assert res.extra["coefficients"]["(Intercept)"] == pytest.approx(1.0, abs=1e-12)
    assert res.extra["r2"] == pytest.approx(1.0, abs=1e-12)
    assert res.extra["residual_se"] == pytest.approx(0.0, abs=1e-12)
    assert res.n == 20
    assert res.model == "linear"
    assert res.effect_label == "beta"
    assert res.effect == res.beta
    assert res.n_cases is None and res.n_controls is None
    assert res.converged is True


def test_linear_matches_scipy_linregress():
    """Noisy simple regression against scipy's independent implementation."""
    rng = np.random.default_rng(20240501)
    x = rng.normal(size=150)
    y = 0.7 * x - 0.3 + rng.normal(scale=0.8, size=150)

    ref = stats.linregress(x, y)
    res = fit_linear(x.reshape(-1, 1), y, 0, ["x"])

    assert res.beta == pytest.approx(ref.slope, rel=1e-12)
    assert res.se == pytest.approx(ref.stderr, rel=1e-12)
    assert res.pvalue == pytest.approx(ref.pvalue, rel=1e-10)
    assert res.extra["r2"] == pytest.approx(ref.rvalue ** 2, rel=1e-12)
    assert res.extra["coefficients"]["(Intercept)"] == pytest.approx(
        ref.intercept, rel=1e-12)

    # 95% CI must be the t-interval, not a normal one.
    tcrit = stats.t.ppf(0.975, 148)
    assert res.ci_low == pytest.approx(ref.slope - tcrit * ref.stderr, rel=1e-12)
    assert res.ci_high == pytest.approx(ref.slope + tcrit * ref.stderr, rel=1e-12)


def test_linear_multivariable_matches_normal_equations():
    """With a covariate, the tested term must equal the OLS solution."""
    rng = np.random.default_rng(7)
    n = 120
    X = np.column_stack([rng.normal(size=n), rng.normal(size=n)])
    y = 1.5 * X[:, 0] - 0.4 * X[:, 1] + 2.0 + rng.normal(scale=0.5, size=n)

    res = fit_linear(X, y, 0, ["exposure", "age"])
    D = np.column_stack([np.ones(n), X])
    beta_ref = np.linalg.solve(D.T @ D, D.T @ y)

    assert res.beta == pytest.approx(beta_ref[1], rel=1e-12)
    assert res.covariates == ["age"]
    assert res.term == "exposure"


def test_linear_drops_incomplete_rows_and_reports_n():
    """§4.1 requires missingness on every result; n is the complete-case count."""
    x = np.arange(30.0)
    y = 3.0 * x + 2.0
    x_missing = x.copy()
    x_missing[[2, 5]] = np.nan
    y_missing = y.copy()
    y_missing[7] = np.nan

    res = fit_linear(x_missing.reshape(-1, 1), y_missing, 0, ["x"])
    assert res.n == 27
    assert res.extra["n_dropped_missing"] == 3
    assert res.beta == pytest.approx(3.0, abs=1e-10)


# ===========================================================================
# 2. Logistic regression
# ===========================================================================

def _neg_loglik(beta: np.ndarray, D: np.ndarray, y: np.ndarray) -> float:
    """Bernoulli negative log-likelihood, written independently of glm.py."""
    eta = D @ beta
    return float(np.sum(np.logaddexp(0.0, eta) - y * eta))


def test_logistic_matches_direct_likelihood_maximisation():
    """IRLS betas must equal the argmax found by a general-purpose optimiser."""
    rng = np.random.default_rng(99)
    n = 400
    X = np.column_stack([rng.binomial(1, 0.3, n).astype(float), rng.normal(size=n)])
    eta = -0.5 + 1.1 * X[:, 0] + 0.6 * X[:, 1]
    y = rng.binomial(1, 1.0 / (1.0 + np.exp(-eta))).astype(float)

    res = fit_logistic(X, y, 0, ["carrier", "pc1"])

    D = np.column_stack([np.ones(n), X])
    ref = optimize.minimize(_neg_loglik, np.zeros(3), args=(D, y),
                            method="BFGS", options={"gtol": 1e-12, "maxiter": 2000})
    assert ref.success

    assert res.beta == pytest.approx(float(ref.x[1]), abs=1e-5)
    assert res.converged is True
    assert res.warnings == []
    assert res.effect_label == "OR"
    assert res.effect == pytest.approx(np.exp(res.beta), rel=1e-12)
    assert res.n_cases == int(y.sum())
    assert res.n_controls == int(n - y.sum())
    assert res.n_cases + res.n_controls == res.n

    # Wald SE from the inverse observed information, recomputed here.
    pi = 1.0 / (1.0 + np.exp(-(D @ ref.x)))
    cov = np.linalg.inv(D.T @ (D * (pi * (1 - pi))[:, None]))
    assert res.se == pytest.approx(float(np.sqrt(cov[1, 1])), rel=1e-5)


def test_logistic_saturated_2x2_reproduces_the_sample_log_odds_ratio():
    """Closed form: with one binary predictor the model is saturated, so the
    MLE is exactly the empirical log OR. No optimiser needed as a reference."""
    # exposed: 12 cases, 18 controls; unexposed: 20 cases, 50 controls.
    a, c, b, d = 12, 18, 20, 50
    x = np.concatenate([np.ones(a + c), np.zeros(b + d)])
    y = np.concatenate([np.ones(a), np.zeros(c), np.ones(b), np.zeros(d)])

    res = fit_logistic(x.reshape(-1, 1), y, 0, ["x"])
    log_or = float(np.log((a * d) / (b * c)))

    assert res.beta == pytest.approx(log_or, abs=1e-8)
    assert res.converged is True
    assert res.warnings == []
    # The Wald SE of a saturated 2x2 log OR is the Woolf standard error.
    woolf = float(np.sqrt(1 / a + 1 / b + 1 / c + 1 / d))
    assert res.se == pytest.approx(woolf, rel=1e-8)


def test_logistic_flags_complete_separation():
    """Perfect prediction: the MLE does not exist and the fit must say so."""
    x = np.array([0.0] * 20 + [1.0] * 20)
    y = np.array([0.0] * 20 + [1.0] * 20)

    res = fit_logistic(x.reshape(-1, 1), y, 0, ["carrier"])

    assert res.converged is False, "separated fit must not be reported as converged"
    assert res.extra["separation"] is True
    assert any("separation" in w for w in res.warnings)
    # The coefficient has run off toward infinity, which is the diagnostic.
    assert abs(res.beta) > glm.SEPARATION_BETA


def test_logistic_flags_quasi_complete_separation_zero_cell():
    """A zero cell in the exposure x outcome table is the case we care about."""
    rng = np.random.default_rng(3)
    n = 200
    y = rng.binomial(1, 0.3, n).astype(float)
    x = np.zeros(n)
    carriers = np.where(y == 1)[0][:6]
    x[carriers] = 1.0  # every carrier is a case -> zero cell

    res = fit_logistic(x.reshape(-1, 1), y, 0, ["carrier"])
    assert any("zero cell" in w or "separation" in w for w in res.warnings)
    assert res.converged is False


def test_logistic_reports_non_convergence_when_iterations_are_starved():
    rng = np.random.default_rng(11)
    n = 300
    X = rng.normal(size=(n, 1))
    y = rng.binomial(1, 1.0 / (1.0 + np.exp(-X[:, 0]))).astype(float)

    res = fit_logistic(X, y, 0, ["x"], max_iter=1, tol=1e-12)
    assert res.converged is False
    assert any("did not converge" in w for w in res.warnings)
    assert res.extra["irls_converged"] is False


# ===========================================================================
# 3. Firth penalised logistic
# ===========================================================================

def _firth_closed_form_2x2(a: float, b: float, c: float, d: float) -> float:
    """Firth log-OR for a saturated 2x2 model = the "+1/2 to every cell" estimate.

    For two independent binomial groups the Jeffreys prior is
    |I|^(1/2) ∝ [π1(1-π1)]^(1/2) [π2(1-π2)]^(1/2), so the penalised likelihood is
    the ordinary likelihood with half an observation added to each of the four
    cells. This makes the Firth estimate available in closed form and gives us a
    reference that does not go through glm.py at all.

    a = exposed cases, c = exposed controls, b = unexposed cases,
    d = unexposed controls.
    """
    return float(np.log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))))


@pytest.mark.parametrize("table", [
    (12.0, 40.0, 8.0, 60.0),    # ordinary table
    (3.0, 55.0, 1.0, 70.0),     # rare exposure, small counts
    (10.0, 0.0, 0.0, 10.0),     # complete separation
    (6.0, 90.0, 0.0, 120.0),    # quasi-complete separation (one zero cell)
])
def test_firth_saturated_2x2_matches_closed_form(table):
    """The strongest known-answer check available for Firth."""
    a, b, c, d = table  # a=exposed cases, b=unexposed cases, c=exposed ctrl, d=unexp ctrl
    x = np.concatenate([np.ones(int(a)), np.ones(int(c)),
                        np.zeros(int(b)), np.zeros(int(d))])
    y = np.concatenate([np.ones(int(a)), np.zeros(int(c)),
                        np.ones(int(b)), np.zeros(int(d))])

    res = fit_firth(x.reshape(-1, 1), y, 0, ["carrier"])
    expected = _firth_closed_form_2x2(a, b, c, d)

    assert res.beta == pytest.approx(expected, abs=1e-6)
    assert res.effect == pytest.approx(np.exp(expected), rel=1e-6)
    assert res.converged is True


def test_firth_is_finite_where_logistic_diverges():
    """The classic separation test. Standard logistic fails; Firth must not."""
    x = np.array([0.0] * 10 + [1.0] * 10)
    y = np.array([0.0] * 10 + [1.0] * 10)

    plain = fit_logistic(x.reshape(-1, 1), y, 0, ["carrier"])
    firth = fit_firth(x.reshape(-1, 1), y, 0, ["carrier"])

    # Standard logistic: no finite MLE.
    assert plain.converged is False
    assert any("separation" in w for w in plain.warnings)
    assert abs(plain.beta) > glm.SEPARATION_BETA

    # Firth: everything finite, and equal to the closed-form +1/2 answer,
    # log(10.5 * 10.5 / (0.5 * 0.5)) = log(441).
    assert np.isfinite(firth.beta)
    assert firth.beta == pytest.approx(np.log(441.0), abs=1e-6)
    assert np.isfinite(firth.pvalue) and 0.0 < firth.pvalue < 1.0
    assert firth.pvalue < 1e-3
    assert np.isfinite(firth.ci_low) and np.isfinite(firth.ci_high)
    assert firth.ci_low > 0.0, "CI must exclude the null under perfect prediction"
    assert firth.extra["test"] == "penalised likelihood ratio"

    # Firth must still *say* the data were separated. Its own fitted
    # probabilities give no hint (0.5/11 and 10.5/11 are unremarkable), so the
    # detector has to interrogate the unpenalised likelihood.
    assert firth.extra["separation"] is True
    assert any("separation" in w for w in firth.warnings)
    assert 0.01 < firth.extra["min_fitted_prob"] < 0.1


def test_firth_does_not_cry_separation_on_well_behaved_data():
    rng = np.random.default_rng(19)
    n = 400
    X = np.column_stack([rng.binomial(1, 0.3, n).astype(float), rng.normal(size=n)])
    eta = -0.5 + 0.8 * X[:, 0] + 0.3 * X[:, 1]
    y = rng.binomial(1, 1.0 / (1.0 + np.exp(-eta))).astype(float)

    res = fit_firth(X, y, 0, ["carrier", "pc1"])
    assert res.extra["separation"] is False
    assert res.warnings == []
    assert res.converged is True


def test_firth_lrt_beats_wald_under_separation():
    """Why the p-value is an LRT: the Wald test collapses near separation.

    This is the Hauck-Donner effect. With near-separated data the SE grows
    faster than the coefficient, so beta/SE shrinks and the Wald test loses the
    signal exactly where it is strongest. The penalised LRT does not.
    """
    # Quasi-complete separation of the realistic kind: 3 carriers in 1,000
    # subjects, all three of them cases. This is the germline situation the
    # spec's worked example describes.
    n_cases, n_controls, n_carriers = 300, 700, 3
    y = np.array([1.0] * n_cases + [0.0] * n_controls)
    x = np.zeros(n_cases + n_controls)
    x[:n_carriers] = 1.0

    res = fit_firth(x.reshape(-1, 1), y, 0, ["carrier"])
    p_wald = float(2.0 * stats.norm.sf(abs(res.beta / res.se)))

    # The two tests reach opposite conclusions at alpha = 0.05, and the LRT is
    # the one that is right: the carriers really are all cases.
    assert res.pvalue < 0.05 < p_wald, (
        "LRT p={:.4g}, Wald p={:.4g}".format(res.pvalue, p_wald))
    assert res.pvalue < p_wald / 5.0

    # Under full separation the gap is even larger.
    xs = np.array([0.0] * 10 + [1.0] * 10)
    ys = np.array([0.0] * 10 + [1.0] * 10)
    sep = fit_firth(xs.reshape(-1, 1), ys, 0, ["carrier"])
    assert sep.pvalue < 2.0 * stats.norm.sf(abs(sep.beta / sep.se))


def test_firth_shrinks_estimates_toward_zero_on_well_behaved_data():
    """Without separation Firth should agree with logistic, slightly shrunk."""
    rng = np.random.default_rng(2718)
    n = 150
    X = np.column_stack([rng.binomial(1, 0.25, n).astype(float), rng.normal(size=n)])
    eta = -1.0 + 1.4 * X[:, 0] + 0.5 * X[:, 1]
    y = rng.binomial(1, 1.0 / (1.0 + np.exp(-eta))).astype(float)

    plain = fit_logistic(X, y, 0, ["carrier", "pc1"])
    firth = fit_firth(X, y, 0, ["carrier", "pc1"])
    assert plain.converged is True

    assert abs(firth.beta) < abs(plain.beta), "Firth must shrink toward the null"
    assert firth.beta == pytest.approx(plain.beta, rel=0.25), (
        "shrinkage should be modest on a non-separated dataset")
    assert np.sign(firth.beta) == np.sign(plain.beta)
    assert firth.n == plain.n == n


def test_firth_with_covariates_matches_direct_penalised_maximisation():
    """The closed form only covers the saturated 2x2 case. With covariates the
    reference is a general-purpose optimiser run on the penalised objective,
    which is itself verified against l + 0.5*logdet by the next test."""
    rng = np.random.default_rng(1234)
    n = 250
    X = np.column_stack([
        rng.binomial(1, 0.08, n).astype(float),   # rare exposure
        rng.normal(size=n),                        # a PC
        rng.uniform(40, 80, n),                    # age
    ])
    eta = -4.0 + 1.3 * X[:, 0] + 0.4 * X[:, 1] + 0.03 * X[:, 2]
    y = rng.binomial(1, 1.0 / (1.0 + np.exp(-eta))).astype(float)

    res = fit_firth(X, y, 0, ["carrier", "pc1", "age"])

    D = np.column_stack([np.ones(n), X])

    def neg_penalised(beta):
        eta_b = D @ beta
        pi = 1.0 / (1.0 + np.exp(-eta_b))
        pi = np.clip(pi, 1e-12, 1 - 1e-12)
        ll = float(np.sum(y * np.log(pi) + (1 - y) * np.log(1 - pi)))
        info = D.T @ (D * (pi * (1 - pi))[:, None])
        return -(ll + 0.5 * float(np.linalg.slogdet(info)[1]))

    ref = optimize.minimize(neg_penalised, np.zeros(4), method="Nelder-Mead",
                            options={"xatol": 1e-10, "fatol": 1e-12,
                                     "maxiter": 200000, "maxfev": 200000})
    assert ref.success
    assert res.beta == pytest.approx(float(ref.x[1]), abs=1e-4)
    assert res.extra["loglik_penalised"] == pytest.approx(-float(ref.fun), abs=1e-6)
    assert res.covariates == ["pc1", "age"]


def test_firth_modified_score_is_zero_at_the_estimate():
    """U*(beta) = X'[y - pi + h(1/2 - pi)] must vanish at the Firth solution."""
    rng = np.random.default_rng(555)
    n = 180
    X = np.column_stack([rng.binomial(1, 0.1, n).astype(float), rng.normal(size=n)])
    y = rng.binomial(1, 0.25, n).astype(float)
    D = np.column_stack([np.ones(n), X])

    beta, _, converged, _ = glm._firth_newton(D, y, 100, 1e-10)
    assert converged

    pi = 1.0 / (1.0 + np.exp(-(D @ beta)))
    w = pi * (1 - pi)
    sw = np.sqrt(w)
    H = (sw[:, None] * D) @ np.linalg.inv(D.T @ (D * w[:, None])) @ (sw[:, None] * D).T
    h = np.diag(H)
    score = D.T @ (y - pi + h * (0.5 - pi))
    assert np.max(np.abs(score)) < 1e-6


def test_firth_penalised_loglik_equals_loglik_plus_half_logdet():
    """l*(beta) = l(beta) + 0.5 * log det(X'WX), verified independently."""
    rng = np.random.default_rng(5150)
    n = 80
    X = np.column_stack([rng.binomial(1, 0.3, n).astype(float), rng.normal(size=n)])
    y = rng.binomial(1, 0.4, n).astype(float)
    D = np.column_stack([np.ones(n), X])
    beta = np.array([-0.4, 0.9, 0.2])

    pi = 1.0 / (1.0 + np.exp(-(D @ beta)))
    plain_ll = float(np.sum(y * np.log(pi) + (1 - y) * np.log(1 - pi)))
    info = D.T @ (D * (pi * (1 - pi))[:, None])
    expected = plain_ll + 0.5 * float(np.log(np.linalg.det(info)))

    assert glm._firth_penalised_loglik(D, y, beta) == pytest.approx(expected, abs=1e-10)


def test_firth_score_uses_the_hat_matrix_diagonal():
    """h_i must be the diagonal of W^(1/2) X (X'WX)^-1 X' W^(1/2)."""
    rng = np.random.default_rng(31415)
    n = 40
    X = np.column_stack([rng.normal(size=n), rng.normal(size=n)])
    D = np.column_stack([np.ones(n), X])
    w = rng.uniform(0.05, 0.25, n)

    h, info = glm._hat_diagonal(D, w)
    sw = np.sqrt(w)
    H_full = (sw[:, None] * D) @ np.linalg.inv(D.T @ (D * w[:, None])) @ (sw[:, None] * D).T

    assert np.allclose(h, np.diag(H_full), atol=1e-12)
    assert np.allclose(info, D.T @ (D * w[:, None]), atol=1e-12)
    # Leverages are a projection diagonal: they must lie in [0, 1] and sum to p.
    assert h.min() >= -1e-12 and h.max() <= 1.0 + 1e-12
    assert h.sum() == pytest.approx(D.shape[1], abs=1e-10)


def test_firth_profile_ci_endpoints_invert_the_reported_test():
    """The CI endpoints must be where the penalised LRT statistic hits 3.841.

    This is what makes the interval and the p-value mutually consistent: if the
    CI excludes 0 then the LRT rejects at 0.05, and vice versa.
    """
    rng = np.random.default_rng(606)
    n = 120
    X = np.column_stack([rng.binomial(1, 0.15, n).astype(float), rng.normal(size=n)])
    eta = -1.2 + 1.0 * X[:, 0] + 0.4 * X[:, 1]
    y = rng.binomial(1, 1.0 / (1.0 + np.exp(-eta))).astype(float)

    res = fit_firth(X, y, 0, ["carrier", "pc1"])
    assert res.extra["ci_method"] == "profile penalised likelihood"

    D = np.column_stack([np.ones(n), X])
    beta_hat, ll_full, _, _ = glm._firth_newton(D, y, 100, 1e-8)

    for endpoint in (res.ci_low, res.ci_high):
        _, ll_c, _, _ = glm._firth_newton(
            D, y, 100, 1e-10, fixed_col=1, fixed_value=endpoint, start=beta_hat)
        assert 2.0 * (ll_full - ll_c) == pytest.approx(glm.CHI2_95_1DF, abs=1e-4)

    assert res.ci_low < res.beta < res.ci_high
    # Interval and test must agree about the null.
    assert (res.ci_low > 0 or res.ci_high < 0) == (res.pvalue < 0.05)


def test_firth_lrt_statistic_matches_the_two_penalised_fits():
    rng = np.random.default_rng(808)
    n = 100
    X = rng.binomial(1, 0.2, (n, 1)).astype(float)
    y = rng.binomial(1, 0.35, n).astype(float)

    res = fit_firth(X, y, 0, ["carrier"])
    lrt = 2.0 * (res.extra["loglik_penalised"] - res.extra["loglik_penalised_null"])

    assert res.extra["lrt_stat"] == pytest.approx(max(lrt, 0.0), abs=1e-8)
    assert res.pvalue == pytest.approx(stats.chi2.sf(res.extra["lrt_stat"], 1), rel=1e-12)


def test_firth_drops_missing_rows():
    rng = np.random.default_rng(123)
    n = 90
    X = rng.binomial(1, 0.3, (n, 1)).astype(float)
    y = rng.binomial(1, 0.4, n).astype(float)
    X[[1, 4, 9], 0] = np.nan

    res = fit_firth(X, y, 0, ["carrier"])
    assert res.n == n - 3
    assert res.n_cases + res.n_controls == res.n
    assert res.extra["n_dropped_missing"] == 3


# ===========================================================================
# 4. Robust linear (HC3)
# ===========================================================================

def _hc3_reference(D: np.ndarray, y: np.ndarray) -> np.ndarray:
    """HC3 sandwich built elementwise, independently of glm.py."""
    xtx_inv = np.linalg.inv(D.T @ D)
    beta = xtx_inv @ D.T @ y
    resid = y - D @ beta
    H = D @ xtx_inv @ D.T
    meat = np.zeros((D.shape[1], D.shape[1]))
    for i in range(D.shape[0]):
        wi = (resid[i] ** 2) / (1.0 - H[i, i]) ** 2
        meat += wi * np.outer(D[i], D[i])
    return xtx_inv @ meat @ xtx_inv


def test_hc3_matches_an_independent_sandwich_computation():
    rng = np.random.default_rng(4242)
    n = 90
    X = np.column_stack([rng.normal(size=n), rng.normal(size=n)])
    y = 0.8 * X[:, 0] + rng.normal(scale=0.5 + np.abs(X[:, 0]), size=n)

    res = fit_robust_linear(X, y, 0, ["exposure", "age"])
    D = np.column_stack([np.ones(n), X])
    cov_ref = _hc3_reference(D, y)

    assert res.se == pytest.approx(float(np.sqrt(cov_ref[1, 1])), rel=1e-10)
    assert res.extra["se_method"] == "HC3"


def test_hc3_se_differs_from_ols_under_deliberate_heteroskedasticity():
    """Residual variance grows with |x|; the classical SE is then wrong."""
    rng = np.random.default_rng(1006)
    n = 400
    x = rng.normal(size=n)
    y = 0.5 * x + rng.normal(scale=0.2 + 2.0 * np.abs(x), size=n)

    ols = fit_linear(x.reshape(-1, 1), y, 0, ["x"])
    rob = fit_robust_linear(x.reshape(-1, 1), y, 0, ["x"])

    # Point estimates are identical: only the variance estimator changed.
    assert rob.beta == pytest.approx(ols.beta, rel=1e-12)
    ratio = rob.se / ols.se
    assert ratio > 1.2, "HC3 should be materially larger here, got ratio {:.3f}".format(ratio)
    assert rob.pvalue > ols.pvalue, "the honest SE must widen the p-value"
    assert rob.extra["se_classical"] == pytest.approx(ols.se, rel=1e-10)


def test_hc3_se_close_to_ols_when_errors_are_homoskedastic():
    """Under the OLS assumptions HC3 should cost little."""
    rng = np.random.default_rng(1007)
    n = 800
    x = rng.normal(size=n)
    y = 0.5 * x + rng.normal(scale=1.0, size=n)

    ols = fit_linear(x.reshape(-1, 1), y, 0, ["x"])
    rob = fit_robust_linear(x.reshape(-1, 1), y, 0, ["x"])
    assert rob.se == pytest.approx(ols.se, rel=0.15)


# ===========================================================================
# 5. recommend_model (§4.1 automatic model selection)
# ===========================================================================

def test_recommend_model_quantitative_is_linear():
    rng = np.random.default_rng(2)
    y = rng.normal(size=500)
    rec = recommend_model(y, "quantitative", exposure=rng.normal(size=500))

    assert rec["model"] == "linear"
    assert "quantitative" in rec["rationale"]
    assert "linear" in " ".join(rec["alternatives"]) or rec["alternatives"]


def test_recommend_model_rare_exposure_few_cases_is_firth():
    """Mirrors the spec's own worked example: 43 cases, 0.4% carrier frequency."""
    rng = np.random.default_rng(43)
    n = 4000
    y = np.zeros(n)
    y[rng.choice(n, 43, replace=False)] = 1.0
    exposure = np.zeros(n)
    exposure[rng.choice(n, 16, replace=False)] = 1.0  # 0.4%

    rec = recommend_model(y, "binary", exposure=exposure)

    assert rec["model"] == "firth"
    assert "43 cases" in rec["rationale"]
    assert "0.4%" in rec["rationale"]
    assert "separate" in rec["rationale"]
    assert "Firth" in rec["rationale"]
    assert "logistic" in rec["alternatives"]


def test_recommend_model_adequate_binary_is_logistic():
    rng = np.random.default_rng(5)
    n = 5000
    y = rng.binomial(1, 0.3, n).astype(float)
    exposure = rng.binomial(1, 0.2, n).astype(float)

    rec = recommend_model(y, "binary", exposure=exposure)
    assert rec["model"] == "logistic"
    assert "firth" in rec["alternatives"]


def test_recommend_model_few_cases_alone_triggers_firth():
    y = np.array([1.0] * 20 + [0.0] * 2000)
    exposure = np.array([1.0] * 400 + [0.0] * 1620)
    rec = recommend_model(y, "binary", exposure=exposure)
    assert rec["model"] == "firth"
    assert "only 20 cases" in rec["rationale"]


def test_recommend_model_mentions_relatedness_and_pruning():
    rng = np.random.default_rng(6)
    y = rng.normal(size=300)
    rec = recommend_model(y, "quantitative", exposure=None, n_related=37)
    assert "37 related sample" in rec["rationale"]
    assert "mixed model" in rec["rationale"]
    assert "prun" in rec["rationale"]
    assert "mixed" in rec["alternatives"]


def test_recommend_model_rejects_unknown_kind():
    with pytest.raises(ValueError):
        recommend_model(np.zeros(10), "categorical")


# ===========================================================================
# 6. Rank-based inverse-normal transform
# ===========================================================================

def test_rank_inverse_normal_is_approximately_standard_normal():
    rng = np.random.default_rng(31)
    v = rng.lognormal(mean=2.0, sigma=1.5, size=4000)  # violently skewed input
    z = rank_inverse_normal(v)

    assert z.mean() == pytest.approx(0.0, abs=1e-6)   # symmetric by construction
    assert z.std(ddof=0) == pytest.approx(1.0, abs=0.02)
    assert abs(stats.skew(z)) < 0.02
    assert stats.kstest(z, "norm").pvalue > 0.05


def test_rank_inverse_normal_is_monotone_and_averages_ties():
    v = np.array([5.0, 1.0, 3.0, 3.0, 9.0])
    z = rank_inverse_normal(v)

    assert z[2] == pytest.approx(z[3]), "tied inputs must map to the same value"
    assert z[1] < z[2] < z[0] < z[4]
    # Blom scores for n=5 with ranks 1, 3.5, 3.5, 2... check rank 1 exactly.
    assert z[1] == pytest.approx(stats.norm.ppf((1 - 0.375) / (5 + 0.25)), abs=1e-12)


def test_rank_inverse_normal_passes_nan_through_and_ranks_only_observed():
    v = np.array([2.0, np.nan, 8.0, 4.0, np.nan])
    z = rank_inverse_normal(v)

    assert np.isnan(z[1]) and np.isnan(z[4])
    assert np.all(np.isfinite(z[[0, 2, 3]]))
    # n = 3 for the ranking, so the extremes are the n=3 Blom scores.
    assert z[0] == pytest.approx(stats.norm.ppf((1 - 0.375) / (3 + 0.25)), abs=1e-12)
    assert z[2] == pytest.approx(stats.norm.ppf((3 - 0.375) / (3 + 0.25)), abs=1e-12)
    assert rank_inverse_normal(np.array([np.nan, np.nan])).shape == (2,)


def test_rank_inverse_normal_then_linear_recovers_direction():
    """The transform destroys units but must preserve the association."""
    rng = np.random.default_rng(77)
    n = 500
    x = rng.binomial(1, 0.3, n).astype(float)
    y = np.exp(0.8 * x + rng.normal(size=n))  # lognormal outcome
    res = fit_linear(x.reshape(-1, 1), rank_inverse_normal(y), 0, ["carrier"])
    assert res.beta > 0
    assert res.pvalue < 1e-6


# ===========================================================================
# 7. Multiple testing (§5.6)
# ===========================================================================

def test_bonferroni_is_p_times_m_capped_at_one():
    p = np.array([0.001, 0.02, 0.3, 0.9])
    adj = multiple_testing(p, "bonferroni")
    assert np.allclose(adj, [0.004, 0.08, 1.0, 1.0])


def test_benjamini_hochberg_hand_worked_example():
    """Ten p-values worked through by hand; matches R p.adjust(method="BH").

    i   raw     p*m/i       running min from the largest p downward
    1   0.001   0.010000 -> 0.010000
    2   0.008   0.040000 -> 0.040000
    3   0.039   0.130000 -> 0.084000   (capped by row 5)
    4   0.041   0.102500 -> 0.084000   (capped by row 5)
    5   0.042   0.084000 -> 0.084000
    6   0.060   0.100000 -> 0.100000   (row 7 is larger, so no cap applies)
    7   0.074   0.105714 -> 0.105714
    8   0.205   0.256250 -> 0.216000   (capped by row 10)
    9   0.212   0.235556 -> 0.216000   (capped by row 10)
    10  0.216   0.216000 -> 0.216000
    """
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074,
                  0.205, 0.212, 0.216])
    expected = np.array([0.01, 0.04, 0.084, 0.084, 0.084, 0.10,
                         0.074 * 10 / 7, 0.216, 0.216, 0.216])
    adj = multiple_testing(p, "fdr_bh")
    assert np.allclose(adj, expected, atol=1e-12)


def test_benjamini_hochberg_preserves_input_order():
    p = np.array([0.216, 0.001, 0.042, 0.008])
    adj = multiple_testing(p, "fdr_bh")
    ref = multiple_testing(np.sort(p), "fdr_bh")
    assert adj[1] == pytest.approx(ref[0])
    assert adj[3] == pytest.approx(ref[1])
    assert adj[2] == pytest.approx(ref[2])
    assert adj[0] == pytest.approx(ref[3])


def test_benjamini_hochberg_is_monotone_and_never_below_raw_p():
    rng = np.random.default_rng(9)
    p = rng.uniform(size=200)
    adj = multiple_testing(p, "fdr_bh")
    assert np.all(adj >= p - 1e-12)
    assert np.all(adj <= 1.0)
    order = np.argsort(p)
    assert np.all(np.diff(adj[order]) >= -1e-12)
    # BH is never more conservative than Bonferroni.
    assert np.all(adj <= multiple_testing(p, "bonferroni") + 1e-12)


def test_multiple_testing_handles_nan_and_rejects_bad_input():
    p = np.array([0.01, np.nan, 0.04])
    adj = multiple_testing(p, "bonferroni")
    assert np.isnan(adj[1])
    assert np.allclose(adj[[0, 2]], [0.02, 0.08])  # m = 2, the NaN is excluded

    with pytest.raises(ValueError):
        multiple_testing([0.1, 1.5], "fdr_bh")
    with pytest.raises(ValueError):
        multiple_testing([0.1, 0.2], "holm")


# ===========================================================================
# Cross-cutting contract
# ===========================================================================

@pytest.mark.parametrize("fit,kind", [
    (fit_linear, "quantitative"),
    (fit_logistic, "binary"),
    (fit_firth, "binary"),
    (fit_robust_linear, "quantitative"),
])
def test_every_fit_populates_the_required_reporting_fields(fit, kind):
    """§4.1: n, cases, controls, effect, SE, CI, p, model, covariates,
    missingness and convergence on every result."""
    rng = np.random.default_rng(4321)
    n = 200
    X = np.column_stack([rng.binomial(1, 0.3, n).astype(float), rng.normal(size=n)])
    if kind == "binary":
        y = rng.binomial(1, 1.0 / (1.0 + np.exp(-(X[:, 0] + 0.3 * X[:, 1])))).astype(float)
    else:
        y = X[:, 0] + 0.3 * X[:, 1] + rng.normal(size=n)

    res = fit(X, y, 0, ["carrier", "pc1"])
    assert isinstance(res, FitResult)
    assert res.term == "carrier"
    assert res.covariates == ["pc1"]
    assert res.n == n
    assert 0.0 <= res.pvalue <= 1.0
    assert res.ci_low <= res.beta <= res.ci_high
    assert res.effect_ci_low <= res.effect <= res.effect_ci_high
    assert res.effect_label in ("beta", "OR", "HR")
    assert res.model in ("linear", "logistic", "firth", "cox")
    assert isinstance(res.warnings, list)
    assert "n_dropped_missing" in res.extra
    # §5.3 diagnostics: the full coefficient vector, covariates included.
    assert list(res.extra["coefficients"]) == ["(Intercept)", "carrier", "pc1"]
    assert res.extra["coefficients"]["carrier"] == pytest.approx(res.beta, rel=1e-12)
    if kind == "binary":
        assert res.n_cases is not None and res.n_controls is not None
        assert res.effect == pytest.approx(np.exp(res.beta), rel=1e-12)
        assert res.effect_ci_low == pytest.approx(np.exp(res.ci_low), rel=1e-12)
    else:
        assert res.n_cases is None and res.n_controls is None
        assert res.effect == res.beta


@pytest.mark.parametrize("fit", [fit_linear, fit_logistic, fit_firth, fit_robust_linear])
def test_fits_reject_an_out_of_range_term_index(fit):
    X = np.zeros((10, 2))
    y = np.array([0.0, 1.0] * 5)
    with pytest.raises(ValueError):
        fit(X, y, 5, ["a", "b"])


@pytest.mark.parametrize("fit", [fit_linear, fit_logistic, fit_firth, fit_robust_linear])
def test_fits_refuse_a_constant_tested_term_rather_than_returning_nan(fit):
    """A monomorphic variant has no estimable effect (§4.4 guardrail).

    Before this guard the Firth path produced beta=0, CI=(0,0) and p=NaN, which
    is exactly the kind of value that gets formatted into a results table and
    read as "no effect" rather than "not testable".
    """
    rng = np.random.default_rng(64)
    n = 120
    X = np.column_stack([np.zeros(n), rng.normal(size=n)])  # nobody carries it
    y = rng.binomial(1, 0.3, n).astype(float)

    with pytest.raises(ValueError, match="constant"):
        fit(X, y, 0, ["carrier", "pc1"])


def test_constant_term_is_judged_after_missing_rows_are_dropped():
    """The only carriers are dropped for a missing covariate -> not estimable."""
    rng = np.random.default_rng(65)
    n = 100
    carrier = np.zeros(n)
    carrier[[3, 17]] = 1.0
    pc = rng.normal(size=n)
    pc[[3, 17]] = np.nan
    y = rng.binomial(1, 0.3, n).astype(float)

    with pytest.raises(ValueError, match="complete cases"):
        fit_firth(np.column_stack([carrier, pc]), y, 0, ["carrier", "pc1"])


def test_binary_fits_reject_a_non_binary_outcome():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(50, 1))
    y = rng.normal(size=50)
    with pytest.raises(ValueError):
        fit_logistic(X, y, 0, ["x"])
    with pytest.raises(ValueError):
        fit_firth(X, y, 0, ["x"])
