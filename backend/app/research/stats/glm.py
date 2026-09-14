"""Regression kernels for Research Mode association analysis (Part II §4.1).

Every association analysis in the platform — single-variant, gene burden,
PheWAS, PRS-vs-outcome — bottoms out in one of these four fits. They are kept
in one module, dependency-free apart from numpy/scipy, because the results they
produce are the numbers a researcher will publish: it has to be possible to
read the estimator end to end without chasing a library's dispatch layers.

Three design decisions worth stating up front, because they are the ones a
reviewer should argue with:

1.  **Firth p-values come from the penalised likelihood ratio test, never from
    a Wald test.** Firth regression exists for the separation regime. In that
    regime the Wald statistic β̂/SE is actively misleading: as separation is
    approached the standard error grows faster than the coefficient, so the
    Wald statistic collapses toward zero and the test loses power exactly where
    the effect is strongest (the Hauck–Donner effect). The penalised LRT does
    not have this failure mode. See `fit_firth`.

2.  **Non-convergence is reported, not smoothed over.** An IRLS loop that hits
    its iteration cap has not found the MLE, and the last iterate is not "close
    enough" — under separation it is an arbitrary point on a ridge that runs to
    infinity. `converged=False` plus a warning string is the honest answer, and
    `recommend_model` exists so the user is steered to Firth *before* getting
    there.

3.  **Everything reports n after complete-case deletion.** §4.1 requires
    missingness on every result. Rows with any NaN in the design or outcome are
    dropped and `n` is the count that actually entered the fit, never the count
    the user selected.

The model literal "cox" is reserved for the survival kernel (§4.7); it is not
implemented here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import optimize, stats

# ---------------------------------------------------------------------------
# Thresholds. Named because they are policy, not arithmetic, and a researcher
# should be able to find and disagree with them.
# ---------------------------------------------------------------------------

# |β| beyond this on the logit scale is not a real effect estimate. exp(20) is
# an odds ratio of 4.8e8; nothing in germline genetics produces that. It means
# the likelihood is flat out to infinity, i.e. separation.
SEPARATION_BETA = 20.0

# A fitted probability this close to 0 or 1 means some subject is perfectly
# classified by the linear predictor — the other diagnostic for separation.
SEPARATION_PROB_EPS = 1e-8

# §4.1 automatic model selection thresholds.
MIN_CASES_FOR_LOGISTIC = 50
MIN_CARRIER_FREQ_FOR_LOGISTIC = 0.01
MIN_EXPECTED_CELL_COUNT = 5.0

Z_975 = float(stats.norm.ppf(0.975))
CHI2_95_1DF = float(stats.chi2.ppf(0.95, 1))


@dataclass
class FitResult:
    """One tested term from one model fit.

    §4.1: "Every result reports: n, cases, controls, effect (OR/beta/HR), SE,
    95% CI, p, adjusted p, model, covariates, missingness, convergence status."
    Adjusted p is added by `multiple_testing` at the family level, not here —
    a single fit does not know how many tests it belongs to.
    """

    model: str                  # "linear" | "logistic" | "firth" | "cox"
    term: str
    beta: float                 # effect on the model's natural scale (logit for binary)
    se: float
    ci_low: float               # 95% CI on the natural scale
    ci_high: float
    pvalue: float
    n: int
    n_cases: Optional[int]      # binary outcomes only
    n_controls: Optional[int]
    converged: bool
    n_iter: int
    effect_label: str           # "beta" | "OR" | "HR"
    effect: float               # beta for linear, exp(beta) for logistic/firth
    effect_ci_low: float
    effect_ci_high: float
    covariates: List[str]
    warnings: List[str] = field(default_factory=list)
    extra: Dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Shared plumbing
# ---------------------------------------------------------------------------

def _prepare(
    X: np.ndarray, y: np.ndarray, term_index: int, names: Optional[Sequence[str]]
) -> Tuple[np.ndarray, np.ndarray, str, List[str], int, int, List[str]]:
    """Complete-case filter, intercept prepend, name resolution.

    Returns (design_with_intercept, y, term_name, covariate_names, design_term_col,
    n_dropped, design_names). The design term column is term_index + 1 because the
    intercept occupies column 0 — callers pass X *without* an intercept so they
    never have to think about whether a helper already added one.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    y = np.asarray(y, dtype=float).ravel()
    if X.shape[0] != y.shape[0]:
        raise ValueError(
            "X has {} rows but y has {}".format(X.shape[0], y.shape[0]))
    n_cols = X.shape[1]
    if not (0 <= term_index < n_cols):
        raise ValueError(
            "term_index {} out of range for {} columns".format(term_index, n_cols))

    if names is None:
        names = ["x{}".format(j) for j in range(n_cols)]
    names = list(names)
    if len(names) != n_cols:
        raise ValueError(
            "names has {} entries but X has {} columns".format(len(names), n_cols))

    keep = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    n_dropped = int((~keep).sum())
    X, y = X[keep], y[keep]

    # A tested term with no variance left after complete-case filtering — a
    # variant that is monomorphic in the analysed set, or whose only carriers
    # were dropped for missing covariates — has no identifiable coefficient.
    # §4.4's guardrail is that this is *reported as such* rather than handed
    # back as an unstable estimate, so it is an error here rather than a silent
    # NaN that a caller might format into a results table.
    if X.shape[0] > 0 and np.ptp(X[:, term_index]) == 0.0:
        raise ValueError(
            "tested term '{}' is constant across all {} complete cases "
            "(value {:g}); its effect is not estimable".format(
                names[term_index], X.shape[0], float(X[0, term_index])))

    design = np.column_stack([np.ones(X.shape[0]), X])
    term = names[term_index]
    covariates = [nm for j, nm in enumerate(names) if j != term_index]
    design_names = ["(Intercept)"] + names
    return design, y, term, covariates, term_index + 1, n_dropped, design_names


def _check_binary(y: np.ndarray) -> Tuple[int, int]:
    uniq = np.unique(y)
    if not np.all(np.isin(uniq, (0.0, 1.0))):
        raise ValueError(
            "binary outcome must be coded 0/1; found {}".format(uniq[:8]))
    n_cases = int(np.sum(y == 1.0))
    n_controls = int(np.sum(y == 0.0))
    return n_cases, n_controls


def _solve_psd(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Solve A x = b, falling back to a pseudo-inverse when A is singular.

    A singular XᵀWX is itself diagnostic (collinear covariates, or weights that
    have collapsed to zero under separation). We do not want the fit to raise —
    we want it to finish so the caller can read the warnings that explain why
    the estimate is untrustworthy.
    """
    try:
        return np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(A) @ b


def _inv_psd(A: np.ndarray) -> np.ndarray:
    try:
        return np.linalg.inv(A)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(A)


def _expit(z: np.ndarray) -> np.ndarray:
    """Numerically stable logistic sigmoid.

    The naive 1/(1+exp(-z)) overflows for z < -700, which happens routinely
    while IRLS is walking away from a separated dataset.
    """
    out = np.empty_like(z)
    pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    ez = np.exp(z[~pos])
    out[~pos] = ez / (1.0 + ez)
    return out


def _safe_exp(z: float) -> float:
    """exp() that returns inf rather than raising an overflow warning.

    Under separation the logit-scale CI runs to +/-inf, so an infinite odds
    ratio is the *correct* answer to report; it just should not look like a
    numerical accident in the logs.
    """
    with np.errstate(over="ignore"):
        return float(np.exp(z))


def _log_likelihood(y: np.ndarray, eta: np.ndarray) -> float:
    """Bernoulli log-likelihood from the linear predictor.

    Written in terms of eta rather than pi so that log(pi) never becomes
    log(0.0) = -inf for a subject the model has classified perfectly:
    y*eta - log1p(exp(eta)) is exact, and the log1p term is evaluated in the
    stable branch.
    """
    # log(1 + e^eta) = max(eta,0) + log1p(e^-|eta|)
    log1pexp = np.maximum(eta, 0.0) + np.log1p(np.exp(-np.abs(eta)))
    return float(np.sum(y * eta - log1pexp))


# ---------------------------------------------------------------------------
# 1. Linear
# ---------------------------------------------------------------------------

def fit_linear(
    X: np.ndarray,
    y: np.ndarray,
    term_index: int,
    names: Optional[Sequence[str]] = None,
) -> FitResult:
    """OLS for a quantitative outcome (§4.1 "Quantitative: linear regression").

    Solved by least squares rather than by forming (XᵀX)⁻¹ explicitly: with
    ancestry PCs and age² in the covariate set the design is routinely
    ill-conditioned, and the normal equations square the condition number.
    """
    design, yy, term, covariates, col, n_dropped, design_names = _prepare(
        X, y, term_index, names)
    n, p = design.shape
    warnings: List[str] = []

    if n <= p:
        raise ValueError(
            "n={} observations for {} parameters — model is not identifiable".format(n, p))

    beta, _, rank, _ = np.linalg.lstsq(design, yy, rcond=None)
    if rank < p:
        warnings.append(
            "design matrix is rank deficient (rank {} of {}); estimates are not "
            "uniquely identified".format(rank, p))

    resid = yy - design @ beta
    df_resid = n - p
    rss = float(resid @ resid)
    sigma2 = rss / df_resid
    xtx_inv = _inv_psd(design.T @ design)
    se_all = np.sqrt(np.maximum(np.diag(xtx_inv) * sigma2, 0.0))

    b = float(beta[col])
    se = float(se_all[col])
    tcrit = float(stats.t.ppf(0.975, df_resid))
    lo, hi = b - tcrit * se, b + tcrit * se
    tstat = b / se if se > 0 else np.inf
    pval = float(2.0 * stats.t.sf(abs(tstat), df_resid))

    tss = float(np.sum((yy - yy.mean()) ** 2))
    r2 = 1.0 - rss / tss if tss > 0 else float("nan")
    adj_r2 = 1.0 - (1.0 - r2) * (n - 1) / df_resid if tss > 0 else float("nan")

    return FitResult(
        model="linear", term=term, beta=b, se=se, ci_low=lo, ci_high=hi,
        pvalue=pval, n=n, n_cases=None, n_controls=None,
        converged=True, n_iter=1,
        effect_label="beta", effect=b, effect_ci_low=lo, effect_ci_high=hi,
        covariates=covariates, warnings=warnings,
        extra={
            "r2": r2,
            "adj_r2": adj_r2,
            "residual_se": float(np.sqrt(sigma2)),
            "df_resid": int(df_resid),
            "rss": rss,
            "n_dropped_missing": n_dropped,
            "coefficients": {nm: float(bv) for nm, bv in zip(design_names, beta)},
            "se_method": "classical (homoskedastic)",
        },
    )


# ---------------------------------------------------------------------------
# 2. Logistic
# ---------------------------------------------------------------------------

def _irls_logistic(
    design: np.ndarray, y: np.ndarray, max_iter: int, tol: float
) -> Tuple[np.ndarray, np.ndarray, bool, int]:
    """Newton–Raphson / IRLS for the unpenalised Bernoulli likelihood.

    Step halving is applied whenever a full Newton step lowers the
    log-likelihood. Without it the loop can overshoot into a region where W is
    numerically zero and the next Hessian is unusable.
    """
    n, p = design.shape
    beta = np.zeros(p)
    # Start the intercept at the sample logit: closer to the answer and it
    # keeps the first Hessian well conditioned when cases are rare.
    ybar = float(np.clip(y.mean(), 1e-6, 1 - 1e-6))
    beta[0] = np.log(ybar / (1.0 - ybar))
    ll = _log_likelihood(y, design @ beta)

    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        eta = design @ beta
        pi = _expit(eta)
        w = np.maximum(pi * (1.0 - pi), 1e-12)
        score = design.T @ (y - pi)
        hess = design.T @ (design * w[:, None])
        step = _solve_psd(hess, score)

        # Step halving on the likelihood, not on the gradient norm.
        t = 1.0
        for _ in range(30):
            cand = beta + t * step
            ll_new = _log_likelihood(y, design @ cand)
            if ll_new >= ll - 1e-12:
                break
            t *= 0.5
        else:
            cand, ll_new = beta + step, _log_likelihood(y, design @ (beta + step))

        delta = float(np.max(np.abs(cand - beta)))
        beta, ll = cand, ll_new
        if delta < tol:
            converged = True
            break
    return beta, _expit(design @ beta), converged, it


def fit_logistic(
    X: np.ndarray,
    y: np.ndarray,
    term_index: int,
    names: Optional[Sequence[str]] = None,
    max_iter: int = 50,
    tol: float = 1e-8,
) -> FitResult:
    """Standard (unpenalised) logistic regression for a binary outcome.

    Separation is reported rather than hidden. What it means for the result:
    under complete or quasi-complete separation the maximum likelihood estimate
    of β does not exist — the likelihood increases monotonically as β → ±∞ — so
    the returned coefficient is wherever the iteration stopped, and its standard
    error, CI and Wald p-value are all meaningless. The correct response is to
    re-run with `fit_firth`, which is what `recommend_model` will say.
    """
    design, yy, term, covariates, col, n_dropped, design_names = _prepare(
        X, y, term_index, names)
    n_cases, n_controls = _check_binary(yy)
    n, p = design.shape
    warnings: List[str] = []

    if n <= p:
        raise ValueError(
            "n={} observations for {} parameters — model is not identifiable".format(n, p))

    beta, pi, converged, n_iter = _irls_logistic(design, yy, max_iter, tol)

    # --- separation diagnostics -------------------------------------------
    max_abs_beta = float(np.max(np.abs(beta)))
    extreme_fitted = bool(
        np.any(pi < SEPARATION_PROB_EPS) or np.any(pi > 1.0 - SEPARATION_PROB_EPS))
    separated = max_abs_beta > SEPARATION_BETA or extreme_fitted
    if separated:
        warnings.append(
            "complete/quasi-complete separation detected (max |beta| = {:.1f}"
            "{}); the maximum likelihood estimate does not exist and this "
            "coefficient, SE, CI and p-value are not interpretable — refit with "
            "Firth penalised logistic regression".format(
                max_abs_beta, ", fitted probabilities at 0/1" if extreme_fitted else ""))
    # A zero cell in the exposure x outcome table is quasi-complete separation
    # in the tested term specifically, which is the case we most care about.
    tcol = design[:, col]
    if np.all(np.isin(np.unique(tcol), (0.0, 1.0))):
        cells = [int(np.sum((tcol == e) & (yy == o))) for e in (0.0, 1.0) for o in (0.0, 1.0)]
        if min(cells) == 0:
            warnings.append(
                "tested term '{}' has a zero cell in its 2x2 table with the "
                "outcome (counts {}) — the odds ratio is unbounded".format(term, cells))
    if not converged:
        warnings.append(
            "IRLS did not converge in {} iterations (max coefficient change "
            "exceeded tol={:g}); estimates are the last iterate and should not "
            "be reported".format(max_iter, tol))

    w = np.maximum(pi * (1.0 - pi), 1e-12)
    cov = _inv_psd(design.T @ (design * w[:, None]))
    se = float(np.sqrt(max(cov[col, col], 0.0)))
    b = float(beta[col])
    lo, hi = b - Z_975 * se, b + Z_975 * se
    z = b / se if se > 0 else np.inf
    pval = float(2.0 * stats.norm.sf(abs(z)))

    ll = _log_likelihood(yy, design @ beta)
    ybar = yy.mean()
    ll_null = float(n * (ybar * np.log(ybar) + (1 - ybar) * np.log(1 - ybar))) \
        if 0 < ybar < 1 else 0.0

    return FitResult(
        model="logistic", term=term, beta=b, se=se, ci_low=lo, ci_high=hi,
        pvalue=pval, n=n, n_cases=n_cases, n_controls=n_controls,
        converged=converged and not separated, n_iter=n_iter,
        effect_label="OR", effect=_safe_exp(b),
        effect_ci_low=_safe_exp(lo), effect_ci_high=_safe_exp(hi),
        covariates=covariates, warnings=warnings,
        extra={
            "loglik": ll,
            "loglik_null": ll_null,
            "pseudo_r2_mcfadden": float(1.0 - ll / ll_null) if ll_null != 0 else float("nan"),
            "max_abs_beta": max_abs_beta,
            "min_fitted_prob": float(pi.min()),
            "max_fitted_prob": float(pi.max()),
            "separation": separated,
            "irls_converged": converged,
            "coefficients": {nm: float(bv) for nm, bv in zip(design_names, beta)},
            "n_dropped_missing": n_dropped,
            "test": "Wald",
        },
    )


# ---------------------------------------------------------------------------
# 3. Firth penalised logistic
# ---------------------------------------------------------------------------

def _firth_penalised_loglik(design: np.ndarray, y: np.ndarray, beta: np.ndarray) -> float:
    """l*(β) = l(β) + ½·log det(XᵀWX), W = diag(π(1−π)).

    The penalty is the log of the Jeffreys prior. It vanishes at the boundary
    (as any π → 0 or 1, W → 0, det → 0, log det → −∞), which is precisely why a
    penalised fit cannot run off to infinity under separation: the penalty
    charges an unbounded price for perfect prediction.
    """
    eta = design @ beta
    pi = _expit(eta)
    w = pi * (1.0 - pi)
    info = design.T @ (design * w[:, None])
    sign, logdet = np.linalg.slogdet(info)
    if sign <= 0:
        # Information matrix numerically singular: the penalty is -inf.
        return -np.inf
    return _log_likelihood(y, eta) + 0.5 * logdet


def _hat_diagonal(design: np.ndarray, w: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Diagonal of H = W^{1/2} X (XᵀWX)⁻¹ Xᵀ W^{1/2}, plus XᵀWX.

    Computed row-wise as h_i = w_i · x_iᵀ (XᵀWX)⁻¹ x_i rather than by forming
    the n x n hat matrix, which would be a biobank-sized allocation.
    """
    info = design.T @ (design * w[:, None])
    info_inv = _inv_psd(info)
    h = w * np.einsum("ij,jk,ik->i", design, info_inv, design)
    return h, info


def _firth_newton(
    design: np.ndarray,
    y: np.ndarray,
    max_iter: int,
    tol: float,
    fixed_col: Optional[int] = None,
    fixed_value: float = 0.0,
    start: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, float, bool, int]:
    """Maximise the Firth-penalised likelihood, optionally with one β fixed.

    The modified score is

        U*(β_j) = Σ_i x_ij · [ y_i − π_i + h_i·(½ − π_i) ]

    i.e. the ordinary score with the working response y_i replaced by
    y_i + h_i(½ − π_i). Each observation is nudged toward ½ in proportion to its
    leverage, which is the finite-sample bias correction.

    When `fixed_col` is given, that coefficient is held at `fixed_value` and only
    the remaining coefficients are updated — but h_i and the penalty are still
    computed from the **full** design. That matters: the profile penalised
    likelihood must use the same Jeffreys penalty as the unconstrained fit, or
    the two log-likelihoods are on different scales and their difference is not
    a likelihood ratio.
    """
    n, p = design.shape
    free = np.array([j for j in range(p) if j != fixed_col], dtype=int)

    beta = np.zeros(p) if start is None else np.array(start, dtype=float)
    if fixed_col is not None:
        beta[fixed_col] = fixed_value
    if start is None:
        ybar = float(np.clip(y.mean(), 1e-6, 1 - 1e-6))
        beta[0] = np.log(ybar / (1.0 - ybar))

    ll = _firth_penalised_loglik(design, y, beta)
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        pi = _expit(design @ beta)
        w = np.maximum(pi * (1.0 - pi), 1e-12)
        h, _ = _hat_diagonal(design, w)
        u = design[:, free].T @ (y - pi + h * (0.5 - pi))
        info_free = design[:, free].T @ (design[:, free] * w[:, None])
        step = _solve_psd(info_free, u)

        # Cap the raw step: early iterations on a separated dataset can propose
        # a jump of many hundreds on the logit scale, from which the line
        # search cannot recover in a reasonable number of halvings.
        step_norm = float(np.max(np.abs(step)))
        if step_norm > 5.0:
            step = step * (5.0 / step_norm)

        t = 1.0
        cand = beta.copy()
        ll_new = ll
        for _ in range(40):
            trial = beta.copy()
            trial[free] = beta[free] + t * step
            ll_trial = _firth_penalised_loglik(design, y, trial)
            if np.isfinite(ll_trial) and ll_trial >= ll - 1e-12:
                cand, ll_new = trial, ll_trial
                break
            t *= 0.5
        else:
            # No uphill step found: we are at (or numerically at) the optimum.
            converged = True
            break

        delta = float(np.max(np.abs(cand[free] - beta[free])))
        beta, ll = cand, ll_new
        if delta < tol:
            converged = True
            break
    return beta, ll, converged, it


def _firth_profile_ci(
    design: np.ndarray,
    y: np.ndarray,
    col: int,
    beta_hat: np.ndarray,
    ll_full: float,
    se: float,
    max_iter: int,
    tol: float,
) -> Tuple[Optional[float], Optional[float], str]:
    """Profile penalised-likelihood 95% CI for one coefficient.

    The endpoints are the β values where the penalised LRT statistic against the
    unconstrained fit equals χ²₀.₉₅,₁ = 3.841. This is the CI that inverts the
    test we actually report, so the CI and the p-value agree by construction —
    a Wald CI can exclude zero while the LRT p-value is > 0.05, which is
    indefensible on a results page.

    Returns (lo, hi, method). Either endpoint may be None if the root could not
    be bracketed, in which case the caller falls back to Wald and warns.
    """
    b_hat = float(beta_hat[col])

    def deviance(c: float) -> float:
        _, ll_c, _, _ = _firth_newton(
            design, y, max_iter, tol, fixed_col=col, fixed_value=c, start=beta_hat)
        return 2.0 * (ll_full - ll_c) - CHI2_95_1DF

    step = max(se, 0.5) if np.isfinite(se) and se > 0 else 0.5
    endpoints: List[Optional[float]] = []
    for direction in (-1.0, 1.0):
        lo_c, hi_c = b_hat, None
        d = step
        for _ in range(40):
            trial = b_hat + direction * d
            if deviance(trial) > 0:
                hi_c = trial
                break
            lo_c = trial
            d *= 1.6
            if abs(d) > 1e4:
                break
        if hi_c is None:
            endpoints.append(None)
            continue
        try:
            root = optimize.brentq(deviance, min(lo_c, hi_c), max(lo_c, hi_c),
                                   xtol=1e-6, rtol=1e-10, maxiter=200)
            endpoints.append(float(root))
        except (ValueError, RuntimeError):
            endpoints.append(None)

    lo, hi = endpoints
    if lo is not None and hi is not None and lo > hi:
        lo, hi = hi, lo
    method = "profile penalised likelihood" if (lo is not None and hi is not None) \
        else "partial (profile where bracketed)"
    return lo, hi, method


def fit_firth(
    X: np.ndarray,
    y: np.ndarray,
    term_index: int,
    names: Optional[Sequence[str]] = None,
    max_iter: int = 100,
    tol: float = 1e-8,
) -> FitResult:
    """Firth penalised logistic regression — required by §4.1.

    Why it is required: germline association tests are dominated by rare
    exposures. A variant carried by 0.4% of a cohort with 43 cases will often
    produce a 2x2 table with a zero cell, and standard logistic regression then
    has no finite MLE. Firth's penalty (the Jeffreys prior) removes the O(1/n)
    bias of the MLE and, as a side effect, guarantees a finite estimate even
    under complete separation.

    p-value: penalised likelihood ratio test, 1 df, comparing the unconstrained
    fit against the fit with β_term ≡ 0. Not Wald — see the module docstring.

    CI: profile penalised likelihood, inverting the same test. If an endpoint
    cannot be bracketed the result falls back to a Wald interval and says so in
    `warnings`, because a Wald CI next to an LRT p-value is a real
    inconsistency and the reader is entitled to know.
    """
    design, yy, term, covariates, col, n_dropped, design_names = _prepare(
        X, y, term_index, names)
    n_cases, n_controls = _check_binary(yy)
    n, p = design.shape
    warnings: List[str] = []

    if n <= p:
        raise ValueError(
            "n={} observations for {} parameters — model is not identifiable".format(n, p))

    beta, ll_full, converged, n_iter = _firth_newton(design, yy, max_iter, tol)
    if not converged:
        warnings.append(
            "Firth Newton iteration did not converge in {} iterations".format(max_iter))

    pi = _expit(design @ beta)
    w = np.maximum(pi * (1.0 - pi), 1e-12)
    info = design.T @ (design * w[:, None])
    cov = _inv_psd(info)
    se = float(np.sqrt(max(cov[col, col], 0.0)))
    b = float(beta[col])

    # --- penalised LRT ----------------------------------------------------
    _, ll_null, conv_null, _ = _firth_newton(
        design, yy, max_iter, tol, fixed_col=col, fixed_value=0.0, start=beta)
    if not conv_null:
        warnings.append(
            "constrained (beta_{}=0) Firth fit did not converge; the LRT "
            "p-value may be inaccurate".format(term))
    if not (np.isfinite(ll_full) and np.isfinite(ll_null)):
        # Both penalised likelihoods are -inf, which means X'WX was singular at
        # the optimum: the design is degenerate (collinear covariates) and the
        # likelihood ratio is not defined. Refusing to produce a number beats
        # producing chi2.sf(nan).
        raise ValueError(
            "penalised likelihood is not finite for term '{}' — X'WX is "
            "singular, which means the design is collinear or degenerate; "
            "check the covariate set".format(term))
    lrt = max(2.0 * (ll_full - ll_null), 0.0)
    pval = float(stats.chi2.sf(lrt, 1))

    # --- profile CI -------------------------------------------------------
    lo, hi, ci_method = _firth_profile_ci(
        design, yy, col, beta, ll_full, se, max_iter, tol)
    if lo is None or hi is None:
        wald_lo, wald_hi = b - Z_975 * se, b + Z_975 * se
        if lo is None:
            lo = wald_lo
        if hi is None:
            hi = wald_hi
        ci_method = "Wald (profile likelihood root could not be bracketed)"
        warnings.append(
            "confidence interval is Wald-based while the p-value is from the "
            "penalised likelihood ratio test; the two may disagree near the "
            "significance boundary")

    # Separation is not an error for Firth — it is the case it was built for —
    # but the reader should still know the data were separated, because the
    # estimate is then driven by the penalty as much as by the data.
    #
    # It cannot be detected from the Firth fit itself: the whole point of the
    # penalty is that the Firth coefficients and fitted probabilities stay well
    # inside the interior even when the data are perfectly separated (for a
    # separated 2x2 the fitted probabilities are 0.5/11 and 10.5/11, nowhere
    # near 0 or 1). So we ask the *unpenalised* likelihood, which is the one
    # that blows up. A short bounded IRLS run is enough to see it happen.
    # A fixed, loose tolerance: this probe only has to distinguish "converges"
    # from "runs away", and a caller who asks for tol=1e-14 should not get a
    # spurious separation warning just because IRLS needed one more iteration.
    _, pi_ml, ml_converged, _ = _irls_logistic(design, yy, 25, max(tol, 1e-8))
    separated = (
        not ml_converged
        or np.any(pi_ml < SEPARATION_PROB_EPS)
        or np.any(pi_ml > 1.0 - SEPARATION_PROB_EPS))
    if separated:
        warnings.append(
            "complete or quasi-complete separation detected: the unpenalised "
            "maximum likelihood estimate does not exist. The Firth estimate is "
            "finite by construction, but it is shaped by the Jeffreys penalty "
            "as much as by the data and should be read as a shrunk estimate, "
            "not as the sample odds ratio")

    return FitResult(
        model="firth", term=term, beta=b, se=se, ci_low=float(lo), ci_high=float(hi),
        pvalue=pval, n=n, n_cases=n_cases, n_controls=n_controls,
        converged=converged and conv_null, n_iter=n_iter,
        effect_label="OR", effect=_safe_exp(b),
        effect_ci_low=_safe_exp(lo), effect_ci_high=_safe_exp(hi),
        covariates=covariates, warnings=warnings,
        extra={
            "loglik_penalised": float(ll_full),
            "loglik_penalised_null": float(ll_null),
            "lrt_stat": float(lrt),
            "lrt_df": 1,
            "test": "penalised likelihood ratio",
            "ci_method": ci_method,
            "se_method": "Wald (inverse Fisher information at the Firth estimate)",
            "penalty": "Jeffreys prior, 0.5*log det(X'WX)",
            "min_fitted_prob": float(pi.min()),
            "max_fitted_prob": float(pi.max()),
            "separation": bool(separated),
            "coefficients": {nm: float(bv) for nm, bv in zip(design_names, beta)},
            "n_dropped_missing": n_dropped,
        },
    )


# ---------------------------------------------------------------------------
# 4. Robust linear (HC3)
# ---------------------------------------------------------------------------

def fit_robust_linear(
    X: np.ndarray,
    y: np.ndarray,
    term_index: int,
    names: Optional[Sequence[str]] = None,
) -> FitResult:
    """OLS point estimates with HC3 heteroskedasticity-consistent SEs.

    Why HC3 and not HC0: the classical SE assumes constant residual variance,
    which quantitative phenotypes routinely violate (biomarker variance scales
    with the mean; age-dependent traits fan out). HC0 fixes the bias
    asymptotically but is anti-conservative in small samples. HC3 divides the
    squared residual by (1−h_i)², which inflates the contribution of
    high-leverage points — exactly the rare-genotype carriers who would
    otherwise drive an association with a deceptively small SE. It is the
    recommended default below n ≈ 250 and is never much worse above it.

    The point estimates are identical to `fit_linear`; only the variance and
    everything derived from it change.
    """
    design, yy, term, covariates, col, n_dropped, design_names = _prepare(
        X, y, term_index, names)
    n, p = design.shape
    warnings: List[str] = []

    if n <= p:
        raise ValueError(
            "n={} observations for {} parameters — model is not identifiable".format(n, p))

    beta, _, rank, _ = np.linalg.lstsq(design, yy, rcond=None)
    if rank < p:
        warnings.append(
            "design matrix is rank deficient (rank {} of {})".format(rank, p))

    resid = yy - design @ beta
    xtx_inv = _inv_psd(design.T @ design)
    h = np.einsum("ij,jk,ik->i", design, xtx_inv, design)
    # Guard leverage against exactly 1: a point with h_i = 1 is fit perfectly
    # and HC3's (1-h)^-2 is undefined there.
    denom = np.maximum((1.0 - h) ** 2, 1e-10)
    if np.any(h > 1.0 - 1e-6):
        warnings.append(
            "{} observation(s) have leverage h_i ~ 1; HC3 weights there are "
            "capped and the robust SE is understated".format(int(np.sum(h > 1.0 - 1e-6))))
    omega = (resid ** 2) / denom
    meat = design.T @ (design * omega[:, None])
    cov_hc3 = xtx_inv @ meat @ xtx_inv

    b = float(beta[col])
    se = float(np.sqrt(max(cov_hc3[col, col], 0.0)))
    df_resid = n - p
    tcrit = float(stats.t.ppf(0.975, df_resid))
    lo, hi = b - tcrit * se, b + tcrit * se
    tstat = b / se if se > 0 else np.inf
    pval = float(2.0 * stats.t.sf(abs(tstat), df_resid))

    rss = float(resid @ resid)
    tss = float(np.sum((yy - yy.mean()) ** 2))
    r2 = 1.0 - rss / tss if tss > 0 else float("nan")
    se_classical = float(np.sqrt(max(np.diag(xtx_inv)[col] * rss / df_resid, 0.0)))

    return FitResult(
        model="linear", term=term, beta=b, se=se, ci_low=lo, ci_high=hi,
        pvalue=pval, n=n, n_cases=None, n_controls=None,
        converged=True, n_iter=1,
        effect_label="beta", effect=b, effect_ci_low=lo, effect_ci_high=hi,
        covariates=covariates, warnings=warnings,
        extra={
            "se_method": "HC3",
            "se_classical": se_classical,
            "r2": r2,
            "df_resid": int(df_resid),
            "max_leverage": float(h.max()),
            "coefficients": {nm: float(bv) for nm, bv in zip(design_names, beta)},
            "n_dropped_missing": n_dropped,
        },
    )


# ---------------------------------------------------------------------------
# 5. Automatic model selection (§4.1)
# ---------------------------------------------------------------------------

def recommend_model(
    y: np.ndarray,
    kind: str,
    exposure: Optional[np.ndarray] = None,
    n_related: int = 0,
) -> Dict[str, Any]:
    """Propose a model, state why, and leave the decision with the user.

    §4.1: "The platform proposes a model from the variable types and data
    structure, states why, and lets the user override." The rationale is written
    as a sentence a researcher can disagree with — it names the counts that
    drove the choice, not just the choice.

    `kind` is "quantitative" | "binary" (anything else is rejected rather than
    guessed at). `exposure` is the tested variable; carrier frequency is the
    fraction of non-missing subjects with a non-zero value, which is the right
    notion for both 0/1 carrier flags and 0/1/2 dosages.
    """
    y = np.asarray(y, dtype=float).ravel()
    y = y[np.isfinite(y)]
    reasons: List[str] = []

    if kind == "quantitative":
        model = "linear"
        rationale = (
            "Outcome is quantitative with {:,} non-missing observations; "
            "ordinary least squares is the default. Check the residual "
            "distribution — if it is skewed, apply the rank-based "
            "inverse-normal transform or use robust (HC3) standard errors."
        ).format(len(y))
        alternatives = ["rank_inverse_normal + linear", "robust_linear", "mixed"]

    elif kind == "binary":
        n_cases = int(np.sum(y == 1.0))
        n_controls = int(np.sum(y == 0.0))

        carrier_freq: Optional[float] = None
        min_expected: Optional[float] = None
        if exposure is not None:
            e = np.asarray(exposure, dtype=float).ravel()
            ok = np.isfinite(e)
            if len(e) == len(np.asarray(y)) or ok.sum() > 0:
                e_ok = e[ok]
                if e_ok.size:
                    carrier_freq = float(np.mean(e_ok != 0.0))
            # Expected cell counts under independence for the collapsed
            # carrier x outcome 2x2 table — the standard "is this table too
            # sparse to trust an asymptotic test" check.
            e_full = np.asarray(exposure, dtype=float).ravel()
            y_full = np.asarray(y, dtype=float).ravel()
            if e_full.shape == y_full.shape:
                m = np.isfinite(e_full) & np.isfinite(y_full)
                if m.sum() > 0:
                    carrier = (e_full[m] != 0.0).astype(float)
                    case = y_full[m]
                    n_t = int(m.sum())
                    row = [float(np.sum(carrier == 0)), float(np.sum(carrier == 1))]
                    colc = [float(np.sum(case == 0)), float(np.sum(case == 1))]
                    min_expected = min(r * c / n_t for r in row for c in colc)

        if n_cases < MIN_CASES_FOR_LOGISTIC:
            reasons.append("only {} cases".format(n_cases))
        if carrier_freq is not None and carrier_freq < MIN_CARRIER_FREQ_FOR_LOGISTIC:
            reasons.append("a carrier frequency of {:.1f}%".format(carrier_freq * 100))
        if min_expected is not None and min_expected < MIN_EXPECTED_CELL_COUNT:
            reasons.append(
                "a minimum expected cell count of {:.1f}".format(min_expected))

        if reasons:
            model = "firth"
            freq_txt = ("" if carrier_freq is None
                        else " and a carrier frequency of {:.1f}%".format(carrier_freq * 100))
            rationale = (
                "Outcome is binary with {} cases and {} controls{}. "
                "Standard logistic regression will likely separate ({}). "
                "Recommending Firth logistic regression."
            ).format(n_cases, n_controls, freq_txt, "; ".join(reasons))
            alternatives = ["logistic", "fisher_exact"]
        else:
            model = "logistic"
            freq_txt = ("" if carrier_freq is None
                        else " and a carrier frequency of {:.1f}%".format(carrier_freq * 100))
            rationale = (
                "Outcome is binary with {} cases and {} controls{}; case count "
                "and exposure frequency are both adequate for the asymptotic "
                "approximation, so standard logistic regression applies."
            ).format(n_cases, n_controls, freq_txt)
            alternatives = ["firth", "fisher_exact"]
    else:
        raise ValueError(
            "kind must be 'quantitative' or 'binary', got {!r}".format(kind))

    if n_related > 0:
        rationale += (
            " {} related sample(s) are present: a mixed model with a kinship "
            "random effect would be preferable, and relatedness is currently "
            "handled by pruning rather than modelled, so the effective sample "
            "size is reduced."
        ).format(n_related)
        alternatives = list(alternatives) + ["mixed"]

    return {"model": model, "rationale": rationale, "alternatives": alternatives}


# ---------------------------------------------------------------------------
# 6. Rank-based inverse-normal transform
# ---------------------------------------------------------------------------

def rank_inverse_normal(v: np.ndarray, c: float = 0.375) -> np.ndarray:
    """Rank-based inverse-normal transform (§4.1 quantitative).

    z_i = Φ⁻¹((r_i − c) / (n − 2c + 1)) with c = 3/8 (Blom). Applied to a
    phenotype it guarantees the outcome is normal, which makes a linear model's
    p-values valid regardless of the original distribution — at the cost of
    destroying the original units, so the resulting beta is in SD units of the
    transformed trait and must be labelled as such.

    Ties take the average rank, which maps tied values to the same z. NaN is
    passed through and excluded from the ranking, so n is the non-missing count.
    """
    v = np.asarray(v, dtype=float).ravel()
    out = np.full(v.shape, np.nan)
    obs = np.isfinite(v)
    n = int(obs.sum())
    if n == 0:
        return out
    ranks = stats.rankdata(v[obs], method="average")
    out[obs] = stats.norm.ppf((ranks - c) / (n - 2.0 * c + 1.0))
    return out


# ---------------------------------------------------------------------------
# 7. Multiple testing (§5.6)
# ---------------------------------------------------------------------------

def multiple_testing(pvalues: Sequence[float], method: str = "fdr_bh") -> np.ndarray:
    """Adjusted p-values, returned in the input order.

    §5.6: "The platform never displays an uncorrected p-value from a multi-test
    analysis without its correction adjacent." This function is the other half
    of that contract — the caller must attach the output to the result rows.

    "bonferroni" controls the family-wise error rate: p·m, capped at 1.
    "fdr_bh" controls the false discovery rate by the Benjamini–Hochberg step-up
    procedure, with the enforced monotonicity (each adjusted value is the
    running minimum from the largest p downward) that makes the output a valid
    q-value sequence. BH assumes independence or positive regression
    dependence — for a PheWAS over correlated phenotypes that assumption is
    approximate, which §4.3 requires be stated alongside the numbers.

    NaN p-values propagate as NaN and are excluded from m.
    """
    p = np.asarray(pvalues, dtype=float).ravel()
    out = np.full(p.shape, np.nan)
    obs = np.isfinite(p)
    m = int(obs.sum())
    if m == 0:
        return out
    if np.any((p[obs] < 0) | (p[obs] > 1)):
        raise ValueError("p-values must lie in [0, 1]")

    pv = p[obs]
    if method == "bonferroni":
        out[obs] = np.minimum(pv * m, 1.0)
    elif method == "fdr_bh":
        order = np.argsort(pv, kind="mergesort")
        ranked = pv[order]
        q = ranked * m / np.arange(1, m + 1)
        # Step-up: walking down from the largest p, an adjusted value can never
        # exceed the one above it.
        q = np.minimum.accumulate(q[::-1])[::-1]
        adj = np.empty(m)
        adj[order] = np.minimum(q, 1.0)
        out[obs] = adj
    else:
        raise ValueError(
            "method must be 'bonferroni' or 'fdr_bh', got {!r}".format(method))
    return out


__all__ = [
    "FitResult",
    "fit_linear",
    "fit_logistic",
    "fit_firth",
    "fit_robust_linear",
    "recommend_model",
    "rank_inverse_normal",
    "multiple_testing",
]
