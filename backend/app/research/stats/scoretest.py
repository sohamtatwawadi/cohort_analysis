"""Vectorised score test — the scan half of a two-stage GWAS.

The GWAS scan used to fit a full GLM per variant. That is the right answer for
one variant and the wrong one for eighty thousand: the covariate block is
identical every time, so the iteration that dominates the cost is repeated
unchanged for every marker. On a 3,000-sample / 82,000-variant scan it ran for
over a quarter of an hour.

The score test evaluates the null model ONCE and then tests each variant against
those residuals in closed form, which reduces to two matrix products over the
whole genotype block. Under the null the score statistic and the likelihood-ratio
statistic are asymptotically equivalent, so the scan ranks variants the same way;
where they differ is at large effects, which is precisely where a scan hands off
anyway.

So this module is a filter, not a replacement. It finds the candidates; the exact
GLM then characterises them, because an odds ratio and a confidence interval have
to come from the fitted model, not from a score statistic that never estimated
one. That two-stage shape is what production GWAS tools do, and for the same
reason.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
from scipy import stats

from ..types import MISSING

# Below this, the variance of the score is numerically indistinguishable from
# zero — a monomorphic variant, or one collinear with a covariate.
MIN_VARIANCE = 1e-10


def _design(covariates: Optional[np.ndarray], n: int) -> np.ndarray:
    """Intercept plus covariates. The intercept is never optional: without it
    the residuals are not centred and every score is biased by the mean."""
    if covariates is None or covariates.size == 0:
        return np.ones((n, 1))
    cov = np.asarray(covariates, dtype=float)
    if cov.ndim == 1:
        cov = cov.reshape(-1, 1)
    return np.column_stack([np.ones(n), cov])


def _null_logistic(y: np.ndarray, X: np.ndarray,
                   max_iter: int = 50, tol: float = 1e-9
                   ) -> Tuple[np.ndarray, np.ndarray]:
    """IRLS on the covariate-only model. Returns fitted probabilities and the
    working weights."""
    n, k = X.shape
    beta = np.zeros(k)
    # A sensible start: the log-odds of the overall rate.
    rate = float(np.clip(y.mean(), 1e-6, 1 - 1e-6))
    beta[0] = np.log(rate / (1 - rate))

    for _ in range(max_iter):
        eta = X @ beta
        p = 1.0 / (1.0 + np.exp(-np.clip(eta, -500, 500)))
        w = np.clip(p * (1 - p), 1e-10, None)
        z = eta + (y - p) / w
        XtW = X.T * w
        try:
            new = np.linalg.solve(XtW @ X, XtW @ z)
        except np.linalg.LinAlgError:
            new = np.linalg.lstsq(XtW @ X, XtW @ z, rcond=None)[0]
        if np.max(np.abs(new - beta)) < tol:
            beta = new
            break
        beta = new

    eta = X @ beta
    p = 1.0 / (1.0 + np.exp(-np.clip(eta, -500, 500)))
    return p, np.clip(p * (1 - p), 1e-10, None)


def _impute_and_centre(G: np.ndarray) -> np.ndarray:
    """Mean-impute missing dosages, per variant.

    Dropping samples per variant would give every variant a different n and make
    the single null fit invalid. Mean imputation keeps the design aligned; it
    biases toward the null (an imputed genotype carries no signal), which is the
    safe direction for a scan.
    """
    G = G.astype(float)
    missing = G == MISSING
    if missing.any():
        G[missing] = np.nan
        means = np.nanmean(G, axis=1)
        means = np.where(np.isfinite(means), means, 0.0)
        G = np.where(np.isnan(G), means[:, None], G)
    return G


# Variants per block. The whole point of the vectorised form is the matrix
# product, but materialising 80,000 x 3,000 as float64 is ~2GB and the working
# copies double it. Blocking keeps the speed and bounds the memory.
CHUNK = 4000


def score_scan(G: np.ndarray, y: np.ndarray,
               covariates: Optional[np.ndarray] = None,
               kind: str = "binary",
               progress=None) -> Dict[str, np.ndarray]:
    """Score test for every row of `G` (variants x samples) against `y`.

    Returns arrays aligned to G's rows: chi2, p, and the score and its variance.
    A variant whose score variance underflows gets p = nan rather than a p-value
    manufactured from a division by ~0.
    """
    G = np.asarray(G)
    y = np.asarray(y, dtype=float)
    m, n = G.shape
    if y.shape[0] != n:
        raise ValueError(
            "genotypes are {} variants x {} samples but y has {} entries — the "
            "matrix is probably transposed".format(m, n, y.shape[0]))

    X = _design(covariates, n)

    if kind == "binary":
        p0, w = _null_logistic(y, X)
        resid = y - p0
    else:
        # Linear: the null fit is a least-squares projection and the weights are
        # constant, so sigma^2 scales out of the ratio and is applied once.
        beta0 = np.linalg.lstsq(X, y, rcond=None)[0]
        fitted = X @ beta0
        resid = y - fitted
        dof = max(n - X.shape[1], 1)
        sigma2 = float(resid @ resid) / dof
        w = np.full(n, 1.0 / sigma2)
        resid = resid / sigma2

    # V = G'WG - (G'WX)(X'WX)^-1(X'WG), the covariate-adjusted variance.
    XtW = X.T * w
    XtWX = XtW @ X
    try:
        XtWX_inv = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        XtWX_inv = np.linalg.pinv(XtWX)

    U = np.empty(m)
    V = np.empty(m)
    for start in range(0, m, CHUNK):
        stop = min(start + CHUNK, m)
        block = _impute_and_centre(G[start:stop])
        U[start:stop] = block @ resid
        BW = block * w                                  # (b x n)
        diag = np.einsum("ij,ij->i", BW, block)         # G'WG per variant
        B = BW @ X                                      # (b x k)
        V[start:stop] = diag - np.einsum("ij,ij->i", B @ XtWX_inv, B)
        if progress:
            progress(stop, m)

    chi2 = np.full(m, np.nan)
    ok = V > MIN_VARIANCE
    chi2[ok] = (U[ok] ** 2) / V[ok]
    pvals = np.full(m, np.nan)
    pvals[ok] = stats.chi2.sf(chi2[ok], df=1)

    return {
        "chi2": chi2,
        "p": pvals,
        "score": U,
        "variance": V,
        "usable": ok,
        "n_used": n,
        "method": (
            "Score test against a single null fit of the covariate model. "
            "Asymptotically equivalent to the likelihood-ratio test under the "
            "null; top hits are refitted exactly to obtain effect sizes."),
    }


def genomic_inflation(pvalues: np.ndarray) -> Optional[float]:
    """Lambda_GC: observed median chi2 over its expected value.

    Around 1.0 means the null is behaving. Materially above means residual
    structure — relatedness or ancestry not fully adjusted for — and the scan
    should not be read as it stands.
    """
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p) & (p > 0)]
    if p.size < 10:
        return None
    observed = float(np.median(stats.chi2.isf(p, df=1)))
    return observed / float(stats.chi2.ppf(0.5, df=1))
