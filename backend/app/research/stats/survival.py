"""Time-to-event analysis (Part II §4.7).

§4.7: "Kaplan-Meier, Cox with diagnostics for the proportional-hazards
assumption, competing risks, cumulative incidence by carrier status with
confidence bands and numbers at risk."

Every clause there is load-bearing, so each is implemented rather than
approximated:

*   **Numbers at risk** are returned alongside every curve, at every time the
    curve is evaluated. A KM curve without them is unreadable: the tail of a
    survival curve is drawn from a handful of subjects and looks exactly like
    the head, which is drawn from hundreds. Readers cannot tell a 20% survival
    estimate based on 200 people from one based on 3 unless you show them.

*   **Log-log confidence bands**, not plain Wald bands on S(t). A symmetric
    Wald interval on a probability runs outside [0, 1] near the ends of the
    curve — exactly where the uncertainty is largest — and has poor coverage
    there even when it does not. The log(-log S) transform is unbounded, so
    back-transforming always lands inside [0, 1] and coverage holds far better
    in small risk sets.

*   **Efron's method for tied event times**, not Breslow. Breslow's
    approximation treats d simultaneous events as if each faced the full risk
    set, which systematically shrinks |beta| toward zero when ties are common.
    Clinical follow-up recorded at year granularity is nothing *but* ties — a
    cohort with ages of onset in whole years will have dozens of subjects
    sharing each event time — so the "approximation" would be biasing the
    primary estimate, not smoothing a rare edge case. Efron's is available via
    ``ties="breslow"`` for comparison only.

*   **Competing risks** via the Aalen-Johansen cumulative incidence function.
    1 - KM computed on cause-specific events is *not* the cumulative incidence
    when a competing event can remove a subject from risk; it overstates the
    absolute risk, sometimes grossly. The CIF weights each cause-specific
    hazard increment by the all-cause survival still standing at that moment.

A §4.7 guardrail this module cannot enforce on its own: penetrance estimated
from an ascertained clinical cohort is biased upward. That detection lives with
the data profile; these functions compute what they are asked to compute, and
the caller is responsible for labelling the output.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import stats as sps

__all__ = [
    "kaplan_meier",
    "logrank_test",
    "cox_ph",
    "ph_assumption_test",
    "cumulative_incidence",
]

_Z = 1.959963984540054  # two-sided 95%


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _clean(time: Any, event: Any) -> Tuple[np.ndarray, np.ndarray]:
    t = np.asarray(time, dtype=float).ravel()
    e = np.asarray(event, dtype=float).ravel()
    if t.shape != e.shape:
        raise ValueError("time and event must be the same length")
    ok = np.isfinite(t) & np.isfinite(e)
    if not ok.all():
        t, e = t[ok], e[ok]
    if np.any(t < 0):
        raise ValueError("negative follow-up times")
    return t, (e != 0).astype(int)


def _risk_table(t: np.ndarray, e: np.ndarray) -> Dict[str, np.ndarray]:
    """Per distinct observed time: at risk, events, censored.

    Distinct *observed* times, not just event times: censoring times appear as
    rows with 0 events so the numbers-at-risk column is complete. The survival
    estimate is flat across those rows, which is correct.
    """
    times = np.unique(t)
    n_risk = np.array([int(np.sum(t >= u)) for u in times], dtype=int)
    n_event = np.array([int(np.sum((t == u) & (e == 1))) for u in times], dtype=int)
    n_cens = np.array([int(np.sum((t == u) & (e == 0))) for u in times], dtype=int)
    return {"time": times, "n_risk": n_risk, "n_event": n_event, "n_censored": n_cens}


def _km_one(t: np.ndarray, e: np.ndarray, alpha: float) -> Dict[str, Any]:
    tab = _risk_table(t, e)
    times, n, d, c = tab["time"], tab["n_risk"], tab["n_event"], tab["n_censored"]
    z = float(sps.norm.isf(alpha / 2.0))

    frac = np.where(n > 0, d / np.maximum(n, 1), 0.0)
    surv = np.cumprod(1.0 - frac)

    # Greenwood: Var(S) = S^2 * sum d / (n (n-d)). The n == d term is infinite
    # (the curve hits zero); carry it as inf so the band opens up rather than
    # pretending the last point is precise.
    with np.errstate(divide="ignore", invalid="ignore"):
        inc = np.where((n > 0) & (n > d), d / (n * (n - d)), np.where(d > 0, np.inf, 0.0))
    cum = np.cumsum(inc)
    # S = 0 with cum = inf is 0 * inf; the variance there is genuinely
    # undefined, so carry it as inf rather than letting NaN leak out.
    with np.errstate(invalid="ignore"):
        var_s = np.where(np.isfinite(cum), (surv ** 2) * cum, np.inf)
    std_err = np.sqrt(np.where(np.isnan(var_s), np.inf, var_s))

    # log-log band: CI on log(-log S), back-transformed. Var(log(-log S)) =
    # cum / (log S)^2.
    lo = np.empty_like(surv)
    hi = np.empty_like(surv)
    for i in range(len(surv)):
        s = surv[i]
        if not np.isfinite(cum[i]) or s <= 0.0 or s >= 1.0:
            lo[i] = 0.0 if s <= 0.0 else (1.0 if s >= 1.0 else np.nan)
            hi[i] = 0.0 if s <= 0.0 else (1.0 if s >= 1.0 else np.nan)
            continue
        se_ll = math.sqrt(cum[i]) / abs(math.log(s))
        lo[i] = s ** math.exp(z * se_ll)
        hi[i] = s ** math.exp(-z * se_ll)

    median = _median_from_curve(times, surv)
    median_ci = _median_ci(times, lo, hi)

    return {
        "time": times.tolist(),
        "survival": surv.tolist(),
        "std_err": [float(x) for x in std_err],
        "ci_lower": [float(x) for x in lo],
        "ci_upper": [float(x) for x in hi],
        "n_risk": n.tolist(),
        "n_event": d.tolist(),
        "n_censored": c.tolist(),
        "n": int(len(t)),
        "n_events": int(e.sum()),
        "median": median,
        "median_ci": median_ci,
    }


def _median_from_curve(times: np.ndarray, surv: np.ndarray) -> Optional[float]:
    idx = np.nonzero(surv <= 0.5)[0]
    if idx.size == 0:
        return None          # curve never reaches 0.5 — median not estimable
    return float(times[idx[0]])


def _median_ci(times: np.ndarray, lo: np.ndarray, hi: np.ndarray
               ) -> Optional[List[Optional[float]]]:
    """Brookmeyer-Crowley style: times where the band still straddles 0.5."""
    def first_at_or_below(arr: np.ndarray) -> Optional[float]:
        a = np.where(np.isnan(arr), 1.0, arr)
        idx = np.nonzero(a <= 0.5)[0]
        return float(times[idx[0]]) if idx.size else None
    return [first_at_or_below(hi), first_at_or_below(lo)]


def _group_labels(groups: Any, n: int) -> np.ndarray:
    g = np.asarray(groups).ravel()
    if len(g) != n:
        raise ValueError("groups has {} entries, expected {}".format(len(g), n))
    return g


# --------------------------------------------------------------------------- #
# Kaplan-Meier
# --------------------------------------------------------------------------- #
def kaplan_meier(
    time: Any,
    event: Any,
    groups: Optional[Any] = None,
    alpha: float = 0.05,
) -> Dict[str, Any]:
    """Kaplan-Meier estimator with Greenwood variance and log-log bands.

    `event` is 1 for the event, 0 for right-censored. `groups` (optional) is a
    per-subject label; one curve is returned per level.

    Returned per curve: `time`, `survival`, `std_err` (Greenwood), `ci_lower`
    / `ci_upper` (log-log transform — see the module docstring for why not
    Wald), `n_risk` / `n_event` / `n_censored` **at every time point** as §4.7
    requires, `median` and `median_ci`.

    Rows exist for censoring-only times too (with `n_event` 0 and a flat
    survival estimate) so that the numbers-at-risk column is complete.

    With `groups=None` the single curve's fields also appear at the top level
    for convenience, and under `groups["all"]`.
    """
    t, e = _clean(time, event)
    out: Dict[str, Any] = {"alpha": float(alpha)}

    if groups is None:
        curve = _km_one(t, e, alpha)
        out["groups"] = {"all": curve}
        out.update(curve)
        out["group_names"] = ["all"]
        return out

    g = _group_labels(groups, len(np.asarray(time, dtype=float).ravel()))
    ok = np.isfinite(np.asarray(time, dtype=float).ravel()) & \
        np.isfinite(np.asarray(event, dtype=float).ravel())
    g = g[ok]
    names = [x for x in _sorted_levels(g)]
    curves = {}
    for name in names:
        m = g == name
        curves[str(name)] = _km_one(t[m], e[m], alpha)
    out["groups"] = curves
    out["group_names"] = [str(x) for x in names]
    return out


def _sorted_levels(g: np.ndarray) -> List[Any]:
    lv = list(np.unique(g))
    try:
        return sorted(lv)
    except TypeError:
        return lv


# --------------------------------------------------------------------------- #
# Log-rank
# --------------------------------------------------------------------------- #
def logrank_test(time: Any, event: Any, groups: Any) -> Dict[str, Any]:
    """Mantel-Haenszel log-rank test, k groups, chi-square with k-1 df.

    At each distinct event time the test compares the observed events in each
    group with the number expected if the groups shared one hazard, pooling
    those 2xk tables across times. The hypergeometric covariance

        V[g,h] = sum_t  d_t (n_t - d_t) / (n_t - 1) * (n_tg/n_t)(delta_gh - n_th/n_t)

    handles tied event times exactly, and the quadratic form on any k-1 of the
    (O - E) components gives the k-group test. Times with n_t == 1 contribute
    nothing (no variance) and are skipped.
    """
    t_raw = np.asarray(time, dtype=float).ravel()
    e_raw = np.asarray(event, dtype=float).ravel()
    g_raw = _group_labels(groups, len(t_raw))
    ok = np.isfinite(t_raw) & np.isfinite(e_raw)
    t, e, g = t_raw[ok], (e_raw[ok] != 0).astype(int), g_raw[ok]

    names = _sorted_levels(g)
    k = len(names)
    if k < 2:
        raise ValueError("logrank_test needs at least 2 groups")

    masks = [g == nm for nm in names]
    obs = np.array([float(np.sum(e[m] == 1)) for m in masks])
    exp = np.zeros(k)
    V = np.zeros((k, k))

    for u in np.unique(t[e == 1]):
        n_t = float(np.sum(t >= u))
        d_t = float(np.sum((t == u) & (e == 1)))
        if n_t <= 1.0 or d_t <= 0.0:
            continue
        n_tg = np.array([float(np.sum(t[m] >= u)) for m in masks])
        exp += n_tg * (d_t / n_t)
        f = d_t * (n_t - d_t) / (n_t - 1.0) / (n_t * n_t)
        V += f * (np.diag(n_tg * n_t) - np.outer(n_tg, n_tg))

    u_vec = (obs - exp)[: k - 1]
    V_sub = V[: k - 1, : k - 1]
    try:
        chi2 = float(u_vec @ np.linalg.solve(V_sub, u_vec))
    except np.linalg.LinAlgError:
        chi2 = float(u_vec @ np.linalg.pinv(V_sub) @ u_vec)
    chi2 = max(chi2, 0.0)
    df = k - 1
    p = float(sps.chi2.sf(chi2, df))

    return {
        "chi2": chi2,
        "df": int(df),
        "p": p,
        "groups": [str(nm) for nm in names],
        "observed": obs.tolist(),
        "expected": exp.tolist(),
        "n": [int(m.sum()) for m in masks],
        "variance": V.tolist(),
        "test": "logrank_mantel_haenszel",
    }


# --------------------------------------------------------------------------- #
# Cox proportional hazards
# --------------------------------------------------------------------------- #
def _cox_loglik(
    t: np.ndarray, e: np.ndarray, X: np.ndarray, beta: np.ndarray, ties: str
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Partial log-likelihood, score and observed information.

    Data must arrive sorted by time ascending. Risk-set sums are reverse
    cumulative sums, so the whole thing is O(n p^2) per iteration rather than
    O(n^2 p^2).
    """
    n, p = X.shape
    eta = X @ beta
    eta = eta - eta.max() if n else eta          # keep exp() in range
    w = np.exp(np.clip(eta, -500.0, 500.0))

    wx = w[:, None] * X
    wxx = wx[:, :, None] * X[:, None, :]
    S0 = np.cumsum(w[::-1])[::-1]
    S1 = np.cumsum(wx[::-1], axis=0)[::-1]
    S2 = np.cumsum(wxx[::-1], axis=0)[::-1]

    loglik = 0.0
    score = np.zeros(p)
    info = np.zeros((p, p))

    ev_times = np.unique(t[e == 1])
    for u in ev_times:
        first = int(np.searchsorted(t, u, side="left"))
        D = np.nonzero((t == u) & (e == 1))[0]
        d = len(D)
        S0R, S1R, S2R = S0[first], S1[first], S2[first]
        S0D = float(np.sum(w[D]))
        S1D = np.sum(wx[D], axis=0)
        S2D = np.sum(wxx[D], axis=0)

        loglik += float(np.sum(eta[D]))
        score += np.sum(X[D], axis=0)

        if ties == "breslow":
            # Breslow: every tied event faces the full risk set. Cruder, kept
            # only so the Efron/Breslow difference can be demonstrated.
            steps = [(0.0, d)]
        else:
            steps = [(l / d, 1.0) for l in range(d)]

        for frac, mult in steps:
            phi = S0R - frac * S0D
            if phi <= 0:
                phi = 1e-300
            a = (S1R - frac * S1D) / phi
            b = (S2R - frac * S2D) / phi
            loglik -= mult * math.log(phi)
            score -= mult * a
            info += mult * (b - np.outer(a, a))

    return loglik, score, info


def cox_ph(
    time: Any,
    event: Any,
    X: Any,
    names: Optional[Sequence[str]] = None,
    ties: str = "efron",
    max_iter: int = 50,
    tol: float = 1e-9,
) -> Dict[str, Any]:
    """Cox proportional-hazards model fitted by Newton-Raphson.

    **Ties are handled by Efron's method** (the default). Breslow's
    approximation is available as ``ties="breslow"`` for comparison but should
    not be used for reporting: it pretends each of d simultaneous events faced
    the undiminished risk set, which biases |beta| toward zero, and clinical
    follow-up recorded to the nearest year produces ties in bulk. Efron
    averages the risk set over the d! orderings the data cannot distinguish,
    which is essentially exact for the tie sizes seen in practice and costs
    nothing here.

    `X` is n x p (a 1-D vector is treated as a single covariate). No intercept:
    the baseline hazard is left unspecified, which is the point of the model.

    Returns per covariate beta, se, hazard_ratio, ci_lower/ci_upper (95% on the
    HR scale), z and Wald p; plus `loglik` (the maximised log partial
    likelihood), `loglik_null`, `converged` and `n_iter`.
    """
    if ties not in ("efron", "breslow"):
        raise ValueError("ties must be 'efron' or 'breslow'")

    t_raw = np.asarray(time, dtype=float).ravel()
    e_raw = np.asarray(event, dtype=float).ravel()
    Xa = np.asarray(X, dtype=float)
    if Xa.ndim == 1:
        Xa = Xa[:, None]
    if Xa.shape[0] != len(t_raw):
        raise ValueError("X has {} rows, time has {}".format(Xa.shape[0], len(t_raw)))

    ok = np.isfinite(t_raw) & np.isfinite(e_raw) & np.all(np.isfinite(Xa), axis=1)
    t_raw, e_raw, Xa = t_raw[ok], e_raw[ok], Xa[ok]
    e_raw = (e_raw != 0).astype(int)

    order = np.argsort(t_raw, kind="mergesort")
    t, e, Xs = t_raw[order], e_raw[order], Xa[order]
    n, p = Xs.shape
    if names is None:
        names = ["x{}".format(j) for j in range(p)]
    names = [str(s) for s in names]
    if len(names) != p:
        raise ValueError("names has {} entries, X has {} columns".format(len(names), p))

    beta = np.zeros(p)
    loglik, score, info = _cox_loglik(t, e, Xs, beta, ties)
    loglik_null = loglik
    converged = False
    it = 0

    for it in range(1, max_iter + 1):
        try:
            step = np.linalg.solve(info, score)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(info) @ score
        new_beta = beta + step
        new_ll, new_score, new_info = _cox_loglik(t, e, Xs, new_beta, ties)
        # Step-halving: Newton on the partial likelihood can overshoot badly
        # when a covariate nearly separates the event ordering.
        halvings = 0
        while (not np.isfinite(new_ll) or new_ll < loglik - 1e-9) and halvings < 25:
            step = step / 2.0
            new_beta = beta + step
            new_ll, new_score, new_info = _cox_loglik(t, e, Xs, new_beta, ties)
            halvings += 1
        delta = abs(new_ll - loglik)
        beta, loglik, score, info = new_beta, new_ll, new_score, new_info
        if delta < tol * (abs(loglik) + tol):
            converged = True
            break

    try:
        cov = np.linalg.inv(info)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(info)
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(se > 0, beta / se, np.nan)
    pv = 2.0 * sps.norm.sf(np.abs(z))

    covariates = []
    for j in range(p):
        covariates.append({
            "name": names[j],
            "beta": float(beta[j]),
            "se": float(se[j]),
            "hazard_ratio": float(math.exp(beta[j])),
            "ci_lower": float(math.exp(beta[j] - _Z * se[j])),
            "ci_upper": float(math.exp(beta[j] + _Z * se[j])),
            "z": float(z[j]),
            "p": float(pv[j]),
        })

    lr_stat = float(2.0 * (loglik - loglik_null))
    return {
        "covariates": covariates,
        "names": names,
        "beta": beta.tolist(),
        "se": se.tolist(),
        "hazard_ratio": [float(math.exp(b)) for b in beta],
        "p": [float(x) for x in pv],
        "cov": cov.tolist(),
        "loglik": float(loglik),
        "loglik_null": float(loglik_null),
        "lr_chi2": lr_stat,
        "lr_p": float(sps.chi2.sf(max(lr_stat, 0.0), p)),
        "converged": bool(converged),
        "n_iter": int(it),
        "ties": ties,
        "n": int(n),
        "n_events": int(e.sum()),
    }


# --------------------------------------------------------------------------- #
# Proportional-hazards diagnostics
# --------------------------------------------------------------------------- #
def ph_assumption_test(
    time: Any,
    event: Any,
    X: Any,
    beta: Any,
    names: Optional[Sequence[str]] = None,
    transform: str = "rank",
) -> Dict[str, Any]:
    """Schoenfeld-residual test of proportional hazards (Grambsch-Therneau).

    §4.7 requires "Cox with diagnostics for the proportional-hazards
    assumption". The assumption is that each covariate's log hazard ratio is
    constant over follow-up; the diagnostic is that the Schoenfeld residual at
    event time t,

        s_k = x_(k) - xbar(t_k),   xbar(t) = sum_R w_j x_j / sum_R w_j

    has expectation zero at every event time under PH, but drifts with time
    when the effect is time-varying. Scaling them by d * I(beta)^-1 (Grambsch &
    Therneau 1994) turns the drift into an estimate of beta(t) - beta, so a
    non-zero slope against time is the violation.

    Per covariate this reports `rho` (correlation of the scaled residual with
    g(t)) and a 1-df chi-square

        T = sum_k (g_k - gbar) * s*_k,   Var(T) = d * I^-1_pp * sum_k (g_k-gbar)^2

    plus a p-df global test on the same quantities. A small p means the PH
    assumption is rejected for that covariate and its hazard ratio is an
    average over follow-up rather than a constant.

    `transform` picks g(t): "rank" (default, robust to the time scale and the
    heavy right tail of follow-up), "identity", or "log".

    Ties: all events at the same time share one risk-set mean (the Breslow
    convention for residuals). With heavy ties this is mildly conservative; it
    is the standard implementation.
    """
    t_raw = np.asarray(time, dtype=float).ravel()
    e_raw = np.asarray(event, dtype=float).ravel()
    Xa = np.asarray(X, dtype=float)
    if Xa.ndim == 1:
        Xa = Xa[:, None]
    b = np.asarray(beta, dtype=float).ravel()

    ok = np.isfinite(t_raw) & np.isfinite(e_raw) & np.all(np.isfinite(Xa), axis=1)
    t_raw, e_raw, Xa = t_raw[ok], e_raw[ok], Xa[ok]
    e_raw = (e_raw != 0).astype(int)

    order = np.argsort(t_raw, kind="mergesort")
    t, e, Xs = t_raw[order], e_raw[order], Xa[order]
    n, p = Xs.shape
    if len(b) != p:
        raise ValueError("beta has {} entries, X has {} columns".format(len(b), p))
    if names is None:
        names = ["x{}".format(j) for j in range(p)]
    names = [str(s) for s in names]

    eta = Xs @ b
    eta = eta - eta.max()
    w = np.exp(np.clip(eta, -500.0, 500.0))
    wx = w[:, None] * Xs
    wxx = wx[:, :, None] * Xs[:, None, :]
    S0 = np.cumsum(w[::-1])[::-1]
    S1 = np.cumsum(wx[::-1], axis=0)[::-1]
    S2 = np.cumsum(wxx[::-1], axis=0)[::-1]

    resid: List[np.ndarray] = []
    ev_time: List[float] = []
    info = np.zeros((p, p))

    for u in np.unique(t[e == 1]):
        first = int(np.searchsorted(t, u, side="left"))
        S0R = S0[first]
        xbar = S1[first] / S0R
        V_t = S2[first] / S0R - np.outer(xbar, xbar)
        D = np.nonzero((t == u) & (e == 1))[0]
        for i in D:
            resid.append(Xs[i] - xbar)
            ev_time.append(float(u))
            info += V_t

    d = len(resid)
    if d < 3:
        return {"status": "insufficient_events", "p_global": None, "n_events": d,
                "covariates": [], "reason": "need at least 3 events for a PH test"}

    R = np.asarray(resid)                     # d x p
    tt = np.asarray(ev_time, dtype=float)

    if transform == "rank":
        g = sps.rankdata(tt)
    elif transform == "identity":
        g = tt.copy()
    elif transform == "log":
        g = np.log(np.maximum(tt, np.finfo(float).tiny))
    else:
        raise ValueError("transform must be 'rank', 'identity' or 'log'")
    g = g - g.mean()
    sg2 = float(np.sum(g ** 2))

    try:
        inv_info = np.linalg.inv(info)
    except np.linalg.LinAlgError:
        inv_info = np.linalg.pinv(info)

    # Scaled Schoenfeld residuals: an estimate of beta(t_k).
    scaled = b[None, :] + d * (R @ inv_info.T)

    T = d * (inv_info @ (R.T @ g))            # length p
    var_diag = d * np.diag(inv_info) * sg2

    covariates = []
    for j in range(p):
        vj = float(var_diag[j])
        if vj <= 0 or not np.isfinite(vj) or sg2 <= 0:
            covariates.append({"name": names[j], "rho": float("nan"),
                               "chi2": float("nan"), "p": None})
            continue
        chi2 = float(T[j] ** 2 / vj)
        sd_s = float(np.std(scaled[:, j]))
        rho = float(np.corrcoef(g, scaled[:, j])[0, 1]) if sd_s > 0 else float("nan")
        covariates.append({
            "name": names[j],
            "rho": rho,
            "chi2": chi2,
            "df": 1,
            "p": float(sps.chi2.sf(chi2, 1)),
        })

    cov_T = d * inv_info * sg2
    try:
        chi2_g = float(T @ np.linalg.solve(cov_T, T))
    except np.linalg.LinAlgError:
        chi2_g = float(T @ np.linalg.pinv(cov_T) @ T)
    chi2_g = max(chi2_g, 0.0)

    return {
        "status": "ok",
        "covariates": covariates,
        "chi2_global": chi2_g,
        "df_global": int(p),
        "p_global": float(sps.chi2.sf(chi2_g, p)),
        "n_events": int(d),
        "transform": transform,
        "scaled_schoenfeld": scaled.tolist(),
        "event_times": tt.tolist(),
        "test": "grambsch_therneau_schoenfeld",
    }


# --------------------------------------------------------------------------- #
# Competing risks
# --------------------------------------------------------------------------- #
def _cif_one(
    t: np.ndarray, e: np.ndarray, ct: np.ndarray, causes: List[Any], alpha: float
) -> Dict[str, Any]:
    """Aalen-Johansen CIF for every cause, on one group."""
    times = np.unique(t)
    k = len(times)
    z = float(sps.norm.isf(alpha / 2.0))

    n_risk = np.array([int(np.sum(t >= u)) for u in times], dtype=float)
    d_all = np.array([int(np.sum((t == u) & (e == 1))) for u in times], dtype=float)
    n_cens = np.array([int(np.sum((t == u) & (e == 0))) for u in times], dtype=int)

    # All-cause KM; S_prev[i] is S(t_(i-1)), the survival still standing when
    # the increment at t_i is taken. This is the factor that 1-KM on
    # cause-specific events gets wrong.
    surv = np.cumprod(1.0 - d_all / np.maximum(n_risk, 1.0))
    s_prev = np.concatenate([[1.0], surv[:-1]])

    with np.errstate(divide="ignore", invalid="ignore"):
        A = np.where(n_risk > d_all, d_all / (n_risk * (n_risk - d_all)), 0.0)

    out_causes: Dict[str, Any] = {}
    for cause in causes:
        d_j = np.array(
            [int(np.sum((t == u) & (e == 1) & (ct == cause))) for u in times],
            dtype=float)
        incr = s_prev * d_j / np.maximum(n_risk, 1.0)
        F = np.cumsum(incr)

        # Delta-method variance (Marubini & Valsecchi), expanded into running
        # sums so it is O(k) rather than O(k^2):
        #   Var(F_k) = sum_i (F_k-F_i)^2 A_i + sum_i B_i - 2 sum_i (F_k-F_i) C_i
        B = (s_prev ** 2) * ((n_risk - d_j) / np.maximum(n_risk, 1.0)) * \
            d_j / np.maximum(n_risk ** 2, 1.0)
        C = s_prev * d_j / np.maximum(n_risk ** 2, 1.0)
        cA = np.cumsum(A)
        cAF = np.cumsum(A * F)
        cAF2 = np.cumsum(A * F * F)
        cB = np.cumsum(B)
        cC = np.cumsum(C)
        cCF = np.cumsum(C * F)
        var = (F ** 2) * cA - 2.0 * F * cAF + cAF2 + cB - 2.0 * F * cC + 2.0 * cCF
        var = np.clip(var, 0.0, None)
        se = np.sqrt(var)

        lo = np.zeros(k)
        hi = np.zeros(k)
        for i in range(k):
            f = F[i]
            if f <= 0.0 or f >= 1.0 or se[i] <= 0.0:
                lo[i] = max(min(f - z * se[i], 1.0), 0.0)
                hi[i] = max(min(f + z * se[i], 1.0), 0.0)
                continue
            # complementary log-log, matching the KM band's logic: the
            # transform is unbounded so back-transforming stays in [0,1].
            theta = math.log(-math.log(1.0 - f))
            var_theta = var[i] / (((1.0 - f) ** 2) * (math.log(1.0 - f) ** 2))
            half = z * math.sqrt(max(var_theta, 0.0))
            # theta increases with F, so exp(+half) raises (1-F) to a power > 1
            # and therefore raises F: that branch is the UPPER limit.
            lo[i] = 1.0 - (1.0 - f) ** math.exp(-half)
            hi[i] = 1.0 - (1.0 - f) ** math.exp(half)
        out_causes[str(cause)] = {
            "cif": F.tolist(),
            "std_err": se.tolist(),
            "ci_lower": lo.tolist(),
            "ci_upper": hi.tolist(),
            "n_event": d_j.astype(int).tolist(),
        }

    return {
        "time": times.tolist(),
        "n_risk": n_risk.astype(int).tolist(),
        "n_event_any": d_all.astype(int).tolist(),
        "n_censored": n_cens.tolist(),
        "overall_survival": surv.tolist(),
        "causes": out_causes,
        "n": int(len(t)),
    }


def cumulative_incidence(
    time: Any,
    event: Any,
    event_type: Any,
    groups: Optional[Any] = None,
    cause: Optional[Any] = None,
    alpha: float = 0.05,
) -> Dict[str, Any]:
    """Aalen-Johansen cumulative incidence under competing risks (§4.7).

    `event` is 1 if any event occurred and 0 if right-censored; `event_type`
    labels which event it was (its value is ignored where `event == 0`).
    `cause` names the event of interest — the others are competing events. If
    omitted, the first cause in sorted order is used and reported in
    `cause_of_interest`.

    The estimator is

        F_j(t) = sum_{t_i <= t}  S(t_i-) * d_ji / n_i

    where S is the **all-cause** Kaplan-Meier. Using 1 - KM on the
    cause-specific events instead — i.e. treating competing events as
    censoring — assumes a subject removed by the competing event could still
    later experience the event of interest, and so overstates absolute risk.
    With a strong competing risk (an older cohort, a lethal competing
    condition) the overstatement is large enough to change clinical
    conclusions, which is why §4.7 asks for competing risks explicitly.

    Confidence bands use a complementary log-log transform of the delta-method
    variance. `n_risk` is returned at every time point.

    A CIF is returned for *every* observed cause under `causes`, keyed by cause
    label, so the competing events can be plotted alongside rather than left
    implicit; `cause_labels` lists them.
    """
    t_raw = np.asarray(time, dtype=float).ravel()
    e_raw = np.asarray(event, dtype=float).ravel()
    ct_raw = np.asarray(event_type).ravel()
    if not (len(t_raw) == len(e_raw) == len(ct_raw)):
        raise ValueError("time, event and event_type must be the same length")

    ok = np.isfinite(t_raw) & np.isfinite(e_raw)
    t, e, ct = t_raw[ok], (e_raw[ok] != 0).astype(int), ct_raw[ok]

    observed = ct[e == 1]
    causes = _sorted_levels(observed) if observed.size else []
    if cause is None:
        cause_of_interest = causes[0] if causes else None
    else:
        cause_of_interest = cause
        # Keep a requested-but-unobserved cause in the output with a flat
        # all-zero CIF, rather than silently dropping the column the caller
        # asked for.
        if not any(c == cause for c in causes):
            causes = list(causes) + [cause]

    out: Dict[str, Any] = {
        "cause_of_interest": None if cause_of_interest is None else str(cause_of_interest),
        # `cause_labels` is the list; `causes` is the per-cause result dict
        # (filled in below, and per group when grouped) -- they must not share
        # a key or the flattened single-group view overwrites the list.
        "cause_labels": [str(c) for c in causes],
        "alpha": float(alpha),
        "estimator": "aalen_johansen",
    }

    if groups is None:
        res = _cif_one(t, e, ct, causes, alpha)
        out["groups"] = {"all": res}
        out.update(res)
        out["group_names"] = ["all"]
        return out

    g = _group_labels(groups, len(t_raw))[ok]
    names = _sorted_levels(g)
    out["groups"] = {
        str(nm): _cif_one(t[g == nm], e[g == nm], ct[g == nm], causes, alpha)
        for nm in names
    }
    out["group_names"] = [str(nm) for nm in names]
    return out
