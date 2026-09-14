"""Power calculation — Part II §5.1.

    "Every analysis shows estimated power before it runs."

    n = 4,200 · cases = 180 · MAF = 0.008 · assumed OR = 2.0 · α = 5×10⁻⁸
    Estimated power: 11%

    This analysis is unlikely to detect an effect of this size.

This is the guardrail that does the most work in practice. A researcher who
sees 11% before running is far more likely to reconsider than one who sees a
null result afterwards and concludes the gene is not involved — which is the
same number read backwards, and the wrong conclusion.

Power is also better evidence than the raw-n thresholds in the capability
matrix (§3.2 says so explicitly), because n alone says nothing about case
balance or exposure frequency.

All formulas here are the standard normal/non-central-chi-square approximations.
They are accurate enough to drive a decision and are labelled as estimates
everywhere they surface.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional

import numpy as np
from scipy import stats

GENOME_WIDE_ALPHA = 5e-8
EXOME_WIDE_ALPHA = 2.5e-6
CANDIDATE_ALPHA = 0.05


def _z(alpha: float) -> float:
    return float(stats.norm.ppf(1.0 - alpha / 2.0))


# ------------------------------------------------------------ case-control ---
def power_case_control(n_cases: int, n_controls: int, maf: float,
                       odds_ratio: float, alpha: float = GENOME_WIDE_ALPHA,
                       model: str = "additive") -> Dict[str, Any]:
    """Power for an allelic/additive case-control association test.

    Works on the log-odds scale: derive the case allele frequency implied by the
    odds ratio, then compare the expected log-OR against its standard error.
    """
    if n_cases <= 0 or n_controls <= 0 or not (0 < maf < 1) or odds_ratio <= 0:
        return _empty(alpha, "inputs out of range")

    p0 = float(maf)
    orr = float(odds_ratio)
    # Frequency in cases implied by the OR, from p1/(1-p1) = OR * p0/(1-p0).
    p1 = (orr * p0) / (1.0 - p0 + orr * p0)

    # Allele counts: each subject contributes 2 alleles under an additive model.
    m = 2 if model == "additive" else 1
    var = (1.0 / (m * n_cases * p1 * (1 - p1))
           + 1.0 / (m * n_controls * p0 * (1 - p0)))
    if var <= 0:
        return _empty(alpha, "degenerate variance")

    ncp = abs(math.log(orr)) / math.sqrt(var)
    crit = _z(alpha)
    power = float(stats.norm.cdf(ncp - crit) + stats.norm.cdf(-ncp - crit))

    return {
        "power": round(power, 4),
        "power_pct": round(100 * power, 1),
        "alpha": alpha,
        "n": n_cases + n_controls,
        "n_cases": n_cases,
        "n_controls": n_controls,
        "maf": p0,
        "odds_ratio": orr,
        "expected_carriers_cases": int(round(2 * n_cases * p1)),
        "expected_carriers_controls": int(round(2 * n_controls * p0)),
        "model": model,
        "assumptions": ("Additive allelic model, no covariate adjustment, "
                        "exact allele frequencies. Estimate only."),
        "verdict": verdict(power),
        "summary": ("n = {:,} · cases = {:,} · MAF = {:.4g} · assumed OR = {:.3g} "
                    "· α = {:.3g}".format(n_cases + n_controls, n_cases, p0, orr, alpha)),
    }


# ------------------------------------------------------------- quantitative --
def power_quantitative(n: int, maf: float, beta_sd: float,
                       alpha: float = GENOME_WIDE_ALPHA,
                       model: str = "additive") -> Dict[str, Any]:
    """Power for a quantitative outcome, effect expressed in outcome SD units
    per allele — the scale a researcher actually reasons in."""
    if n <= 2 or not (0 < maf < 1):
        return _empty(alpha, "inputs out of range")

    p = float(maf)
    var_g = 2 * p * (1 - p) if model == "additive" else p * (1 - p)
    r2 = (beta_sd ** 2) * var_g
    r2 = min(max(r2, 0.0), 0.999999)
    if r2 <= 0:
        return _empty(alpha, "zero effect size")

    ncp = n * r2 / (1 - r2)
    crit = stats.chi2.ppf(1 - alpha, df=1)
    power = float(stats.ncx2.sf(crit, df=1, nc=ncp))

    return {
        "power": round(power, 4),
        "power_pct": round(100 * power, 1),
        "alpha": alpha,
        "n": n,
        "maf": p,
        "beta_sd": beta_sd,
        "variance_explained": round(r2, 6),
        "expected_carriers": int(round(n * (1 - (1 - p) ** 2))),
        "model": model,
        "assumptions": ("Additive model, outcome standardised, no covariate "
                        "adjustment. Estimate only."),
        "verdict": verdict(power),
        "summary": ("n = {:,} · MAF = {:.4g} · assumed β = {:.3g} SD · α = {:.3g}"
                    .format(n, p, beta_sd, alpha)),
    }


# ------------------------------------------------- two-group / carrier tests --
def power_two_group_mean(n_exposed: int, n_unexposed: int, delta_sd: float,
                         alpha: float = 0.05) -> Dict[str, Any]:
    """Power to detect a mean difference between carriers and non-carriers —
    the shape of the spec's §7 worked example (312 carriers vs 41,868)."""
    if n_exposed <= 1 or n_unexposed <= 1:
        return _empty(alpha, "group too small")
    se = math.sqrt(1.0 / n_exposed + 1.0 / n_unexposed)
    ncp = abs(delta_sd) / se
    crit = _z(alpha)
    power = float(stats.norm.cdf(ncp - crit) + stats.norm.cdf(-ncp - crit))
    return {
        "power": round(power, 4), "power_pct": round(100 * power, 1),
        "alpha": alpha, "n": n_exposed + n_unexposed,
        "n_exposed": n_exposed, "n_unexposed": n_unexposed,
        "delta_sd": delta_sd, "verdict": verdict(power),
        "assumptions": "Equal variances, two-sided test. Estimate only.",
        "summary": ("carriers = {:,} · non-carriers = {:,} · assumed β = {:.3g} SD "
                    "· α = {:.3g}".format(n_exposed, n_unexposed, delta_sd, alpha)),
    }


# ---------------------------------------------------- minimum detectable ----
def minimum_detectable_or(n_cases: int, n_controls: int, maf: float,
                          alpha: float = GENOME_WIDE_ALPHA,
                          target_power: float = 0.8) -> Optional[float]:
    """The OR this design could detect at the target power.

    More actionable than a bare power number: it answers "what would I have to
    believe for this study to be worth running?"
    """
    lo, hi = 1.0001, 100.0
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        p = power_case_control(n_cases, n_controls, maf, mid, alpha)["power"]
        if p < target_power:
            lo = mid
        else:
            hi = mid
    return round(hi, 3) if hi < 99 else None


def minimum_detectable_beta(n: int, maf: float, alpha: float = GENOME_WIDE_ALPHA,
                            target_power: float = 0.8) -> Optional[float]:
    lo, hi = 1e-4, 10.0
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        p = power_quantitative(n, maf, mid, alpha)["power"]
        if p < target_power:
            lo = mid
        else:
            hi = mid
    return round(hi, 4) if hi < 9 else None


# ------------------------------------------------------------------ verdict --
def verdict(power: float) -> Dict[str, Any]:
    """The sentence the researcher reads, plus the options the spec lists."""
    if power >= 0.8:
        return {"level": "adequate",
                "text": "This analysis is adequately powered for the assumed effect.",
                "options": []}
    if power >= 0.5:
        return {"level": "marginal",
                "text": "This analysis is marginally powered. A null result will not "
                        "distinguish 'no effect' from 'not enough data'.",
                "options": ["aggregate to gene level",
                            "restrict to a candidate-gene analysis with relaxed α",
                            "report descriptively"]}
    return {"level": "underpowered",
            "text": "This analysis is unlikely to detect an effect of this size.",
            "options": ["relax α for a candidate-gene analysis",
                        "aggregate to gene level",
                        "increase sample size",
                        "report descriptively"]}


def _empty(alpha: float, why: str) -> Dict[str, Any]:
    return {"power": None, "power_pct": None, "alpha": alpha, "error": why,
            "verdict": {"level": "unknown",
                        "text": "Power could not be estimated: {}.".format(why),
                        "options": []},
            "summary": ""}


def suggested_alpha(n_variants: int, scope: str = "genome") -> Dict[str, Any]:
    """Multiple-testing threshold appropriate to the analysis scope (§5.6)."""
    if scope == "genome":
        return {"alpha": GENOME_WIDE_ALPHA,
                "label": "genome-wide (5×10⁻⁸)",
                "rationale": "Conventional genome-wide threshold."}
    if scope == "exome":
        return {"alpha": EXOME_WIDE_ALPHA, "label": "exome-wide (2.5×10⁻⁶)",
                "rationale": "Bonferroni over ~20,000 genes."}
    bonf = 0.05 / max(1, n_variants)
    return {"alpha": bonf,
            "label": "Bonferroni over {:,} tests ({:.3g})".format(n_variants, bonf),
            "rationale": "0.05 corrected for the number of tests actually run."}
