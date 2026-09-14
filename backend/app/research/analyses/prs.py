"""Polygenic risk score — Part II §4.5.

    "Mandatory: ancestry-stratified performance reporting. A score derived in
     European cohorts and applied to South Asian subjects will be miscalibrated,
     and the platform must show that rather than let it pass silently."

Mandatory is implemented as: the stratified table is always computed and always
returned. There is no flag that turns it off, because the failure mode is a
score that looks fine in aggregate and is badly calibrated in the subgroup the
tenant actually serves.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
from scipy import stats

from ..jobs import register
from ..stats import glm
from ..types import MISSING, GenotypeMatrix, as_float


def score_samples(gm: GenotypeMatrix, weights: Dict[str, float],
                  effect_alleles: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    """Weighted allele-dosage sum.

    Allele orientation is the classic silent failure: if the score's effect
    allele is the REF in our data, the dosage must be flipped (2 - d) or the
    score comes out anti-correlated with risk while looking perfectly
    plausible. Variants whose alleles cannot be matched are DROPPED and counted,
    never guessed.
    """
    effect_alleles = effect_alleles or {}
    n = gm.n_samples
    score = np.zeros(n, dtype=float)
    n_used = 0
    n_flipped = 0
    missing_per_sample = np.zeros(n, dtype=float)
    unmatched: List[str] = []

    index = {}
    for i, v in enumerate(gm.variants):
        index.setdefault(v.key, i)
        index.setdefault("{}:{}".format(v.chrom, v.pos), i)
        if v.vid:
            index.setdefault(v.vid, i)

    for key, w in weights.items():
        i = index.get(key)
        if i is None:
            unmatched.append(key)
            continue
        v = gm.variants[i]
        dose = gm.dosages[i, :].astype(float)
        observed = gm.dosages[i, :] != MISSING

        ea = effect_alleles.get(key)
        if ea is not None:
            if ea == v.alt:
                pass
            elif ea == v.ref:
                dose = 2.0 - dose
                n_flipped += 1
            else:
                unmatched.append(key)
                continue

        # Mean-impute missing dosage, which is the standard PRS convention —
        # dropping the sample entirely would lose it from the whole score.
        af = float(np.mean(dose[observed]) / 2.0) if observed.any() else 0.0
        dose = np.where(observed, dose, 2.0 * af)
        missing_per_sample += (~observed).astype(float)

        score += w * dose
        n_used += 1

    return {
        "score": score,
        "n_variants_used": n_used,
        "n_variants_requested": len(weights),
        "n_flipped": n_flipped,
        "n_unmatched": len(unmatched),
        "unmatched_examples": unmatched[:20],
        "mean_missing_variants_per_sample": float(np.mean(missing_per_sample))
        if n else 0.0,
    }


def standardise(score: np.ndarray) -> np.ndarray:
    sd = float(np.std(score))
    return (score - float(np.mean(score))) / sd if sd > 0 else score * 0.0


def percentiles(score: np.ndarray) -> np.ndarray:
    order = stats.rankdata(score, method="average")
    return 100.0 * order / (len(score) + 1)


def auc(y: np.ndarray, score: np.ndarray) -> Optional[float]:
    """AUC via the Mann-Whitney U identity — no sklearn needed."""
    finite = np.isfinite(y) & np.isfinite(score)
    yy, ss = y[finite], score[finite]
    pos, neg = ss[yy == 1], ss[yy == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(u / (len(pos) * len(neg)))


def calibration(y: np.ndarray, z: np.ndarray, kind: str) -> Dict[str, Any]:
    """Performance metrics per §4.5: AUC, OR per SD, PPV/NPV in the top decile."""
    out: Dict[str, Any] = {}
    finite = np.isfinite(y) & np.isfinite(z)
    yy, zz = y[finite], z[finite]
    out["n"] = int(finite.sum())
    if out["n"] < 20:
        out["note"] = "too few samples for stable calibration metrics"
        return out

    if kind == "binary":
        out["auc"] = auc(yy, zz)
        fit = glm.fit_logistic(zz.reshape(-1, 1), yy, 0, ["prs"])
        out["or_per_sd"] = fit.effect
        out["or_per_sd_ci"] = [fit.effect_ci_low, fit.effect_ci_high]
        out["p"] = fit.pvalue
        out["prevalence"] = float(np.mean(yy))

        pct = percentiles(zz)
        top = pct >= 90
        out["top_decile"] = {
            "n": int(top.sum()),
            "cases": int(np.sum(yy[top] == 1)),
            "ppv": float(np.mean(yy[top])) if top.any() else None,
            "npv": float(1 - np.mean(yy[~top])) if (~top).any() else None,
            "or_vs_rest": _or_2x2(yy, top),
        }
    else:
        fit = glm.fit_linear(zz.reshape(-1, 1), yy, 0, ["prs"])
        out["beta_per_sd"] = fit.beta
        out["beta_per_sd_ci"] = [fit.ci_low, fit.ci_high]
        out["p"] = fit.pvalue
        out["r2"] = fit.extra.get("r_squared")
    return out


def _or_2x2(y: np.ndarray, top: np.ndarray) -> Optional[float]:
    a = float(np.sum((y == 1) & top))
    b = float(np.sum((y == 0) & top))
    c = float(np.sum((y == 1) & ~top))
    d = float(np.sum((y == 0) & ~top))
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5   # Haldane correction
    return round((a * d) / (b * c), 4) if b * c else None


@register("prs")
def prs_job(spec: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    from .. import store as research_store

    log = context.get("log", lambda m: None)
    dataset_id = context["dataset_id"]
    gm = research_store.load_genotypes(dataset_id)
    ph = research_store.load_phenotypes(dataset_id)

    weights = {str(k): float(v) for k, v in (spec.get("weights") or {}).items()}
    if not weights:
        raise ValueError("A PRS needs score weights (PGS Catalog import or upload).")
    effect_alleles = spec.get("effect_alleles") or {}

    scored = score_samples(gm, weights, effect_alleles)
    log("scored {} of {} weighted variants ({} allele-flipped)".format(
        scored["n_variants_used"], scored["n_variants_requested"], scored["n_flipped"]))

    if scored["n_variants_used"] == 0:
        raise ValueError(
            "None of the {} score variants matched this dataset. Check the genome "
            "build and variant identifiers.".format(len(weights)))

    raw = scored["score"]
    z = standardise(raw)
    pct = percentiles(raw)

    result: Dict[str, Any] = {
        "analysis": "prs",
        "score_name": spec.get("score_name", "custom score"),
        "score_source": spec.get("score_source", "uploaded weights"),
        "n_variants_used": scored["n_variants_used"],
        "n_variants_requested": scored["n_variants_requested"],
        "n_unmatched": scored["n_unmatched"],
        "unmatched_examples": scored["unmatched_examples"],
        "n_flipped": scored["n_flipped"],
        "coverage_note": (
            "{} of {} score variants were found in this dataset ({:.1%}). A score "
            "computed on a subset of its variants is not the published score and "
            "its calibration will differ.".format(
                scored["n_variants_used"], scored["n_variants_requested"],
                scored["n_variants_used"] / max(1, scored["n_variants_requested"]))),
        "distribution": {
            "mean": float(np.mean(raw)), "sd": float(np.std(raw)),
            "min": float(np.min(raw)), "max": float(np.max(raw)),
            "histogram": _histogram(z),
        },
        "samples": [{"sample_id": s, "score": float(r), "z": float(zz),
                     "percentile": round(float(p), 2)}
                    for s, r, zz, p in zip(gm.sample_ids, raw, z, pct)][:2000],
    }

    outcome = spec.get("outcome")
    if outcome and ph is not None and outcome in ph.columns:
        pheno_index = {s: i for i, s in enumerate(ph.sample_ids)}
        keep = [i for i, s in enumerate(gm.sample_ids) if s in pheno_index]
        rows = [pheno_index[gm.sample_ids[i]] for i in keep]
        y = as_float(ph.get(outcome))[rows]
        zk = z[keep]
        kind = ph.kind(outcome)
        result["outcome"] = outcome
        result["calibration"] = calibration(y, zk, kind)

        # §4.5 mandatory: never optional, never suppressed.
        result["ancestry_stratified"] = _stratified(
            y, zk, kind, spec, ph, rows, gm, keep)
        result["ancestry_note"] = (
            "Polygenic scores are calibrated in the ancestry they were derived in. "
            "Applying a European-derived score to another population produces "
            "systematically shifted absolute risk even when the ranking is "
            "preserved. Per-group performance is shown above; if a group is "
            "absent, the score has not been validated in it here.")

    result["diagnostics"] = {
        "variant_coverage": scored["n_variants_used"] / max(1, scored["n_variants_requested"]),
        "mean_missing_variants_per_sample": scored["mean_missing_variants_per_sample"],
        "n_flipped": scored["n_flipped"],
    }
    return result


def _stratified(y, z, kind, spec, ph, rows, gm, keep) -> Dict[str, Any]:
    """Per-ancestry performance. Absence of a group is itself the finding."""
    group_col = spec.get("ancestry_column", "ancestry")
    if ph is None or group_col not in ph.columns:
        return {"available": False,
                "reason": ("No ancestry column ('{}') is present, so stratified "
                           "performance cannot be reported. The overall metrics "
                           "may not transfer across ancestries.".format(group_col))}

    labels = np.asarray(ph.columns[group_col])[rows]
    out: List[Dict[str, Any]] = []
    for g in sorted(set(str(x) for x in labels)):
        m = np.array([str(x) == g for x in labels])
        if m.sum() < 20:
            out.append({"group": g, "n": int(m.sum()),
                        "note": "too few samples for stable metrics"})
            continue
        metrics = calibration(y[m], z[m], kind)
        metrics["group"] = g
        metrics["mean_z"] = float(np.mean(z[m]))
        out.append(metrics)
    return {"available": True, "groups": out,
            "shift_note": ("Differences in mean standardised score between groups "
                           "indicate the score is not centred in those populations; "
                           "absolute risk from such a score will be biased.")}


def _histogram(z: np.ndarray, bins: int = 30) -> List[Dict[str, Any]]:
    counts, edges = np.histogram(z[np.isfinite(z)], bins=bins)
    return [{"lo": round(float(edges[i]), 3), "hi": round(float(edges[i + 1]), 3),
             "count": int(counts[i])} for i in range(len(counts))]
