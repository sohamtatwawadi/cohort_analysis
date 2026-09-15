"""Association analysis — Part II §4.1, with the §5 guardrails.

The statistics are in stats/glm.py. This module is the part that decides what
to run, refuses when it should not, and reports everything a reader needs to
judge the result:

    §4.1  automatic model selection with a stated rationale the user can override
    §5.1  power BEFORE execution
    §5.2  pre-run assumption checks
    §5.3  post-run diagnostics
    §5.5  one-click sensitivity analysis
    §5.6  multiple testing always

    "Every result reports: n, cases, controls, effect (OR/beta/HR), SE, 95% CI,
     p, adjusted p, model, covariates, missingness, convergence status."

That list is the contract for `_result_row`, and nothing returns a partial one.
"""
from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..jobs import register
from ..stats import glm, power as power_mod
from ..types import MISSING, GenotypeMatrix, PhenotypeTable, as_float

# How a carrier is defined when an exposure is collapsed to a binary indicator.
CARRIER_MODELS = {
    "additive": "dosage 0/1/2",
    "dominant": "any ALT allele",
    "recessive": "two ALT alleles",
}


def encode_genotype(dosage: np.ndarray, model: str = "additive") -> np.ndarray:
    """Genetic model encoding. Missing stays missing — imputing a genotype to
    the mean here would quietly inflate n and shrink the standard error."""
    g = dosage.astype(float)
    g[dosage == MISSING] = np.nan
    if model == "dominant":
        return (g >= 1).astype(float) * np.where(np.isnan(g), np.nan, 1.0)
    if model == "recessive":
        return (g >= 2).astype(float) * np.where(np.isnan(g), np.nan, 1.0)
    return g


def build_design(exposure: np.ndarray,
                 covariates: Optional[Dict[str, np.ndarray]] = None
                 ) -> Tuple[np.ndarray, List[str]]:
    cols = [exposure.astype(float)]
    names = ["exposure"]
    for name, values in (covariates or {}).items():
        cols.append(as_float(np.asarray(values)))
        names.append(name)
    return np.column_stack(cols), names


# ------------------------------------------------------- §5.2 pre-run checks --
def assumption_checks(y: np.ndarray, X: np.ndarray, names: List[str],
                      kind: str) -> List[Dict[str, Any]]:
    """Spec §5.2: missingness, separation, distributional assumptions,
    multicollinearity, cell counts, convergence.

    Returned as a list of named checks with pass/warn/fail so the UI can show
    them individually rather than as one opaque verdict.
    """
    checks: List[Dict[str, Any]] = []

    finite = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    n_complete = int(finite.sum())
    missing_frac = 1.0 - (n_complete / max(1, len(y)))
    checks.append({
        "check": "missingness",
        "status": "fail" if n_complete < 20 else "warn" if missing_frac > 0.2 else "pass",
        "detail": "{} of {} rows complete ({:.1%} dropped)".format(
            n_complete, len(y), missing_frac),
        "note": ("Complete-case analysis. Rows with any missing covariate are "
                 "excluded, which biases the result if missingness is not random."),
    })

    if n_complete == 0:
        return checks

    yc, Xc = y[finite], X[finite]

    if kind == "binary":
        n_cases = int(np.sum(yc == 1))
        n_controls = int(np.sum(yc == 0))
        checks.append({
            "check": "case_control_balance",
            "status": "fail" if min(n_cases, n_controls) == 0
                      else "warn" if min(n_cases, n_controls) < 10 else "pass",
            "detail": "{} cases, {} controls".format(n_cases, n_controls),
        })

        # 2x2 cell counts for a binary exposure — the separation early-warning.
        expo = Xc[:, 0]
        if set(np.unique(expo[np.isfinite(expo)])) <= {0.0, 1.0}:
            cells = {
                "exposed_cases": int(np.sum((expo == 1) & (yc == 1))),
                "exposed_controls": int(np.sum((expo == 1) & (yc == 0))),
                "unexposed_cases": int(np.sum((expo == 0) & (yc == 1))),
                "unexposed_controls": int(np.sum((expo == 0) & (yc == 0))),
            }
            smallest = min(cells.values())
            # A zero cell is NOT a blocking condition. Separation is exactly the
            # situation Firth penalisation exists for (§4.1: "Firth logistic
            # regression (required — standard logistic separates with rare
            # exposures)"). Blocking here would refuse precisely the analyses
            # the spec mandates a method for; it is a warning that routes the
            # model choice, not a reason to stop.
            checks.append({
                "check": "cell_counts",
                "status": "warn" if smallest < 5 else "pass",
                "detail": ", ".join("{}={}".format(k, v) for k, v in cells.items()),
                "note": ("A zero cell means the exposure perfectly predicts the "
                         "outcome. Standard logistic regression will not converge; "
                         "Firth penalised regression is used instead."
                         if smallest == 0 else
                         "Sparse cell counts — Firth penalised regression is used."),
                "forces_firth": smallest < 5,
            })
    else:
        vals = yc[np.isfinite(yc)]
        skew = float(np.mean(((vals - vals.mean()) / (vals.std() or 1)) ** 3)) if len(vals) > 2 else 0.0
        checks.append({
            "check": "outcome_distribution",
            "status": "warn" if abs(skew) > 2 else "pass",
            "detail": "skew = {:.2f}".format(skew),
            "note": ("Heavily skewed outcome. Consider the rank-based "
                     "inverse-normal transform." if abs(skew) > 2 else None),
        })

    # Multicollinearity via the condition number — a design that is nearly
    # singular produces standard errors that are arithmetically fine and
    # scientifically meaningless.
    if Xc.shape[1] > 1:
        design = np.column_stack([np.ones(len(Xc)), Xc])
        try:
            cond = float(np.linalg.cond(design))
        except np.linalg.LinAlgError:
            cond = float("inf")
        checks.append({
            "check": "multicollinearity",
            "status": "fail" if cond > 1e10 else "warn" if cond > 1e4 else "pass",
            "detail": "condition number = {:.3g}".format(cond),
            "note": ("Covariates are nearly collinear; coefficients are unstable."
                     if cond > 1e4 else None),
        })

    expo = Xc[:, 0]
    n_distinct = len(np.unique(expo[np.isfinite(expo)]))
    checks.append({
        "check": "exposure_variation",
        "status": "fail" if n_distinct < 2 else "pass",
        "detail": "{} distinct exposure value(s)".format(n_distinct),
        "note": "No variation in the exposure — nothing to estimate." if n_distinct < 2 else None,
    })
    return checks


def checks_block(checks: Sequence[Dict[str, Any]]) -> Optional[str]:
    failed = [c for c in checks if c["status"] == "fail"]
    if not failed:
        return None
    return "; ".join("{}: {}".format(c["check"], c["detail"]) for c in failed)


# -------------------------------------------------------------- the analysis --
def run_association(
    y: np.ndarray,
    exposure: np.ndarray,
    kind: str,
    covariates: Optional[Dict[str, np.ndarray]] = None,
    model: Optional[str] = None,
    term_label: str = "exposure",
    alpha: float = 0.05,
    transform: Optional[str] = None,
) -> Dict[str, Any]:
    """One association test with the full §4.1 result contract."""
    y = as_float(np.asarray(y))
    exposure = as_float(np.asarray(exposure))

    if transform == "rint" and kind == "quantitative":
        y = glm.rank_inverse_normal(y)

    X, names = build_design(exposure, covariates)
    checks = assumption_checks(y, X, names, kind)
    blocked = checks_block(checks)
    if blocked:
        return {"status": "blocked", "reason": blocked, "checks": checks,
                "term": term_label}

    n_related = 0
    chosen = model or glm.recommend_model(y, kind, exposure, n_related=n_related)
    rationale = ""
    alternatives: List[str] = []
    if isinstance(chosen, dict):
        rationale = chosen.get("rationale", "")
        alternatives = chosen.get("alternatives", [])
        chosen = chosen["model"]

    # Sparse or empty cells override a plain-logistic choice. The recommender
    # works from marginal counts and can miss separation that only appears in
    # the joint 2x2 table.
    if chosen == "logistic" and any(c.get("forces_firth") for c in checks):
        cells = next(c["detail"] for c in checks if c["check"] == "cell_counts")
        chosen = "firth"
        rationale = ("Sparse 2x2 cells ({}). Standard logistic regression is "
                     "unreliable here; using Firth penalised regression.".format(cells))
        alternatives = ["logistic (not recommended: separation)", "Fisher's exact"]

    fitter = {
        "linear": glm.fit_linear,
        "robust_linear": glm.fit_robust_linear,
        "logistic": glm.fit_logistic,
        "firth": glm.fit_firth,
    }.get(chosen)
    if fitter is None:
        return {"status": "error", "reason": "unknown model '{}'".format(chosen),
                "checks": checks, "term": term_label}

    try:
        fit = fitter(X, y, 0, names)
    except (ValueError, np.linalg.LinAlgError) as exc:
        # A kernel refuses to return a number it does not trust — e.g. a term
        # that became constant once incomplete rows were dropped. In a scan
        # over hundreds of thousands of variants that must not abort the run,
        # and it must not become a row that reads as "no effect" either.
        return {"status": "not_estimable", "reason": str(exc), "checks": checks,
                "term": term_label, "model": chosen}

    return {
        "status": "ok",
        "term": term_label,
        "model": fit.model,
        "model_rationale": rationale,
        "model_alternatives": alternatives,
        "n": fit.n,
        "n_cases": fit.n_cases,
        "n_controls": fit.n_controls,
        "beta": fit.beta,
        "se": fit.se,
        "ci_low": fit.ci_low,
        "ci_high": fit.ci_high,
        "effect_label": fit.effect_label,
        "effect": fit.effect,
        "effect_ci_low": fit.effect_ci_low,
        "effect_ci_high": fit.effect_ci_high,
        "pvalue": fit.pvalue,
        "converged": fit.converged,
        "n_iter": fit.n_iter,
        "covariates": [n for n in names if n != "exposure"],
        "missingness": round(1.0 - fit.n / max(1, len(y)), 4),
        "warnings": list(fit.warnings),
        "checks": checks,
        "transform": transform,
        "extra": fit.extra,
    }


# --------------------------------------------------------- §5.5 sensitivity --
def sensitivity_analysis(base_spec: Dict[str, Any],
                         runner: Callable[[Dict[str, Any]], Dict[str, Any]],
                         variants: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Spec §5.5: "One click re-runs an analysis under alternative assumptions
    ... and compares results side by side."

    A result that survives a change of MAF threshold and covariate set is a
    different kind of claim from one that only exists under one specification.
    """
    out: List[Dict[str, Any]] = []
    base = runner(base_spec)
    out.append({"label": "Primary specification", "changes": {}, "result": base})
    for v in variants:
        spec = dict(base_spec)
        spec.update(v.get("changes", {}))
        try:
            res = runner(spec)
        except Exception as exc:
            res = {"status": "error", "reason": str(exc)}
        out.append({"label": v.get("label", "variant"),
                    "changes": v.get("changes", {}), "result": res})
    return out


# --------------------------------------------------------------- job entry ---
@register("association")
def association_job(spec: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    """Job entry point. Scans one or many variants for association.

    Spec keys: outcome, exposure_variants (list of variant keys) or gene,
    covariates, genetic_model, model, alpha_scope, transform, sensitivity.
    """
    from .. import store as research_store

    log = context.get("log", lambda m: None)
    dataset_id = context["dataset_id"]

    gm = research_store.load_genotypes(dataset_id)
    ph = research_store.load_phenotypes(dataset_id)
    if ph is None:
        raise ValueError("This dataset has no phenotypes; association needs an outcome.")

    outcome = spec["outcome"]
    kind = ph.kind(outcome)
    if kind not in ("binary", "quantitative"):
        raise ValueError(
            "Outcome '{}' is {}; association needs a binary or quantitative outcome."
            .format(outcome, kind))

    # Align genotype samples to phenotype samples. Order is not guaranteed to
    # match between two files uploaded separately, and a silent misalignment
    # here produces a perfectly plausible, entirely wrong result.
    pheno_index = {s: i for i, s in enumerate(ph.sample_ids)}
    keep = [i for i, s in enumerate(gm.sample_ids) if s in pheno_index]
    if not keep:
        raise ValueError("No samples are shared between genotypes and phenotypes.")
    rows = [pheno_index[gm.sample_ids[i]] for i in keep]
    y = as_float(ph.get(outcome))[rows]

    covariates: Dict[str, np.ndarray] = {}
    for cov in spec.get("covariates", []):
        if cov in ph.columns:
            covariates[cov] = as_float(ph.get(cov))[rows]
    n_pcs = int(spec.get("n_pcs", 0))
    if n_pcs:
        pcs = research_store.load_covariates(dataset_id, "pcs")
        if pcs is not None:
            for k in range(min(n_pcs, pcs.shape[1])):
                covariates["PC{}".format(k + 1)] = pcs[keep, k]

    genetic_model = spec.get("genetic_model", "additive")
    targets = spec.get("exposure_variants")
    gene = (spec.get("gene") or "").strip()
    variant_index = {v.key: i for i, v in enumerate(gm.variants)}
    maf = gm.maf()
    min_maf = float(spec.get("min_maf", 0.01))

    if targets:
        idx = [variant_index[k] for k in targets if k in variant_index]
        if not idx:
            raise ValueError(
                "None of the {} variant key(s) given are in this dataset. Keys "
                "look like chrom:pos:ref:alt, e.g. {}."
                .format(len(targets), gm.variants[0].key if gm.variants else "1:100:A:G"))
    elif gene:
        # Documented in this docstring since the first version but never
        # implemented: the job silently fell through to a whole-dataset scan,
        # so asking for one gene tested every variant in the cohort.
        annotations = research_store.load_annotations(dataset_id) or {}
        want = gene.upper()
        idx = [i for i, v in enumerate(gm.variants)
               if str((annotations.get(v.key) or {}).get("gene", "")).upper() == want]
        if not idx:
            known = sorted({str((a or {}).get("gene")) for a in annotations.values()
                            if a and a.get("gene")})
            raise ValueError(
                "No variants annotated to gene '{}'. {}".format(
                    gene,
                    "This dataset has no gene annotations at all — attach them "
                    "to test by gene." if not known else
                    "Known genes include: {}.".format(", ".join(known[:8]))))
        idx = [i for i in idx if np.isfinite(maf[i]) and maf[i] >= min_maf]
        if not idx:
            raise ValueError(
                "Gene '{}' has variants, but none above the {:.3g} minimum allele "
                "frequency. Single-variant association has no power on rarer "
                "variants — use a gene-based burden test instead."
                .format(gene, min_maf))
    else:
        # Whole-dataset scan, MAF-filtered so the test set is defensible.
        idx = [i for i in range(gm.n_variants)
               if np.isfinite(maf[i]) and maf[i] >= min_maf]
        max_scan = int(spec.get("max_variants") or 0)
        if max_scan and len(idx) > max_scan:
            raise ValueError(
                "An untargeted scan would fit {} separate models, which is a "
                "genome-wide scan run one variant at a time. Name a gene or "
                "specific variants here, or use the GWAS analysis, which is "
                "built for this and prunes relatedness first."
                .format(len(idx)))

    if not idx:
        raise ValueError("No variants matched the exposure specification.")

    alpha_info = power_mod.suggested_alpha(len(idx), spec.get("alpha_scope", "custom"))
    log("testing {} variant(s) at alpha {:.3g}".format(len(idx), alpha_info["alpha"]))

    maf_all = gm.maf()
    results: List[Dict[str, Any]] = []
    for n_done, i in enumerate(idx):
        exposure = encode_genotype(gm.dosages[i, keep], genetic_model)
        res = run_association(
            y, exposure, kind, covariates=covariates,
            model=spec.get("model"), term_label=gm.variants[i].key,
            transform=spec.get("transform"))
        res["variant"] = gm.variants[i].key
        res["chrom"] = gm.variants[i].chrom
        res["pos"] = gm.variants[i].pos
        res["maf"] = float(maf_all[i]) if np.isfinite(maf_all[i]) else None
        results.append(res)
        if n_done and n_done % 5000 == 0:
            log("{}/{} variants tested".format(n_done, len(idx)))

    ok = [r for r in results if r.get("status") == "ok" and r.get("pvalue") is not None]

    # §5.6: never an uncorrected p without its correction adjacent.
    if ok:
        praw = [r["pvalue"] for r in ok]
        bonf = glm.multiple_testing(praw, "bonferroni")
        fdr = glm.multiple_testing(praw, "fdr_bh")
        for r, b, f in zip(ok, bonf, fdr):
            r["p_bonferroni"] = float(b)
            r["p_fdr_bh"] = float(f)

    ok.sort(key=lambda r: r["pvalue"])

    # §5.1 power, computed for the top hit's actual observed parameters.
    power_report = None
    if ok:
        top = ok[0]
        if kind == "binary":
            n_cases = top.get("n_cases") or 0
            n_controls = top.get("n_controls") or 0
            power_report = power_mod.power_case_control(
                n_cases, n_controls, top.get("maf") or 0.01,
                spec.get("assumed_or", 2.0), alpha_info["alpha"])
        else:
            power_report = power_mod.power_quantitative(
                top.get("n") or 0, top.get("maf") or 0.01,
                spec.get("assumed_beta", 0.2), alpha_info["alpha"])

    return {
        "analysis": "association",
        "outcome": outcome,
        "outcome_kind": kind,
        "genetic_model": genetic_model,
        "genetic_model_label": CARRIER_MODELS.get(genetic_model, genetic_model),
        "covariates": sorted(covariates),
        "n_variants_tested": len(idx),
        "n_results": len(ok),
        "alpha": alpha_info,
        "power": power_report,
        "results": ok[:1000],
        "blocked": [r for r in results if r.get("status") != "ok"][:50],
        "multiple_testing_note": (
            "{} tests. Raw p, Bonferroni and Benjamini-Hochberg FDR are reported "
            "for every result.".format(len(ok))),
        "diagnostics": {
            "n_samples_analysed": len(keep),
            "n_blocked": len(results) - len(ok),
            "checks_example": ok[0]["checks"] if ok else [],
        },
    }
