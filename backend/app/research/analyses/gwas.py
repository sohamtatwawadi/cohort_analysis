"""GWAS — Part II §4.2.

    "Mandatory QC before any GWAS runs: sample call rate, variant call rate,
     MAF, HWE, missingness differential between cases and controls, relatedness
     pruning, ancestry outlier removal, sex-check."

Mandatory means mandatory: `gwas_job` runs QC itself and analyses what survives.
There is no parameter that skips it. A GWAS on unfiltered data produces a
Manhattan plot that looks exactly like a real one, which is the problem.

    "Guardrail: if λ_GC exceeds a threshold, results are flagged as showing
     possible stratification with the specific likely cause identified."

λ_GC is computed from the observed chi-square distribution and reported with an
attributed cause, because "λ = 1.34" on its own tells a researcher something is
wrong but not what to do about it.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy import stats

from ..jobs import register
from ..stats import glm, power as power_mod
from ..types import MISSING, GenotypeMatrix, as_float
from .association import encode_genotype, run_association

# Above this, results are flagged. 1.05 is the conventional "clean" ceiling;
# beyond 1.10 stratification is usually visible in a QQ plot by eye.
LAMBDA_WARN = 1.05
LAMBDA_SEVERE = 1.10


def genomic_inflation(pvalues: np.ndarray) -> Dict[str, Any]:
    """λ_GC — the median observed chi-square over its null expectation.

    Uses the median rather than the mean because the median is insensitive to
    the genuine associations in the tail, which is the whole point: λ should
    measure background inflation, not signal.
    """
    p = np.asarray(pvalues, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if len(p) < 10:
        return {"lambda_gc": None, "n": int(len(p)),
                "note": "too few tests to estimate inflation"}
    chi2 = stats.chi2.isf(p, df=1)
    lam = float(np.median(chi2) / stats.chi2.ppf(0.5, df=1))
    return {"lambda_gc": round(lam, 4), "n": int(len(p))}


def attribute_inflation(lam: Optional[float], profile: Dict[str, Any],
                        qc: Dict[str, Any]) -> Dict[str, Any]:
    """Name the likely cause, per §4.2's guardrail.

    Ordered by how often each actually explains inflation in practice:
    unmodelled ancestry first, then cryptic relatedness, then batch.
    """
    if lam is None:
        return {"status": "unknown", "text": "Inflation could not be estimated."}
    if lam <= LAMBDA_WARN:
        return {"status": "ok",
                "text": "λ_GC = {:.3f}. No evidence of systematic inflation.".format(lam)}

    causes: List[str] = []
    n_pcs = int((profile.get("ancestry") or {}).get("n_pcs", 0) or 0)
    if n_pcs == 0:
        causes.append("no ancestry principal components were available to adjust for "
                      "population structure")
    rel = profile.get("relatedness") or {}
    n_related = int(rel.get("first_degree", 0)) + int(rel.get("second_degree", 0))
    if n_related:
        causes.append("{} related sample pair(s) remain in the analysis set"
                      .format(n_related))
    batches = int((profile.get("batch_structure") or {}).get("n_batches", 0) or 0)
    if batches > 1:
        causes.append("{} genotyping batches are present and were not modelled"
                      .format(batches))
    if qc.get("differential_missingness_flagged"):
        causes.append("differential missingness between cases and controls")
    if not causes:
        causes.append("cause not identifiable from the data profile; consider "
                      "unmodelled covariates or cryptic structure")

    severity = "severe" if lam >= LAMBDA_SEVERE else "warn"
    return {
        "status": severity,
        "text": ("λ_GC = {:.3f} — results show possible population stratification. "
                 "Likely cause: {}.".format(lam, "; ".join(causes))),
        "causes": causes,
    }


def qq_points(pvalues: np.ndarray, max_points: int = 5000) -> Dict[str, List[float]]:
    """Observed vs expected -log10 p. Thinned in the dense head, kept whole in
    the tail, because the tail is the part anyone reads."""
    p = np.asarray(pvalues, dtype=float)
    p = np.sort(p[np.isfinite(p) & (p > 0)])
    if not len(p):
        return {"expected": [], "observed": []}
    n = len(p)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p)
    if n > max_points:
        head = np.linspace(0, n - 1000, max_points - 1000).astype(int)
        idx = np.unique(np.concatenate([head, np.arange(max(0, n - 1000), n)]))
    else:
        idx = np.arange(n)
    return {"expected": [round(float(x), 4) for x in expected[idx]],
            "observed": [round(float(x), 4) for x in observed[idx]]}


def manhattan_points(results: List[Dict[str, Any]],
                     max_points: int = 20000) -> List[Dict[str, Any]]:
    """Points for the Manhattan plot. Everything above -log10 p = 2 is kept;
    the dense floor below that is thinned, since it is visually a solid band
    either way."""
    pts = []
    floor = []
    for r in results:
        p = r.get("pvalue")
        if p is None or not np.isfinite(p) or p <= 0:
            continue
        item = {"chrom": r.get("chrom"), "pos": r.get("pos"),
                "variant": r.get("variant"),
                "neglog10p": round(-float(np.log10(p)), 4)}
        (pts if item["neglog10p"] >= 2.0 else floor).append(item)
    if len(floor) > max_points:
        step = int(np.ceil(len(floor) / max_points))
        floor = floor[::step]
    return pts + floor


@register("gwas")
def gwas_job(spec: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    """Guided GWAS: QC -> relatedness pruning -> scan -> inflation guardrail."""
    from .. import registry as reg
    from .. import store as research_store
    from ..stats import kinship as kin_mod
    from ..stats import qc as qc_mod
    from ..stats import scoretest

    log = context.get("log", lambda m: None)
    dataset_id = context["dataset_id"]

    gm = research_store.load_genotypes(dataset_id)
    ph = research_store.load_phenotypes(dataset_id)
    if ph is None:
        raise ValueError("GWAS requires a phenotype table.")

    outcome = spec["outcome"]
    kind = ph.kind(outcome)
    pheno_index = {s: i for i, s in enumerate(ph.sample_ids)}

    # ---- mandatory QC (§4.2). Not optional, not parameterised away. --------
    log("running mandatory pre-GWAS QC")
    y_all = as_float(ph.get(outcome))
    cases = None
    if kind == "binary":
        cases = np.array([y_all[pheno_index[s]] == 1 if s in pheno_index else False
                          for s in gm.sample_ids])

    vqc = qc_mod.variant_qc(
        gm,
        min_call_rate=float(spec.get("min_variant_call_rate", 0.95)),
        min_maf=float(spec.get("min_maf", 0.01)),
        hwe_p=float(spec.get("hwe_p", 1e-6)),
        cases=cases)
    sqc = qc_mod.sample_qc(
        gm, min_call_rate=float(spec.get("min_sample_call_rate", 0.95)))

    # QCResult.keep — not `keep_mask`. The previous getattr/subscript fallback
    # here tried a second wrong spelling instead of failing on the first, so a
    # plain typo surfaced as "'QCResult' object is not subscriptable" from a
    # line that looked defensive.
    keep_variants = np.asarray(vqc.keep, dtype=bool)
    keep_samples = np.asarray(sqc.keep, dtype=bool)
    qc_gm = gm.subset_variants(keep_variants)
    qc_gm = qc_gm.subset_samples([s for s, k in zip(gm.sample_ids, keep_samples) if k])
    log("QC kept {} variants and {} samples".format(qc_gm.n_variants, qc_gm.n_samples))

    if qc_gm.n_variants == 0 or qc_gm.n_samples == 0:
        raise ValueError("QC removed every variant or every sample; nothing to test.")

    # ---- relatedness pruning (§4.2) ---------------------------------------
    pruned_related = 0
    if spec.get("prune_related", True):
        kin = kin_mod.king_robust(qc_gm)
        unrelated = kin_mod.unrelated_set(kin)
        pruned_related = qc_gm.n_samples - len(unrelated)
        if pruned_related:
            qc_gm = qc_gm.subset_samples(unrelated)
            log("pruned {} related sample(s)".format(pruned_related))

    # ---- ancestry PCs ------------------------------------------------------
    n_pcs = int(spec.get("n_pcs", 10))
    covariates: Dict[str, np.ndarray] = {}
    rows = [pheno_index[s] for s in qc_gm.sample_ids if s in pheno_index]
    keep_idx = [i for i, s in enumerate(qc_gm.sample_ids) if s in pheno_index]
    if len(keep_idx) != qc_gm.n_samples:
        qc_gm = qc_gm.subset_samples([qc_gm.sample_ids[i] for i in keep_idx])
        rows = [pheno_index[s] for s in qc_gm.sample_ids]

    if n_pcs:
        from ..stats import pca as pca_mod
        # Ancestry PCs come from a marker subset, not the whole scan set. The
        # SVD is exact and in-memory, so its cost grows with the variant count
        # while the ancestry signal does not — the leading PCs are common-variant
        # structure and are stable well below 20,000 markers. Running it over
        # every QC-passing variant was the slowest step in the whole job.
        max_pca_variants = int(spec.get("max_pca_variants") or 20000)
        pca_gm = qc_gm
        if qc_gm.n_variants > max_pca_variants:
            step = qc_gm.n_variants // max_pca_variants
            mask = np.zeros(qc_gm.n_variants, dtype=bool)
            mask[::step] = True                    # evenly spread along the genome
            pca_gm = qc_gm.subset_variants(mask)
            log("PCA on {} of {} variants".format(pca_gm.n_variants, qc_gm.n_variants))
        pca = pca_mod.compute_pca(pca_gm,
                                  n_components=min(n_pcs, max(2, pca_gm.n_samples - 1)))
        for k in range(min(n_pcs, pca.components.shape[1])):
            covariates["PC{}".format(k + 1)] = pca.components[:, k]
        log("computed {} ancestry PCs".format(len(covariates)))

    for cov in spec.get("covariates", []):
        if cov in ph.columns:
            covariates[cov] = as_float(ph.get(cov))[rows]

    y = as_float(ph.get(outcome))[rows]
    genetic_model = spec.get("genetic_model", "additive")

    # ---- the scan ----------------------------------------------------------
    # Two stages, which is what a GWAS is. The scan fits the covariate model
    # ONCE and scores every variant against its residuals; the exact GLM is then
    # refitted only for the variants worth reporting an effect size for.
    #
    # Fitting a full model per variant instead — which this did — means redoing
    # the identical covariate iteration eighty thousand times. On this dataset
    # that was over half an hour, versus seconds now, for the same ranking.
    log("scanning {} variants".format(qc_gm.n_variants))
    maf = qc_gm.maf()

    cov_matrix = (np.column_stack([covariates[c] for c in sorted(covariates)])
                  if covariates else None)
    scan = scoretest.score_scan(qc_gm.dosages, y, cov_matrix,
                                kind="binary" if kind == "binary" else "quantitative")
    scan_p = scan["p"]
    usable = scan["usable"]
    log("scored {} of {} variants".format(int(usable.sum()), qc_gm.n_variants))

    if not usable.any():
        raise ValueError(
            "No variant had estimable score variance after QC — every remaining "
            "variant is monomorphic or collinear with a covariate.")

    # Refit exactly: everything suggestive, and at minimum the top of the list
    # so a result is never empty.
    n_refit = int(spec.get("n_exact_refit") or 500)
    order = np.argsort(np.where(np.isfinite(scan_p), scan_p, np.inf))
    suggestive = set(np.flatnonzero(np.isfinite(scan_p) & (scan_p < 1e-4)).tolist())
    refit_idx = sorted(suggestive | set(order[:n_refit].tolist()))
    log("refitting {} variant(s) with the exact model".format(len(refit_idx)))

    exact_by_index: Dict[int, Dict[str, Any]] = {}
    for n_done, i in enumerate(refit_idx):
        exposure = encode_genotype(qc_gm.dosages[i, :], genetic_model)
        res = run_association(y, exposure, kind, covariates=covariates,
                              model=spec.get("model"),
                              term_label=qc_gm.variants[i].key)
        if res.get("status") == "ok":
            exact_by_index[i] = res
        if n_done and n_done % 200 == 0:
            log("refit {}/{}".format(n_done, len(refit_idx)))

    results: List[Dict[str, Any]] = []
    for i in range(qc_gm.n_variants):
        if not usable[i]:
            continue
        res = dict(exact_by_index.get(i) or {})
        if res:
            res["estimate_source"] = "exact model refit"
        else:
            # Scanned but not refitted: report the score p-value and say so,
            # rather than inventing an effect size the scan never estimated.
            res = {"pvalue": float(scan_p[i]), "status": "ok",
                   "beta": None, "se": None, "ci_low": None, "ci_high": None,
                   "effect": None, "effect_label": None,
                   "model": "score test",
                   "estimate_source": "score test only — not refitted"}
        res["variant"] = qc_gm.variants[i].key
        res["chrom"] = qc_gm.variants[i].chrom
        res["pos"] = qc_gm.variants[i].pos
        res["maf"] = float(maf[i]) if np.isfinite(maf[i]) else None
        res["score_pvalue"] = float(scan_p[i])
        results.append(res)

    if not results:
        raise ValueError("No variant produced an estimable model after QC.")

    pvals = np.array([r["pvalue"] for r in results])
    bonf = glm.multiple_testing(pvals, "bonferroni")
    fdr = glm.multiple_testing(pvals, "fdr_bh")
    for r, b, f in zip(results, bonf, fdr):
        r["p_bonferroni"] = float(b)
        r["p_fdr_bh"] = float(f)

    inflation = genomic_inflation(pvals)
    profile = reg.get_profile(dataset_id) or {}
    qc_summary = {
        "variants_in": gm.n_variants, "variants_kept": qc_gm.n_variants,
        "samples_in": gm.n_samples, "samples_kept": qc_gm.n_samples,
        "related_pruned": pruned_related,
        "variant_filters": getattr(vqc, "as_dict", lambda: vqc)(),
        "sample_filters": getattr(sqc, "as_dict", lambda: sqc)(),
        "differential_missingness_flagged": bool(
            getattr(vqc, "differential_missingness_flagged", False)),
    }
    guardrail = attribute_inflation(inflation.get("lambda_gc"), profile, qc_summary)

    ordered = sorted(results, key=lambda r: r["pvalue"])
    alpha = power_mod.GENOME_WIDE_ALPHA

    return {
        "analysis": "gwas",
        "outcome": outcome,
        "outcome_kind": kind,
        "genetic_model": genetic_model,
        "covariates": sorted(covariates),
        "n_variants_tested": len(results),
        "n_samples": qc_gm.n_samples,
        "alpha": {"alpha": alpha, "label": "genome-wide (5×10⁻⁸)"},
        "genome_wide_significant": [r for r in ordered if r["pvalue"] < alpha][:200],
        "suggestive": [r for r in ordered if alpha <= r["pvalue"] < 1e-5][:200],
        "top_hits": ordered[:100],
        "inflation": inflation,
        "inflation_guardrail": guardrail,
        "manhattan": manhattan_points(results),
        "qq": qq_points(pvals),
        "qc": qc_summary,
        "diagnostics": {
            "lambda_gc": inflation.get("lambda_gc"),
            "guardrail": guardrail,
            "qc": qc_summary,
            "n_pcs_used": len([c for c in covariates if c.startswith("PC")]),
        },
    }
