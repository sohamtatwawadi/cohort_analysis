"""Ancestry principal components.

Part II §3.1 profiles "ancestry composition (PCA vs reference)"; §4.2 makes
"ancestry outlier removal" mandatory QC before a GWAS runs and lists PCs among
the required covariates; §3.2 will not unlock GWAS or single-variant
association without them.

Why PCs at all: unadjusted association testing in a structured cohort finds
allele-frequency differences between ancestry groups, not disease biology.
That is the classic false-positive mechanism, and it shows up as λ_GC
inflation (§4.2 guardrail). Projecting samples onto the top eigenvectors of the
genotype covariance gives continuous axes that absorb that structure as
regression covariates.

Why LD-prune first: PCA finds directions of maximum variance. A block of
variants in tight LD is one signal duplicated a thousand times, so an
uncorrected PCA happily spends its leading eigenvectors on a single large
inversion or an HLA-style region and never reaches ancestry. Pruning to a
roughly independent marker set is the standard fix, and it is the same
pre-processing kinship estimation needs (see kinship.py).

Why §4.8 provenance: ancestry PCs are a genetic measurement, not an ethnicity
label, and they are meaningless without the marker set that produced them. The
result carries that provenance so it cannot be silently mixed with a
self-reported category.

Scale note: this is an exact in-memory SVD, appropriate for prototype and
mid-size cohorts. Biobank-scale work wants a randomised/truncated SVD over a
memory-mapped matrix, or delegation to an established tool (§4.2 takes the same
position about GWAS engines).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from ..types import MISSING, GenotypeMatrix

# Below this many markers, PCs are noise-dominated and do not describe ancestry.
# 10k LD-pruned common variants is the usual working floor for genome-wide data.
MIN_VARIANTS_FOR_ANCESTRY = 10000
SOFT_MIN_VARIANTS = 1000


@dataclass
class PCAResult:
    components: np.ndarray               # (n_samples, n_components) PC scores
    explained_variance_ratio: np.ndarray  # (n_components,)
    n_variants_used: int
    sample_ids: List[str]
    warnings: List[str] = field(default_factory=list)
    # Provenance (§4.8): these PCs are only interpretable alongside the marker
    # set and standardisation that produced them.
    variant_mask: Optional[np.ndarray] = None   # aligned to the input gm.variants
    provenance: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_components(self) -> int:
        return int(self.components.shape[1]) if self.components.ndim == 2 else 0


# --------------------------------------------------------------- internals --
def _row_zscores(dosages: np.ndarray) -> np.ndarray:
    """Mean-impute missing calls, then z-score each variant (row).

    Imputation is to the per-variant mean dosage, which is the same choice the
    PCA standardisation makes (2p): a missing call contributes nothing to any
    covariance rather than pretending to be a reference call, which would look
    like a real allele-frequency difference.
    """
    d = dosages.astype(np.float64, copy=True)
    miss = dosages == MISSING
    if miss.any():
        d[miss] = 0.0
        n_called = (~miss).sum(axis=1)
        totals = d.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            means = np.where(n_called > 0, totals / np.maximum(n_called, 1), 0.0)
        d = d - means[:, None]
        d[miss] = 0.0          # imputed == mean == 0 after centring
    else:
        d = d - d.mean(axis=1)[:, None]
    sd = np.sqrt((d * d).mean(axis=1))
    # sd == 0 is a monomorphic (or fully missing) variant: correlation with it
    # is undefined, so zero it out and let it correlate with nothing.
    safe = sd > 0
    out = np.zeros_like(d)
    out[safe] = d[safe] / sd[safe, None]
    return out


def _standardise(dosages: np.ndarray, af: np.ndarray) -> np.ndarray:
    """Centre by 2p and scale by sqrt(2p(1-p)); missing calls impute to 2p.

    This is the Patterson/Price standardisation: dividing by the binomial SD
    under HWE, rather than the empirical SD, up-weights rare variants in
    proportion to the information they carry about drift. Callers must have
    already dropped p == 0 / p == 1 or this divides by zero.
    """
    p = af[:, None]
    miss = dosages == MISSING
    d = dosages.astype(np.float64, copy=True)
    d[miss] = 0.0
    d = d - 2.0 * p
    d[miss] = 0.0                       # imputed to 2p, i.e. centred to zero
    return d / np.sqrt(2.0 * p * (1.0 - p))


# ------------------------------------------------------------- LD pruning --
def ld_prune(gm: GenotypeMatrix,
             window: int = 50,
             step: int = 5,
             r2_threshold: float = 0.2) -> np.ndarray:
    """Sliding-window pairwise-r² pruning; returns a keep mask over variants.

    Same semantics as `plink --indep-pairwise <window> <step> <r2>`: window and
    step are counts of *variants*, not base pairs. Inside each window, any
    still-kept pair with r² above the threshold loses its later member.

    Windowing is what keeps this tractable — LD decays with distance, so the
    only pairs worth testing are nearby ones. Cost is
    O(n_variants/step · window² · n_samples) rather than the O(n_variants²) an
    all-pairs r² matrix would need.

    Keeping the earlier variant of a correlated pair (rather than, say, the
    higher-MAF one) makes the output a deterministic function of input order,
    which §5.5/§4.2 reproducibility requires.
    """
    m = gm.n_variants
    keep = np.ones(m, dtype=bool)
    if m < 2 or gm.n_samples < 3:
        return keep
    if window < 2 or step < 1:
        raise ValueError("window must be >= 2 and step >= 1")

    z = _row_zscores(gm.dosages)
    n = float(gm.n_samples)

    start = 0
    while start < m:
        end = min(start + window, m)
        idx = np.nonzero(keep[start:end])[0] + start
        if idx.size > 1:
            zz = z[idx]
            # rows are z-scored (mean 0, unit second moment), so ZZ'/n is r.
            r2 = (zz @ zz.T / n) ** 2
            alive = np.ones(idx.size, dtype=bool)
            for a in range(idx.size - 1):
                if not alive[a]:
                    continue
                alive[a + 1:] &= ~(r2[a, a + 1:] > r2_threshold)
            keep[idx[~alive]] = False
        if end >= m:
            break
        start += step
    return keep


# --------------------------------------------------------------------- PCA --
def compute_pca(gm: GenotypeMatrix,
                n_components: int = 10,
                prune: bool = True,
                window: int = 50,
                step: int = 5,
                r2_threshold: float = 0.2) -> PCAResult:
    """Ancestry PCs from a genotype matrix.

    Monomorphic variants are dropped rather than special-cased: 2p(1-p) is zero
    for them, so they carry no information and would only produce a division by
    zero.
    """
    warnings: List[str] = []
    n_samples = gm.n_samples
    if n_samples == 0 or gm.n_variants == 0:
        raise ValueError("cannot run PCA on an empty genotype matrix")

    mask = ld_prune(gm, window=window, step=step,
                    r2_threshold=r2_threshold) if prune else np.ones(
                        gm.n_variants, dtype=bool)
    n_after_prune = int(mask.sum())
    if prune:
        n_pruned = gm.n_variants - n_after_prune
        if n_pruned:
            warnings.append(
                "LD pruning removed {} of {} variants (r2 > {})".format(
                    n_pruned, gm.n_variants, r2_threshold))

    dos = gm.dosages[mask, :]
    sub = gm.subset_variants(mask)
    af = sub.allele_frequency()

    polymorphic = np.isfinite(af) & (af > 0.0) & (af < 1.0)
    n_dropped = int((~polymorphic).sum())
    if n_dropped:
        warnings.append(
            "dropped {} monomorphic or uncallable variants "
            "(2p(1-p) == 0)".format(n_dropped))
    dos = dos[polymorphic, :]
    af = af[polymorphic]
    n_used = int(dos.shape[0])
    if n_used == 0:
        raise ValueError("no polymorphic variants left after pruning/filtering")

    # Composite mask back onto the original variant list, for provenance.
    used_mask = np.zeros(gm.n_variants, dtype=bool)
    used_mask[np.nonzero(mask)[0][polymorphic]] = True

    missing_rate = float((dos == MISSING).mean())
    if missing_rate > 0.05:
        warnings.append(
            "{:.1%} of genotypes were missing and mean-imputed to 2p; "
            "imputation shrinks samples toward the origin of the PC "
            "space".format(missing_rate))

    if n_used < SOFT_MIN_VARIANTS:
        warnings.append(
            "PCs computed from only {} variants: ancestry axes are normally "
            "estimated from >= {} LD-pruned common variants, and PCs on this "
            "few markers are noise-dominated and do not reliably capture "
            "ancestry".format(n_used, MIN_VARIANTS_FOR_ANCESTRY))
    elif n_used < MIN_VARIANTS_FOR_ANCESTRY:
        warnings.append(
            "PCs computed from {} variants, below the {} usually wanted for "
            "stable ancestry axes; treat trailing PCs with "
            "caution".format(n_used, MIN_VARIANTS_FOR_ANCESTRY))
    if n_used < n_samples:
        warnings.append(
            "fewer variants ({}) than samples ({}): the genotype matrix is "
            "rank-deficient and PCs beyond rank {} are arbitrary".format(
                n_used, n_samples, n_used))
    if n_samples < 10:
        warnings.append(
            "only {} samples: PC structure is not interpretable at this "
            "size".format(n_samples))

    # samples x variants, so the left singular vectors are the sample loadings.
    x = _standardise(dos, af).T
    u, s, _vt = np.linalg.svd(x, full_matrices=False)

    max_k = int(min(x.shape))
    k = int(max(1, min(n_components, max_k)))
    if k < n_components:
        warnings.append(
            "requested {} components but the data supports at most {}".format(
                n_components, max_k))

    total_var = float((s ** 2).sum())
    evr = (s ** 2) / total_var if total_var > 0 else np.zeros_like(s)

    return PCAResult(
        components=u[:, :k] * s[:k],
        explained_variance_ratio=evr[:k],
        n_variants_used=n_used,
        sample_ids=list(gm.sample_ids),
        warnings=warnings,
        variant_mask=used_mask,
        provenance={
            "method": "SVD of 2p-centred, sqrt(2p(1-p))-scaled dosages",
            "ld_pruned": bool(prune),
            "ld_window": window,
            "ld_step": step,
            "ld_r2_threshold": r2_threshold,
            "n_variants_input": int(gm.n_variants),
            "n_variants_after_prune": n_after_prune,
            "n_variants_used": n_used,
            "n_samples": n_samples,
            "build": gm.build,
            "source_format": gm.source_format,
            # §4.8: these axes describe genetic similarity within THIS marker
            # set and cohort. They are not transferable ethnicity labels and
            # must not be compared across datasets or to self-report.
            "basis": "cohort-internal (not projected onto a reference panel)",
        },
    )


def ancestry_outliers(pcs: Any, n_sd: float = 6.0, n_pcs: int = 4) -> np.ndarray:
    """Boolean mask of samples lying beyond n_sd on any of the first n_pcs.

    §4.2 mandates "ancestry outlier removal": a handful of samples far from the
    cohort's ancestry cloud cannot be adjusted for by PC covariates (they are
    leverage points, not a population), so they are dropped instead.

    Mean/SD per PC, the smartpca convention. It is deliberately not a robust
    (median/MAD) centre: with the conventional n_sd = 6 the inflation of SD by
    the outliers themselves is what keeps this conservative — it flags the
    samples nobody disputes and leaves borderline ones in.
    """
    arr = getattr(pcs, "components", pcs)
    arr = np.asarray(arr, dtype=float)
    if arr.ndim != 2:
        raise ValueError("pcs must be 2-D (n_samples, n_components)")
    n_samples = arr.shape[0]
    out = np.zeros(n_samples, dtype=bool)
    if n_samples < 3:
        return out
    k = int(min(n_pcs, arr.shape[1]))
    for j in range(k):
        col = arr[:, j]
        sd = float(col.std(ddof=1))
        if sd <= 0 or not np.isfinite(sd):
            continue
        out |= np.abs(col - col.mean()) > n_sd * sd
    return out
