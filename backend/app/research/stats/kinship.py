"""Relatedness — KING-robust kinship.

Part II §3.1 profiles "relatedness structure (kinship distribution)"; §3.2
gates GWAS on an "unrelated set"; §4.2 makes "relatedness pruning" mandatory
QC. Cryptic relatedness is a false-positive generator with the same signature
as population stratification: related samples are not independent
observations, so the effective n is smaller than the nominal n and the null
distribution of the test statistic is too wide.

Why KING-robust and not a GRM. The obvious estimator — correlate two samples'
standardised genotypes against a cohort allele frequency, i.e. a GRM entry —
assumes every sample is drawn from one population with one allele frequency
vector. In an admixed or multi-ancestry cohort (§3.1 explicitly profiles
ancestry composition, so we must assume one) that assumption fails and the GRM
confounds ancestry sharing with family sharing: two unrelated South Asian
samples in a mostly European cohort look like cousins. KING-robust conditions
on the *pair's own* genotypes — it only counts heterozygote and opposite-
homozygote configurations between i and j — so no population frequency enters
the estimate at all, and it stays approximately unbiased under structure.

The price is that KING-robust is biased toward zero for pairs from genuinely
different ancestries, and can go negative. We report the raw value rather than
clipping, because a strongly negative φ is itself diagnostic (it means the two
samples are ancestrally divergent, which §4.2's stratification guardrail wants
to know about).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List

import numpy as np

from ..types import MISSING, GenotypeMatrix
from .pca import ld_prune

# Conventional kinship bands (Manichaikul et al. 2010). The cut points are the
# midpoints on a log2 scale between the expected values 0.5, 0.25, 0.125,
# 0.0625, 0 — which is why they look like odd decimals.
DUPLICATE_MIN = 0.354      # phi > 0.354            -> duplicate / MZ twin
FIRST_DEGREE_MIN = 0.177   # 0.177  < phi <= 0.354  -> parent-offspring, full sib
SECOND_DEGREE_MIN = 0.0884  # 0.0884 < phi <= 0.177  -> half sib, avuncular, grandparent
THIRD_DEGREE_MIN = 0.0442  # 0.0442 < phi <= 0.0884 -> first cousin

DEGREES = ("duplicate", "first-degree", "second-degree", "third-degree", "unrelated")


@dataclass
class KinshipResult:
    sample_ids: List[str]
    matrix: np.ndarray        # (n, n) kinship coefficients, diagonal 0.5
    pairs: List[Dict]         # related pairs above threshold, with degree
    warnings: List[str] = field(default_factory=list)
    n_variants_used: int = 0

    def sample_index(self) -> Dict[str, int]:
        return {s: i for i, s in enumerate(self.sample_ids)}


def classify_degree(phi: float) -> str:
    """Relationship band for a kinship coefficient."""
    if not np.isfinite(phi):
        return "unrelated"          # undefined pair; treated as unrelated, flagged in warnings
    if phi > DUPLICATE_MIN:
        return "duplicate"
    if phi > FIRST_DEGREE_MIN:
        return "first-degree"
    if phi > SECOND_DEGREE_MIN:
        return "second-degree"
    if phi > THIRD_DEGREE_MIN:
        return "third-degree"
    return "unrelated"


def king_robust(gm: GenotypeMatrix,
                prune: bool = True,
                pair_threshold: float = THIRD_DEGREE_MIN,
                window: int = 50,
                step: int = 5,
                r2_threshold: float = 0.2) -> KinshipResult:
    """KING-robust kinship for every sample pair.

    The estimator, over the variants called in BOTH members of the pair:

        phi_ij = (N_AaAa - 2 * N_AAaa) / (2 * min(N_Aa_i, N_Aa_j))

      N_AaAa  = # variants where both i and j are heterozygous
      N_AAaa  = # variants where one is hom-ref and the other hom-alt
      N_Aa_i  = # variants where i is heterozygous
      N_Aa_j  = # variants where j is heterozygous

    Opposite homozygotes are the signal: two people who share recent ancestry
    essentially cannot be AA and aa at the same site, so N_AAaa is a direct
    count of IBD violations. min() in the denominator rather than a sum or
    average is what makes it robust — it normalises by the less heterozygous of
    the two, so an ancestry difference in overall heterozygosity cannot inflate
    the ratio.

    LD pruning first for the same reason as PCA: correlated markers are not
    independent draws, and repeated counting of the same haplotype block
    distorts both numerator and denominator.
    """
    warnings: List[str] = []
    n = gm.n_samples
    if n == 0:
        raise ValueError("cannot compute kinship on zero samples")

    if prune:
        mask = ld_prune(gm, window=window, step=step, r2_threshold=r2_threshold)
        n_pruned = int(gm.n_variants - mask.sum())
        if n_pruned:
            warnings.append("LD pruning removed {} of {} variants".format(
                n_pruned, gm.n_variants))
        dos = gm.dosages[mask, :]
    else:
        dos = gm.dosages

    m = int(dos.shape[0])
    if m == 0:
        raise ValueError("no variants left for kinship estimation")
    if m < 1000:
        warnings.append(
            "kinship estimated from only {} variants; KING needs a few "
            "thousand independent markers before the degree bands "
            "separate cleanly".format(m))

    het = (dos == 1)
    hom_ref = (dos == 0)
    hom_alt = (dos == 2)
    called = (dos != MISSING)

    h = het.astype(np.float64)
    c = called.astype(np.float64)
    # Every product below is implicitly restricted to variants called in both
    # samples, because a missing call is neither het nor hom in either operand.
    n_het_het = h.T @ h                                   # N_AaAa
    n_opp_hom = (hom_ref.astype(np.float64).T @ hom_alt.astype(np.float64))
    n_opp_hom = n_opp_hom + n_opp_hom.T                   # N_AAaa, symmetric
    het_i = h.T @ c                                       # [i, j] = N_Aa_i over shared calls

    denom = 2.0 * np.minimum(het_i, het_i.T)
    with np.errstate(invalid="ignore", divide="ignore"):
        phi = np.where(denom > 0,
                       (n_het_het - 2.0 * n_opp_hom) / np.where(denom > 0, denom, 1.0),
                       np.nan)
    phi = np.asarray(phi, dtype=float)
    phi = 0.5 * (phi + phi.T)          # kill float asymmetry from the matmuls
    np.fill_diagonal(phi, 0.5)         # a sample's kinship with itself is 0.5 by definition

    off = ~np.eye(n, dtype=bool)
    n_undefined = int(np.isnan(phi[off]).sum() // 2)
    if n_undefined:
        warnings.append(
            "{} pairs have no heterozygous site in common; KING is undefined "
            "for them and they are reported as NaN".format(n_undefined))

    n_shared = (c.T @ c)
    pairs: List[Dict] = []
    iu, ju = np.triu_indices(n, k=1)
    for i, j in zip(iu, ju):
        val = phi[i, j]
        if np.isfinite(val) and val > pair_threshold:
            pairs.append({
                "sample1": gm.sample_ids[i],
                "sample2": gm.sample_ids[j],
                "kinship": float(val),
                "degree": classify_degree(float(val)),
                "n_variants": int(n_shared[i, j]),
            })
    # Deterministic ordering: strongest first, then by id (§5.5 reproducibility).
    pairs.sort(key=lambda d: (-d["kinship"], d["sample1"], d["sample2"]))

    return KinshipResult(sample_ids=list(gm.sample_ids), matrix=phi,
                         pairs=pairs, warnings=warnings, n_variants_used=m)


def unrelated_set(kin: KinshipResult, threshold: float = SECOND_DEGREE_MIN) -> List[str]:
    """Greedy maximal unrelated set: sample ids with no pair above threshold.

    Repeatedly drop the sample involved in the most remaining relationships.
    This is the standard greedy approximation to maximum independent set (the
    exact problem is NP-hard); it is not guaranteed optimal but it removes a
    hub — a proband with both parents and three sibs in the cohort — in one
    step instead of four.

    Default threshold is second-degree: GWAS practice keeps third-degree and
    more distant pairs, whose non-independence a PC-adjusted model tolerates.

    Ties are broken by dropping the lexicographically last sample id, so the
    returned set is a pure function of the input and reruns reproduce it
    exactly (§5.5).
    """
    ids = list(kin.sample_ids)
    n = len(ids)
    if n == 0:
        return []
    phi = np.array(kin.matrix, dtype=float, copy=True)
    phi[~np.isfinite(phi)] = 0.0       # undefined pair -> not evidence of relatedness
    related = phi > threshold
    np.fill_diagonal(related, False)

    alive = np.ones(n, dtype=bool)
    while True:
        active = related & alive[:, None] & alive[None, :]
        degree = active.sum(axis=1)
        degree[~alive] = 0
        top = int(degree.max()) if n else 0
        if top == 0:
            break
        candidates = np.nonzero(degree == top)[0]
        drop = max(candidates, key=lambda i: ids[i])
        alive[drop] = False
    return [ids[i] for i in range(n) if alive[i]]


def kinship_summary(kin: KinshipResult) -> Dict[str, Any]:
    """Kinship distribution for the §3.1 data profile.

    Counts per degree band over all pairs, not just the ones stored in
    `pairs` — the profile has to be able to say "no relatedness detected" as a
    positive finding, which requires the denominator.
    """
    n = len(kin.sample_ids)
    phi = np.asarray(kin.matrix, dtype=float)
    iu, ju = np.triu_indices(n, k=1)
    vals = phi[iu, ju]
    finite = np.isfinite(vals)
    v = vals[finite]

    counts = {
        "duplicate": int((v > DUPLICATE_MIN).sum()),
        "first-degree": int(((v > FIRST_DEGREE_MIN) & (v <= DUPLICATE_MIN)).sum()),
        "second-degree": int(((v > SECOND_DEGREE_MIN) & (v <= FIRST_DEGREE_MIN)).sum()),
        "third-degree": int(((v > THIRD_DEGREE_MIN) & (v <= SECOND_DEGREE_MIN)).sum()),
        "unrelated": int((v <= THIRD_DEGREE_MIN).sum()),
    }

    related_mask = np.zeros((n, n), dtype=bool)
    with np.errstate(invalid="ignore"):
        related_mask[iu, ju] = np.where(finite, vals > THIRD_DEGREE_MIN, False)
    related_mask |= related_mask.T
    n_with_relative = int(related_mask.any(axis=1).sum())

    return {
        "n_samples": n,
        "n_pairs": int(vals.size),
        "n_pairs_undefined": int((~finite).sum()),
        "n_variants_used": kin.n_variants_used,
        "counts": counts,
        "n_related_pairs": int(vals.size - counts["unrelated"] - (~finite).sum()),
        "n_samples_with_relative": n_with_relative,
        "max_kinship": float(v.max()) if v.size else float("nan"),
        "median_kinship": float(np.median(v)) if v.size else float("nan"),
        "has_family_structure": counts["first-degree"] > 0 or counts["second-degree"] > 0,
        "has_duplicates": counts["duplicate"] > 0,
        "thresholds": {
            "duplicate": DUPLICATE_MIN,
            "first-degree": FIRST_DEGREE_MIN,
            "second-degree": SECOND_DEGREE_MIN,
            "third-degree": THIRD_DEGREE_MIN,
        },
    }
