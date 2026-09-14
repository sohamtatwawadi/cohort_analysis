"""Mandatory pre-GWAS QC — Part II §4.2.

    "Mandatory QC before any GWAS runs: sample call rate, variant call rate,
     MAF, HWE, missingness differential between cases and controls,
     relatedness pruning, ancestry outlier removal, sex-check."

Relatedness pruning lives in kinship.py and ancestry outlier removal in
pca.py; this module covers the rest, plus the §2.3 ingestion check "Reported
vs genetic sex concordance".

Every filter reports what it removed and why. A QC step that silently drops
90% of a dataset is indistinguishable from a bug, and §5.3 requires results to
carry inspectable diagnostics, so `QCResult.filters` is an ordered log rather
than a single final mask.

Filters are applied sequentially and each one reports removals among the
variants still standing. That makes the counts additive and readable ("MAF
removed 40k of the 900k that survived call rate") instead of overlapping.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
from scipy import stats

from ..types import MISSING, GenotypeMatrix, as_float

# Pseudo-autosomal regions: males are diploid here, so these sites are
# genuinely heterozygous in males and would wreck an X-heterozygosity sex
# check. Coordinates are 1-based inclusive.
PAR_REGIONS = {
    "GRCh37": ((60001, 2699520), (154931044, 155260560)),
    "GRCh38": ((10001, 2781479), (155701383, 156030895)),
}

MAX_LISTED_IDS = 1000


@dataclass
class FilterStep:
    name: str
    description: str
    threshold: Optional[float]
    n_input: int                 # how many were still standing when this ran
    n_removed: int
    removed_ids: List[str] = field(default_factory=list)
    detail: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "threshold": self.threshold,
            "n_input": self.n_input,
            "n_removed": self.n_removed,
            "n_kept": self.n_input - self.n_removed,
            "removed_ids": list(self.removed_ids),
            "detail": dict(self.detail),
        }


@dataclass
class QCResult:
    level: str                       # "variant" | "sample"
    ids: List[str]                   # variant keys or sample ids, input order
    keep: np.ndarray                 # boolean mask aligned to `ids`
    filters: List[FilterStep] = field(default_factory=list)
    metrics: Dict[str, np.ndarray] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    params: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_input(self) -> int:
        return len(self.ids)

    @property
    def n_kept(self) -> int:
        return int(self.keep.sum())

    @property
    def n_removed(self) -> int:
        return self.n_input - self.n_kept

    def kept_ids(self) -> List[str]:
        return [i for i, k in zip(self.ids, self.keep) if k]

    def removed_ids(self) -> List[str]:
        return [i for i, k in zip(self.ids, self.keep) if not k]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "level": self.level,
            "n_input": self.n_input,
            "n_kept": self.n_kept,
            "n_removed": self.n_removed,
            "filters": [f.to_dict() for f in self.filters],
            "warnings": list(self.warnings),
            "params": dict(self.params),
        }


# ------------------------------------------------------ Hardy-Weinberg test --
@lru_cache(maxsize=200000)
def hwe_exact_p(n_het: int, n_hom1: int, n_hom2: int, midp: bool = True) -> float:
    """Hardy-Weinberg exact test (Wigginton, Cutler & Abecasis, AJHG 2005).

    Why exact and not chi-square: the chi-square goodness-of-fit approximation
    needs expected genotype counts of roughly 5 or more per cell. For a variant
    at MAF 0.5% in 2,000 samples the expected homozygote count is 0.05, and the
    approximation returns wildly anti-conservative p-values — so the test fails
    exactly on the rare variants where a genotyping artefact is most likely and
    where HWE is the only handle we have on it.

    The exact test conditions on the observed minor-allele count and enumerates
    the full distribution of heterozygote counts consistent with it, summing
    the probabilities of every configuration no more likely than the observed
    one. The recurrence below walks that distribution outward from its mode, so
    it is O(minor allele count) with no factorials to overflow.

    midp (default, as in PLINK2 --hwe midp) counts only half the probability of
    the observed configuration. The plain exact test is conservative because
    the discrete atom at the observed value is counted in full; the mid-p
    adjustment restores approximately nominal type-I error, which matters when
    a genome-wide HWE filter is applied a million times.

    Returns 1.0 for monomorphic sites — HWE is vacuous there, and they are the
    MAF filter's business.
    """
    n_het = int(n_het)
    n_hom1 = int(n_hom1)
    n_hom2 = int(n_hom2)
    if n_het < 0 or n_hom1 < 0 or n_hom2 < 0:
        raise ValueError("genotype counts must be non-negative")

    obs_hom_rare = min(n_hom1, n_hom2)
    obs_hom_common = max(n_hom1, n_hom2)
    rare = 2 * obs_hom_rare + n_het          # minor allele count
    n = n_het + obs_hom_rare + obs_hom_common
    if n == 0 or rare == 0 or rare == 2 * n:
        return 1.0

    probs = np.zeros(rare + 1, dtype=np.float64)
    # Start at the mode of the conditional distribution: E[het] under HWE,
    # nudged to the parity the minor allele count forces.
    mid = int(rare * (2 * n - rare) // (2 * n))
    if mid % 2 != rare % 2:
        mid += 1
    probs[mid] = 1.0
    total = 1.0

    curr_het = mid
    curr_hom_rare = (rare - mid) // 2
    curr_hom_common = n - curr_het - curr_hom_rare
    while curr_het >= 2:
        probs[curr_het - 2] = (probs[curr_het] * curr_het * (curr_het - 1.0)
                               / (4.0 * (curr_hom_rare + 1.0) * (curr_hom_common + 1.0)))
        total += probs[curr_het - 2]
        curr_het -= 2
        curr_hom_rare += 1
        curr_hom_common += 1

    curr_het = mid
    curr_hom_rare = (rare - mid) // 2
    curr_hom_common = n - curr_het - curr_hom_rare
    while curr_het <= rare - 2:
        probs[curr_het + 2] = (probs[curr_het] * 4.0 * curr_hom_rare * curr_hom_common
                               / ((curr_het + 2.0) * (curr_het + 1.0)))
        total += probs[curr_het + 2]
        curr_het += 2
        curr_hom_rare -= 1
        curr_hom_common -= 1

    probs /= total
    obs_prob = probs[n_het]
    if midp:
        p = float(probs[probs < obs_prob * (1.0 - 1e-7)].sum() + 0.5 * obs_prob)
    else:
        p = float(probs[probs <= obs_prob * (1.0 + 1e-7)].sum())
    return float(min(1.0, max(0.0, p)))


def genotype_counts(dosages: np.ndarray) -> np.ndarray:
    """Per-variant (n_hom_ref, n_het, n_hom_alt), missing calls excluded."""
    d = np.atleast_2d(dosages)
    return np.stack([(d == 0).sum(axis=1),
                     (d == 1).sum(axis=1),
                     (d == 2).sum(axis=1)], axis=1).astype(np.int64)


# ------------------------------------------------------------------ helpers --
def _norm_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    c = c.upper()
    return {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT"}.get(c, c)


def _is_autosome(chrom: str) -> bool:
    return _norm_chrom(chrom).isdigit()


def _in_par(build: Optional[str], pos: int) -> bool:
    regions = PAR_REGIONS.get(build or "")
    if not regions:
        return False
    return any(lo <= pos <= hi for lo, hi in regions)


def _variant_keys(gm: GenotypeMatrix) -> List[str]:
    return [v.vid or v.key for v in gm.variants]


def _case_masks(cases: Optional[Sequence[Any]], n: int):
    """Normalise a case/control vector to (is_case, is_control) boolean masks.

    Accepts booleans, 0/1, or anything coercible with NaN for unknown; unknown
    status makes a sample neither a case nor a control rather than silently a
    control.
    """
    if cases is None:
        return None, None
    arr = np.asarray(cases)
    if arr.shape[0] != n:
        raise ValueError("cases has length {} but there are {} samples".format(
            arr.shape[0], n))
    if arr.dtype == bool:
        return arr.copy(), ~arr
    v = as_float(arr)
    return (v == 1) & np.isfinite(v), (v == 0) & np.isfinite(v)


def _chi2_2x2(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Vectorised 2x2 chi-square (1 df, no continuity correction)."""
    a = a.astype(np.float64)
    b = b.astype(np.float64)
    c = c.astype(np.float64)
    d = d.astype(np.float64)
    n = a + b + c + d
    num = n * (a * d - b * c) ** 2
    den = (a + b) * (c + d) * (a + c) * (b + d)
    with np.errstate(invalid="ignore", divide="ignore"):
        chi2 = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)
    return stats.chi2.sf(chi2, 1)


# ------------------------------------------------------------- variant QC --
def variant_qc(gm: GenotypeMatrix,
               min_call_rate: float = 0.95,
               min_maf: float = 0.01,
               hwe_p: float = 1e-6,
               cases: Optional[Sequence[Any]] = None,
               diff_miss_p: float = 1e-5,
               hwe_midp: bool = True) -> QCResult:
    """Variant-level GWAS QC (§4.2). Returns a keep mask plus a per-filter log.

    `cases` is an optional per-sample case/control vector. Supplying it changes
    two things: HWE is computed in controls only, and the differential
    missingness test becomes available.
    """
    ids = _variant_keys(gm)
    m = gm.n_variants
    keep = np.ones(m, dtype=bool)
    warnings: List[str] = []
    filters: List[FilterStep] = []
    metrics: Dict[str, np.ndarray] = {}

    is_case, is_control = _case_masks(cases, gm.n_samples)

    def _record(name, description, threshold, failed, detail=None):
        failed = failed & keep
        step = FilterStep(
            name=name, description=description, threshold=threshold,
            n_input=int(keep.sum()), n_removed=int(failed.sum()),
            removed_ids=[ids[i] for i in np.nonzero(failed)[0][:MAX_LISTED_IDS]],
            detail=detail or {})
        if step.n_removed > MAX_LISTED_IDS:
            step.detail["removed_ids_truncated"] = True
        filters.append(step)
        keep[failed] = False

    # 1. Call rate ----------------------------------------------------------
    call_rate = gm.call_rate_per_variant()
    metrics["call_rate"] = call_rate
    _record("variant_call_rate",
            "fraction of samples with a non-missing call",
            min_call_rate, call_rate < min_call_rate)

    # 2. MAF ----------------------------------------------------------------
    # Low-MAF variants are not wrong, they are underpowered: a handful of minor
    # allele carriers gives an unstable effect estimate and a p-value driven by
    # one or two people. Rare variants belong in the gene-based tests (§4.4).
    maf = gm.maf()
    metrics["maf"] = maf
    _record("maf", "minor allele frequency", min_maf,
            ~np.isfinite(maf) | (maf < min_maf))

    # 3. HWE ----------------------------------------------------------------
    # In controls only when status is known: a variant that genuinely raises
    # disease risk is *expected* to deviate from HWE in cases (the cases are an
    # ascertained, non-random sample of genotypes), so testing the combined set
    # throws away true positives. HWE here is a genotyping-artefact detector —
    # cluster-plot failures and copy-number artefacts show up as heterozygote
    # excess or deficit — not a biological finding.
    #
    # Autosomes only: males are hemizygous on non-PAR X, so an X site has no
    # Hardy-Weinberg expectation across a mixed-sex cohort.
    autosomal = np.array([_is_autosome(v.chrom) for v in gm.variants], dtype=bool)
    if is_control is not None and is_control.sum() > 0:
        hwe_dos = gm.dosages[:, is_control]
        hwe_group = "controls"
        n_hwe_samples = int(is_control.sum())
    else:
        hwe_dos = gm.dosages
        hwe_group = "all samples"
        n_hwe_samples = gm.n_samples
        if cases is not None:
            warnings.append(
                "case/control status supplied but no controls found; HWE fell "
                "back to all samples, which can remove true associations")
        else:
            warnings.append(
                "HWE computed over all samples because no case/control status "
                "was supplied; disease-associated variants legitimately "
                "deviate in cases and may be removed")

    counts = genotype_counts(hwe_dos)
    hwe_pvals = np.full(m, np.nan)
    for i in range(m):
        if not autosomal[i] or not keep[i]:
            continue
        hwe_pvals[i] = hwe_exact_p(int(counts[i, 1]), int(counts[i, 0]),
                                   int(counts[i, 2]), hwe_midp)
    metrics["hwe_p"] = hwe_pvals
    _record("hwe",
            "Hardy-Weinberg exact test ({}, mid-p={}) on autosomes".format(
                hwe_group, hwe_midp),
            hwe_p,
            np.isfinite(hwe_pvals) & (hwe_pvals < hwe_p),
            {"group": hwe_group, "n_samples_tested": n_hwe_samples,
             "n_skipped_non_autosomal": int((~autosomal).sum()),
             "test": "Wigginton et al. 2005 exact"})

    # 4. Differential missingness ------------------------------------------
    # Named explicitly in §4.2. If a variant fails to call more often in cases
    # than controls — a plate effect, a DNA-quality difference, cases and
    # controls genotyped in different batches — the missingness itself is
    # associated with phenotype, and whatever calls do survive carry that
    # association. This is a classic source of genome-wide-significant
    # nonsense that survives every other filter.
    if is_case is not None and is_case.sum() > 0 and is_control.sum() > 0:
        miss = gm.dosages == MISSING
        miss_case = miss[:, is_case].sum(axis=1)
        call_case = int(is_case.sum()) - miss_case
        miss_ctrl = miss[:, is_control].sum(axis=1)
        call_ctrl = int(is_control.sum()) - miss_ctrl

        pvals = _chi2_2x2(miss_case, call_case, miss_ctrl, call_ctrl)
        # Fisher where the chi-square approximation does not hold. Missingness
        # tables are usually sparse in exactly this way (a handful of no-calls
        # against thousands of successful ones), so the fallback is the norm,
        # not the exception — but it is slow, so it is applied only where
        # needed.
        # Smallest expected cell = min(row marginal) * min(col marginal) / n.
        n_tot = float(int(is_case.sum()) + int(is_control.sum()))
        n_missing = (miss_case + miss_ctrl).astype(float)
        n_called = (call_case + call_ctrl).astype(float)
        min_status = float(min(int(is_case.sum()), int(is_control.sum())))
        e_min = np.minimum(n_missing, n_called) * min_status / max(n_tot, 1.0)
        need_fisher = np.nonzero(keep & (e_min < 5.0) & (n_missing > 0))[0]
        for i in need_fisher:
            table = [[int(miss_case[i]), int(call_case[i])],
                     [int(miss_ctrl[i]), int(call_ctrl[i])]]
            pvals[i] = float(stats.fisher_exact(table)[1])
        pvals[n_missing == 0] = 1.0     # nothing missing anywhere: nothing to test
        metrics["diff_missingness_p"] = pvals
        _record("differential_missingness",
                "case/control difference in call rate (chi-square, Fisher "
                "exact for sparse tables)",
                diff_miss_p,
                np.isfinite(pvals) & (pvals < diff_miss_p),
                {"n_cases": int(is_case.sum()), "n_controls": int(is_control.sum()),
                 "n_fisher": int(need_fisher.size)})
    else:
        if cases is None:
            warnings.append(
                "differential missingness between cases and controls was not "
                "tested because no case/control status was supplied (§4.2 "
                "requires it before a GWAS)")
        else:
            warnings.append(
                "differential missingness not tested: need both cases and "
                "controls")

    if keep.sum() == 0:
        warnings.append("every variant was removed by QC")

    return QCResult(level="variant", ids=ids, keep=keep, filters=filters,
                    metrics=metrics, warnings=warnings,
                    params={"min_call_rate": min_call_rate, "min_maf": min_maf,
                            "hwe_p": hwe_p, "hwe_midp": hwe_midp,
                            "diff_miss_p": diff_miss_p,
                            "case_control_supplied": cases is not None})


# -------------------------------------------------------------- sex check --
def _x_variant_mask(gm: GenotypeMatrix) -> np.ndarray:
    """Non-PAR X variants — the only ones where male hemizygosity is expected."""
    out = np.zeros(gm.n_variants, dtype=bool)
    for i, v in enumerate(gm.variants):
        if _norm_chrom(v.chrom) != "X":
            continue
        if _in_par(gm.build, int(v.pos)):
            continue
        out[i] = True
    return out


def normalise_sex(value: Any) -> str:
    """Map whatever the phenotype file used onto male/female/unknown."""
    if value is None:
        return "unknown"
    if isinstance(value, (bool, np.bool_)):
        return "unknown"
    if isinstance(value, (int, float, np.integer, np.floating)):
        if not np.isfinite(float(value)):
            return "unknown"
        return {1: "male", 2: "female"}.get(int(value), "unknown")
    s = str(value).strip().lower()
    if s in ("m", "male", "1"):
        return "male"
    if s in ("f", "female", "2"):
        return "female"
    return "unknown"


def sex_check(gm: GenotypeMatrix,
              reported_sex: Optional[Sequence[Any]] = None,
              male_max_x_het: float = 0.10,
              female_min_x_het: float = 0.20,
              min_x_variants: int = 20) -> Dict[str, Any]:
    """Infer genetic sex from non-PAR X heterozygosity and compare to report.

    §2.3 makes reported-vs-genetic sex concordance a mandatory ingestion check
    and §4.2 repeats it as GWAS QC. A discordance is rarely an interesting
    biological finding — it is a sample swap, a mislabelled plate well, or a
    pedigree transcription error, and it means at least one other sample in the
    dataset is also wrong.

    Males have one X, so a non-PAR X call can only be hemizygous; any
    heterozygous call there is a genotyping error. Females sit at whatever the
    marker set's heterozygosity is (typically 0.2-0.4). The gap between the two
    thresholds is deliberate: samples landing in it are reported "ambiguous"
    rather than forced into a call, because that band is where sample
    contamination and sex-chromosome aneuploidy live.
    """
    n = gm.n_samples
    xmask = _x_variant_mask(gm)
    n_x = int(xmask.sum())
    het_rate = np.full(n, np.nan)
    n_called = np.zeros(n, dtype=int)
    inferred = np.array(["unknown"] * n, dtype=object)
    warnings: List[str] = []

    if n_x == 0:
        warnings.append(
            "no non-PAR X variants present; genetic sex cannot be inferred "
            "and the §2.3 sex-concordance check is unmet")
    else:
        if gm.build is None:
            warnings.append(
                "genome build unknown, so pseudo-autosomal X regions could not "
                "be excluded; male X heterozygosity may be overstated")
        xd = gm.dosages[xmask, :]
        called = xd != MISSING
        n_called = called.sum(axis=0).astype(int)
        het = (xd == 1).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            het_rate = np.where(n_called > 0, het / np.maximum(n_called, 1), np.nan)
        if n_x < min_x_variants:
            warnings.append(
                "only {} non-PAR X variants; sex inference from this few "
                "markers is unreliable".format(n_x))
        for i in range(n):
            if n_called[i] < min_x_variants or not np.isfinite(het_rate[i]):
                inferred[i] = "unknown"
            elif het_rate[i] <= male_max_x_het:
                inferred[i] = "male"
            elif het_rate[i] >= female_min_x_het:
                inferred[i] = "female"
            else:
                inferred[i] = "ambiguous"

    reported = np.array(["unknown"] * n, dtype=object)
    discordant: List[Dict[str, Any]] = []
    if reported_sex is not None:
        rs = list(reported_sex)
        if len(rs) != n:
            raise ValueError("reported_sex has length {} but there are {} "
                             "samples".format(len(rs), n))
        reported = np.array([normalise_sex(x) for x in rs], dtype=object)
        for i in range(n):
            if reported[i] in ("male", "female") and inferred[i] in ("male", "female") \
                    and reported[i] != inferred[i]:
                discordant.append({
                    "sample_id": gm.sample_ids[i],
                    "reported": reported[i],
                    "inferred": inferred[i],
                    "x_het_rate": float(het_rate[i]),
                    "n_x_called": int(n_called[i]),
                })
    else:
        warnings.append(
            "no reported sex supplied; genetic sex was inferred but the "
            "concordance check (§2.3) could not run")

    return {
        "sample_ids": list(gm.sample_ids),
        "inferred_sex": inferred,
        "reported_sex": reported,
        "x_het_rate": het_rate,
        "n_x_called": n_called,
        "n_x_variants": n_x,
        "discordant": discordant,
        "discordant_ids": [d["sample_id"] for d in discordant],
        "n_ambiguous": int((inferred == "ambiguous").sum()),
        "warnings": warnings,
        "thresholds": {"male_max_x_het": male_max_x_het,
                       "female_min_x_het": female_min_x_het,
                       "min_x_variants": min_x_variants},
    }


# --------------------------------------------------------------- sample QC --
def sample_qc(gm: GenotypeMatrix,
              min_call_rate: float = 0.95,
              reported_sex: Optional[Sequence[Any]] = None,
              drop_sex_discordant: bool = True,
              male_max_x_het: float = 0.10,
              female_min_x_het: float = 0.20) -> QCResult:
    """Sample-level GWAS QC (§4.2): call rate and sex concordance.

    A low sample call rate is a DNA-quality signal that correlates with
    everything else that went wrong for that sample, so these are dropped
    before variant-level statistics are computed on the survivors.

    Ancestry outliers (pca.ancestry_outliers) and relatedness
    (kinship.unrelated_set) are the other two mandatory sample-level steps;
    they need their own inputs and are applied by the caller.
    """
    ids = list(gm.sample_ids)
    n = gm.n_samples
    keep = np.ones(n, dtype=bool)
    filters: List[FilterStep] = []
    warnings: List[str] = []
    metrics: Dict[str, np.ndarray] = {}

    def _record(name, description, threshold, failed, detail=None):
        failed = failed & keep
        step = FilterStep(
            name=name, description=description, threshold=threshold,
            n_input=int(keep.sum()), n_removed=int(failed.sum()),
            removed_ids=[ids[i] for i in np.nonzero(failed)[0][:MAX_LISTED_IDS]],
            detail=detail or {})
        filters.append(step)
        keep[failed] = False

    call_rate = gm.call_rate_per_sample()
    metrics["call_rate"] = call_rate
    _record("sample_call_rate", "fraction of variants with a non-missing call",
            min_call_rate, call_rate < min_call_rate)

    # Autosomal heterozygosity is reported but not filtered on: extreme values
    # flag contamination (high) or inbreeding/consanguinity (low), and the
    # second of those is a real population feature, not a defect, so the
    # threshold is a judgement call left to the analyst.
    autosomal = np.array([_is_autosome(v.chrom) for v in gm.variants], dtype=bool)
    if autosomal.any():
        ad = gm.dosages[autosomal, :]
        a_called = (ad != MISSING).sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            metrics["autosomal_het_rate"] = np.where(
                a_called > 0, (ad == 1).sum(axis=0) / np.maximum(a_called, 1), np.nan)
    else:
        metrics["autosomal_het_rate"] = np.full(n, np.nan)

    sx = sex_check(gm, reported_sex=reported_sex,
                   male_max_x_het=male_max_x_het,
                   female_min_x_het=female_min_x_het)
    metrics["x_het_rate"] = sx["x_het_rate"]
    metrics["inferred_sex"] = sx["inferred_sex"]
    metrics["reported_sex"] = sx["reported_sex"]
    warnings.extend(sx["warnings"])

    if reported_sex is not None and sx["n_x_variants"] > 0:
        disc_ids = set(sx["discordant_ids"])
        failed = np.array([s in disc_ids for s in ids], dtype=bool)
        if not drop_sex_discordant:
            failed = np.zeros(n, dtype=bool)
        _record("sex_check",
                "reported vs genetic sex (non-PAR X heterozygosity)",
                None, failed,
                {"discordant": sx["discordant"],
                 "n_discordant": len(sx["discordant"]),
                 "n_ambiguous": sx["n_ambiguous"],
                 "n_x_variants": sx["n_x_variants"],
                 "dropped": bool(drop_sex_discordant)})
        if sx["discordant"] and not drop_sex_discordant:
            warnings.append(
                "{} sex-discordant samples were flagged but kept".format(
                    len(sx["discordant"])))
        if sx["n_ambiguous"]:
            warnings.append(
                "{} samples have ambiguous X heterozygosity (possible "
                "contamination or sex-chromosome aneuploidy) and were left "
                "in".format(sx["n_ambiguous"]))

    if keep.sum() == 0:
        warnings.append("every sample was removed by QC")

    return QCResult(level="sample", ids=ids, keep=keep, filters=filters,
                    metrics=metrics, warnings=warnings,
                    params={"min_call_rate": min_call_rate,
                            "reported_sex_supplied": reported_sex is not None,
                            "drop_sex_discordant": drop_sex_discordant,
                            "male_max_x_het": male_max_x_het,
                            "female_min_x_het": female_min_x_het})
