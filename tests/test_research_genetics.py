"""Research-Mode population genetics: PCA, KING kinship, GWAS QC.

Every assertion here is against constructed ground truth — two populations
simulated with known allele-frequency divergence, a parent-offspring pair built
by transmitting one haplotype from each parent, an exact duplicate sample, a
hand-enumerable Hardy-Weinberg null, males made hemizygous on X by
construction. Nothing is checked against the implementation's own output.
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial
from typing import List, Optional

import numpy as np
import pytest

from backend.app.research.types import GenotypeMatrix, Variant
from backend.app.research.stats import kinship as kin_mod
from backend.app.research.stats import pca as pca_mod
from backend.app.research.stats import qc as qc_mod


# ------------------------------------------------------------- construction --
def make_gm(dosages: np.ndarray,
            chrom: str = "1",
            sample_ids: Optional[List[str]] = None,
            build: str = "GRCh38",
            start_pos: int = 1000000,
            spacing: int = 5000) -> GenotypeMatrix:
    dosages = np.asarray(dosages, dtype=np.int8)
    m, n = dosages.shape
    variants = [Variant(chrom=chrom, pos=start_pos + i * spacing, ref="A", alt="G",
                        vid="{}:{}".format(chrom, start_pos + i * spacing))
                for i in range(m)]
    if sample_ids is None:
        sample_ids = ["S{:03d}".format(i) for i in range(n)]
    return GenotypeMatrix(sample_ids=list(sample_ids), variants=variants,
                          dosages=dosages, build=build, source_format="test")


def stack_gm(a: GenotypeMatrix, b: GenotypeMatrix) -> GenotypeMatrix:
    """Concatenate two matrices over variants (same samples)."""
    return GenotypeMatrix(
        sample_ids=list(a.sample_ids),
        variants=list(a.variants) + list(b.variants),
        dosages=np.vstack([a.dosages, b.dosages]).astype(np.int8),
        build=a.build, source_format=a.source_format)


# ====================================================================== PCA ==
def simulate_two_populations(rng, n_variants=800, n_per_pop=40, fst=0.15):
    """Balding-Nichols: an ancestral frequency drifted independently in two
    populations. Fst 0.15 is roughly a continental-scale split."""
    p_anc = rng.uniform(0.15, 0.85, n_variants)
    a = p_anc * (1.0 - fst) / fst
    b = (1.0 - p_anc) * (1.0 - fst) / fst
    p1 = rng.beta(a, b)
    p2 = rng.beta(a, b)
    g1 = rng.binomial(2, np.repeat(p1[:, None], n_per_pop, axis=1))
    g2 = rng.binomial(2, np.repeat(p2[:, None], n_per_pop, axis=1))
    dos = np.concatenate([g1, g2], axis=1).astype(np.int8)
    labels = np.array([0] * n_per_pop + [1] * n_per_pop)
    return make_gm(dos), labels


def test_pca_separates_two_populations():
    rng = np.random.default_rng(20260912)
    gm, labels = simulate_two_populations(rng)
    res = pca_mod.compute_pca(gm, n_components=5, prune=False)

    assert res.components.shape == (gm.n_samples, 5)
    assert res.sample_ids == gm.sample_ids
    assert res.n_variants_used == gm.n_variants

    pc1 = res.components[:, 0]
    a, b = pc1[labels == 0], pc1[labels == 1]
    # The standardisation centres every variant on the pooled mean, so the two
    # group means must straddle zero if PC1 is the ancestry axis.
    assert np.sign(a.mean()) != np.sign(b.mean())
    gap = abs(a.mean() - b.mean())
    assert gap > 3.0 * (a.std() + b.std()), (
        "PC1 gap {:.2f} did not exceed within-group spread".format(gap))
    # No other PC should separate the populations anywhere near as well.
    for j in range(1, 5):
        other = res.components[:, j]
        sep = abs(other[labels == 0].mean() - other[labels == 1].mean())
        assert sep < gap / 3.0


def test_pca_explained_variance_is_descending_and_bounded():
    rng = np.random.default_rng(7)
    gm, _ = simulate_two_populations(rng, n_variants=400, n_per_pop=25)
    res = pca_mod.compute_pca(gm, n_components=6, prune=False)
    evr = res.explained_variance_ratio
    assert evr.shape == (6,)
    assert np.all(np.diff(evr) <= 1e-12)
    assert np.all(evr >= 0) and evr.sum() <= 1.0 + 1e-9
    # Population structure means PC1 carries appreciably more than an even share.
    assert evr[0] > evr[1] * 2


def test_pca_drops_monomorphic_and_warns_on_few_variants():
    rng = np.random.default_rng(3)
    gm, _ = simulate_two_populations(rng, n_variants=100, n_per_pop=20)
    mono = np.zeros((5, gm.n_samples), dtype=np.int8)      # p == 0 everywhere
    gm2 = stack_gm(gm, make_gm(mono, start_pos=90000000))
    res = pca_mod.compute_pca(gm2, n_components=3, prune=False)
    assert res.n_variants_used == 100                      # the 5 monomorphic are gone
    assert any("monomorphic" in w for w in res.warnings)
    assert any("noise-dominated" in w or "only 100 variants" in w
               for w in res.warnings)
    assert res.variant_mask is not None and res.variant_mask.sum() == 100


def test_pca_handles_missing_calls():
    rng = np.random.default_rng(11)
    gm, labels = simulate_two_populations(rng, n_variants=600, n_per_pop=30)
    dos = gm.dosages.copy()
    holes = rng.random(dos.shape) < 0.03
    dos[holes] = -1
    gm_missing = make_gm(dos)
    res = pca_mod.compute_pca(gm_missing, n_components=3, prune=False)
    pc1 = res.components[:, 0]
    assert np.sign(pc1[labels == 0].mean()) != np.sign(pc1[labels == 1].mean())
    assert np.isfinite(res.components).all()


def test_ancestry_outliers_flags_extreme_samples():
    rng = np.random.default_rng(5)
    pcs = rng.normal(size=(200, 4))
    pcs[17, 0] = 60.0                    # unmistakable outlier on PC1
    pcs[42, 2] = -55.0                   # and one on PC3
    mask = pca_mod.ancestry_outliers(pcs, n_sd=6.0)
    assert mask[17] and mask[42]
    assert mask.sum() == 2
    # Only the first n_pcs axes are considered.
    assert not pca_mod.ancestry_outliers(pcs, n_sd=6.0, n_pcs=1)[42]


def test_ancestry_outliers_accepts_a_pca_result():
    rng = np.random.default_rng(6)
    gm, _ = simulate_two_populations(rng, n_variants=300, n_per_pop=20)
    res = pca_mod.compute_pca(gm, n_components=3, prune=False)
    mask = pca_mod.ancestry_outliers(res, n_sd=6.0)
    assert mask.shape == (gm.n_samples,)


# =============================================================== LD pruning ==
def test_ld_prune_removes_a_perfectly_correlated_duplicate():
    rng = np.random.default_rng(99)
    m, n = 30, 200
    dos = rng.binomial(2, 0.3, size=(m, n)).astype(np.int8)
    dos[4] = dos[3]                                  # exact duplicate of variant 3
    gm = make_gm(dos)
    keep = pca_mod.ld_prune(gm, window=10, step=5, r2_threshold=0.2)
    assert keep[3] and not keep[4], "the later member of the r2 == 1 pair must go"
    # Independent random variants should overwhelmingly survive at n = 200.
    assert keep.sum() >= m - 3


def test_ld_prune_window_limits_comparisons():
    """A duplicate further apart than the window is out of reach by design."""
    rng = np.random.default_rng(100)
    dos = rng.binomial(2, 0.3, size=(60, 300)).astype(np.int8)
    dos[50] = dos[0]
    gm = make_gm(dos)
    keep = pca_mod.ld_prune(gm, window=10, step=10, r2_threshold=0.2)
    assert keep[0] and keep[50]
    keep_wide = pca_mod.ld_prune(gm, window=60, step=5, r2_threshold=0.2)
    assert keep_wide[0] and not keep_wide[50]


def test_ld_prune_is_deterministic():
    rng = np.random.default_rng(101)
    gm = make_gm(rng.binomial(2, 0.25, size=(80, 150)).astype(np.int8))
    a = pca_mod.ld_prune(gm)
    b = pca_mod.ld_prune(gm)
    assert np.array_equal(a, b)


# ================================================================== kinship ==
def simulate_family(rng, n_variants=4000, n_founders=10):
    """Unrelated founders under HWE, plus one offspring of founders 0 and 1
    built by transmitting one haplotype from each parent, plus an exact
    duplicate of founder 5."""
    p = rng.uniform(0.1, 0.5, n_variants)
    # haplotypes[v, f, c] -> allele count 0/1 for founder f, chromosome copy c
    haps = (rng.random((n_variants, n_founders, 2)) < p[:, None, None]).astype(np.int8)
    founders = haps.sum(axis=2)                              # (m, n_founders)

    pick_mum = rng.integers(0, 2, n_variants)
    pick_dad = rng.integers(0, 2, n_variants)
    rows = np.arange(n_variants)
    child = (haps[rows, 0, pick_mum] + haps[rows, 1, pick_dad]).astype(np.int8)

    dos = np.concatenate([founders, child[:, None], founders[:, 5][:, None]], axis=1)
    ids = ["F{:02d}".format(i) for i in range(n_founders)] + ["CHILD", "F05_DUP"]
    return make_gm(dos.astype(np.int8), sample_ids=ids), p


def test_king_recovers_parent_offspring_duplicate_and_unrelated():
    rng = np.random.default_rng(424242)
    gm, _ = simulate_family(rng)
    # prune=False: LD pruning on independently simulated markers in a 12-sample
    # cohort would drop markers on sampling noise alone.
    kin = kin_mod.king_robust(gm, prune=False)

    idx = kin.sample_index()
    assert kin.matrix.shape == (12, 12)
    assert np.allclose(np.diag(kin.matrix), 0.5)
    assert np.allclose(kin.matrix, kin.matrix.T)

    for parent in ("F00", "F01"):
        phi = kin.matrix[idx[parent], idx["CHILD"]]
        assert 0.177 < phi <= 0.354, "{}-CHILD phi = {:.4f}".format(parent, phi)
        assert abs(phi - 0.25) < 0.04
        assert kin_mod.classify_degree(phi) == "first-degree"

    # An identical sample has N_AAaa == 0 and identical het counts, so the
    # estimator is exactly 0.5 by construction.
    dup = kin.matrix[idx["F05"], idx["F05_DUP"]]
    assert abs(dup - 0.5) < 1e-12
    assert kin_mod.classify_degree(dup) == "duplicate"

    # Every founder-founder pair is unrelated by construction.
    founders = ["F{:02d}".format(i) for i in range(10)]
    for a in range(10):
        for b in range(a + 1, 10):
            phi = kin.matrix[idx[founders[a]], idx[founders[b]]]
            assert abs(phi) < 0.0442, "{}-{} phi = {:.4f}".format(
                founders[a], founders[b], phi)
            assert kin_mod.classify_degree(phi) == "unrelated"


def test_king_pairs_and_summary():
    rng = np.random.default_rng(424242)
    gm, _ = simulate_family(rng)
    kin = kin_mod.king_robust(gm, prune=False)

    reported = {(p["sample1"], p["sample2"]): p["degree"] for p in kin.pairs}
    assert reported == {
        ("F00", "CHILD"): "first-degree",
        ("F01", "CHILD"): "first-degree",
        ("F05", "F05_DUP"): "duplicate",
    }
    # Sorted strongest-first for a reproducible report.
    assert [p["kinship"] for p in kin.pairs] == sorted(
        [p["kinship"] for p in kin.pairs], reverse=True)
    assert kin.pairs[0]["n_variants"] == gm.n_variants

    summary = kin_mod.kinship_summary(kin)
    assert summary["n_samples"] == 12
    assert summary["n_pairs"] == 66
    assert summary["counts"]["duplicate"] == 1
    assert summary["counts"]["first-degree"] == 2
    assert summary["counts"]["second-degree"] == 0
    assert summary["counts"]["third-degree"] == 0
    assert summary["counts"]["unrelated"] == 63
    assert summary["n_related_pairs"] == 3
    assert summary["n_samples_with_relative"] == 5     # F00 F01 CHILD F05 DUP
    assert summary["has_family_structure"] and summary["has_duplicates"]
    assert abs(summary["max_kinship"] - 0.5) < 1e-9


def test_king_handles_missing_calls():
    rng = np.random.default_rng(424242)
    gm, _ = simulate_family(rng)
    dos = gm.dosages.copy()
    dos[rng.random(dos.shape) < 0.05] = -1
    gm2 = make_gm(dos, sample_ids=gm.sample_ids)
    kin = kin_mod.king_robust(gm2, prune=False)
    idx = kin.sample_index()
    assert abs(kin.matrix[idx["F00"], idx["CHILD"]] - 0.25) < 0.05
    assert kin.pairs[0]["n_variants"] < gm.n_variants     # shared-call denominator


def test_unrelated_set_drops_one_of_each_related_pair_and_is_deterministic():
    rng = np.random.default_rng(424242)
    gm, _ = simulate_family(rng)
    kin = kin_mod.king_robust(gm, prune=False)

    keep = kin_mod.unrelated_set(kin, threshold=0.0884)
    assert kin_mod.unrelated_set(kin, threshold=0.0884) == keep      # reproducible

    idx = kin.sample_index()
    kept = set(keep)
    # The child is the hub (related to two founders), so greedy removes it
    # first rather than removing both parents.
    assert "CHILD" not in kept
    assert {"F00", "F01"} <= kept
    # Exactly one of the duplicate pair survives; the tie-break drops the
    # lexicographically last id.
    assert ("F05" in kept) and ("F05_DUP" not in kept)
    assert len(keep) == 10

    # Nothing above threshold remains among the survivors.
    rows = [idx[s] for s in keep]
    sub = kin.matrix[np.ix_(rows, rows)]
    np.fill_diagonal(sub, 0.0)
    assert sub.max() <= 0.0884


def test_unrelated_set_on_a_fully_unrelated_cohort_keeps_everyone():
    rng = np.random.default_rng(8)
    gm = make_gm(rng.binomial(2, 0.3, size=(3000, 15)).astype(np.int8))
    kin = kin_mod.king_robust(gm, prune=False)
    assert kin_mod.unrelated_set(kin) == list(gm.sample_ids)
    assert kin.pairs == []


# ============================================================ HWE exact test ==
def brute_force_hwe_p(n_het, n_hom1, n_hom2):
    """Independent reference: enumerate the exact conditional distribution.

        P(n_Aa = h | n, n_A) = 2^h * n! / (nAA! nAa! naa!)
                               * nA! nB! / (2n)!

    Computed in exact rationals, so this is ground truth, not a second
    approximation. Two-sided p = total probability of every configuration no
    more likely than the observed one.
    """
    n = n_het + n_hom1 + n_hom2
    n_a = 2 * n_hom1 + n_het          # copies of allele A
    n_b = 2 * n - n_a
    const = Fraction(factorial(n_a) * factorial(n_b), factorial(2 * n))

    probs = {}
    for h in range(0, min(n_a, n_b) + 1):
        if (n_a - h) % 2:
            continue
        hom_a = (n_a - h) // 2
        hom_b = (n_b - h) // 2
        if hom_a < 0 or hom_b < 0 or hom_a + hom_b + h != n:
            continue
        probs[h] = (Fraction(2 ** h) * Fraction(
            factorial(n), factorial(hom_a) * factorial(h) * factorial(hom_b)) * const)
    assert abs(float(sum(probs.values())) - 1.0) < 1e-12
    obs = probs[n_het]
    return float(sum(v for v in probs.values() if v <= obs))


@pytest.mark.parametrize("counts", [
    (10, 5, 5), (0, 10, 10), (20, 0, 0), (3, 14, 3), (1, 18, 1),
    (7, 2, 11), (12, 4, 4), (5, 1, 14),
])
def test_hwe_exact_matches_exhaustive_enumeration(counts):
    n_het, n_hom1, n_hom2 = counts
    got = qc_mod.hwe_exact_p(n_het, n_hom1, n_hom2, midp=False)
    assert got == pytest.approx(brute_force_hwe_p(n_het, n_hom1, n_hom2), abs=1e-10)


def test_hwe_exact_flags_heterozygote_excess_and_passes_a_perfect_fit():
    # p = 0.5, 1000 samples, textbook expectation 250/500/250.
    assert qc_mod.hwe_exact_p(500, 250, 250) > 0.9
    # Every single sample heterozygous: impossible under HWE.
    assert qc_mod.hwe_exact_p(1000, 0, 0) < 1e-100
    # Heterozygote deficit (the null-allele / probe-dropout signature).
    assert qc_mod.hwe_exact_p(10, 495, 495) < 1e-100
    # Mild, tolerable departure.
    assert qc_mod.hwe_exact_p(480, 260, 260) > 0.05


def test_hwe_exact_is_symmetric_and_handles_edges():
    assert qc_mod.hwe_exact_p(30, 10, 60) == qc_mod.hwe_exact_p(30, 60, 10)
    assert qc_mod.hwe_exact_p(0, 100, 0) == 1.0          # monomorphic
    assert qc_mod.hwe_exact_p(0, 0, 0) == 1.0            # nothing called
    # mid-p is less conservative than the plain exact test.
    assert qc_mod.hwe_exact_p(57, 14, 929, midp=True) < \
        qc_mod.hwe_exact_p(57, 14, 929, midp=False)


def test_hwe_exact_beats_chi_square_on_a_rare_variant():
    """The reason the exact test is mandatory: chi-square on a rare variant
    with a single heterozygote-excess artefact is anti-conservative."""
    from scipy import stats as sp
    n_hom_ref, n_het, n_hom_alt = 1990, 5, 5
    n = n_hom_ref + n_het + n_hom_alt
    p = (2 * n_hom_alt + n_het) / (2.0 * n)
    exp = np.array([(1 - p) ** 2, 2 * p * (1 - p), p ** 2]) * n
    obs = np.array([n_hom_ref, n_het, n_hom_alt], dtype=float)
    chi2 = float(((obs - exp) ** 2 / exp).sum())
    chi2_p = float(sp.chi2.sf(chi2, 1))
    exact_p = qc_mod.hwe_exact_p(n_het, n_hom_ref, n_hom_alt, midp=False)
    assert exp.min() < 5, "this case is supposed to violate the chi-square assumption"
    assert chi2_p < exact_p          # chi-square overstates the evidence
    assert exact_p == pytest.approx(brute_force_hwe_p(n_het, n_hom_ref, n_hom_alt),
                                    rel=1e-9)


# ================================================================ variant QC ==
def build_qc_matrix():
    """Six variants, each engineered to fail exactly one filter.

    0 healthy · 1 low call rate · 2 low MAF · 3 HWE violation in controls
    4 differential missingness · 5 healthy
    """
    rng = np.random.default_rng(31415)
    n_case, n_ctrl = 200, 200
    n = n_case + n_ctrl
    cases = np.array([1] * n_case + [0] * n_ctrl)

    def hwe_row(p):
        return rng.binomial(2, p, size=n).astype(np.int8)

    rows = [hwe_row(0.3)]                               # 0: healthy

    bad_call = hwe_row(0.3)                             # 1: 20% missing, evenly
    miss_idx = rng.permutation(n)[:int(0.2 * n)]
    bad_call[miss_idx] = -1
    rows.append(bad_call)

    rare = np.zeros(n, dtype=np.int8)                   # 2: MAF = 2/800
    rare[[3, 300]] = 1
    rows.append(rare)

    het_all = np.ones(n, dtype=np.int8)                 # 3: everyone het
    rows.append(het_all)

    # 4: 19 no-calls, all in cases. Overall call rate 95.25% so it sails
    # through the call-rate filter — differential missingness is the only
    # filter that can catch it, which is exactly why §4.2 names it.
    diff_miss = hwe_row(0.3)
    diff_miss[rng.permutation(n_case)[:19]] = -1
    rows.append(diff_miss)

    rows.append(hwe_row(0.4))                           # 5: healthy

    dos = np.vstack(rows).astype(np.int8)
    return make_gm(dos), cases


def test_variant_qc_attributes_each_removal_to_the_right_filter():
    gm, cases = build_qc_matrix()
    res = qc_mod.variant_qc(gm, min_call_rate=0.95, min_maf=0.01,
                            hwe_p=1e-6, cases=cases, diff_miss_p=1e-5)

    ids = res.ids
    by_name = {f.name: f for f in res.filters}
    assert list(by_name) == ["variant_call_rate", "maf", "hwe",
                             "differential_missingness"]

    assert by_name["variant_call_rate"].n_removed == 1
    assert by_name["variant_call_rate"].removed_ids == [ids[1]]
    assert by_name["maf"].n_removed == 1
    assert by_name["maf"].removed_ids == [ids[2]]
    assert by_name["hwe"].n_removed == 1
    assert by_name["hwe"].removed_ids == [ids[3]]
    assert by_name["differential_missingness"].n_removed == 1
    assert by_name["differential_missingness"].removed_ids == [ids[4]]

    assert res.n_input == 6
    assert res.n_kept == 2
    assert res.kept_ids() == [ids[0], ids[5]]
    # Sequential application: each filter reports the pool it actually saw.
    assert by_name["variant_call_rate"].n_input == 6
    assert by_name["maf"].n_input == 5
    assert by_name["hwe"].n_input == 4
    assert by_name["differential_missingness"].n_input == 3

    d = res.to_dict()
    assert d["n_removed"] == 4
    assert sum(f["n_removed"] for f in d["filters"]) == 4


def test_variant_qc_computes_hwe_in_controls_only():
    """A variant in perfect HWE among controls but wildly out among cases —
    a real risk locus under a recessive-ish model — must survive."""
    rng = np.random.default_rng(2718)
    n_case, n_ctrl = 300, 300
    controls = rng.binomial(2, 0.3, size=n_ctrl)           # clean HWE
    affected = np.full(n_case, 1)                          # all heterozygous
    dos = np.concatenate([affected, controls])[None, :].astype(np.int8)
    gm = make_gm(dos)
    cases = np.array([1] * n_case + [0] * n_ctrl)

    with_status = qc_mod.variant_qc(gm, cases=cases, min_maf=0.0, hwe_p=1e-6)
    assert with_status.keep[0], "HWE must be judged on controls only"

    # Without case/control status the same variant is thrown away, and the
    # result says so rather than hiding it.
    blind = qc_mod.variant_qc(gm, min_maf=0.0, hwe_p=1e-6)
    assert not blind.keep[0]
    assert any("no case/control status" in w for w in blind.warnings)
    hwe_step = [f for f in with_status.filters if f.name == "hwe"][0]
    assert hwe_step.detail["group"] == "controls"
    assert hwe_step.detail["n_samples_tested"] == n_ctrl


def test_variant_qc_skips_hwe_on_non_autosomes():
    """Non-PAR X has no Hardy-Weinberg expectation in a mixed-sex cohort, so
    male hemizygosity must not be read as a heterozygote deficit."""
    rng = np.random.default_rng(17)
    n = 400
    males = np.where(rng.random(n // 2) < 0.3, 2, 0)
    females = rng.binomial(2, 0.3, size=n // 2)
    dos = np.concatenate([males, females])[None, :].astype(np.int8)
    gm_x = make_gm(dos, chrom="X", start_pos=50000000)
    res = qc_mod.variant_qc(gm_x, min_maf=0.0, hwe_p=1e-6)
    assert res.keep[0]
    assert np.isnan(res.metrics["hwe_p"][0])
    hwe_step = [f for f in res.filters if f.name == "hwe"][0]
    assert hwe_step.detail["n_skipped_non_autosomal"] == 1
    # The same genotypes on an autosome would be a gross HWE failure.
    gm_auto = make_gm(dos, chrom="1")
    assert not qc_mod.variant_qc(gm_auto, min_maf=0.0, hwe_p=1e-6).keep[0]


def test_variant_qc_differential_missingness_ignores_balanced_missingness():
    rng = np.random.default_rng(1234)
    n_case = n_ctrl = 250
    n = n_case + n_ctrl
    dos = rng.binomial(2, 0.3, size=(2, n)).astype(np.int8)
    # Variant 0: 10% missing in both arms. Variant 1: 10% missing in cases only.
    dos[0, rng.permutation(n_case)[:25]] = -1
    dos[0, n_case + rng.permutation(n_ctrl)[:25]] = -1
    dos[1, rng.permutation(n_case)[:25]] = -1
    gm = make_gm(dos)
    cases = np.array([1] * n_case + [0] * n_ctrl)
    res = qc_mod.variant_qc(gm, min_call_rate=0.8, cases=cases, diff_miss_p=1e-3)
    p = res.metrics["diff_missingness_p"]
    assert p[0] > 0.05 and p[1] < 1e-3
    assert res.keep[0] and not res.keep[1]


# ================================================================= sample QC ==
def build_sex_matrix(rng, n_male=25, n_female=25, n_x=300, n_auto=200):
    """Males hemizygous on X (never heterozygous), females diploid."""
    p = rng.uniform(0.2, 0.5, n_x)
    male_x = np.where(rng.random((n_x, n_male)) < p[:, None], 2, 0)
    female_x = rng.binomial(2, np.repeat(p[:, None], n_female, axis=1))
    x = np.concatenate([male_x, female_x], axis=1).astype(np.int8)

    n = n_male + n_female
    auto = rng.binomial(2, 0.3, size=(n_auto, n)).astype(np.int8)

    ids = (["M{:02d}".format(i) for i in range(n_male)]
           + ["F{:02d}".format(i) for i in range(n_female)])
    gm_x = make_gm(x, chrom="X", sample_ids=ids, start_pos=50000000)
    gm_auto = make_gm(auto, chrom="1", sample_ids=ids)
    return stack_gm(gm_auto, gm_x), ids, n_male


def test_sex_check_infers_sex_and_flags_a_mislabelled_sample():
    rng = np.random.default_rng(606)
    gm, ids, n_male = build_sex_matrix(rng)
    truth = ["M"] * n_male + ["F"] * (len(ids) - n_male)

    reported = list(truth)
    reported[3] = "F"                      # M03 mislabelled female
    reported[30] = "M"                     # F05 mislabelled male

    sx = qc_mod.sex_check(gm, reported_sex=reported)
    inferred = list(sx["inferred_sex"])
    assert inferred[:n_male] == ["male"] * n_male
    assert inferred[n_male:] == ["female"] * (len(ids) - n_male)
    assert sx["n_x_variants"] == 300
    assert sorted(sx["discordant_ids"]) == ["F05", "M03"]
    assert sx["discordant"][0]["reported"] != sx["discordant"][0]["inferred"]
    assert sx["n_ambiguous"] == 0

    res = qc_mod.sample_qc(gm, reported_sex=reported)
    step = [f for f in res.filters if f.name == "sex_check"][0]
    assert step.n_removed == 2
    assert sorted(step.removed_ids) == ["F05", "M03"]
    assert sorted(res.removed_ids()) == ["F05", "M03"]

    kept = qc_mod.sample_qc(gm, reported_sex=reported, drop_sex_discordant=False)
    assert kept.n_removed == 0
    assert any("flagged but kept" in w for w in kept.warnings)


def test_sex_check_excludes_pseudoautosomal_variants():
    """Males are diploid in the PAR, so PAR heterozygosity must not count
    against a male's sex inference."""
    rng = np.random.default_rng(77)
    n_male = 10
    x = np.where(rng.random((100, n_male)) < 0.4, 2, 0).astype(np.int8)
    gm_nonpar = make_gm(x, chrom="X", start_pos=50000000,
                        sample_ids=["M{}".format(i) for i in range(n_male)])
    par = np.ones((100, n_male), dtype=np.int8)          # all het, inside PAR1
    gm_par = make_gm(par, chrom="chrX", start_pos=100000, spacing=10000,
                     sample_ids=gm_nonpar.sample_ids)
    gm = stack_gm(gm_nonpar, gm_par)
    assert gm.build == "GRCh38"

    sx = qc_mod.sex_check(gm)
    assert sx["n_x_variants"] == 100                      # PAR sites excluded
    assert list(sx["inferred_sex"]) == ["male"] * n_male

    # Without a build we cannot locate the PAR, and the het calls there drag
    # the inference to ambiguous — the result must say why.
    gm_nobuild = GenotypeMatrix(sample_ids=list(gm.sample_ids),
                                variants=list(gm.variants),
                                dosages=gm.dosages, build=None)
    sx2 = qc_mod.sex_check(gm_nobuild)
    assert sx2["n_x_variants"] == 200
    assert any("build unknown" in w for w in sx2["warnings"])
    assert "male" not in list(sx2["inferred_sex"])


def test_sample_qc_call_rate():
    rng = np.random.default_rng(909)
    gm, ids, n_male = build_sex_matrix(rng)
    dos = gm.dosages.copy()
    dos[:, 7] = np.where(rng.random(gm.n_variants) < 0.5, -1, dos[:, 7])
    gm2 = GenotypeMatrix(sample_ids=list(ids), variants=list(gm.variants),
                         dosages=dos.astype(np.int8), build="GRCh38")
    res = qc_mod.sample_qc(gm2, min_call_rate=0.95)
    step = [f for f in res.filters if f.name == "sample_call_rate"][0]
    assert step.n_removed == 1 and step.removed_ids == [ids[7]]
    assert not res.keep[7]
    assert res.metrics["call_rate"][7] < 0.95
    assert np.isfinite(res.metrics["autosomal_het_rate"]).all()


def test_sample_qc_without_x_data_says_the_check_could_not_run():
    rng = np.random.default_rng(4)
    gm = make_gm(rng.binomial(2, 0.3, size=(100, 20)).astype(np.int8))
    res = qc_mod.sample_qc(gm, reported_sex=["M"] * 20)
    assert res.n_removed == 0
    assert any("no non-PAR X variants" in w for w in res.warnings)
    assert not any(f.name == "sex_check" for f in res.filters)


def test_normalise_sex_accepts_the_usual_encodings():
    assert qc_mod.normalise_sex("M") == "male"
    assert qc_mod.normalise_sex("female") == "female"
    assert qc_mod.normalise_sex(1) == "male"          # PLINK .fam convention
    assert qc_mod.normalise_sex(2) == "female"
    assert qc_mod.normalise_sex(0) == "unknown"
    assert qc_mod.normalise_sex(None) == "unknown"
    assert qc_mod.normalise_sex(float("nan")) == "unknown"


# ================================================= end-to-end QC composition ==
def test_full_qc_pipeline_composes():
    """The §4.2 mandatory sequence: sample QC, variant QC, relatedness,
    ancestry outliers — each step inspectable on its own."""
    rng = np.random.default_rng(2026)
    gm, labels = simulate_two_populations(rng, n_variants=1200, n_per_pop=40)
    cases = np.array([1, 0] * (gm.n_samples // 2))

    s_res = qc_mod.sample_qc(gm, min_call_rate=0.95)
    gm1 = gm.subset_samples(s_res.kept_ids())
    assert gm1.n_samples == gm.n_samples

    v_res = qc_mod.variant_qc(gm1, cases=cases, min_maf=0.01, hwe_p=1e-6)
    gm2 = gm1.subset_variants(v_res.keep)
    assert gm2.n_variants == v_res.n_kept and gm2.n_variants > 1000

    kin = kin_mod.king_robust(gm2, prune=True)
    unrelated = kin_mod.unrelated_set(kin)
    assert len(unrelated) == gm2.n_samples          # nothing related by construction

    pcs = pca_mod.compute_pca(gm2.subset_samples(unrelated), n_components=4,
                              prune=True)
    outliers = pca_mod.ancestry_outliers(pcs)
    assert outliers.sum() == 0                      # two clean clusters, no strays
    assert pcs.provenance["ld_pruned"] is True
    assert pcs.provenance["n_variants_used"] == pcs.n_variants_used
