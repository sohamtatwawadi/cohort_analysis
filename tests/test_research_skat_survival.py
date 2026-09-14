"""Verification for the gene-based rare-variant tests (§4.4) and the
time-to-event suite (§4.7).

The bar here is ground truth, not smoke:

*   Davies' quadratic-form tail is checked against cases where the mixture
    collapses to a plain chi-square and the answer is known in closed form, and
    against Monte Carlo where it does not.
*   SKAT is checked for **null calibration over many replicates**. This is the
    only test that matters for a variance-component statistic: a miscalibrated
    SKAT produces a perfectly plausible-looking p-value on any single dataset
    and is worthless across a genome.
*   Kaplan-Meier is checked against a curve computed by hand.
*   Cox is checked against a simulation with a known hazard ratio.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats as sps

from backend.app.research.stats import skat as sk
from backend.app.research.stats import survival as sv
from backend.app.research.types import MISSING, GenotypeMatrix, Variant


# =========================================================================== #
# helpers
# =========================================================================== #
def _rare_genotypes(rng, n, mafs):
    """Samples x variants dosage matrix drawn under HWE."""
    m = len(mafs)
    G = np.zeros((n, m), dtype=float)
    for j, f in enumerate(mafs):
        G[:, j] = rng.binomial(2, f, size=n)
    return G


# =========================================================================== #
# §4.4 -- quadratic-form tail probabilities
# =========================================================================== #
class TestQuadraticForm:
    def test_davies_matches_plain_chisquare(self):
        """lambda = (1,1,1) makes Q exactly chi-square with 3 df."""
        for q in (0.5, 1.0, 3.0, 7.8147, 16.27, 30.0):
            p, info = sk.davies_pvalue(q, [1.0, 1.0, 1.0], acc=1e-9)
            assert info["fault"] in (0, 2)
            # `acc` is an absolute tolerance, so assert against it as one.
            assert p == pytest.approx(float(sps.chi2.sf(q, 3)), rel=1e-6, abs=1e-9)

    def test_davies_matches_scaled_chisquare(self):
        """lambda = (2,2) makes Q = 2 * chi2_2, so P(Q>q) = P(chi2_2 > q/2)."""
        for q in (1.0, 4.0, 12.0, 25.0, 40.0):
            p, info = sk.davies_pvalue(q, [2.0, 2.0])
            assert info["fault"] in (0, 2)
            assert p == pytest.approx(
                float(sps.chi2.sf(q / 2.0, 2)), rel=1e-6, abs=1e-9)

    def test_davies_far_tail_stays_accurate(self):
        """Davies' `acc` is an absolute tolerance, so the relative error grows
        as p shrinks -- but it is still correct to several digits where a
        gene-based scan reports its hits."""
        for q, tol in ((30.0, 1e-3), (45.0, 1e-2)):
            p, info = sk.davies_pvalue(q, [1.0, 1.0, 1.0])
            truth = float(sps.chi2.sf(q, 3))
            assert info["fault"] in (0, 2)
            assert truth < 1e-5
            assert p == pytest.approx(truth, rel=tol), (q, p, truth)

    def test_chisquare_approximation_is_wrong_where_davies_is_right(self):
        """The justification for implementing Davies at all.

        A Satterthwaite-style scaled chi-square matches the first two moments
        of the mixture exactly and is perfectly respectable in the body of the
        distribution. In the tail -- the only place a p-value is read -- it is
        off by orders of magnitude, in the anti-conservative direction.
        """
        lam = 2.0 ** np.arange(10)[::-1].astype(float)   # 512, 256, ..., 1
        c1, c2 = lam.sum(), (lam ** 2).sum()
        a, dof = c2 / c1, c1 ** 2 / c2          # Satterthwaite match

        # Body of the distribution: the approximation is fine, which is exactly
        # why it survives casual checking.
        q_body = float(lam.sum() * 2)
        p_d_body, _ = sk.davies_pvalue(q_body, lam)
        assert sps.chi2.sf(q_body / a, dof) == pytest.approx(p_d_body, rel=0.15)

        # Tail, where p-values are actually read. Davies is validated against
        # exact chi-square and Monte Carlo in the tests above; the moment match
        # is out by more than two orders of magnitude, and in the dangerous
        # direction (it calls the gene far more significant than it is).
        q_tail = float(lam.sum() * 16)
        p_davies, info = sk.davies_pvalue(q_tail, lam)
        p_chisq = float(sps.chi2.sf(q_tail / a, dof))
        assert info["fault"] in (0, 2)
        assert 1e-9 < p_davies < 1e-3
        assert p_davies / p_chisq > 50.0, (p_davies, p_chisq)

    def test_davies_matches_monte_carlo_for_mixed_weights(self):
        """No closed form for an unequal mixture -- check against simulation."""
        lam = np.array([5.0, 2.0, 1.0, 0.5, 0.1])
        rng = np.random.default_rng(11)
        draws = (rng.chisquare(1, size=(400000, len(lam))) * lam).sum(axis=1)
        for q in np.percentile(draws, [50, 75, 90, 95, 99]):
            mc = float(np.mean(draws > q))
            p, info = sk.davies_pvalue(float(q), lam)
            assert info["fault"] == 0
            assert abs(p - mc) < 4.0 * math.sqrt(mc * (1 - mc) / len(draws)) + 1e-4

    def test_davies_and_liu_agree_within_an_order_of_magnitude(self):
        """Mid-range case: the fallback must be in the same ballpark."""
        lam = np.array([4.0, 3.0, 2.0, 1.0, 0.5, 0.2])
        for q in (25.0, 40.0, 60.0):
            p_d, info = sk.davies_pvalue(q, lam)
            p_l = sk.liu_pvalue(q, lam)
            assert info["fault"] == 0
            assert 0.0 < p_d < 1.0 and 0.0 < p_l < 1.0
            ratio = p_d / p_l
            assert 0.1 < ratio < 10.0, (q, p_d, p_l)

    def test_liu_is_also_close_to_chisquare_truth_in_the_body(self):
        for q in (2.0, 4.0, 7.0):
            assert sk.liu_pvalue(q, [1.0, 1.0, 1.0]) == pytest.approx(
                float(sps.chi2.sf(q, 3)), rel=0.05)

    def test_davies_monotone_and_bounded(self):
        lam = [3.0, 1.0, 0.25]
        ps = [sk.davies_pvalue(q, lam)[0] for q in (0.5, 2, 5, 10, 20, 40)]
        assert all(0.0 <= p <= 1.0 for p in ps)
        assert all(ps[i] > ps[i + 1] for i in range(len(ps) - 1))


# =========================================================================== #
# §4.4 -- weights
# =========================================================================== #
class TestBetaWeights:
    def test_beta_1_25_upweights_rare(self):
        w = sk.beta_weights(np.array([0.0005, 0.005, 0.01, 0.05, 0.2]))
        assert w[0] > w[1] > w[2] > w[3] > w[4]
        assert w[0] == pytest.approx(25.0, rel=0.05)     # density at MAF ~ 0
        assert w[0] / w[3] > 3.0                         # near-singleton vs 5%
        assert w[0] / w[4] > 100.0                       # near-singleton vs 20%

    def test_matches_scipy_beta_density(self):
        maf = np.array([0.001, 0.01, 0.1, 0.3])
        assert np.allclose(sk.beta_weights(maf), sps.beta.pdf(maf, 1, 25))

    def test_alternate_parameters(self):
        w = sk.beta_weights(np.array([0.01, 0.2]), a=0.5, b=0.5)
        assert np.all(np.isfinite(w)) and w[0] > w[1]


# =========================================================================== #
# §4.4 -- burden
# =========================================================================== #
class TestBurden:
    def test_enriched_carriers_give_small_p_and_positive_beta(self):
        rng = np.random.default_rng(7)
        n, mafs = 1500, [0.01, 0.015, 0.02, 0.008, 0.012, 0.02, 0.01, 0.015]
        G = _rare_genotypes(rng, n, mafs)
        carrier = (G > 0).any(axis=1)
        logit = -1.2 + 1.6 * carrier
        y = rng.binomial(1, 1.0 / (1.0 + np.exp(-logit))).astype(float)

        res = sk.burden_test(G, y, binary=True)
        assert res["status"] == "ok"
        assert res["p"] < 1e-3
        assert res["beta"] > 0                       # enriched in cases
        assert res["odds_ratio"] > 1.0
        assert res["n_carriers_cases"] + res["n_carriers_controls"] == res["n_carriers"]

    def test_protective_direction_is_recovered(self):
        rng = np.random.default_rng(8)
        G = _rare_genotypes(rng, 1500, [0.02] * 8)
        carrier = (G > 0).any(axis=1)
        logit = 0.0 - 1.6 * carrier
        y = rng.binomial(1, 1.0 / (1.0 + np.exp(-logit))).astype(float)
        res = sk.burden_test(G, y, binary=True)
        assert res["status"] == "ok"
        assert res["p"] < 1e-3
        assert res["beta"] < 0

    def test_seeded_null_is_not_significant(self):
        rng = np.random.default_rng(2024)
        G = _rare_genotypes(rng, 1200, [0.02] * 10)
        y = rng.binomial(1, 0.5, size=1200).astype(float)
        res = sk.burden_test(G, y, binary=True)
        assert res["status"] == "ok"
        assert res["p"] > 0.05

    def test_null_p_values_are_uniform_ish(self):
        rng = np.random.default_rng(99)
        ps = []
        for _ in range(200):
            G = _rare_genotypes(rng, 500, [0.03] * 8)
            y = rng.binomial(1, 0.5, size=500).astype(float)
            r = sk.burden_test(G, y, binary=True)
            if r["status"] == "ok":
                ps.append(r["p"])
        frac = float(np.mean(np.asarray(ps) < 0.05))
        assert 0.01 <= frac <= 0.12, frac

    def test_quantitative_trait(self):
        rng = np.random.default_rng(5)
        G = _rare_genotypes(rng, 1000, [0.02] * 10)
        burden = G.sum(axis=1)
        y = 0.9 * burden + rng.normal(size=1000)
        res = sk.burden_test(G, y, binary=False)
        assert res["status"] == "ok"
        assert res["model"] == "linear"
        assert res["beta"] > 0 and res["p"] < 1e-6
        assert res["n_carriers_cases"] is None       # undefined for a QT

    def test_covariate_adjustment_removes_a_confounded_signal(self):
        """Ancestry-like confounder drives both genotype and outcome."""
        rng = np.random.default_rng(31)
        n = 2000
        pop = rng.binomial(1, 0.5, size=n).astype(float)
        G = np.zeros((n, 6))
        for j in range(6):
            f = np.where(pop > 0, 0.06, 0.005)
            G[:, j] = rng.binomial(2, f)
        y = rng.binomial(1, 1.0 / (1.0 + np.exp(-(-0.5 + 1.5 * pop)))).astype(float)

        naive = sk.burden_test(G, y, binary=True)
        adj = sk.burden_test(G, y, X=pop[:, None], binary=True)
        assert naive["p"] < 1e-4                     # spurious
        assert adj["p"] > 0.01                       # confounder absorbed

    def test_accepts_genotype_matrix_and_handles_missing(self):
        rng = np.random.default_rng(3)
        n, m = 600, 6
        d = rng.binomial(2, 0.03, size=(m, n)).astype(np.int8)
        d[0, :20] = MISSING                          # a no-call block
        gm = GenotypeMatrix(
            sample_ids=["s{}".format(i) for i in range(n)],
            variants=[Variant("1", 100 + j, "A", "T") for j in range(m)],
            dosages=d)
        y = rng.binomial(1, 0.4, size=n).astype(float)
        res = sk.burden_test(gm, y, binary=True)
        assert res["status"] == "ok"
        assert res["n_samples"] == n and res["n_variants"] == m


# =========================================================================== #
# §4.4 -- SKAT
# =========================================================================== #
class TestSkat:
    def test_null_calibration_binary(self):
        """THE test for a variance-component statistic.

        Random genotypes, random outcome, no relationship whatsoever. The
        p-values must be approximately uniform; in particular the type-I error
        at the nominal 0.05 must actually be near 0.05. A SKAT that is
        anti-conservative here looks completely fine on any single gene and
        produces a genome full of false positives.
        """
        rng = np.random.default_rng(20240501)
        reps = 250
        n, mafs = 400, [0.05, 0.04, 0.03, 0.06, 0.02, 0.05, 0.03, 0.04, 0.02, 0.05]
        ps = []
        for _ in range(reps):
            G = _rare_genotypes(rng, n, mafs)
            y = rng.binomial(1, 0.5, size=n).astype(float)
            r = sk.skat_test(G, y, binary=True)
            if r["status"] == "ok":
                ps.append(r["p"])
        ps = np.asarray(ps)
        assert len(ps) > 0.95 * reps
        frac05 = float(np.mean(ps < 0.05))
        assert 0.01 <= frac05 <= 0.12, "type-I error at 0.05 was {}".format(frac05)
        # Mean of a Uniform(0,1) is 0.5; a badly-scaled Q shifts this hard.
        assert 0.40 <= float(ps.mean()) <= 0.60, ps.mean()
        # And the whole distribution should look uniform.
        ks = float(sps.kstest(ps, "uniform").pvalue)
        assert ks > 0.001, "KS vs uniform p = {}".format(ks)

    def test_null_calibration_quantitative(self):
        rng = np.random.default_rng(777)
        reps = 200
        ps = []
        for _ in range(reps):
            G = _rare_genotypes(rng, 300, [0.05] * 8)
            y = rng.normal(size=300)
            r = sk.skat_test(G, y, binary=False)
            if r["status"] == "ok":
                ps.append(r["p"])
        ps = np.asarray(ps)
        frac05 = float(np.mean(ps < 0.05))
        assert 0.01 <= frac05 <= 0.12, frac05
        assert float(sps.kstest(ps, "uniform").pvalue) > 0.001

    def test_detects_mixed_direction_effects_where_burden_cannot(self):
        """The reason SKAT exists: half the variants protective."""
        rng = np.random.default_rng(42)
        n, m = 2500, 10
        G = _rare_genotypes(rng, n, [0.05] * m)
        sign = np.array([1.0, -1.0] * (m // 2))
        eff = (G * (sign * 1.1)[None, :]).sum(axis=1)
        y = rng.binomial(1, 1.0 / (1.0 + np.exp(-(0.0 + eff)))).astype(float)

        s = sk.skat_test(G, y, binary=True)
        b = sk.burden_test(G, y, binary=True)
        assert s["p"] < 1e-4
        assert b["p"] > s["p"] * 100          # burden cancels out, SKAT does not

    def test_records_which_method_produced_the_p_value(self):
        rng = np.random.default_rng(13)
        G = _rare_genotypes(rng, 500, [0.03] * 8)
        y = rng.binomial(1, 0.5, size=500).astype(float)
        r = sk.skat_test(G, y, binary=True)
        assert r["p_method"] in ("davies", "liu", "exact")
        assert "p_liu" in r and "davies_fallback_used" in r
        # Davies should be doing the work on an ordinary gene.
        assert r["p_method"] == "davies"

    def test_min_carriers_guardrail(self):
        """§4.4: report the gene as under-powered, do not invent a p-value."""
        rng = np.random.default_rng(1)
        n = 400
        G = np.zeros((n, 3))
        G[0, 0] = 1                                  # exactly one carrier
        y = rng.binomial(1, 0.5, size=n).astype(float)

        for fn in (sk.burden_test, sk.skat_test, sk.skat_o_test):
            r = fn(G, y, binary=True, min_carriers=2)
            assert r["status"] == "insufficient_carriers", fn.__name__
            assert r["p"] is None
            assert r["n_carriers"] == 1
            assert r["min_carriers"] == 2
            assert "reason" in r and r["reason"]

    def test_min_carriers_threshold_is_respected(self):
        rng = np.random.default_rng(17)
        n = 400
        G = np.zeros((n, 3))
        G[:5, 0] = 1                                 # five carriers
        y = rng.binomial(1, 0.5, size=n).astype(float)
        assert sk.skat_test(G, y, min_carriers=2)["status"] == "ok"
        assert sk.skat_test(G, y, min_carriers=5)["status"] == "ok"
        assert sk.skat_test(G, y, min_carriers=6)["status"] == "insufficient_carriers"

    def test_carrier_counts_reported_even_when_gated(self):
        rng = np.random.default_rng(19)
        n = 200
        G = np.zeros((n, 2))
        G[0, 0] = 2
        y = np.zeros(n)
        y[:100] = 1.0
        r = sk.burden_test(G, y, binary=True, min_carriers=3)
        assert r["n_carriers_cases"] + r["n_carriers_controls"] == 1


# =========================================================================== #
# §4.4 -- SKAT-O
# =========================================================================== #
class TestSkatO:
    def test_grid_and_approximation_are_declared(self):
        rng = np.random.default_rng(4)
        G = _rare_genotypes(rng, 800, [0.04] * 8)
        y = rng.binomial(1, 0.5, size=800).astype(float)
        r = sk.skat_o_test(G, y, binary=True)
        assert r["status"] == "ok"
        assert r["rho_grid"] == list(sk.DEFAULT_RHOS)
        assert r["correction"] == "sidak_galwey_meff"
        assert r["correction_is_approximate"] is True
        assert "approximation" in r["correction_note"].lower() or \
               "APPROXIMATION" in r["correction_note"]
        assert 1.0 <= r["m_eff"] <= len(r["rho_grid"])
        assert r["p_min"] <= r["p"] <= r["p_bonferroni"] + 1e-12

    def test_tracks_burden_when_effects_are_concordant(self):
        rng = np.random.default_rng(21)
        n, m = 1500, 10
        G = _rare_genotypes(rng, n, [0.03] * m)
        eff = 0.35 * G.sum(axis=1)
        y = rng.binomial(1, 1.0 / (1.0 + np.exp(-(-0.5 + eff)))).astype(float)
        r = sk.skat_o_test(G, y, binary=True)
        assert r["status"] == "ok"
        assert r["p"] < 1e-3
        assert r["p_grid_method"] == "davies"
        assert r["rho_optimal"] >= 0.5          # leans toward burden
        # Burden beats plain SKAT here, and SKAT-O should sit near burden.
        assert sk.burden_test(G, y, binary=True)["p"] < sk.skat_test(G, y, binary=True)["p"]

    def test_tracks_skat_when_effects_are_mixed(self):
        rng = np.random.default_rng(22)
        n, m = 2000, 10
        G = _rare_genotypes(rng, n, [0.05] * m)
        sign = np.array([1.0, -1.0] * (m // 2))
        eff = (G * (sign * 1.1)[None, :]).sum(axis=1)
        y = rng.binomial(1, 1.0 / (1.0 + np.exp(-eff))).astype(float)
        r = sk.skat_o_test(G, y, binary=True)
        assert r["status"] == "ok"
        assert r["p"] < 1e-3
        assert r["rho_optimal"] <= 0.3          # leans toward SKAT

    def test_q_rho_is_the_stated_linear_combination(self):
        rng = np.random.default_rng(23)
        G = _rare_genotypes(rng, 600, [0.04] * 6)
        y = rng.binomial(1, 0.5, size=600).astype(float)
        r = sk.skat_o_test(G, y, binary=True)
        for rho, q in zip(r["rho_grid"], r["q_by_rho"]):
            expected = (1 - rho) * r["Q_skat"] + rho * r["Q_burden"]
            assert q == pytest.approx(expected, rel=1e-10)

    def test_null_calibration(self):
        """The Sidak/M_eff correction is approximate, so it gets the same
        treatment as SKAT itself: it must hold the nominal level under the
        null. An uncorrected min-over-grid would sit near 0.08 here."""
        rng = np.random.default_rng(20240502)
        ps = []
        p_mins = []
        for _ in range(200):
            G = _rare_genotypes(rng, 400, [0.05] * 8)
            y = rng.binomial(1, 0.5, size=400).astype(float)
            r = sk.skat_o_test(G, y, binary=True)
            if r["status"] == "ok":
                ps.append(r["p"])
                p_mins.append(r["p_min"])
        ps = np.asarray(ps)
        frac = float(np.mean(ps < 0.05))
        assert 0.01 <= frac <= 0.12, frac
        # The correction has to actually be doing something: the raw minimum
        # over the grid is anti-conservative.
        assert float(np.mean(np.asarray(p_mins) < 0.05)) > frac


# =========================================================================== #
# §4.7 -- Kaplan-Meier
# =========================================================================== #
class TestKaplanMeier:
    def test_hand_computed_six_subject_curve(self):
        """Six subjects, one censored at t=3. Worked by hand:

            t=1  n=6 d=1  S = 5/6
            t=2  n=5 d=1  S = 5/6 * 4/5 = 2/3
            t=3  n=4 d=0  S = 2/3            (censoring, curve flat)
            t=4  n=3 d=1  S = 2/3 * 2/3 = 4/9
            t=5  n=2 d=1  S = 4/9 * 1/2 = 2/9
            t=6  n=1 d=1  S = 0
        """
        time = [1, 2, 3, 4, 5, 6]
        event = [1, 1, 0, 1, 1, 1]
        km = sv.kaplan_meier(time, event)

        assert km["time"] == [1, 2, 3, 4, 5, 6]
        expected = [5 / 6, 2 / 3, 2 / 3, 4 / 9, 2 / 9, 0.0]
        assert km["survival"] == pytest.approx(expected, abs=1e-12)
        assert km["n_risk"] == [6, 5, 4, 3, 2, 1]
        assert km["n_event"] == [1, 1, 0, 1, 1, 1]
        assert km["n_censored"] == [0, 0, 1, 0, 0, 0]
        # Median = first time S(t) <= 0.5, i.e. t = 4 where S = 4/9.
        assert km["median"] == 4.0

    def test_greenwood_variance_hand_checked(self):
        """At t=2, S=2/3 and Greenwood gives
        Var = S^2 * (1/(6*5) + 1/(5*4)) = (4/9) * (1/30 + 1/20)."""
        km = sv.kaplan_meier([1, 2, 3, 4, 5, 6], [1, 1, 0, 1, 1, 1])
        var = (2 / 3) ** 2 * (1 / 30 + 1 / 20)
        assert km["std_err"][1] == pytest.approx(math.sqrt(var), rel=1e-12)

    def test_no_censoring_reduces_to_empirical_survival(self):
        km = sv.kaplan_meier([1, 2, 3, 4], [1, 1, 1, 1])
        assert km["survival"] == pytest.approx([0.75, 0.5, 0.25, 0.0])
        assert km["n_risk"] == [4, 3, 2, 1]

    def test_tied_event_times(self):
        """Three events at t=1 out of 5 at risk -> S = 2/5."""
        km = sv.kaplan_meier([1, 1, 1, 2, 3], [1, 1, 1, 1, 0])
        assert km["time"] == [1, 2, 3]
        assert km["survival"] == pytest.approx([2 / 5, 1 / 5, 1 / 5])
        assert km["n_risk"] == [5, 2, 1]
        assert km["n_event"] == [3, 1, 0]

    def test_confidence_band_stays_inside_unit_interval(self):
        rng = np.random.default_rng(6)
        t = rng.exponential(5.0, size=150)
        e = rng.binomial(1, 0.7, size=150)
        km = sv.kaplan_meier(t, e)
        lo = np.asarray(km["ci_lower"], dtype=float)
        hi = np.asarray(km["ci_upper"], dtype=float)
        fin = np.isfinite(lo) & np.isfinite(hi)
        assert np.all(lo[fin] >= 0.0) and np.all(hi[fin] <= 1.0)
        s = np.asarray(km["survival"])
        assert np.all(lo[fin] <= s[fin] + 1e-12)
        assert np.all(hi[fin] >= s[fin] - 1e-12)

    def test_grouped_curves(self):
        time = [1, 2, 3, 4, 1, 2, 3, 4]
        event = [1, 1, 0, 1, 1, 1, 0, 1]
        groups = ["a"] * 4 + ["b"] * 4
        km = sv.kaplan_meier(time, event, groups=groups)
        assert sorted(km["groups"]) == ["a", "b"]
        assert km["groups"]["a"]["n_risk"] == [4, 3, 2, 1]
        assert km["groups"]["a"]["survival"] == pytest.approx(
            km["groups"]["b"]["survival"])

    def test_median_none_when_curve_never_reaches_half(self):
        km = sv.kaplan_meier([1, 2, 3, 4, 5], [1, 0, 0, 0, 0])
        assert km["median"] is None


# =========================================================================== #
# §4.7 -- log-rank
# =========================================================================== #
class TestLogRank:
    def test_identical_groups_give_p_of_one(self):
        base_t = list(range(1, 11))
        base_e = [1] * 9 + [0]
        res = sv.logrank_test(base_t * 2, base_e * 2, ["a"] * 10 + ["b"] * 10)
        assert res["chi2"] == pytest.approx(0.0, abs=1e-12)
        assert res["p"] == pytest.approx(1.0)
        assert res["df"] == 1

    def test_strongly_separated_groups_give_tiny_p(self):
        rng = np.random.default_rng(10)
        n = 150
        t_a = rng.exponential(1.0, n)
        t_b = rng.exponential(6.0, n)
        t = np.concatenate([t_a, t_b])
        e = np.ones(2 * n)
        g = np.array(["a"] * n + ["b"] * n)
        res = sv.logrank_test(t, e, g)
        assert res["p"] < 1e-15
        assert res["observed"][0] == res["observed"][1] == n
        # Group a dies early -> fewer expected events than observed.
        assert res["observed"][0] > res["expected"][0]

    def test_null_p_values_are_uniform_ish(self):
        rng = np.random.default_rng(55)
        ps = []
        for _ in range(300):
            t = rng.exponential(1.0, 120)
            e = rng.binomial(1, 0.8, 120)
            g = rng.integers(0, 2, 120)
            ps.append(sv.logrank_test(t, e, g)["p"])
        frac = float(np.mean(np.asarray(ps) < 0.05))
        assert 0.01 <= frac <= 0.12, frac

    def test_three_groups_uses_two_df(self):
        rng = np.random.default_rng(12)
        t = np.concatenate([rng.exponential(s, 100) for s in (1.0, 2.0, 6.0)])
        e = np.ones(300)
        g = np.array([0] * 100 + [1] * 100 + [2] * 100)
        res = sv.logrank_test(t, e, g)
        assert res["df"] == 2
        assert len(res["groups"]) == 3
        assert res["p"] < 1e-10

    def test_expected_counts_sum_to_total_events(self):
        rng = np.random.default_rng(14)
        t = rng.exponential(2.0, 200)
        e = rng.binomial(1, 0.8, 200)
        g = rng.integers(0, 3, 200)
        # Make the last observation censored so no event time has a risk set of
        # one (those contribute no variance and are skipped by construction).
        t[int(np.argmax(t))] += 1.0
        e[int(np.argmax(t))] = 0
        res = sv.logrank_test(t, e, g)
        assert sum(res["expected"]) == pytest.approx(float(e.sum()), rel=1e-9)
        assert sum(res["observed"]) == pytest.approx(float(e.sum()))


# =========================================================================== #
# §4.7 -- Cox
# =========================================================================== #
def _simulate_exponential_ph(rng, n, beta, base_rate=0.10, cens_rate=0.04):
    x = rng.binomial(1, 0.5, size=n).astype(float)
    t_event = rng.exponential(1.0 / (base_rate * np.exp(beta * x)))
    t_cens = rng.exponential(1.0 / cens_rate, size=n)
    t = np.minimum(t_event, t_cens)
    e = (t_event <= t_cens).astype(int)
    return t, e, x


class TestCoxPH:
    def test_recovers_a_known_hazard_ratio(self):
        rng = np.random.default_rng(2718)
        t, e, x = _simulate_exponential_ph(rng, 4000, math.log(2.0))
        fit = sv.cox_ph(t, e, x[:, None], ["carrier"])
        assert fit["converged"]
        hr = fit["covariates"][0]["hazard_ratio"]
        assert abs(hr - 2.0) / 2.0 < 0.15, hr
        assert fit["covariates"][0]["ci_lower"] < 2.0 < fit["covariates"][0]["ci_upper"]
        assert fit["covariates"][0]["p"] < 1e-10
        assert fit["ties"] == "efron"
        assert fit["n_iter"] >= 1

    def test_recovers_a_protective_hazard_ratio(self):
        rng = np.random.default_rng(314)
        t, e, x = _simulate_exponential_ph(rng, 4000, math.log(0.5))
        fit = sv.cox_ph(t, e, x[:, None], ["carrier"])
        hr = fit["covariates"][0]["hazard_ratio"]
        assert abs(hr - 0.5) / 0.5 < 0.15, hr

    def test_two_covariates_recovered_jointly(self):
        rng = np.random.default_rng(8080)
        n = 5000
        x1 = rng.binomial(1, 0.5, n).astype(float)
        x2 = rng.normal(size=n)
        lin = math.log(2.0) * x1 + 0.5 * x2
        t_ev = rng.exponential(1.0 / (0.1 * np.exp(lin)))
        t_c = rng.exponential(25.0, n)
        t = np.minimum(t_ev, t_c)
        e = (t_ev <= t_c).astype(int)
        fit = sv.cox_ph(t, e, np.column_stack([x1, x2]), ["x1", "x2"])
        assert abs(fit["covariates"][0]["beta"] - math.log(2.0)) < 0.10
        assert abs(fit["covariates"][1]["beta"] - 0.5) < 0.06
        assert fit["names"] == ["x1", "x2"]

    def test_null_covariate_is_not_significant(self):
        rng = np.random.default_rng(4242)
        n = 2000
        x = rng.normal(size=n)
        t = rng.exponential(10.0, n)
        e = rng.binomial(1, 0.8, n)
        fit = sv.cox_ph(t, e, x[:, None], ["noise"])
        assert fit["covariates"][0]["p"] > 0.05

    def test_efron_differs_from_breslow_when_ties_are_heavy(self):
        """Year-granularity follow-up: ties everywhere. Breslow attenuates."""
        rng = np.random.default_rng(1234)
        t, e, x = _simulate_exponential_ph(rng, 3000, math.log(2.5))
        t_tied = np.ceil(t / 2.0)                    # coarse, heavy ties
        assert len(np.unique(t_tied)) < 0.1 * len(t_tied)

        ef = sv.cox_ph(t_tied, e, x[:, None], ["carrier"], ties="efron")
        br = sv.cox_ph(t_tied, e, x[:, None], ["carrier"], ties="breslow")
        b_ef = ef["covariates"][0]["beta"]
        b_br = br["covariates"][0]["beta"]

        assert b_ef != b_br
        assert abs(b_ef - b_br) > 1e-3, (b_ef, b_br)
        # Breslow pretends each tied event faced the full risk set, which
        # shrinks the coefficient toward the null.
        assert abs(b_ef) > abs(b_br)
        assert abs(b_ef - math.log(2.5)) < abs(b_br - math.log(2.5))

    def test_efron_and_breslow_agree_without_ties(self):
        rng = np.random.default_rng(600)
        t, e, x = _simulate_exponential_ph(rng, 800, math.log(2.0))
        assert len(np.unique(t)) == len(t)           # continuous, no ties
        ef = sv.cox_ph(t, e, x[:, None], ["c"], ties="efron")
        br = sv.cox_ph(t, e, x[:, None], ["c"], ties="breslow")
        assert ef["covariates"][0]["beta"] == pytest.approx(
            br["covariates"][0]["beta"], rel=1e-8)

    def test_log_partial_likelihood_is_maximised_at_the_fit(self):
        rng = np.random.default_rng(501)
        t, e, x = _simulate_exponential_ph(rng, 1200, math.log(2.0))
        fit = sv.cox_ph(t, e, x[:, None], ["c"])
        b = fit["covariates"][0]["beta"]
        for delta in (-0.2, -0.05, 0.05, 0.2):
            order = np.argsort(t, kind="mergesort")
            ll, _, _ = sv._cox_loglik(
                t[order], e[order], x[order][:, None], np.array([b + delta]), "efron")
            assert ll < fit["loglik"] + 1e-9
        assert fit["loglik"] > fit["loglik_null"]
        assert fit["lr_p"] < 1e-8


# =========================================================================== #
# §4.7 -- proportional-hazards diagnostics
# =========================================================================== #
class TestPHAssumption:
    def test_does_not_flag_a_genuinely_proportional_effect(self):
        rng = np.random.default_rng(2000)
        t, e, x = _simulate_exponential_ph(rng, 3000, math.log(2.0))
        fit = sv.cox_ph(t, e, x[:, None], ["carrier"])
        ph = sv.ph_assumption_test(t, e, x[:, None], fit["beta"], ["carrier"])
        assert ph["status"] == "ok"
        assert ph["p_global"] > 0.05, ph["p_global"]
        assert ph["covariates"][0]["p"] > 0.05

    def test_flags_a_time_varying_effect(self):
        """Crossing hazards: Weibull shape 0.5 vs 3.0. No constant HR exists."""
        rng = np.random.default_rng(3000)
        n = 1200
        x = rng.binomial(1, 0.5, size=n).astype(float)
        shape = np.where(x > 0, 3.0, 0.5)
        t = rng.weibull(shape) * 1.0
        e = np.ones(n, dtype=int)
        fit = sv.cox_ph(t, e, x[:, None], ["carrier"])
        ph = sv.ph_assumption_test(t, e, x[:, None], fit["beta"], ["carrier"])
        assert ph["status"] == "ok"
        assert ph["p_global"] < 1e-6, ph["p_global"]
        assert ph["covariates"][0]["p"] < 1e-6
        assert abs(ph["covariates"][0]["rho"]) > 0.2

    def test_flags_only_the_offending_covariate(self):
        rng = np.random.default_rng(3100)
        n = 2000
        x_bad = rng.binomial(1, 0.5, size=n).astype(float)
        x_ok = rng.normal(size=n)
        shape = np.where(x_bad > 0, 3.0, 0.5)
        t = rng.weibull(shape) * np.exp(-0.3 * x_ok)
        e = np.ones(n, dtype=int)
        X = np.column_stack([x_bad, x_ok])
        fit = sv.cox_ph(t, e, X, ["bad", "ok"])
        ph = sv.ph_assumption_test(t, e, X, fit["beta"], ["bad", "ok"])
        by_name = {c["name"]: c for c in ph["covariates"]}
        assert by_name["bad"]["p"] < 1e-6
        assert by_name["ok"]["p"] > 0.01

    def test_null_calibration_of_the_ph_test(self):
        rng = np.random.default_rng(3200)
        ps = []
        for _ in range(150):
            t, e, x = _simulate_exponential_ph(rng, 300, math.log(2.0))
            fit = sv.cox_ph(t, e, x[:, None], ["c"])
            ps.append(sv.ph_assumption_test(t, e, x[:, None], fit["beta"], ["c"])["p_global"])
        frac = float(np.mean(np.asarray(ps) < 0.05))
        assert 0.005 <= frac <= 0.15, frac

    def test_too_few_events_reports_rather_than_guesses(self):
        ph = sv.ph_assumption_test([1, 2, 3], [1, 0, 0], np.zeros((3, 1)), [0.0], ["c"])
        assert ph["status"] == "insufficient_events"
        assert ph["p_global"] is None


# =========================================================================== #
# §4.7 -- competing risks
# =========================================================================== #
class TestCumulativeIncidence:
    def test_hand_computed_four_subject_cif(self):
        """Four subjects; cause 1 at t=1 and t=3, competing cause 2 at t=2,
        censored at t=4.

            t=1  n=4  S(0-)=1     F1 += 1 * 1/4   = 0.25   S = 0.75
            t=2  n=3  S(1)=0.75   F2 += 0.75/3    = 0.25   S = 0.50
            t=3  n=2  S(2)=0.50   F1 += 0.50/2    = 0.25   S = 0.25
            t=4  n=1  censored                             S = 0.25

        F1(4) = 0.50, F2(4) = 0.25, S(4) = 0.25 -- and they sum to 1.
        """
        ci = sv.cumulative_incidence(
            time=[1, 2, 3, 4], event=[1, 1, 1, 0], event_type=[1, 2, 1, 0])
        assert ci["time"] == [1, 2, 3, 4]
        assert ci["n_risk"] == [4, 3, 2, 1]
        f1 = ci["causes"]["1"]["cif"]
        f2 = ci["causes"]["2"]["cif"]
        assert f1 == pytest.approx([0.25, 0.25, 0.5, 0.5])
        assert f2 == pytest.approx([0.0, 0.25, 0.25, 0.25])
        assert ci["overall_survival"] == pytest.approx([0.75, 0.5, 0.25, 0.25])
        assert f1[-1] + f2[-1] + ci["overall_survival"][-1] == pytest.approx(1.0)
        assert ci["cause_of_interest"] == "1"

    def test_naive_one_minus_km_overstates_absolute_risk(self):
        """The whole reason §4.7 asks for competing risks."""
        ci = sv.cumulative_incidence([1, 2, 3, 4], [1, 1, 1, 0], [1, 2, 1, 0])
        cif = ci["causes"]["1"]["cif"][-1]
        # Naive: treat the competing event as censoring, take 1 - KM.
        km = sv.kaplan_meier([1, 2, 3, 4], [1, 0, 1, 0])
        naive = 1.0 - km["survival"][-1]
        assert naive == pytest.approx(0.625)
        assert cif == pytest.approx(0.5)
        assert cif < naive

    def test_cifs_and_survival_always_sum_to_one(self):
        rng = np.random.default_rng(909)
        n = 400
        t = rng.exponential(3.0, n)
        e = rng.binomial(1, 0.75, n)
        ctype = np.where(rng.random(n) < 0.6, 1, 2)
        ctype = np.where(e == 1, ctype, 0)
        ci = sv.cumulative_incidence(t, e, ctype)
        total = np.asarray(ci["overall_survival"], dtype=float)
        for c in ("1", "2"):
            total = total + np.asarray(ci["causes"][c]["cif"], dtype=float)
        assert np.allclose(total, 1.0, atol=1e-10)

    def test_bands_are_within_unit_interval_and_bracket_the_estimate(self):
        rng = np.random.default_rng(910)
        n = 300
        t = rng.exponential(3.0, n)
        e = rng.binomial(1, 0.8, n)
        ctype = np.where(e == 1, np.where(rng.random(n) < 0.5, 1, 2), 0)
        ci = sv.cumulative_incidence(t, e, ctype)
        f = np.asarray(ci["causes"]["1"]["cif"])
        lo = np.asarray(ci["causes"]["1"]["ci_lower"])
        hi = np.asarray(ci["causes"]["1"]["ci_upper"])
        assert np.all(lo >= -1e-12) and np.all(hi <= 1.0 + 1e-12)
        assert np.all(lo <= f + 1e-9) and np.all(hi >= f - 1e-9)

    def test_grouped_by_carrier_status(self):
        rng = np.random.default_rng(911)
        n = 400
        carrier = rng.binomial(1, 0.5, n)
        t = rng.exponential(np.where(carrier == 1, 1.5, 5.0))
        e = rng.binomial(1, 0.85, n)
        ctype = np.where(e == 1, np.where(rng.random(n) < 0.7, 1, 2), 0)
        ci = sv.cumulative_incidence(t, e, ctype,
                                     groups=np.where(carrier == 1, "carrier", "non"))
        assert sorted(ci["groups"]) == ["carrier", "non"]
        for g in ci["groups"].values():
            assert len(g["n_risk"]) == len(g["time"])
            assert "1" in g["causes"]
        # Carriers reach a higher cumulative incidence of the event of interest.
        assert ci["groups"]["carrier"]["causes"]["1"]["cif"][-1] > \
            ci["groups"]["non"]["causes"]["1"]["cif"][-1]

    def test_explicit_cause_of_interest(self):
        ci = sv.cumulative_incidence([1, 2, 3, 4], [1, 1, 1, 0], [1, 2, 1, 0], cause=2)
        assert ci["cause_of_interest"] == "2"
        assert "1" in ci["causes"] and "2" in ci["causes"]

    def test_no_competing_event_reduces_to_one_minus_km(self):
        t = [1, 2, 3, 4, 5]
        e = [1, 1, 0, 1, 1]
        ci = sv.cumulative_incidence(t, e, [1, 1, 0, 1, 1])
        km = sv.kaplan_meier(t, e)
        cif = np.asarray(ci["causes"]["1"]["cif"])
        surv = np.asarray(km["survival"])
        assert np.allclose(cif, 1.0 - surv, atol=1e-12)
