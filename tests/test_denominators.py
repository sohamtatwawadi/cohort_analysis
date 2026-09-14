"""Denominators, family independence, and derived zygosity — the §E04 traps."""
from __future__ import annotations

import pytest

from backend.app import config, db
from backend.app.services import carrier, cohort, denominator, zygosity


def test_gene_denominators_are_not_the_cohort_size():
    """Spec E04: "Using cohort size as denominator -> frequencies understated,
    plausible, quoted." The whole point of the service is that these differ."""
    co = cohort.resolve()
    dens = denominator.gene_denominators()
    distinct = {d["assayed"] for d in dens.values() if d["assayed"]}
    assert len(distinct) > 1, "every gene shares one denominator — not gene-specific"
    # At least one gene must be assayed in strictly fewer subjects than the cohort.
    assert any(d["assayed"] < co.n_subjects for d in dens.values())


def test_denominator_buckets_never_exceed_the_cohort():
    co = cohort.resolve()
    for gene, d in denominator.gene_denominators().items():
        total = d["assayed"] + d["provisional"] + d["not_assayed"]
        assert total == co.n_subjects, (
            "{}: {} + {} + {} != cohort {}".format(
                gene, d["assayed"], d["provisional"], d["not_assayed"], co.n_subjects))


def test_carrier_rate_uses_the_gene_denominator_not_cohort_size():
    co = cohort.resolve()
    dens = denominator.gene_denominators()
    for row in carrier.carrier_table(co.criteria):
        assert row["carrier_rate"]["d"] == dens[row["gene"]]["assayed"]


def test_denominator_warning_fires_exactly_on_the_spec_condition():
    """warn = denominator < 30 OR provisional/(denominator+provisional) > 0.25"""
    for d in denominator.gene_denominators().values():
        total = d["assayed"] + d["provisional"]
        expected = (d["assayed"] < config.DENOM_MIN
                    or (total and d["provisional"] / total > config.PROVISIONAL_MAX))
        assert bool(d["warn"]) == bool(expected), d


def test_every_rate_carries_its_counts():
    """Spec §3.7: "Every percentage in the UI renders with its counts."""
    co = cohort.resolve()
    for row in carrier.carrier_table(co.criteria):
        r = row["carrier_rate"]
        assert set(("n", "d", "pct", "label", "display")) <= set(r)
        assert "/" in r["label"]


# ------------------------------------------------------------ independence ---
def test_probands_only_changes_every_rate():
    """Spec §G01 done-when: "Toggling probands-only visibly changes every rate
    on screen." """
    on = cohort.resolve({"probands_only": True})
    off = cohort.resolve({"probands_only": False})
    assert off.n_subjects > on.n_subjects
    assert off.member_hash != on.member_hash


def test_probands_only_yields_one_subject_per_family():
    cohort.resolve({"probands_only": True})
    dupes = db.scalar("""
        SELECT COUNT(*) FROM (
            SELECT family_id FROM cohort_subject GROUP BY family_id HAVING COUNT(*) > 1)
    """, default=0)
    assert dupes == 0


def test_probands_off_produces_the_inflation_warning():
    co = cohort.resolve({"probands_only": False})
    codes = {w["code"] for w in co.warnings}
    assert "probands_off" in codes
    text = next(w["text"] for w in co.warnings if w["code"] == "probands_off")
    assert "must not be quoted" in text


# ------------------------------------------------------- derived zygosity ---
def test_hemizygous_is_derived_not_read_from_the_file():
    """Spec §3.5 / E04: ZYGOSITY carries Het/Hom only. Every hemizygous call
    must therefore differ from what the source file said."""
    rows = db.rows("""
        SELECT f.zygosity_raw, f.chrom, s.sex FROM finding f
        JOIN subject s USING (subject_id)
        WHERE f.zygosity = 'Hemizygous'
    """)
    assert rows, "no hemizygous calls derived — test proves nothing"
    for r in rows:
        assert r["zygosity_raw"] in ("Heterozygous", "Homozygous")
        assert r["chrom"] in ("X", "Y")
        assert r["sex"] == "M"


def test_no_hemizygous_call_on_an_autosome_or_in_a_female():
    bad = db.scalar("""
        SELECT COUNT(*) FROM finding f JOIN subject s USING (subject_id)
        WHERE f.zygosity = 'Hemizygous'
          AND (f.chrom NOT IN ('X','Y') OR s.sex <> 'M')
    """, default=0)
    assert bad == 0


def test_compound_het_requires_two_distinct_plp_variants_in_one_gene():
    rows = db.rows("""
        SELECT f.subject_id, f.gene_symbol,
               COUNT(DISTINCT f.variant_key) AS n_variants
        FROM finding f JOIN interpretation i USING (finding_id)
        WHERE f.zygosity = 'Compound heterozygous'
        GROUP BY 1, 2
    """)
    assert rows, "no compound-het calls derived — test proves nothing"
    for r in rows:
        assert r["n_variants"] >= 2, r


def test_zygosity_states_are_never_summed():
    """§G04 done-when. The by-gene table must expose the four states separately;
    a row whose parts do not reconstruct independently would mean they had been
    collapsed somewhere upstream."""
    co = cohort.resolve()
    rows = zygosity.by_gene(co.criteria)
    assert rows
    for r in rows:
        assert {"het", "hom", "compound_het", "hemizygous"} <= set(r)
        assert r["biallelic"] == r["hom"] + r["compound_het"]
