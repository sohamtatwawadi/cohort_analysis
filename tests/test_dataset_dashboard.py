"""The dataset dashboard.

Two things here are worth a test rather than a glance, because both were wrong
in the first version and neither was visible from the JSON:

  the denominator   a gene's carrier rate divides by the samples CALLED at that
                    gene, never by cohort size
  the threshold     ranking on ALL variants put every gene at 100%, because with
                    common variants included every sample carries something
                    everywhere — a panel that says nothing
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from backend.app import config, db  # noqa: E402
from backend.app.research import dashboard  # noqa: E402


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """A real registered dataset, built through the fixture generator.

    config.DB_PATH is process-global and conftest builds ONE lab store for the
    whole session, so pointing it elsewhere and walking away leaves every module
    that runs afterwards looking at the wrong database. This one holds research
    data and no lab data, so the symptom was eight unrelated failures reporting
    an empty store — which is why the original is saved and put back.
    """
    original_db = config.DB_PATH
    config.DB_PATH = tmp_path_factory.mktemp("dash") / "d.duckdb"
    db.close()
    try:
        from backend.tools.make_research_fixture import main
        rc = main(["--samples", "300", "--variants", "4000", "--seed", "3",
                   "--name", "Dashboard fixture"])
        assert rc == 0, "fixture generation failed"
        row = db.row("SELECT dataset_id FROM dataset WHERE deleted_at IS NULL")
        yield row["dataset_id"]
    finally:
        db.close()
        config.DB_PATH = original_db


@pytest.fixture(scope="module")
def dash(built):
    return dashboard.get(built, refresh=True)


# ------------------------------------------------------------ denominators --
def test_gene_rate_divides_by_samples_called_not_cohort_size(dash, built):
    """The rule the whole product turns on. A sample with no call at a gene is
    not evidence of absence, so it cannot sit in the denominator."""
    from backend.app.research import store
    gm = store.load_genotypes(built)
    assert dash["top_genes"], "no genes ranked"
    for g in dash["top_genes"]:
        assert g["n_called"] <= gm.n_samples
        assert g["carriers"] <= g["n_called"], (
            "{}: {} carriers among {} called".format(
                g["gene"], g["carriers"], g["n_called"]))
        assert g["pct"] == pytest.approx(g["carriers"] / g["n_called"])


def test_every_reported_gene_clears_the_minimum_called(dash):
    for g in dash["top_genes"]:
        assert g["n_called"] >= dashboard.MIN_CALLED


def test_small_carrier_counts_are_suppressed(dash):
    """Small-cell suppression, same rule as the lab side."""
    for g in dash["top_genes"]:
        if 0 < g["carriers"] < dashboard.MIN_CELL:
            assert g["suppressed"], g


# ------------------------------------------------------------- the ranking --
def test_ranking_is_not_saturated(dash):
    """The defect that made the first version useless: counting common variants
    too put every gene at 100%, so the panel conveyed nothing."""
    pcts = [g["pct"] for g in dash["top_genes"]]
    assert pcts, "no genes ranked"
    assert max(pcts) < 1.0, (
        "every gene is at 100% — the frequency threshold is not being applied")
    assert pcts == sorted(pcts, reverse=True), "not ordered by rate"
    # A panel where the top and bottom are indistinguishable is also useless.
    if len(pcts) > 3:
        assert max(pcts) > min(pcts) * 1.5, "no spread across the ranked genes"


def test_only_variants_under_the_threshold_are_counted(dash, built):
    """A gene's carrier count must not exceed what its qualifying variants can
    produce — the check that the frequency filter is actually applied."""
    from backend.app.research import store
    from backend.app.research.types import MISSING
    gm = store.load_genotypes(built)
    ann = store.load_annotations(built) or {}

    d = gm.dosages
    called = d != MISSING
    n_called_v = called.sum(axis=1)
    alt = np.where(called, np.maximum(d, 0), 0).sum(axis=1)
    af = np.where(n_called_v > 0, alt / (2.0 * np.maximum(n_called_v, 1)), np.nan)

    by_gene = {}
    for i, v in enumerate(gm.variants):
        gene = (ann.get(v.key) or {}).get("gene")
        if gene:
            by_gene.setdefault(str(gene), []).append(i)

    for g in dash["top_genes"]:
        idx = [i for i in by_gene[g["gene"]]
               if np.isfinite(af[i]) and 0 < af[i] <= dashboard.GENE_MAX_AF]
        block = gm.dosages[np.array(idx)]
        expected = int(np.any(block >= 1, axis=0).sum())
        assert g["carriers"] == expected, (
            "{}: reported {} carriers, qualifying variants give {}".format(
                g["gene"], g["carriers"], expected))


# ----------------------------------------------------------------- labels --
def test_codes_are_rendered_as_words(dash):
    """A composition panel reading "1.0 / 2.0" is data the reader has to
    decode. The phenotype file stores numbers; the display must not."""
    panels = {p["label"]: p for p in dash["composition"]}
    assert panels, "no composition panels"
    for p in panels.values():
        for g in p["groups"]:
            assert g["name"] not in ("1.0", "0.0", "2.0", "1", "0", "2"), (
                "{} still shows a raw code: {}".format(p["label"], g["name"]))

    if "Sex" in panels:
        names = {g["name"] for g in panels["Sex"]["groups"]}
        assert names <= {"Male", "Female", "Unknown"}, names


def test_acronyms_are_not_sentence_cased():
    assert dashboard._title("ldl") == "LDL"
    assert dashboard._title("followup_years") == "Followup Years"
    assert dashboard._title("affected") == "Affected"


def test_composition_percentages_sum_to_one(dash):
    for p in dash["composition"]:
        total = sum(g["count"] for g in p["groups"])
        assert total == p["total"]
        assert sum(g["pct"] for g in p["groups"]) == pytest.approx(1.0)


# ------------------------------------------------------------------ caching --
def test_the_result_is_cached(built):
    """Top genes walks the whole genotype matrix; recomputing per page view
    would make this the most expensive screen in the product."""
    first = dashboard.get(built, refresh=True)
    again = dashboard.get(built)
    assert again["computed_at"] == first["computed_at"], "not served from cache"

    dashboard.invalidate(built)
    fresh = dashboard.get(built)
    assert fresh["top_genes"] == first["top_genes"], "recompute changed the result"


def test_unknown_dataset_raises():
    with pytest.raises(ValueError, match="unknown dataset"):
        dashboard.compute("ds-nope")
