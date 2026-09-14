"""The gates that must not be bypassable — spec §G01, §G09, E04.

These are the tests that matter most. Every other number in the tool can be
recomputed and argued about; a consent leak cannot be undone.
"""
from __future__ import annotations

import itertools

import pytest

from backend.app import db
from backend.app.reference.tests import CONSENT_SF_OK, CONSENT_WITHDRAWN
from backend.app.services import analysis, carrier, cohort


def _withdrawn_ids():
    return {r["subject_id"] for r in db.rows(
        "SELECT subject_id FROM subject WHERE consent_class = ?", [CONSENT_WITHDRAWN])}


def test_store_actually_contains_withdrawn_subjects():
    """Guard the guard: if no subject is withdrawn, the consent tests below pass
    vacuously and prove nothing."""
    assert len(_withdrawn_ids()) > 0


def test_withdrawn_subjects_unreachable_by_any_filter_combination():
    """Spec §G01: "A withdrawn-consent subject must not be reachable by any
    filter combination." Enumerated rather than asserted once, because the
    claim is about the whole criteria space, not one path."""
    withdrawn = _withdrawn_ids()

    # Criteria chosen to pull in as many directions as possible, including the
    # ones that try to select withdrawn subjects explicitly.
    probes = [
        {},
        {"probands_only": False},
        {"probands_only": False, "reportable_only": False},
        {"consent_class": [CONSENT_WITHDRAWN]},
        {"consent_class": [CONSENT_WITHDRAWN], "probands_only": False},
        {"consent_class": [CONSENT_WITHDRAWN], "reportable_only": False,
         "probands_only": False},
        {"relation": ["Proband", "Mother", "Father", "Sibling", "Child"],
         "probands_only": False},
        {"classification": ["Pathogenic", "Likely pathogenic", "Uncertain significance",
                            "Likely benign", "Benign"], "probands_only": False,
         "reportable_only": False},
        {"sf_only": True, "probands_only": False},
    ]
    for crit in probes:
        co = cohort.resolve(crit)
        members = {r["subject_id"] for r in db.rows(
            "SELECT subject_id FROM cohort_subject")}
        leaked = members & withdrawn
        assert not leaked, "criteria {} surfaced withdrawn subjects {}".format(
            crit, sorted(leaked)[:5])


def test_consent_gate_is_unconditional_in_the_funnel():
    """The gate must appear as a mandatory funnel step regardless of criteria —
    it is a property of the query, not a filter the user chose."""
    for crit in ({}, {"probands_only": False}, {"indication": ["Cardiac"]}):
        co = cohort.resolve(crit)
        consent_steps = [s for s in co.funnel if s["group"] == "consent"
                         and s.get("mandatory")]
        assert consent_steps, "no mandatory consent step for {}".format(crit)


def test_sf_cohort_cannot_include_a_declining_subject():
    """Spec §G09 done-when: "Attempting to construct an SF cohort including a
    declining subject is impossible at the query level, verified by test." """
    co = cohort.resolve({"sf_only": True})
    consents = {r["consent_class"] for r in db.rows("""
        SELECT DISTINCT s.consent_class FROM cohort_subject cs
        JOIN subject s USING (subject_id)
    """)}
    assert consents, "SF cohort is empty — test proves nothing"
    assert consents <= set(CONSENT_SF_OK), (
        "SF cohort contains non-consenting classes: {}".format(consents - set(CONSENT_SF_OK)))


def test_sf_rate_numerator_is_subset_of_denominator():
    """A rate above 100% means the numerator population is not contained in the
    denominator population. Cheap invariant, catches a whole class of error."""
    co = cohort.resolve()
    sf = analysis.secondary_findings(co.criteria)
    k = sf["kpis"]
    assert k["sf_subjects"] <= k["sf_capable"], (
        "{} subjects with SF but only {} SF-capable".format(
            k["sf_subjects"], k["sf_capable"]))
    assert 0 <= k["rate"]["pct"] <= 100


def test_secondary_findings_agree_between_dashboard_and_module():
    """§G02's SF KPI and §G09's count are the same quantity and must not drift."""
    co = cohort.resolve()
    dash = carrier.dashboard(co.criteria)["kpis"]["sf_subjects"]
    module = analysis.secondary_findings(co.criteria)["kpis"]["sf_subjects"]
    assert dash == module


def test_small_cell_suppression_applied_before_display():
    from backend.app import config as cfg
    from backend.app.services.denominator import small_cell

    assert small_cell(0) == "0"
    assert small_cell(1) == "<{}".format(cfg.SMALL_CELL_THRESHOLD)
    assert small_cell(cfg.SMALL_CELL_THRESHOLD - 1) == "<{}".format(cfg.SMALL_CELL_THRESHOLD)
    assert small_cell(cfg.SMALL_CELL_THRESHOLD) == str(cfg.SMALL_CELL_THRESHOLD)


def test_sf_gene_table_suppresses_small_counts():
    co = cohort.resolve()
    rows = analysis.secondary_findings(co.criteria)["by_gene"]
    for r in rows:
        if int(r["subjects"]) < 5:
            assert r["subjects_display"].startswith("<"), (
                "gene {} exposes a count of {}".format(r["gene"], r["subjects"]))
