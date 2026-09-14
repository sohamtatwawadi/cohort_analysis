"""Independent hand count — spec E03.9.

    "Someone counts one cohort by hand, independently, and the numbers match.
     Item 9 is the real gate. The other eight can pass while denominators are
     quietly wrong."

These tests recompute the headline figures from the raw tables in plain Python,
without touching the service layer's SQL. If the service and this file agree,
the aggregation is right; if they disagree, one of them is wrong and the
difference says which.

The recomputation deliberately does NOT import any helper from the services
package. Sharing a helper would make both sides wrong in the same way, which is
exactly the failure the spec is warning about.
"""
from __future__ import annotations

from typing import Dict, Set

import pytest

from backend.app import db
from backend.app.services import carrier, cohort, denominator

PLP = {"Pathogenic", "Likely pathogenic"}
ESTABLISHED = {"Definitive", "Strong"}


# ------------------------------------------------------- raw table readers ---
def _raw():
    """Pull the raw tables once, as plain Python structures."""
    subjects = {r["subject_id"]: r for r in db.rows("SELECT * FROM subject")}
    runs = {r["subject_id"]: r for r in db.rows("SELECT * FROM run")}
    scope = db.rows("SELECT run_id, gene_symbol, reportable, confidence FROM run_scope")
    findings = db.rows("SELECT * FROM finding")
    interps = {r["finding_id"]: r for r in db.rows("SELECT * FROM interpretation")}
    genes = {r["gene_symbol"]: r for r in db.rows("SELECT * FROM gene_disease")}
    relevance = {}
    for r in db.rows("SELECT indication, gene_symbol FROM indication_gene"):
        relevance.setdefault(r["indication"], set()).add(r["gene_symbol"])
    return subjects, runs, scope, findings, interps, genes, relevance


def _hand_cohort_members(subjects) -> Set[str]:
    """Rebuild the default cohort by hand: consent gate, then probands only."""
    return {sid for sid, s in subjects.items()
            if s["consent_class"] != "Withdrawn" and s["is_proband"]}


def test_cohort_membership_matches_a_hand_count():
    subjects, *_ = _raw()
    expected = _hand_cohort_members(subjects)

    co = cohort.resolve()
    actual = {r["subject_id"] for r in db.rows("SELECT subject_id FROM cohort_subject")}

    assert actual == expected, (
        "hand count {} vs tool {} (missing {}, extra {})".format(
            len(expected), len(actual),
            sorted(expected - actual)[:5], sorted(actual - expected)[:5]))


def test_gene_denominator_matches_a_hand_count():
    subjects, runs, scope, *_ = _raw()
    members = _hand_cohort_members(subjects)

    # By hand: for each gene, the members whose run has a reportable scope row
    # of known confidence.
    by_run: Dict[str, Dict[str, str]] = {}
    for s in scope:
        if s["reportable"]:
            by_run.setdefault(s["run_id"], {})[s["gene_symbol"]] = s["confidence"]

    expected: Dict[str, int] = {}
    expected_prov: Dict[str, int] = {}
    for sid in members:
        run = runs.get(sid)
        if not run:
            continue
        for gene, conf in by_run.get(run["run_id"], {}).items():
            if conf in ("declared", "inferred"):
                expected[gene] = expected.get(gene, 0) + 1
            elif conf == "unknown":
                expected_prov[gene] = expected_prov.get(gene, 0) + 1

    cohort.resolve()
    actual = denominator.gene_denominators()

    for gene, d in actual.items():
        assert d["assayed"] == expected.get(gene, 0), (
            "{}: hand {} vs tool {}".format(gene, expected.get(gene, 0), d["assayed"]))
        assert d["provisional"] == expected_prov.get(gene, 0), (
            "{}: provisional hand {} vs tool {}".format(
                gene, expected_prov.get(gene, 0), d["provisional"]))


def test_carrier_rate_matches_a_hand_count():
    subjects, runs, scope, findings, interps, genes, _ = _raw()
    members = _hand_cohort_members(subjects)

    # Numerator by hand: distinct member subjects with a reportable P/LP per gene.
    plp_subjects: Dict[str, Set[str]] = {}
    for f in findings:
        if f["subject_id"] not in members:
            continue
        i = interps.get(f["finding_id"])
        if not i or not i["reportable"] or i["classification"] not in PLP:
            continue
        plp_subjects.setdefault(f["gene_symbol"], set()).add(f["subject_id"])

    co = cohort.resolve()
    table = {r["gene"]: r for r in carrier.carrier_table(co.criteria)}

    for gene, subs in plp_subjects.items():
        assert gene in table, "tool omitted gene {} with {} carriers".format(gene, len(subs))
        assert table[gene]["plp_subjects"] == len(subs), (
            "{}: hand {} vs tool {}".format(gene, len(subs), table[gene]["plp_subjects"]))

    # And nothing in the tool's table that the hand count did not find.
    for gene, row in table.items():
        assert row["plp_subjects"] == len(plp_subjects.get(gene, set()))


def test_diagnostic_yield_matches_a_hand_count():
    """The number clients scrutinise hardest, recomputed from scratch.

    solved = EXISTS observation where classification in P/LP
             AND validity in (Definitive, Strong)
             AND reportable
             AND gene is phenotype-relevant for the subject's indication
    """
    subjects, runs, scope, findings, interps, genes, relevance = _raw()
    members = _hand_cohort_members(subjects)

    solved: Set[str] = set()
    for f in findings:
        sid = f["subject_id"]
        if sid not in members:
            continue
        i = interps.get(f["finding_id"])
        if not i or not i["reportable"] or i["classification"] not in PLP:
            continue
        g = genes.get(f["gene_symbol"])
        if not g or g["validity"] not in ESTABLISHED:
            continue
        indication = subjects[sid]["indication"]
        if f["gene_symbol"] not in relevance.get(indication, set()):
            continue
        solved.add(sid)

    co = cohort.resolve()
    tool = carrier.diagnostic_yield(co.criteria)

    assert tool["total"] == len(members)
    assert tool["solved"] == len(solved), (
        "yield hand count {} vs tool {}".format(len(solved), tool["solved"]))


def test_subject_result_classification_matches_the_yield_computation():
    """Spec §G12 done-when: "The result classification matches the yield
    computation in G02 exactly." """
    co = cohort.resolve()
    y = carrier.diagnostic_yield(co.criteria)
    results = carrier.subject_results(co.criteria)

    tally: Dict[str, int] = {}
    for v in results.values():
        tally[v] = tally.get(v, 0) + 1

    assert tally.get("P/LP", 0) == y["solved"]
    assert tally.get("VUS-only", 0) == y["vus_only"]
    assert tally.get("Negative", 0) == y["negative"]
    assert sum(tally.values()) == co.n_subjects


def test_relevance_gate_actually_suppresses_something():
    """If every P/LP were phenotype-relevant, the gate would be inert and the
    yield test above would pass without testing the gate at all."""
    subjects, runs, scope, findings, interps, genes, relevance = _raw()
    members = _hand_cohort_members(subjects)

    any_established_plp: Set[str] = set()
    relevant_plp: Set[str] = set()
    for f in findings:
        sid = f["subject_id"]
        if sid not in members:
            continue
        i = interps.get(f["finding_id"])
        if not i or not i["reportable"] or i["classification"] not in PLP:
            continue
        g = genes.get(f["gene_symbol"])
        if not g or g["validity"] not in ESTABLISHED:
            continue
        any_established_plp.add(sid)
        if f["gene_symbol"] in relevance.get(subjects[sid]["indication"], set()):
            relevant_plp.add(sid)

    assert len(relevant_plp) < len(any_established_plp), (
        "phenotype relevance excluded nobody — the gate is inert, and an "
        "unfiltered yield would be indistinguishable from a gated one")
