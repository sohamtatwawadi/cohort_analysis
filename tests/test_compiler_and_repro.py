"""Query compiler (D01) and reproducibility (C02, E03.6)."""
from __future__ import annotations

import json

import pytest

from backend.app import db
from backend.app.services import cohort, compiler, governance, library


# ------------------------------------------------------------- D01 compiler --
def test_vocabulary_contains_no_patient_data():
    """Spec D01: "Model never sees patient data — only field names and allowed
    values." """
    vocab = compiler.vocabulary()
    blob = json.dumps(vocab)
    # No identifier prefixes from any fact table may appear.
    for prefix in ("SJ-", "GR-", "GO-", "FAM-", "SP-", "MRN", "MG5"):
        assert prefix not in blob, "vocabulary leaked {} identifiers".format(prefix)
    assert set(vocab) >= {"fields", "scalars", "gene_sets", "target_modules", "rules"}


def test_unmapped_field_is_a_hard_failure():
    """"Unmapped terms fail loudly and stop." Not a warning, not a best guess."""
    r = compiler.validate({"profile": "germline",
                           "criteria": {"tumour_purity": [0.4]},
                           "target_module": "g_dashboard"})
    assert r["validation"] == "FAILED"
    assert "tumour_purity" in r["unmapped"]


def test_unmapped_value_within_a_known_field_also_fails():
    r = compiler.validate({"profile": "germline",
                           "criteria": {"classification": ["Probably bad"]},
                           "target_module": "g_dashboard"})
    assert r["validation"] == "FAILED"
    assert any("Probably bad" in u for u in r["unmapped"])


def test_somatic_profile_is_rejected():
    r = compiler.validate({"profile": "somatic", "criteria": {}})
    assert r["validation"] == "FAILED"
    assert any("germline only" in e for e in r["errors"])


def test_payload_containing_a_computed_figure_is_rejected():
    """"Model emits criteria only, never a number or percentage." A payload
    carrying an answer means the model tried to compute instead of compile."""
    r = compiler.validate({
        "profile": "germline", "criteria": {"gene": ["BRCA1"]},
        "target_module": "g_carrier", "answer": "12.4% of subjects",
    })
    assert r["validation"] == "FAILED"
    assert any("criteria only" in e for e in r["errors"])


def test_compile_does_not_execute_anything():
    """"Nothing executes until the user approves the compiled criteria." """
    r = compiler.compile_and_validate("How many BRCA1 P/LP carriers are there?")
    assert r["status"] == "Nothing has run yet."
    assert "compiled_json" in r
    # No figure of any kind in the response.
    blob = json.dumps(r["compiled"])
    assert "count" not in blob and "rate" not in blob


def test_every_compilation_is_logged_with_its_question():
    """"Every compilation logged with its source question" — builds the
    evaluation corpus the acceptance gate is measured on."""
    before = int(db.scalar("SELECT COUNT(*) FROM compilation_log", default=0))
    compiler.compile_and_validate("Carrier rate for HBB in South Asian subjects")
    after = int(db.scalar("SELECT COUNT(*) FROM compilation_log", default=0))
    assert after == before + 1
    row = db.row("SELECT * FROM compilation_log ORDER BY asked_at DESC LIMIT 1")
    assert "HBB" in row["question"]
    assert row["executed"] is False


@pytest.mark.parametrize("question,field,value", [
    ("Show me all VUS", "classification", "Uncertain significance"),
    ("HBOC panel carriers", "gene_set", "HBOC (breast/ovarian)"),
    ("Lynch syndrome P/LP", "gene_set", "Lynch / MMR"),
    ("compound het subjects", "zygosity", "Compound heterozygous"),
    ("hemizygous males", "zygosity", "Hemizygous"),
    ("BRCA1 carriers", "gene", "BRCA1"),
])
def test_rule_compiler_maps_known_phrases(question, field, value):
    r = compiler.compile_question(question)
    got = r["criteria"].get(field)
    assert (value in got) if isinstance(got, list) else (got == value), r


def test_compiled_criteria_are_executable():
    """A PASSED compilation must produce criteria the engine actually accepts."""
    r = compiler.compile_and_validate("HBOC panel P/LP carriers")
    assert r["validation"] == "PASSED"
    co = cohort.resolve(r["compiled"]["criteria"])
    assert co.n_subjects >= 0
    assert co.criteria["gene_set"] == "HBOC (breast/ovarian)"


# -------------------------------------------------- C02 / E03.6 reproducible --
def test_same_criteria_reproduce_the_same_member_hash():
    """Spec E03.6: "Any saved cohort re-runs to an identical number." """
    crit = {"gene_set": "Lynch / MMR", "classification": ["Pathogenic"]}
    a = cohort.resolve(crit)
    b = cohort.resolve(crit)
    assert a.member_hash == b.member_hash
    assert a.criteria_hash == b.criteria_hash
    assert a.n_subjects == b.n_subjects


def test_criteria_hash_is_order_independent():
    """The hash identifies the cohort definition, not the order it was typed."""
    a = cohort.resolve({"gene": ["BRCA1", "BRCA2"]})
    b = cohort.resolve({"gene": ["BRCA2", "BRCA1"]})
    assert a.criteria_hash == b.criteria_hash
    assert a.member_hash == b.member_hash


def test_different_criteria_produce_different_hashes():
    a = cohort.resolve({"gene": ["BRCA1"]})
    b = cohort.resolve({"gene": ["BRCA2"]})
    assert a.criteria_hash != b.criteria_hash


def test_every_builtin_cohort_is_reproducible():
    for c in library.list_cohorts():
        r = library.verify_reproducible(c["cohort_id"])
        assert r["ok"], "cohort {} drifted".format(c["name"])


def test_manifest_carries_every_required_field():
    """Spec C02. Each field earns its place; a missing one breaks reproducibility
    or hides provenance mixing."""
    co = cohort.resolve()
    mf = governance.manifest(co, module="g_carrier")
    required = {
        "cohort_name", "profile", "tenant_id", "resolved_at", "criteria",
        "member_count", "member_hash", "criteria_hash", "kb_snapshot_id",
        "pipeline_versions", "reference_builds", "output_class",
        "small_cell_suppression", "suppression_threshold",
    }
    assert required <= set(mf)
    assert mf["profile"] == "germline"
    assert mf["small_cell_suppression"] is True
    assert len(mf["member_hash"]) == 8
    assert len(mf["criteria_hash"]) == 8


# ------------------------------------------------------- §2.5 output classes --
def test_every_module_declares_an_output_class():
    """"Store it as module metadata so a screen cannot ship without one." """
    for module in governance.MODULE_TITLES:
        assert governance.output_class(module) in ("RESEARCH", "OPERATIONAL")


def test_unknown_module_has_no_class_and_raises():
    with pytest.raises(KeyError):
        governance.output_class("g_not_a_module")


def test_no_module_is_ever_classified_clinical():
    """Spec §2.5: CLINICAL is "never produced by this module"."""
    assert "CLINICAL" not in set(governance.MODULE_CLASS.values())


def test_every_render_writes_an_audit_entry():
    """Spec §2.7."""
    before = int(db.scalar("SELECT COUNT(*) FROM analysis_run", default=0))
    co = cohort.resolve()
    governance.log_run("g_carrier", co)
    after = int(db.scalar("SELECT COUNT(*) FROM analysis_run", default=0))
    assert after == before + 1
    row = db.row("SELECT * FROM analysis_run ORDER BY executed_at DESC LIMIT 1")
    assert row["output_class"] == "RESEARCH"
    assert row["criteria_hash"] == co.criteria_hash


# --------------------------------------------------------------- C04 export --
def test_research_export_is_watermarked_and_embeds_the_manifest(tmp_path):
    from backend.app.services import carrier

    co = cohort.resolve()
    rows = carrier.carrier_table(co.criteria)[:5]
    flat = [{k: v for k, v in r.items() if not isinstance(v, (dict, list))} for r in rows]
    out = governance.export(co, "g_carrier", flat, fmt="csv")

    assert out["output_class"] == "RESEARCH"
    assert out["watermarked"] is True
    text = open(out["path"]).read()
    assert "RESEARCH USE ONLY" in text
    assert "criteria_hash" in text
    assert "kb_snapshot_id" in text
