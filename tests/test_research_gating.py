"""Capability gating, validation and data governance — Part II §2.3, §3, §8.

§3 calls capability gating "the most important safety feature in the product",
so these tests are about what the platform REFUSES to do. The statistical
kernels are tested elsewhere; this file tests the gate in front of them.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from backend.app import config, db  # noqa: E402
from backend.app.research import capability, jobs, registry  # noqa: E402
from backend.app.research import store as rstore  # noqa: E402
from backend.app.research import analyses  # noqa: E402,F401
from backend.app.research.profile import DataProfile, profile_dataset  # noqa: E402
from backend.app.research.types import GenotypeMatrix, PhenotypeTable, Variant  # noqa: E402
from backend.app.research.validate import validate_upload  # noqa: E402


# ------------------------------------------------------------------ helpers --
def make_gm(n_samples=200, n_variants=400, n_chroms=3, build="GRCh38", seed=5):
    rng = np.random.default_rng(seed)
    chroms = [str(1 + (i % n_chroms)) for i in range(n_variants)]
    variants = [Variant(chrom=chroms[i], pos=1000 + i, ref="A", alt="G")
                for i in range(n_variants)]
    return GenotypeMatrix(
        sample_ids=["S{:04d}".format(i) for i in range(n_samples)],
        variants=variants,
        dosages=rng.integers(0, 3, (n_variants, n_samples)).astype(np.int8),
        build=build)


def make_ph(gm, seed=6, with_controls=True):
    rng = np.random.default_rng(seed)
    y = rng.integers(0, 2, gm.n_samples).astype(float)
    if not with_controls:
        y[:] = 1.0
    return PhenotypeTable(
        sample_ids=list(gm.sample_ids),
        columns={"affected": y, "ldl": rng.normal(3, 1, gm.n_samples),
                 "age": rng.normal(50, 10, gm.n_samples)},
        kinds={"affected": "binary", "ldl": "quantitative", "age": "quantitative"})


@pytest.fixture(scope="module", autouse=True)
def store(tmp_path_factory):
    config.DB_PATH = tmp_path_factory.mktemp("research") / "r.duckdb"
    db.close()
    db.connect()
    import backend.app.research.store as rs
    rs.RESEARCH_DIR = tmp_path_factory.mktemp("payload")
    yield
    db.close()


@pytest.fixture
def project():
    return registry.create_project("Test project", "tester")["project_id"]


# =========================================================== §2.3 validation ==
def test_upload_without_a_declared_build_is_rejected():
    """"Genome build determined and recorded — never inferred from position
    alone." GRCh37 and GRCh38 share a coordinate range; guessing is a coin flip
    that silently corrupts every downstream annotation."""
    gm = make_gm(build=None)
    report, _ = validate_upload(gm, declared_build=None)
    assert not report.ok
    assert any(i.code == "build_unknown" for i in report.failures)


def test_conflicting_build_is_rejected_rather_than_resolved():
    gm = make_gm(build="GRCh38")
    report, _ = validate_upload(gm, declared_build="GRCh37")
    assert not report.ok
    assert any(i.code == "build_conflict" for i in report.failures)


def test_duplicate_sample_ids_fail_the_upload():
    gm = make_gm(n_samples=10)
    gm.sample_ids[3] = gm.sample_ids[0]
    report, _ = validate_upload(gm, declared_build="GRCh38")
    assert not report.ok
    assert any(i.code == "duplicate_sample_ids" for i in report.failures)


def test_non_overlapping_sample_ids_fail_the_upload():
    """Reported in BOTH directions — only reporting one side hides the case
    where the phenotype file is a different cohort entirely."""
    gm = make_gm(n_samples=20)
    ph = PhenotypeTable(sample_ids=["OTHER-{}".format(i) for i in range(20)],
                        columns={"affected": np.zeros(20)},
                        kinds={"affected": "binary"})
    report, facts = validate_upload(gm, declared_build="GRCh38", phenotypes=ph)
    assert not report.ok
    issue = next(i for i in report.failures if i.code == "no_sample_overlap")
    assert issue.detail["genotype_only"]
    assert issue.detail["phenotype_only"]


def test_chromosome_naming_is_harmonised():
    from backend.app.research.validate import harmonise_chrom
    assert harmonise_chrom("chr1") == harmonise_chrom("1") == "1"
    assert harmonise_chrom("chrX") == "X"
    assert harmonise_chrom("chrM") == "MT"


def test_indels_are_left_aligned():
    from backend.app.research.validate import left_align
    # Shared suffix then shared prefix, position shifts with prefix trimming.
    assert left_align("AT", "ATT", 100) == ("A", "AT", 100)
    assert left_align("CTT", "CT", 100) == ("CT", "C", 100)
    assert left_align("A", "G", 100) == ("A", "G", 100)


def test_a_valid_upload_passes_every_mandatory_check():
    gm = make_gm()
    ph = make_ph(gm)
    report, facts = validate_upload(gm, declared_build="GRCh38", phenotypes=ph)
    assert report.ok, report.as_dict()
    for check in ("genome_build", "duplicate_sample_ids", "variant_normalisation",
                  "sample_reconciliation"):
        assert check in report.checks_run


# ============================================================ §3 capability ===
def test_small_targeted_panel_locks_the_large_n_analyses():
    """The spec's own cautionary example: a GWAS on 180 ascertained samples
    must not be runnable."""
    gm = make_gm(n_samples=180, n_variants=400, n_chroms=3)
    prof = profile_dataset(gm, make_ph(gm), compute_genetics=False)
    caps = {c.analysis: c for c in capability.assess(prof)}
    assert not caps["gwas"].available
    assert not caps["phewas"].available
    assert not caps["burden"].available
    # And it explains why, with the observed value.
    unmet = {r.label: r.observed for r in caps["gwas"].unmet}
    assert any("Genome-wide" in k for k in unmet)
    assert any("180" in v for v in unmet.values())


def test_descriptive_analyses_stay_available_on_a_small_panel():
    """Gating must not become feature removal by another name — §3's whole
    argument is against shipping nothing."""
    gm = make_gm(n_samples=180, n_variants=400)
    prof = profile_dataset(gm, make_ph(gm), compute_genetics=False)
    caps = {c.analysis: c for c in capability.assess(prof)}
    assert caps["carrier_frequency"].available
    assert caps["zygosity"].available


def test_every_unmet_requirement_states_the_observed_value():
    """§3.3 — a locked card must say what you have, not just what you need."""
    gm = make_gm(n_samples=180, n_variants=400)
    prof = profile_dataset(gm, make_ph(gm), compute_genetics=False)
    for cap in capability.assess(prof):
        for req in cap.unmet:
            assert req.observed, "{}/{} has no observed value".format(
                cap.analysis, req.label)


def test_missing_data_requirements_are_not_overridable():
    """"Missing data cannot be overridden — it has to be supplied." You cannot
    override your way to genotypes you did not upload."""
    gm = make_gm(n_samples=180, n_variants=400)
    prof = profile_dataset(gm, make_ph(gm), compute_genetics=False)
    caps = {c.analysis: c for c in capability.assess(prof)}
    assert not capability.overridable(caps["gwas"])["can_override"]
    assert not capability.overridable(caps["segregation"])["can_override"]


# ============================================== §3.4 override + §8 governance ==
def test_registration_requires_consent_attestation(project):
    gm = make_gm()
    prof = profile_dataset(gm, make_ph(gm), compute_genetics=False)
    with pytest.raises(PermissionError):
        registry.register_dataset(project, "no consent", prof, "VCF", [], "GRCh38",
                                  "tester", consent_attested=False)


def _register(project, n_samples=180, n_variants=400):
    gm = make_gm(n_samples=n_samples, n_variants=n_variants)
    ph = make_ph(gm)
    prof = profile_dataset(gm, ph, compute_genetics=False)
    ds = registry.register_dataset(project, "ds", prof, "VCF", ["x.vcf"], "GRCh38",
                                   "tester", consent_attested=True,
                                   phenotype_kinds=ph.kinds)
    rstore.save_genotypes(ds["dataset_id"], gm)
    rstore.save_phenotypes(ds["dataset_id"], ph)
    return ds["dataset_id"], gm, ph


def test_locked_analysis_cannot_be_submitted(project):
    did, _, _ = _register(project)
    with pytest.raises(PermissionError) as exc:
        jobs.submit(project, did, "gwas", {"outcome": "affected"})
    # The refusal names what is missing, per §3.3.
    assert "Genome-wide" in str(exc.value) or "variants" in str(exc.value)


def test_override_requires_a_real_justification(project):
    did, _, _ = _register(project)
    with pytest.raises(ValueError):
        registry.create_override(did, "phewas", "ok", "tester")


def test_override_of_missing_data_is_refused(project):
    did, _, _ = _register(project)
    with pytest.raises(PermissionError):
        registry.create_override(
            did, "gwas", "We accept the limits of this exome for exploratory work.",
            "tester")


def test_override_unlocks_and_stamps(project):
    """§3.4 — the escape hatch exists, is recorded, and marks the output."""
    did, _, _ = _register(project)
    before = registry.can_run(did, "phewas")
    assert not before["allowed"]

    registry.create_override(
        did, "phewas", "Exploratory only; reduced power accepted and not quoted.",
        "tester")
    after = registry.can_run(did, "phewas")
    assert after["allowed"]
    assert after["underpowered"] is True
    assert after["justification"]


def test_underpowered_stamp_reaches_the_result(project):
    """The stamp is derived from the override record and attached by the job
    runner, so it cannot be edited off a result payload."""
    did, _, _ = _register(project, n_samples=180)
    registry.create_override(
        did, "association", "Pilot; n below threshold, reporting descriptively only.",
        "tester")
    job = jobs.submit(project, did, "association",
                      {"outcome": "affected", "min_maf": 0.1})
    final = jobs.wait(job["job_id"], timeout=300)
    assert final["status"] == "complete", final.get("error")
    assert final["underpowered"] is True

    res = jobs.get_result(job["job_id"])
    assert res["underpowered"] is True
    assert res["result"]["underpowered_stamp"] == capability.UNDERPOWERED_STAMP
    assert res["result"]["override_justification"]


def test_project_isolation_hides_other_projects_datasets(project):
    """§8 — isolation is a column on every row, not an application-layer
    filter a caller can forget."""
    did, _, _ = _register(project)
    other = registry.create_project("Other project", "someone-else")["project_id"]
    assert registry.get_dataset(did, project) is not None
    assert registry.get_dataset(did, other) is None
    assert did not in [d["dataset_id"] for d in registry.list_datasets(other)]


def test_deletion_propagates_to_every_derived_table(project):
    """§8: "Deletion must actually propagate to derived tables, caches and
    result stores." A flag on `dataset` alone would leave the profile, the
    capability assessment and every stored result in place."""
    did, _, _ = _register(project)
    registry.create_override(
        did, "association", "Pilot; accepted reduced power for this exploratory run.",
        "tester")
    job = jobs.submit(project, did, "association",
                      {"outcome": "affected", "min_maf": 0.2})
    jobs.wait(job["job_id"], timeout=300)

    assert db.scalar("SELECT COUNT(*) FROM dataset_profile WHERE dataset_id = ?",
                     [did], default=0) == 1
    assert db.scalar("SELECT COUNT(*) FROM analysis_result WHERE job_id = ?",
                     [job["job_id"]], default=0) == 1

    out = registry.delete_dataset(did, project)
    assert out["deleted"]

    for table in ("dataset_profile", "dataset_capability", "dataset_phenotype",
                  "capability_override"):
        assert db.scalar("SELECT COUNT(*) FROM {} WHERE dataset_id = ?".format(table),
                         [did], default=0) == 0, table
    assert db.scalar("SELECT COUNT(*) FROM analysis_job WHERE dataset_id = ?",
                     [did], default=0) == 0
    assert db.scalar("SELECT COUNT(*) FROM analysis_result WHERE job_id = ?",
                     [job["job_id"]], default=0) == 0
    assert registry.get_dataset(did, project) is None
    # And the payload itself, not merely the metadata describing it.
    assert not (rstore.dataset_dir(did) / "genotypes.npz").exists()


# ================================================================ §6 compute ==
def test_job_records_reproducibility_information(project):
    """§6 — software versions, reference genome, parameters, seed, dataset."""
    did, _, _ = _register(project)
    registry.create_override(
        did, "association", "Pilot run; below the n threshold, accepted knowingly.",
        "tester")
    job = jobs.submit(project, did, "association",
                      {"outcome": "affected", "min_maf": 0.2})
    import json
    repro = json.loads(jobs.get_job(job["job_id"])["reproducibility_json"])
    assert repro["software"]["numpy"]
    assert repro["reference_genome"] == "GRCh38"
    assert repro["random_seed"] is not None
    assert repro["dataset"]["dataset_id"] == did


def test_cost_estimate_scales_with_the_work(project):
    small = jobs.estimate_cost("association", 500, 1000)
    big = jobs.estimate_cost("association", 50000, 500000)
    assert big["estimated_seconds"] > small["estimated_seconds"] * 100
    assert small["compute_profile"] == "small"
