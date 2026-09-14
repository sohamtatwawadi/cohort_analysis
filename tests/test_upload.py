"""Browser upload.

Until now the only way to get data in was to put it on the server yourself and
type the path. These tests cover the multipart route the UI actually uses.

The bug worth pinning: the endpoint took a single `UploadFile`, so posting a
PLINK triple kept only the LAST part and discarded the other two — no error,
just a fileset that failed to read a moment later for a reason that pointed at
the wrong thing.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from fastapi.testclient import TestClient  # noqa: E402

from backend.app import config, db  # noqa: E402
from backend.app.ingest import synthetic  # noqa: E402
from test_research_formats import write_plink  # noqa: E402


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    config.DB_PATH = tmp_path_factory.mktemp("upload") / "u.duckdb"
    db.close()
    synthetic.rebuild(n_families=20, seed=7)
    from backend.app.api import research_routes
    research_routes.UPLOAD_DIR = tmp_path_factory.mktemp("uploads")
    from backend.app.main import app
    with TestClient(app) as c:
        yield c
    db.close()


@pytest.fixture(scope="module")
def project(client):
    r = client.post("/api/research/projects", json={"name": "Upload tests"})
    assert r.status_code < 300, r.text
    return r.json()["project_id"]


def _vcf(tmp_path: Path) -> Path:
    p = tmp_path / "t.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
        "1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\n"
        "1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1\t1/1\n"
        "1\t300\t.\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0\t0/0\n")
    return p


def _upload(client, project, files, **form):
    body = {"name": "ds", "consent_attested": "true"}
    body.update(form)
    return client.post(
        "/api/research/projects/{}/datasets/upload".format(project),
        data=body, files=files)


def test_vcf_uploads_from_the_browser(client, project, tmp_path):
    v = _vcf(tmp_path)
    r = _upload(client, project,
                [("genotypes", ("t.vcf", v.read_bytes(), "application/octet-stream"))],
                genome_build="GRCh38")
    assert r.status_code < 300, r.text
    ds = r.json()["dataset"]
    assert ds["source_format"] == "VCF"
    assert ds["n_samples"] == 4 and ds["n_variants"] == 3


def test_plink_triple_survives_the_upload(client, project, tmp_path):
    """All three parts must be written. With a single UploadFile the first two
    were silently dropped and only the last survived."""
    write_plink(tmp_path, "array")
    files = [("genotypes", (n, (tmp_path / n).read_bytes(), "application/octet-stream"))
             for n in ("array.fam", "array.bim", "array.bed")]
    r = _upload(client, project, files, genome_build="GRCh37")
    assert r.status_code < 300, r.text
    ds = r.json()["dataset"]
    assert ds["source_format"] == "PLINK1"
    assert ds["n_samples"] == 6 and ds["n_variants"] == 3


def test_the_bed_is_the_entry_point_whatever_order_it_arrives_in(client, project,
                                                                tmp_path):
    """A file picker hands over files in whatever order the OS dialog chose, so
    the endpoint cannot rely on .bed being first."""
    write_plink(tmp_path, "arr2")
    names = ["arr2.bim", "arr2.fam", "arr2.bed"]
    for order in (names, list(reversed(names))):
        files = [("genotypes", (n, (tmp_path / n).read_bytes(),
                                "application/octet-stream")) for n in order]
        r = _upload(client, project, files, genome_build="GRCh37")
        assert r.status_code < 300, "order {}: {}".format(order, r.text)
        assert r.json()["dataset"]["source_format"] == "PLINK1"


def test_an_orphan_bim_is_rejected_rather_than_half_read(client, project, tmp_path):
    write_plink(tmp_path, "lone")
    r = _upload(client, project,
                [("genotypes", ("lone.bim", (tmp_path / "lone.bim").read_bytes(),
                                "application/octet-stream"))])
    assert r.status_code == 400
    assert "incomplete" in r.json()["detail"].lower()


def test_consent_is_required_before_anything_is_read(client, project, tmp_path):
    v = _vcf(tmp_path)
    r = client.post(
        "/api/research/projects/{}/datasets/upload".format(project),
        data={"name": "no consent"},
        files=[("genotypes", ("t.vcf", v.read_bytes(), "application/octet-stream"))])
    assert r.status_code == 400
    assert "consent" in r.json()["detail"].lower()


def test_a_path_traversing_filename_cannot_escape_the_upload_directory(
        client, project, tmp_path):
    """Filenames come from the client and are not trustworthy. `Path(...).name`
    is what stops `../../` from landing somewhere it should not."""
    from backend.app.api import research_routes
    v = _vcf(tmp_path)
    r = _upload(client, project,
                [("genotypes", ("../../escaped.vcf", v.read_bytes(),
                                "application/octet-stream"))],
                genome_build="GRCh38")
    assert r.status_code < 300, r.text
    assert (research_routes.UPLOAD_DIR / "escaped.vcf").exists()
    assert not (research_routes.UPLOAD_DIR.parent.parent / "escaped.vcf").exists()


def test_phenotypes_ride_along_and_unlock_more(client, project, tmp_path):
    v = _vcf(tmp_path)
    ph = tmp_path / "p.csv"
    ph.write_text("sample_id,status,age\nS1,1,50\nS2,0,61\nS3,1,44\nS4,0,58\n")
    r = _upload(client, project,
                [("genotypes", ("t.vcf", v.read_bytes(), "application/octet-stream"))],
                genome_build="GRCh38")
    without = r.json()["dataset"]["dataset_id"]

    r2 = client.post(
        "/api/research/projects/{}/datasets/upload".format(project),
        data={"name": "with pheno", "consent_attested": "true",
              "genome_build": "GRCh38"},
        files=[("genotypes", ("t.vcf", v.read_bytes(), "application/octet-stream")),
               ("phenotypes", ("p.csv", ph.read_bytes(), "text/csv"))])
    assert r2.status_code < 300, r2.text
    assert r2.json()["dataset"]["dataset_id"] != without
