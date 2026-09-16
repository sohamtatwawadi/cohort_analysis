"""Package a registered dataset — with its completed analyses — as one file.

Why this exists: the demo instance is a t3.medium, and generating a genome-wide
cohort and then running a GWAS over it needs more memory than it has. But the
*results* are small. So do the expensive work on a machine that can afford it,
and move the finished artifacts.

    # where there is RAM
    python -m backend.tools.dataset_bundle export ds-abc123 --out demo.tar.gz

    # on the instance
    python -m backend.tools.dataset_bundle import-bundle demo.tar.gz

What travels: the dataset row, its profile, capability matrix, phenotype
summary, sample list, any overrides, the cached dashboard, every analysis job
and its result — and the genotype payload from disk, without which the dashboard
and every drill-down would fail on arrival.

The table list is taken from registry.DERIVED_TABLES rather than written out
again here, because that list is what deletion walks. If a dataset grows a new
derived table, deletion and this both learn about it in the same place — the
alternative is a bundle that silently omits whatever was added last.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import tarfile
import tempfile
from datetime import date, datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

MANIFEST = "manifest.json"
ROWS = "rows.json"
PAYLOAD = "payload"

# Bumped if the on-disk shape changes, so an old bundle is refused with a clear
# message rather than half-imported.
FORMAT_VERSION = 1


def _json_default(o: Any) -> Any:
    if isinstance(o, (datetime, date)):
        return o.isoformat()
    return str(o)


def _table_columns(table: str) -> List[str]:
    from ..app import db
    return [r["column_name"] for r in db.rows(
        "SELECT column_name FROM information_schema.columns "
        "WHERE table_name = ? ORDER BY ordinal_position", [table])]


def _dump(table: str, where: str, params: List[Any]) -> List[Dict[str, Any]]:
    from ..app import db
    if not _table_columns(table):
        return []
    return [dict(r) for r in db.rows(
        "SELECT * FROM {} WHERE {}".format(table, where), params)]


# ------------------------------------------------------------------ export --
def export_dataset(dataset_id: str, out_path: Path) -> Path:
    from ..app import db
    from ..app.research import registry, store as research_store

    db.connect()
    ds = registry.get_dataset(dataset_id)
    if not ds:
        raise SystemExit("unknown dataset: {}".format(dataset_id))

    tables: Dict[str, List[Dict[str, Any]]] = {}
    tables["dataset"] = _dump("dataset", "dataset_id = ?", [dataset_id])

    # The project row travels too. Importing a dataset whose project does not
    # exist would leave it invisible: project_id is the isolation boundary, and
    # every listing filters on it.
    tables["project"] = _dump("project", "project_id = ?", [ds["project_id"]])

    for table, col in registry.DERIVED_TABLES:
        tables[table] = _dump(table, "{} = ?".format(col), [dataset_id])

    tables["analysis_job"] = _dump("analysis_job", "dataset_id = ?", [dataset_id])
    job_ids = [r["job_id"] for r in tables["analysis_job"]]
    if job_ids:
        marks = ",".join("?" * len(job_ids))
        tables["analysis_result"] = _dump(
            "analysis_result", "job_id IN ({})".format(marks), job_ids)
    else:
        tables["analysis_result"] = []

    # Cached dashboard, if one has been built. Carrying it means the instance
    # never has to walk the genotype matrix to render the landing screen.
    try:
        tables["dataset_dashboard"] = _dump(
            "dataset_dashboard", "dataset_id = ?", [dataset_id])
    except Exception:                                    # noqa: BLE001
        tables["dataset_dashboard"] = []

    complete = sum(1 for j in tables["analysis_job"] if j.get("status") == "complete")
    manifest = {
        "format_version": FORMAT_VERSION,
        "dataset_id": dataset_id,
        "dataset_name": ds.get("name"),
        "project_id": ds["project_id"],
        "exported_at": datetime.utcnow().isoformat(timespec="seconds"),
        "n_samples": ds.get("n_samples"),
        "n_variants": ds.get("n_variants"),
        "genome_build": ds.get("genome_build"),
        "row_counts": {t: len(rs) for t, rs in tables.items()},
        "analyses_complete": complete,
    }

    payload_dir = research_store.dataset_dir(dataset_id)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / MANIFEST).write_text(json.dumps(manifest, indent=2,
                                               default=_json_default))
        (tmp / ROWS).write_text(json.dumps(tables, default=_json_default))
        if payload_dir.is_dir():
            shutil.copytree(payload_dir, tmp / PAYLOAD)
        else:
            print("  WARNING: no genotype payload at {} — the dashboard and "
                  "every drill-down will fail after import".format(payload_dir))

        with tarfile.open(out_path, "w:gz") as tar:
            for item in sorted(tmp.iterdir()):
                tar.add(item, arcname=item.name)

    size = out_path.stat().st_size / 1048576
    print("  {}".format(ds.get("name")))
    for t, n in manifest["row_counts"].items():
        if n:
            print("    {:<22} {}".format(t, n))
    print("    {} of {} analyses complete".format(complete,
                                                  len(tables["analysis_job"])))
    print("\n  wrote {} ({:.0f} MB)".format(out_path, size))
    db.close()
    return out_path


# ------------------------------------------------------------------ import --
def import_bundle(bundle: Path, replace: bool = False) -> str:
    from ..app import db
    from ..app.research import registry, store as research_store

    bundle = Path(bundle)
    if not bundle.is_file():
        raise SystemExit("no such bundle: {}".format(bundle))

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        with tarfile.open(bundle, "r:gz") as tar:
            # Refuse absolute paths and traversal. A bundle is a file that
            # arrived from somewhere; extracting it must not be able to write
            # outside the directory we chose.
            for m in tar.getmembers():
                target = (tmp / m.name).resolve()
                if not str(target).startswith(str(tmp.resolve())):
                    raise SystemExit(
                        "refusing to extract {} — escapes the target directory"
                        .format(m.name))
            tar.extractall(tmp)

        manifest = json.loads((tmp / MANIFEST).read_text())
        if manifest.get("format_version") != FORMAT_VERSION:
            raise SystemExit(
                "bundle format {}, this build reads {}".format(
                    manifest.get("format_version"), FORMAT_VERSION))
        tables = json.loads((tmp / ROWS).read_text())
        dataset_id = manifest["dataset_id"]

        db.connect()
        existing = registry.get_dataset(dataset_id)
        if existing and not replace:
            raise SystemExit(
                "{} is already registered here. Pass --replace to overwrite it."
                .format(dataset_id))
        if existing:
            print("  replacing the existing {}".format(dataset_id))
            registry.delete_dataset(dataset_id, existing["project_id"])

        # Created lazily on first use, so a fresh database has no such table
        # and the cached dashboard would be dropped on the way in.
        try:
            from ..app.research import dashboard as dash_mod
            dash_mod.ensure_table()
        except Exception as exc:                         # noqa: BLE001
            print("    could not prepare the dashboard cache table: {}".format(exc))

        print("  {}".format(manifest.get("dataset_name")))
        for table, rows in tables.items():
            if not rows:
                continue
            cols = _table_columns(table)
            if not cols:
                print("    skipping {} — no such table in this build".format(table))
                continue
            keep = [c for c in cols if c in rows[0]]
            # A project may already exist; inserting it twice would violate the
            # primary key and abort an otherwise good import.
            if table == "project":
                pid = manifest["project_id"]
                if db.row("SELECT project_id FROM project WHERE project_id = ?", [pid]):
                    print("    project already present, keeping it")
                    continue
            db.insert_rows(table, keep,
                           [[r.get(c) for c in keep] for r in rows])
            print("    {:<22} {}".format(table, len(rows)))

        src = tmp / PAYLOAD
        if src.is_dir():
            dest = research_store.dataset_dir(dataset_id)
            dest.parent.mkdir(parents=True, exist_ok=True)
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(src, dest)
            n = sum(1 for _ in dest.rglob("*") if _.is_file())
            total = sum(f.stat().st_size for f in dest.rglob("*") if f.is_file())
            print("    payload               {} file(s), {:.0f} MB".format(
                n, total / 1048576))
        else:
            print("    WARNING: bundle carried no payload")

        print("\n  imported {} — {} analyses already complete".format(
            dataset_id, manifest.get("analyses_complete", 0)))
        db.close()
        return dataset_id


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("export", help="package a dataset and its results")
    e.add_argument("dataset_id")
    e.add_argument("--out", default=None)

    i = sub.add_parser("import-bundle", help="restore a packaged dataset")
    i.add_argument("bundle")
    i.add_argument("--replace", action="store_true",
                   help="overwrite a dataset already registered under this id")

    args = p.parse_args(argv)
    if args.cmd == "export":
        out = Path(args.out or "{}.tar.gz".format(args.dataset_id))
        export_dataset(args.dataset_id, out)
    else:
        import_bundle(Path(args.bundle), replace=args.replace)
    return 0


if __name__ == "__main__":
    sys.exit(main())
