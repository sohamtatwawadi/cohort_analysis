"""Command-line interface for building and inspecting the analytics store.

    python -m backend.cli synthetic [--families N] [--seed S]
    python -m backend.cli ingest <varimat-dir> [--clinical <dir>]
    python -m backend.cli inspect <varimat-file>       # dedup funnel, no write
    python -m backend.cli template <path.csv>          # clinical sidecar header
    python -m backend.cli status
    python -m backend.cli verify                       # reproducibility check
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .app import config, db
from .app.ingest import store, synthetic, varimat
from .app.services import cohort as cohort_svc
from .app.services import carrier, library


def cmd_synthetic(args) -> int:
    result = synthetic.rebuild(n_families=args.families, seed=args.seed)
    library.seed_builtins()
    print(json.dumps({k: v for k, v in result.items() if k != "derived"}, indent=2))
    print("derived:", json.dumps(result["derived"]["zygosity"]))
    print("implied panels:", len(result["derived"]["scope"]["clusters"]))
    return 0


def cmd_ingest(args) -> int:
    try:
        result = store.rebuild_from_varimat(
            directory=Path(args.directory),
            clinical_dir=Path(args.clinical) if args.clinical else None)
    except (FileNotFoundError, ValueError) as e:
        print("error: {}".format(e), file=sys.stderr)
        return 1
    library.seed_builtins()
    print("Ingested {} files -> {} subjects, {} runs, {} findings".format(
        result["files"], result["subjects"], result["runs"], result["findings"]))
    print("\nPer-file dedup funnel (spec §3.2):")
    for f in result["funnels"]:
        print("  {}: {} rows -> {} distinct ({}x) -> {} PASS -> {} coding "
              "-> {} rare -> {} protein-altering | {} genes on target".format(
                  f["source_file"], f["raw_rows"], f["distinct_variants"],
                  f["rows_per_variant"], f["pass_filter"], f["coding"],
                  f["rare"], f["protein_altering"], f["on_target_genes"]))
        if f["unresolved_symbols"]:
            print("    unresolved symbols: {}".format(", ".join(f["unresolved_symbols"][:10])))
    print("\nderived:", json.dumps(result["derived"]["zygosity"]))
    return 0


def cmd_inspect(args) -> int:
    """Parse one file and print the funnel without writing anything."""
    sample = varimat.parse_file(Path(args.file))
    print(json.dumps(sample.funnel.as_dict(), indent=2))
    print("kept (reviewable) variants: {}".format(len(sample.variants)))
    print("filename metadata: pipeline={} caller={} build={}".format(
        sample.pipeline_version, sample.caller, sample.reference_build))
    return 0


def cmd_template(args) -> int:
    p = store and varimat  # keep imports used
    from .app.ingest.clinical import write_sidecar_template
    out = write_sidecar_template(Path(args.path))
    print("wrote clinical sidecar template: {}".format(out))
    return 0


def cmd_status(args) -> int:
    db.connect()
    print("store:  {}".format(config.DB_PATH))
    print("source: {}".format(db.meta_get("source", "none")))
    print("built:  {}".format(db.meta_get("built_at", "—")))
    print("kb:     {}".format(db.meta_get("kb_snapshot_id", "—")))
    print("\ntable counts:")
    for t, n in db.table_counts().items():
        print("  {:18s} {:>8,}".format(t, n))
    return 0


def cmd_verify(args) -> int:
    """Spec E03.6 + E03.9 in one command: every saved cohort must re-resolve to
    an identical member hash."""
    db.connect()
    library.seed_builtins()
    failures = 0
    print("{:44s} {:>7s}  {:>9s}  {}".format("cohort", "n", "members", "criteria"))
    for c in library.list_cohorts():
        r = library.verify_reproducible(c["cohort_id"])
        ok = r.get("ok")
        failures += 0 if ok else 1
        print("{:44s} {:>7,}  {:>9s}  {}  {}".format(
            c["name"][:44], r.get("member_count", 0), r.get("member_hash", "—"),
            r.get("criteria_hash", "—"), "OK" if ok else "DRIFT"))
    print("\n{} cohorts, {} reproducibility failures".format(
        len(library.list_cohorts()), failures))
    return 1 if failures else 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="backend.cli", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("synthetic", help="rebuild the store from the seeded generator")
    s.add_argument("--families", type=int, default=640)
    s.add_argument("--seed", type=int, default=20260813)
    s.set_defaults(fn=cmd_synthetic)

    s = sub.add_parser("ingest", help="rebuild the store from VariMAT files")
    s.add_argument("directory")
    s.add_argument("--clinical", help="clinical/LIMS sidecar directory or file")
    s.set_defaults(fn=cmd_ingest)

    s = sub.add_parser("inspect", help="print one file's dedup funnel, write nothing")
    s.add_argument("file")
    s.set_defaults(fn=cmd_inspect)

    s = sub.add_parser("template", help="write the clinical sidecar CSV header")
    s.add_argument("path")
    s.set_defaults(fn=cmd_template)

    s = sub.add_parser("status", help="store summary")
    s.set_defaults(fn=cmd_status)

    s = sub.add_parser("verify", help="re-resolve every saved cohort, check hashes")
    s.set_defaults(fn=cmd_verify)

    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
