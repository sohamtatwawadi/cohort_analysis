"""Register a 1000 Genomes chromosome as a research dataset.

Real public data, from the AWS Registry of Open Data — no credentials, no
account. It earns its place next to the synthetic cohort because it brings the
one thing simulation cannot: genuine population structure. 2,504 samples across
26 populations and 5 super-populations, so the composition panels, population
frequency and PCA show real continental clusters rather than two invented blobs.

    python -m backend.tools.fetch_1000g --chrom 22 --target-variants 150000

What it deliberately does NOT do
--------------------------------
1000 Genomes ships no phenotypes. It has population, super-population and sex —
all real, all used here — and no disease status, no age, no follow-up. So
association, GWAS, burden, survival and PRS stay LOCKED on this dataset, and
that is the correct outcome rather than a gap to paper over: the capability
matrix is telling the truth about what the data can answer.

Nothing here invents a case/control column. These are real, identifiable
consented samples, and attaching a fabricated disease status to them to make a
demo look fuller would be inventing clinical data about real people.

Sizing
------
The full chr22 is 1.1M variants x 2,504 samples = 2.6 GB as int8, which does not
fit a small instance. Thinning keeps one record in N across the whole
chromosome — not the first N, which would be the short arm and nothing else.
At 150,000 variants the matrix is 358 MB.
"""
from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

BUCKET = "https://1000genomes.s3.amazonaws.com"
RELEASE = "release/20130502"
PANEL = "{}/{}/integrated_call_samples_v3.20130502.ALL.panel".format(BUCKET, RELEASE)

# Phase 3 is aligned to GRCh37. Declaring this correctly is not a detail: the
# app refuses to guess a build because GRCh37 and GRCh38 coordinates overlap,
# so a wrong declaration silently corrupts every annotation downstream.
BUILD = "GRCh37"

# chr22 is v5a; the sex chromosomes use a different suffix and are excluded —
# they are haploid in half the cohort and need handling this loader does not do.
VCF_TEMPLATE = ("{bucket}/{release}/ALL.chr{chrom}"
                ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

SUPER_POP_NAMES = {
    "AFR": "African", "AMR": "Admixed American", "EAS": "East Asian",
    "EUR": "European", "SAS": "South Asian",
}


def _ssl_context():
    """A context with a usable CA bundle.

    Python from python.org on macOS ships without wiring up the system trust
    store, so urllib fails with CERTIFICATE_VERIFY_FAILED on a host curl
    handles fine. certifi is present as a transitive dependency and carries the
    bundle; falling back to the default context keeps Linux, where the system
    certificates are already correct, on its normal path.
    """
    import ssl
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:                                    # noqa: BLE001
        return ssl.create_default_context()


def _download(url: str, dest: Path, label: str) -> Path:
    """Fetch to a local cache, skipping if already present and non-empty."""
    if dest.exists() and dest.stat().st_size > 0:
        print("  {} already cached ({:.0f} MB)".format(
            label, dest.stat().st_size / 1048576))
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print("  downloading {} …".format(label))

    # Streamed in chunks rather than urlretrieve, so the CA bundle above can be
    # applied and a 205 MB body never lands in memory.
    ctx = _ssl_context()
    try:
        with urllib.request.urlopen(url, context=ctx, timeout=60) as resp:
            total = int(resp.headers.get("Content-Length") or 0)
            done = 0
            with tmp.open("wb") as out:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
                    done += len(chunk)
                    if total and done % (1 << 24) < (1 << 20):
                        sys.stdout.write("\r    {:.0f} / {:.0f} MB".format(
                            done / 1048576, total / 1048576))
                        sys.stdout.flush()
    except Exception as exc:                             # noqa: BLE001
        tmp.unlink(missing_ok=True)
        raise SystemExit(
            "Could not fetch {}\n  {}: {}\n\nIf this is a certificate error on "
            "macOS, run the 'Install Certificates.command' that ships with your "
            "Python, or pip install certifi.".format(url, type(exc).__name__, exc))
    sys.stdout.write("\r" + " " * 40 + "\r")
    # Rename only once complete, so an interrupted download is never mistaken
    # for a cached file on the next run.
    tmp.rename(dest)
    print("  {}: {:.0f} MB".format(label, dest.stat().st_size / 1048576))
    return dest


def load_panel(path: Path) -> Dict[str, Dict[str, str]]:
    """sample -> {pop, super_pop, gender}. Tab-separated, one header line."""
    out: Dict[str, Dict[str, str]] = {}
    with path.open() as fh:
        header = fh.readline().split()
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.split()
            if len(f) < 4:
                continue
            out[f[idx["sample"]]] = {
                "pop": f[idx["pop"]],
                "super_pop": f[idx["super_pop"]],
                "gender": f[idx["gender"]],
            }
    return out


def count_records(path: Path) -> int:
    """Data lines in the VCF, without parsing any of them."""
    import gzip
    n = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                n += 1
    return n


def build_phenotypes(sample_ids: List[str],
                     panel: Dict[str, Dict[str, str]]):
    """A PhenotypeTable from the panel's real columns, and nothing invented."""
    from ..app.research.types import PhenotypeTable

    pop, spop, sex = [], [], []
    for s in sample_ids:
        rec = panel.get(s, {})
        pop.append(rec.get("pop") or "unknown")
        code = rec.get("super_pop") or "unknown"
        spop.append(SUPER_POP_NAMES.get(code, code))
        g = (rec.get("gender") or "").lower()
        # PLINK's coding, which is what the rest of the app expects.
        sex.append(1.0 if g == "male" else 2.0 if g == "female" else 0.0)

    return PhenotypeTable(
        sample_ids=list(sample_ids),
        columns={
            "population": np.asarray(pop, dtype=object),
            "super_population": np.asarray(spop, dtype=object),
            "sex": np.asarray(sex, dtype=float),
        },
        kinds={
            "population": "categorical",
            "super_population": "categorical",
            "sex": "categorical",
        },
        labels={
            "population": "Population (26 groups)",
            "super_population": "Super-population",
            "sex": "Reported sex",
        },
    )


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chrom", default="22",
                   help="autosome to load; 22 is the smallest and the default")
    p.add_argument("--target-variants", type=int, default=150_000,
                   help="thin to about this many, spread across the chromosome")
    p.add_argument("--cache", default="data/1000g",
                   help="where the downloaded files are kept")
    p.add_argument("--name", default=None)
    p.add_argument("--project", default=None)
    p.add_argument("--max-pca-variants", type=int, default=6000)
    args = p.parse_args(argv)

    if args.chrom.upper() in ("X", "Y", "MT", "M"):
        print("Only autosomes are supported: the sex chromosomes are haploid in "
              "half the cohort and this loader does not handle that.")
        return 2

    from ..app import db
    from ..app.research import registry, store
    from ..app.research.formats import vcf as vcf_reader
    from ..app.research.profile import profile_dataset
    from ..app.research.validate import validate_upload

    cache = Path(args.cache)
    print("1000 Genomes phase 3, chromosome {} ({})".format(args.chrom, BUILD))

    panel_path = _download(PANEL, cache / "samples.panel", "sample panel")
    panel = load_panel(panel_path)
    print("  panel: {} samples, {} populations, {} super-populations".format(
        len(panel),
        len({v["pop"] for v in panel.values()}),
        len({v["super_pop"] for v in panel.values()})))

    url = VCF_TEMPLATE.format(bucket=BUCKET, release=RELEASE, chrom=args.chrom)
    vcf_path = _download(url, cache / "chr{}.vcf.gz".format(args.chrom),
                         "chr{} genotypes".format(args.chrom))

    print("  counting records (no parsing) …")
    total = count_records(vcf_path)
    step = max(1, total // max(1, args.target_variants))
    print("  {} records; keeping 1 in {} → about {}".format(
        format(total, ","), step, format(total // step, ",")))

    print("  reading …")
    gm = vcf_reader.read_vcf(vcf_path, every_nth=step)
    print("  parsed {} variants x {} samples ({:.0f} MB)".format(
        format(gm.n_variants, ","), gm.n_samples, gm.dosages.nbytes / 1048576))

    missing = [s for s in gm.sample_ids if s not in panel]
    if missing:
        print("  {} genotyped samples are absent from the panel; their "
              "metadata reads 'unknown'".format(len(missing)))

    ph = build_phenotypes(gm.sample_ids, panel)
    reported_sex = {s: ("M" if (panel.get(s, {}).get("gender") == "male")
                        else "F" if (panel.get(s, {}).get("gender") == "female")
                        else "U")
                    for s in gm.sample_ids}

    db.connect()
    report, _facts = validate_upload(gm, declared_build=BUILD, phenotypes=ph,
                                     reported_sex=reported_sex,
                                     deep_duplicate_check=False)
    print("  validation ok:", report.ok)
    for issue in report.issues:
        print("    [{}] {}".format(issue.severity, issue.message[:120]))
    if not report.ok:
        return 1

    print("  profiling …")
    prof = profile_dataset(gm, ph, annotations=None,
                           coverage_confidence="declared",
                           compute_genetics=True,
                           max_pca_variants=args.max_pca_variants)
    print("    density: {} ({})".format(prof.density_class, prof.density_rationale))
    print("    ancestry PCs: {}".format((prof.ancestry or {}).get("n_pcs")))
    print("    unrelated: {} of {}".format(prof.n_unrelated, prof.n_samples))

    name = args.name or "1000 Genomes phase 3 · chr{} · {} samples".format(
        args.chrom, format(gm.n_samples, ","))
    project_id = args.project or registry.ensure_default_project()
    ds = registry.register_dataset(
        project_id, name, prof, "VCF", [url], BUILD,
        "1000genomes-loader", consent_attested=True,
        phenotype_kinds=ph.kinds)
    did = ds["dataset_id"]
    store.save_genotypes(did, gm)
    store.save_phenotypes(did, ph)

    print("\nregistered {} in project {}".format(did, project_id))

    from ..app.research import capability
    print("\ncapabilities:")
    for cap in capability.assess(prof):
        mark = "OK  " if cap.available else "LOCK"
        why = "" if cap.available else next(
            ("{} — {}".format(r.label, r.observed)
             for r in cap.requirements if not r.met), "")
        print("  {} {:<22} {}".format(mark, cap.analysis, why))

    print("\nThe locked analyses are correct, not a failure: 1000 Genomes has no "
          "disease phenotype, so nothing case/control can run on it. Population "
          "structure, carrier and allele frequency, and zygosity all can.")
    db.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
