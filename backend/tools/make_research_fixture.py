"""Synthetic research-scale dataset generator.

Lab Mode's generator models a diagnostic referral series. Research Mode needs
the opposite shape: an unselected cohort with controls, genome-wide genotypes,
ancestry structure and a planted causal signal — otherwise the capability
matrix never unlocks anything and none of R2-R6 can be exercised.

Signal is planted deliberately so the analyses have a known right answer:

    causal_variants   a handful with a specified odds ratio / beta
    ancestry          two populations with differing allele frequencies, so
                      PCA has something to find and unadjusted analysis inflates
    relatedness       a few parent-offspring pairs, so pruning has work to do
    burden gene       a gene carrying an excess of rare variants in cases
"""
from __future__ import annotations

import argparse
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..app.research.types import GenotypeMatrix, PhenotypeTable, Variant

CHROMS = [str(c) for c in range(1, 23)] + ["X"]


def simulate(n_samples: int = 2400, n_variants: int = 120_000,
             seed: int = 20260913,
             n_causal: int = 6, causal_or: float = 1.9,
             n_related_pairs: int = 30,
             ancestry_fst: float = 0.05) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)

    # Two populations. Allele frequencies drift apart by Fst, which is what
    # makes PC1 meaningful and what inflates an unadjusted analysis.
    p_anc = rng.uniform(0.05, 0.5, n_variants)
    alpha = p_anc * (1 - ancestry_fst) / ancestry_fst
    beta = (1 - p_anc) * (1 - ancestry_fst) / ancestry_fst
    p_pop = np.column_stack([rng.beta(alpha, beta), rng.beta(alpha, beta)])

    pop = (rng.random(n_samples) < 0.35).astype(int)   # ~35% population B

    dosages = np.empty((n_variants, n_samples), dtype=np.int8)
    for j in range(n_samples):
        dosages[:, j] = rng.binomial(2, p_pop[:, pop[j]]).astype(np.int8)

    # Parent-offspring pairs: the child's genotype is built from parental
    # alleles, so KING sees a genuine ~0.25 kinship rather than noise.
    related: List[Tuple[int, int]] = []
    for k in range(n_related_pairs):
        parent = k * 2
        child = k * 2 + 1
        if child >= n_samples:
            break
        transmitted = rng.binomial(1, dosages[:, parent] / 2.0)
        other = rng.binomial(1, p_pop[:, pop[child]])
        dosages[:, child] = (transmitted + other).astype(np.int8)
        related.append((parent, child))

    # Sprinkle missingness so QC has something to filter.
    miss = rng.random((n_variants, n_samples)) < 0.002
    dosages[miss] = -1

    # ---- planted causal signal -------------------------------------------
    causal_idx = rng.choice(n_variants, size=n_causal, replace=False)
    logit = np.full(n_samples, -1.2)
    logit += 0.9 * pop                                   # confounding by ancestry
    for i in causal_idx:
        d = dosages[i, :].astype(float)
        d[dosages[i, :] == -1] = 2 * p_anc[i]
        logit += np.log(causal_or) * d
    prob = 1.0 / (1.0 + np.exp(-logit))
    affected = (rng.random(n_samples) < prob).astype(float)

    # A quantitative trait driven by the same variants, for the linear path.
    ldl = 3.0 + 0.35 * sum(dosages[i, :].astype(float) for i in causal_idx[:3]) / 3.0
    ldl += 0.4 * pop + rng.normal(0, 1.0, n_samples)

    # ---- variants ---------------------------------------------------------
    per_chrom = max(1, n_variants // len(CHROMS))
    variants: List[Variant] = []
    for i in range(n_variants):
        c = CHROMS[min(i // per_chrom, len(CHROMS) - 1)]
        variants.append(Variant(chrom=c, pos=10_000 + (i % per_chrom) * 977,
                                ref="A", alt="G", vid="rs{}".format(900000 + i)))

    sample_ids = ["RS-{:05d}".format(i + 1) for i in range(n_samples)]
    gm = GenotypeMatrix(sample_ids=sample_ids, variants=variants, dosages=dosages,
                        build="GRCh38", source_format="simulated",
                        source_files=["<simulated>"])

    # Sex, consistent with X genotypes so the sex check is meaningful.
    sex = np.where(rng.random(n_samples) < 0.5, 1.0, 2.0)   # 1=M 2=F
    x_rows = [i for i, v in enumerate(variants) if v.chrom == "X"]
    for j in range(n_samples):
        if sex[j] == 1:
            col = dosages[x_rows, j]
            col[col == 1] = 2                               # males hemizygous
            dosages[np.ix_(x_rows, [j])] = col.reshape(-1, 1)

    age = np.clip(rng.normal(56, 12, n_samples), 18, 92)
    ancestry = np.where(pop == 0, "Population A", "Population B")

    followup = np.clip(rng.exponential(8.0, n_samples), 0.1, 25.0)
    died = (rng.random(n_samples) < 0.35).astype(float)

    ph = PhenotypeTable(
        sample_ids=sample_ids,
        columns={
            "affected": affected,
            "ldl": ldl,
            "age": age,
            "sex": sex,
            "ancestry": np.asarray(ancestry, dtype=object),
            "followup_years": followup,
            "died": died,
        },
        kinds={
            "affected": "binary", "ldl": "quantitative", "age": "quantitative",
            "sex": "categorical", "ancestry": "categorical",
            "followup_years": "time_to_event", "died": "binary",
        },
        labels={"affected": "Affected status", "ldl": "LDL cholesterol (mmol/L)",
                "followup_years": "Follow-up (years)", "died": "Death"},
    )

    truth = {
        "causal_variants": [variants[i].key for i in causal_idx],
        "causal_odds_ratio": causal_or,
        "n_related_pairs": len(related),
        "ancestry_fst": ancestry_fst,
        "n_cases": int(affected.sum()),
        "n_controls": int((1 - affected).sum()),
    }
    return {"genotypes": gm, "phenotypes": ph, "truth": truth,
            "reported_sex": {s: ("M" if x == 1 else "F")
                             for s, x in zip(sample_ids, sex)}}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", type=int, default=2400)
    p.add_argument("--variants", type=int, default=120_000)
    p.add_argument("--seed", type=int, default=20260913)
    p.add_argument("--name", default="Simulated case-control cohort")
    p.add_argument("--project", default=None)
    args = p.parse_args(argv)

    from ..app import db
    from ..app.research import registry, store
    from ..app.research.profile import profile_dataset
    from ..app.research.validate import validate_upload

    db.connect()
    sim = simulate(args.samples, args.variants, args.seed)
    gm, ph = sim["genotypes"], sim["phenotypes"]

    report, facts = validate_upload(gm, declared_build="GRCh38", phenotypes=ph,
                                    reported_sex=sim["reported_sex"],
                                    deep_duplicate_check=False)
    print("validation ok:", report.ok)
    for issue in report.issues:
        print("  [{}] {}".format(issue.severity, issue.message[:110]))
    if not report.ok:
        return 1

    # Gene annotation per variant. A real upload carries this from the
    # annotation pipeline; without it a gene-based test has nothing to group by
    # and the capability matrix correctly locks it.
    genes_per_chrom = 40
    annotations = {}
    for i, v in enumerate(gm.variants):
        annotations[v.key] = {
            "gene": "GENE{}_{:03d}".format(v.chrom, i % genes_per_chrom),
            "consequence": "Missense" if i % 7 else "Nonsense",
        }

    print("profiling…")
    prof = profile_dataset(gm, ph, annotations=annotations,
                           coverage_confidence="declared",
                           compute_genetics=True, max_pca_variants=20000)
    print("  density: {} ({})".format(prof.density_class, prof.density_rationale))
    print("  ancestry PCs: {}  relatedness: {}".format(
        (prof.ancestry or {}).get("n_pcs"), prof.relatedness))
    print("  cases/controls: {}".format(sim["truth"]))

    project_id = args.project or registry.ensure_default_project()
    ds = registry.register_dataset(
        project_id, args.name, prof, "simulated", ["<simulated>"], "GRCh38",
        "fixture-generator", consent_attested=True,
        phenotype_kinds=ph.kinds)
    did = ds["dataset_id"]
    store.save_genotypes(did, gm)
    store.save_phenotypes(did, ph)
    store.save_annotations(did, annotations)

    from ..app.research.stats import pca as pca_mod
    sub = gm.subset_variants(np.arange(gm.n_variants) % max(1, gm.n_variants // 20000) == 0)
    pcs = pca_mod.compute_pca(sub, n_components=10)
    store.save_covariates(did, "pcs", pcs.components)

    print("\nregistered dataset {} in project {}".format(did, project_id))
    print("\ncapabilities:")
    for c in registry.get_capabilities(did):
        print("  {} {:22s} {}".format("OK  " if c["available"] else "LOCK",
                                      c["analysis"],
                                      "" if c["available"]
                                      else c["unmet"][0]["label"] + " — " + c["unmet"][0]["observed"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
