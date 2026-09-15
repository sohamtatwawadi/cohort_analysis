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
             n_related_pairs: int = 30, n_trios: int = 40,
             ancestry_fst: float = 0.05,
             rare_fraction: float = 0.45, variants_per_gene: int = 60,
             n_burden_genes: int = 3, burden_log_or: float = 1.3,
             burden_gene_span: int = 6, target_case_fraction: float = 0.5
             ) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)

    # Allele-frequency spectrum. A cohort of only common variants cannot
    # exercise a gene-based test at all: every rare-variant preset filters on
    # MAF <= 0.001 and removes the entire dataset. Real sequence data is mostly
    # rare, so most variants here are too.
    p_anc = rng.uniform(0.05, 0.5, n_variants)
    rare_mask = rng.random(n_variants) < rare_fraction
    n_rare = int(rare_mask.sum())
    # Aim for 2-8 copies in the cohort: rare enough for every ultra-rare preset
    # (MAF <= 0.001) but actually *observed*. Drawing lower than this produces
    # monomorphic sites — present in the file, absent from every sample — which
    # pass no filter and make a gene-based test look broken.
    rare_p = rng.uniform(0.0004, 0.0012, n_rare)

    # Genes, assigned in contiguous blocks so a gene's variants sit together the
    # way they do on a chromosome. Burden collapses within these.
    gene_of: List[str] = []
    for i in range(n_variants):
        gene_of.append("GENE{:04d}".format(i // variants_per_gene))
    alpha = p_anc * (1 - ancestry_fst) / ancestry_fst
    beta = (1 - p_anc) * (1 - ancestry_fst) / ancestry_fst
    p_pop = np.column_stack([rng.beta(alpha, beta), rng.beta(alpha, beta)])

    # Ultra-rare variants bypass the ancestry drift model. A Beta with these
    # shape parameters is so U-shaped that almost every draw collapses to zero,
    # and population structure in a variant seen four times is not meaningful
    # anyway — so both populations take the intended frequency directly.
    p_pop[rare_mask, 0] = rare_p
    p_pop[rare_mask, 1] = rare_p
    # Keep p_anc in step — it is the imputation value for missing dosages, and
    # imputing a 20%-frequency allele into a variant seen four times would
    # manufacture carriers that do not exist.
    p_anc[rare_mask] = rare_p

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

    # Complete trios, taken from the far end of the cohort so they do not
    # overlap the parent-offspring pairs above. Segregation needs all three
    # members genotyped: a child with only one parent present is not a trio and
    # Mendelian consistency cannot be evaluated for it.
    father_id = np.array(["0"] * n_samples, dtype=object)
    mother_id = np.array(["0"] * n_samples, dtype=object)
    trios: List[Tuple[int, int, int]] = []
    cursor = n_samples - 1
    for _ in range(n_trios):
        if cursor - 2 < 2 * n_related_pairs:
            break
        child, fa, mo = cursor, cursor - 1, cursor - 2
        cursor -= 3
        # One allele transmitted from each parent — the definition of a trio,
        # and what makes an inconsistency meaningful rather than arithmetic.
        from_fa = rng.binomial(1, dosages[:, fa] / 2.0)
        from_mo = rng.binomial(1, dosages[:, mo] / 2.0)
        dosages[:, child] = (from_fa + from_mo).astype(np.int8)
        trios.append((child, fa, mo))

    # Sprinkle missingness so QC has something to filter.
    miss = rng.random((n_variants, n_samples)) < 0.002
    dosages[miss] = -1

    # ---- planted causal signal -------------------------------------------
    # Single-variant causal effects go on COMMON variants. A variant seen four
    # times carries no power at this sample size whatever its true odds ratio,
    # so planting a signal there produces a scan that correctly finds nothing —
    # indistinguishable from a broken scan. Rare-variant signal is planted in
    # the burden genes below, which is the test that can actually detect it.
    common_idx = np.flatnonzero(~rare_mask)
    causal_idx = rng.choice(common_idx, size=min(n_causal, common_idx.size),
                            replace=False)
    logit = np.full(n_samples, -1.2)
    logit += 0.9 * pop                                   # confounding by ancestry
    for i in causal_idx:
        d = dosages[i, :].astype(float)
        d[dosages[i, :] == -1] = 2 * p_anc[i]
        logit += np.log(causal_or) * d

    # A burden signal: a few genes where rare variants collectively raise risk.
    # Individually each is far too rare to detect — which is the entire reason
    # gene-based tests exist — so the signal only appears once collapsed.
    all_genes = sorted(set(gene_of))
    burden_genes = list(rng.choice(all_genes, size=min(n_burden_genes,
                                                       len(all_genes)),
                                   replace=False))

    # Make the burden genes LARGE, the way the genes burden tests actually find
    # are large. Every qualifying variant is capped at MAF 0.001 — six alleles
    # in this cohort — so carriers can only accumulate through variant count.
    # At the default gene size each burden gene reached ~30 carriers, which is a
    # true effect with no power against 2,500 tests: the planted genes ranked
    # 9th and 20th and nothing survived FDR. Absorbing neighbouring blocks is
    # how a real TTN or BRCA2 gets its carrier count.
    gene_start = {g: i * variants_per_gene for i, g in enumerate(all_genes)}
    for bg in burden_genes:
        start = gene_start[bg]
        for i in range(start, min(start + variants_per_gene * burden_gene_span,
                                  n_variants)):
            gene_of[i] = bg

    burden_rows = [i for i, g in enumerate(gene_of)
                   if g in set(burden_genes) and rare_mask[i]]

    # The effect applies ONCE PER GENE, not once per variant — that is what a
    # burden model says, and it is what the test collapses to. Adding it per
    # variant instead made a subject carrying five qualifying variants five
    # times as affected, which drove the cohort to 80% cases and destroyed the
    # case/control balance every other analysis depends on.
    burden_by_gene: Dict[str, List[int]] = {}
    for i in burden_rows:
        burden_by_gene.setdefault(gene_of[i], []).append(i)
    for gene, rows_in_gene in burden_by_gene.items():
        block = dosages[np.array(rows_in_gene), :]
        carries = np.any(block >= 1, axis=0)
        logit += burden_log_or * carries

    # Recentre the intercept on the target prevalence. Every planted effect
    # shifts the mean of the linear predictor, so a fixed intercept makes the
    # case fraction a side effect of how much signal was planted — moving the
    # causal variants onto common sites alone took the cohort to 80% cases,
    # which starves every control-based analysis. Recentring preserves all the
    # relative effects (the analyses still recover them) and puts the baseline
    # risk back under control, which is what an intercept is for.
    target_logit = np.log(target_case_fraction / (1 - target_case_fraction))
    logit += target_logit - np.median(logit)

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
    for child, fa, mo in trios:
        father_id[child] = sample_ids[fa]
        mother_id[child] = sample_ids[mo]
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
            "father_id": father_id,
            "mother_id": mother_id,
        },
        kinds={
            "affected": "binary", "ldl": "quantitative", "age": "quantitative",
            "sex": "categorical", "ancestry": "categorical",
            "followup_years": "time_to_event", "died": "binary",
            "father_id": "categorical", "mother_id": "categorical",
        },
        labels={"affected": "Affected status", "ldl": "LDL cholesterol (mmol/L)",
                "followup_years": "Follow-up (years)", "died": "Death",
                "father_id": "Father (pedigree)", "mother_id": "Mother (pedigree)"},
    )

    truth = {
        "causal_variants": [variants[i].key for i in causal_idx],
        "causal_odds_ratio": causal_or,
        "n_related_pairs": len(related),
        "ancestry_fst": ancestry_fst,
        "n_cases": int(affected.sum()),
        "n_controls": int((1 - affected).sum()),
        "n_trios": len(trios),
        "n_rare_variants": n_rare,
        "burden_genes": burden_genes,
        "burden_variants": len(burden_rows),
        "burden_log_or": burden_log_or,
    }
    pedigree = [{"iid": sample_ids[c], "pat": sample_ids[f], "mat": sample_ids[m]}
                for c, f, m in trios]
    return {"genotypes": gm, "phenotypes": ph, "truth": truth,
            "pedigree": pedigree, "gene_of": gene_of,
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
    gene_of = sim["gene_of"]
    annotations = {}
    for i, v in enumerate(gm.variants):
        # A rare-variant preset filters on consequence as well as frequency, so
        # a share of these must be loss-of-function or the gene-based tests
        # qualify nothing.
        # These strings must match variantset.LOF_CONSEQUENCES exactly — the
        # preset filters on the term, so "Splice" instead of "Splice site"
        # silently qualifies nothing.
        consequence = ("Nonsense" if i % 5 == 0 else
                       "Frameshift" if i % 5 == 1 else
                       "Splice site" if i % 5 == 2 else "Missense")
        is_lof = consequence in ("Nonsense", "Frameshift", "Splice site")
        annotations[v.key] = {
            "gene": gene_of[i],
            "consequence": consequence,
            "lof": is_lof,
            # ultra_rare_lof additionally requires high-confidence LoF; a real
            # annotation pipeline emits this from LOFTEE.
            "lof_confidence": "HC" if is_lof else None,
            # rare_damaging_missense filters on REVEL.
            "REVEL": round(0.2 + 0.75 * ((i * 37) % 100) / 100.0, 3)
                     if consequence == "Missense" else None,
        }

    print("profiling…")
    prof = profile_dataset(gm, ph, annotations=annotations,
                           pedigree=sim.get("pedigree"),
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
