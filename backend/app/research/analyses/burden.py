"""Gene-based rare variant analysis — Part II §4.4.

    "Output: per-gene carrier counts in cases and controls, effect, p, FDR, and
     a gene-ranked plot."

    "Guardrail: genes with too few carriers are reported as such rather than
     given an unstable p-value."

The guardrail is enforced in stats/skat.py (min_carriers) and surfaced here as
a distinct result status. A gene with 1 carrier is not a null result — it is an
untested gene, and collapsing the two is how a rare-disease study concludes a
gene is uninvolved when it simply had no data.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from ..jobs import register
from ..stats import glm, skat as skat_mod
from ..types import MISSING, as_float
from ..variantset import PRESETS, apply_definition, group_by_gene

TESTS = {
    "burden": skat_mod.burden_test,
    "skat": skat_mod.skat_test,
    "skat_o": skat_mod.skat_o_test,
}


def carrier_counts(dosages: np.ndarray, y: np.ndarray) -> Dict[str, int]:
    """Carriers in cases and controls — the numbers a reader checks first."""
    carrier = (dosages > 0).any(axis=0)
    finite = np.isfinite(y)
    return {
        "carriers": int(np.sum(carrier & finite)),
        "carriers_cases": int(np.sum(carrier & finite & (y == 1))),
        "carriers_controls": int(np.sum(carrier & finite & (y == 0))),
        "noncarriers_cases": int(np.sum(~carrier & finite & (y == 1))),
        "noncarriers_controls": int(np.sum(~carrier & finite & (y == 0))),
    }


@register("burden")
def burden_job(spec: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    """Run a gene-based test across every gene with qualifying variants."""
    from .. import store as research_store
    from ..variantset import get_set

    log = context.get("log", lambda m: None)
    dataset_id = context["dataset_id"]

    gm = research_store.load_genotypes(dataset_id)
    ph = research_store.load_phenotypes(dataset_id)
    if ph is None:
        raise ValueError("A gene-based test requires a case/control outcome.")

    outcome = spec["outcome"]
    kind = ph.kind(outcome)
    binary = kind == "binary"

    pheno_index = {s: i for i, s in enumerate(ph.sample_ids)}
    keep = [i for i, s in enumerate(gm.sample_ids) if s in pheno_index]
    if not keep:
        raise ValueError("No samples are shared between genotypes and phenotypes.")
    rows = [pheno_index[gm.sample_ids[i]] for i in keep]
    y = as_float(ph.get(outcome))[rows]

    covariates = []
    for cov in spec.get("covariates", []):
        if cov in ph.columns:
            covariates.append(as_float(ph.get(cov))[rows])
    n_pcs = int(spec.get("n_pcs", 0))
    if n_pcs:
        pcs = research_store.load_covariates(dataset_id, "pcs")
        if pcs is not None:
            for k in range(min(n_pcs, pcs.shape[1])):
                covariates.append(pcs[keep, k])
    X = np.column_stack(covariates) if covariates else None

    # ---- qualifying variants ----------------------------------------------
    if spec.get("variant_set_id"):
        stored = get_set(spec["variant_set_id"])
        if not stored:
            raise KeyError("unknown variant set '{}'".format(spec["variant_set_id"]))
        definition = stored["definition"]
        definition_label = "{} v{}".format(stored["name"], stored["version"])
    elif spec.get("preset"):
        preset = PRESETS.get(spec["preset"])
        if not preset:
            raise KeyError("unknown preset '{}'".format(spec["preset"]))
        definition = preset
        definition_label = preset["label"]
    else:
        definition = spec.get("definition") or PRESETS["ultra_rare_lof"]
        definition_label = definition.get("label", "custom definition")

    # Annotations live with the dataset. Falling back to the spec keeps the
    # scripted path working, but the UI must not have to post them.
    annotations = spec.get("annotations") or research_store.load_annotations(dataset_id)
    selection = apply_definition(gm, definition, annotations)
    mask = selection["mask"]
    log("{} qualifying variants of {}".format(selection["n_qualifying"], gm.n_variants))

    # Two different failures used to share one message, and it named the wrong
    # cause: a definition that filters everything out is not the same as a
    # dataset with no gene annotation, and telling someone to re-annotate when
    # the real problem is a REVEL threshold they cannot meet sends them a long
    # way in the wrong direction.
    if selection["n_qualifying"] == 0:
        culprit = max((f for f in selection["funnel"] if f["dropped"]),
                      key=lambda f: f["dropped"], default=None)
        raise ValueError(
            "No variants qualified under '{}'. {}Relax the definition, or choose a "
            "preset whose annotations this dataset has.".format(
                definition_label,
                "The '{}' filter removed {} of them. ".format(
                    culprit["filter"], culprit["dropped"]) if culprit else ""))

    genes = group_by_gene(gm, mask, annotations)
    if not genes:
        raise ValueError(
            "{} variants qualified, but none carry a gene annotation, so they "
            "cannot be grouped into genes. Re-upload with an annotated VCF, or "
            "annotate the dataset.".format(selection["n_qualifying"]))

    test_name = spec.get("test", "skat_o")
    test_fn = TESTS.get(test_name)
    if test_fn is None:
        raise KeyError("unknown test '{}'".format(test_name))

    min_carriers = int(spec.get("min_carriers", 2))
    maf_all = gm.maf()

    results: List[Dict[str, Any]] = []
    skipped: List[Dict[str, Any]] = []

    for n_done, (gene, idx) in enumerate(sorted(genes.items())):
        G = gm.dosages[np.ix_(idx, keep)]          # variants x samples
        weights = skat_mod.beta_weights(maf_all[idx])
        counts = carrier_counts(G, y)

        # The kernels take samples x variants — the regression layout — while
        # GenotypeMatrix stores variants x samples. Passing the wrong way round
        # is caught by a shape check only when the two happen to differ; for a
        # gene whose variant count equals the sample count it would silently
        # transpose the analysis.
        out = test_fn(G.T, y, X=X, weights=weights, binary=binary,
                      min_carriers=min_carriers)
        row = {
            "gene": gene,
            "n_variants": len(idx),
            "test": test_name,
            **counts,
            **{k: v for k, v in out.items() if k != "status"},
            "status": out.get("status", "ok"),
        }
        if out.get("status") == "insufficient_carriers" or out.get("p") is None:
            row["reason"] = ("Only {} carrier(s) — too few for a stable p-value. "
                             "This gene was not tested; it is not a null result."
                             .format(counts["carriers"]))
            skipped.append(row)
        else:
            results.append(row)
        if n_done and n_done % 200 == 0:
            log("{}/{} genes".format(n_done, len(genes)))

    # §5.6 multiple testing across the genes actually tested — not across all
    # genes, which would over-correct by counting untested ones.
    if results:
        pvals = [r["p"] for r in results]
        fdr = glm.multiple_testing(pvals, "fdr_bh")
        bonf = glm.multiple_testing(pvals, "bonferroni")
        for r, f, b in zip(results, fdr, bonf):
            r["p_fdr_bh"] = float(f)
            r["p_bonferroni"] = float(b)
        results.sort(key=lambda r: r["p"])

    return {
        "analysis": "burden",
        "test": test_name,
        "outcome": outcome,
        "outcome_kind": kind,
        "variant_definition": definition_label,
        "definition": definition,
        "qualifying_funnel": selection["funnel"],
        "n_qualifying_variants": selection["n_qualifying"],
        "n_genes_tested": len(results),
        "n_genes_skipped": len(skipped),
        "min_carriers": min_carriers,
        "covariates": spec.get("covariates", []),
        "results": results[:1000],
        "skipped": skipped[:500],
        "guardrail_note": (
            "{} gene(s) had fewer than {} carriers and were not tested. They are "
            "listed separately — an untested gene is not evidence of no effect."
            .format(len(skipped), min_carriers)),
        "multiple_testing_note": (
            "FDR and Bonferroni computed across the {} genes actually tested."
            .format(len(results))),
        "diagnostics": {
            "n_samples": len(keep),
            "n_genes_with_qualifying_variants": len(genes),
            "p_method": results[0].get("p_method") if results else None,
        },
    }
