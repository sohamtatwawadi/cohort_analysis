"""The R1 descriptive analyses — Part II §4.0.

These were advertised by the capability matrix long before they existed, so the
five cards the matrix marks available under phase R1 did nothing when clicked.
They are descriptive rather than inferential: no model, no p-value, no multiple
testing. What they owe the reader instead is a denominator, because every one of
them is a rate and a rate without its denominator is not a finding.

The rule the lab side already enforces (§G01) applies here too: a frequency
divides by the samples that were actually *called* at that variant, never by the
cohort size. A variant genotyped in 40% of samples and a variant genotyped in
all of them cannot share a denominator.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from ..jobs import register
from ..types import MISSING, GenotypeMatrix, PhenotypeTable

# A frequency computed from a handful of calls is noise presented as a number.
MIN_CALLS = 20

# Small-cell suppression, matching the lab side's governance rule.
MIN_CELL = 5


def _load(context: Dict[str, Any]):
    from .. import store as research_store
    gm = research_store.load_genotypes(context["dataset_id"])
    try:
        ph = research_store.load_phenotypes(context["dataset_id"])
    except Exception:                                    # noqa: BLE001
        ph = None
    return gm, ph


def _called_mask(dosages: np.ndarray) -> np.ndarray:
    """A call is present when it is not the missing sentinel."""
    return dosages != MISSING


def _counts(row: np.ndarray) -> Dict[str, int]:
    """Genotype counts for one variant, over called samples only."""
    called = row[row != MISSING]
    return {
        "n_called": int(called.size),
        "hom_ref": int(np.sum(called == 0)),
        "het": int(np.sum(called == 1)),
        "hom_alt": int(np.sum(called == 2)),
    }


def _af(c: Dict[str, int]) -> Optional[float]:
    """Alternate allele frequency: (het + 2*hom_alt) / (2 * called)."""
    if c["n_called"] < 1:
        return None
    return (c["het"] + 2 * c["hom_alt"]) / (2.0 * c["n_called"])


def _variant_label(v) -> str:
    return v.vid or v.key


def _gene_of(vid: str, annotations: Dict[str, Dict[str, Any]]) -> Optional[str]:
    rec = annotations.get(vid) or {}
    return rec.get("gene") or rec.get("gene_symbol") or None


# ---------------------------------------------------- carrier / allele freq --
@register("carrier_frequency")
def carrier_frequency_job(spec: Dict[str, Any],
                          context: Dict[str, Any]) -> Dict[str, Any]:
    """Per-variant carrier and allele frequency, with the called-sample
    denominator stated on every row."""
    from .. import store as research_store

    log = context.get("log", lambda m: None)
    gm, _ = _load(context)
    annotations = research_store.load_annotations(context["dataset_id"]) or {}

    min_af = float(spec.get("min_af") or 0.0)
    max_af = float(spec.get("max_af") or 1.0)
    model = spec.get("carrier_model") or "dominant"
    top_n = int(spec.get("top_n") or 200)

    log("scanning {} variants across {} samples".format(gm.n_variants, gm.n_samples))

    rows: List[Dict[str, Any]] = []
    skipped_low_call = 0
    for i, v in enumerate(gm.variants):
        c = _counts(gm.dosages[i])
        if c["n_called"] < MIN_CALLS:
            skipped_low_call += 1
            continue
        af = _af(c)
        if af is None or af < min_af or af > max_af:
            continue

        # "Carrier" is a choice, not a fact — a recessive condition's carrier is
        # the heterozygote, a dominant one's is anyone with an ALT allele.
        if model == "recessive":
            carriers = c["hom_alt"]
        elif model == "het_only":
            carriers = c["het"]
        else:
            carriers = c["het"] + c["hom_alt"]

        rows.append({
            "variant": _variant_label(v),
            "chrom": v.chrom, "pos": v.pos, "ref": v.ref, "alt": v.alt,
            "gene": _gene_of(_variant_label(v), annotations),
            "n_called": c["n_called"],
            "carriers": carriers,
            "carrier_rate": carriers / c["n_called"],
            "allele_freq": af,
            "het": c["het"], "hom_alt": c["hom_alt"], "hom_ref": c["hom_ref"],
            "call_rate": c["n_called"] / float(gm.n_samples),
            "suppressed": carriers < MIN_CELL and carriers > 0,
        })

    rows.sort(key=lambda r: r["allele_freq"], reverse=True)
    return {
        "analysis": "carrier_frequency",
        "title": "Carrier / allele frequency",
        "n_samples": gm.n_samples,
        "n_variants_scanned": gm.n_variants,
        "n_variants_reported": len(rows),
        "skipped_low_call_rate": skipped_low_call,
        "carrier_model": model,
        "min_calls_required": MIN_CALLS,
        "rows": rows[:top_n],
        "truncated": len(rows) > top_n,
        "method": (
            "Allele frequency is (het + 2x hom_alt) / (2 x called samples). The "
            "denominator is the samples actually genotyped at that variant, not "
            "the cohort size, so variants with different call rates are not "
            "silently pooled. Variants called in fewer than {} samples are "
            "excluded rather than reported with a wide, unstated interval."
            .format(MIN_CALLS)),
    }


# ------------------------------------------------------------- zygosity ------
@register("zygosity")
def zygosity_job(spec: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    """Zygosity distribution per variant, plus a per-sample burden count."""
    from .. import store as research_store

    gm, ph = _load(context)
    annotations = research_store.load_annotations(context["dataset_id"]) or {}
    top_n = int(spec.get("top_n") or 200)

    called = _called_mask(gm.dosages)
    d = gm.dosages

    # Per-sample: how many variants is this sample het / hom-alt at?
    het_per_sample = np.sum(d == 1, axis=0)
    hom_per_sample = np.sum(d == 2, axis=0)
    called_per_sample = np.sum(called, axis=0)

    samples = []
    for j, sid in enumerate(gm.sample_ids):
        n_called = int(called_per_sample[j])
        het, hom = int(het_per_sample[j]), int(hom_per_sample[j])
        samples.append({
            "sample_id": sid,
            "n_called": n_called,
            "het": het, "hom_alt": hom,
            "het_hom_ratio": (het / hom) if hom else None,
            "call_rate": n_called / float(gm.n_variants) if gm.n_variants else 0.0,
        })

    variants = []
    for i, v in enumerate(gm.variants):
        c = _counts(d[i])
        if c["n_called"] < MIN_CALLS:
            continue
        # A variant where every carrier is homozygous is either recessive,
        # on a haploid chromosome, or a genotyping artefact — all worth seeing.
        carriers = c["het"] + c["hom_alt"]
        variants.append({
            "variant": _variant_label(v),
            "gene": _gene_of(_variant_label(v), annotations),
            "chrom": v.chrom,
            "hom_ref": c["hom_ref"], "het": c["het"], "hom_alt": c["hom_alt"],
            "n_called": c["n_called"],
            "pct_carriers_homozygous": (c["hom_alt"] / carriers) if carriers else None,
        })
    variants.sort(key=lambda r: (r["hom_alt"], r["het"]), reverse=True)

    het_ratios = [s["het_hom_ratio"] for s in samples if s["het_hom_ratio"] is not None]
    return {
        "analysis": "zygosity",
        "title": "Zygosity / inheritance",
        "n_samples": gm.n_samples,
        "n_variants": gm.n_variants,
        "median_het_hom_ratio": float(np.median(het_ratios)) if het_ratios else None,
        "samples": samples[:top_n],
        "variants": variants[:top_n],
        "sex_available": bool(ph and "sex" in (ph.columns or {})),
        "method": (
            "Counts are over called genotypes only. Het/hom ratio is a per-sample "
            "quality signal as much as a biological one: a sample far above the "
            "cohort median is a candidate for contamination, one far below for "
            "a homozygosity-rich or consanguineous background. Hemizygous calls "
            "on X and Y are reported as homozygous by the genotype encoding and "
            "are not separated here."),
    }


# -------------------------------------------------- population frequency -----
@register("population_frequency")
def population_frequency_job(spec: Dict[str, Any],
                             context: Dict[str, Any]) -> Dict[str, Any]:
    """Allele frequency stratified by a grouping column — the internal answer to
    'is this variant actually rare *in this population*?'"""
    from .. import store as research_store

    gm, ph = _load(context)
    annotations = research_store.load_annotations(context["dataset_id"]) or {}
    group_col = spec.get("group") or None
    top_n = int(spec.get("top_n") or 200)

    groups: Dict[str, np.ndarray] = {}
    if group_col and ph is not None and group_col in (ph.columns or {}):
        idx = {s: i for i, s in enumerate(ph.sample_ids)}
        cols = np.array([idx.get(s, -1) for s in gm.sample_ids])
        values = ph.get(group_col)
        labels = np.array([str(values[c]) if c >= 0 else "unknown" for c in cols])
        for lab in sorted(set(labels.tolist())):
            groups[lab] = labels == lab
    else:
        groups["all samples"] = np.ones(gm.n_samples, dtype=bool)

    rows = []
    for i, v in enumerate(gm.variants):
        row = gm.dosages[i]
        overall = _counts(row)
        if overall["n_called"] < MIN_CALLS:
            continue
        per_group = {}
        for lab, mask in groups.items():
            c = _counts(row[mask])
            per_group[lab] = {
                "n_called": c["n_called"],
                "allele_freq": _af(c),
                "carriers": c["het"] + c["hom_alt"],
            }
        freqs = [g["allele_freq"] for g in per_group.values()
                 if g["allele_freq"] is not None and g["n_called"] >= MIN_CALLS]
        spread = (max(freqs) - min(freqs)) if len(freqs) > 1 else None
        rows.append({
            "variant": _variant_label(v),
            "gene": _gene_of(_variant_label(v), annotations),
            "chrom": v.chrom, "pos": v.pos,
            "overall_allele_freq": _af(overall),
            "n_called": overall["n_called"],
            "by_group": per_group,
            "max_group_difference": spread,
        })

    # The interesting rows are the ones that differ between groups.
    rows.sort(key=lambda r: (r["max_group_difference"] is not None,
                             r["max_group_difference"] or 0.0), reverse=True)
    return {
        "analysis": "population_frequency",
        "title": "Population frequency",
        "grouped_by": group_col,
        "groups": {k: int(v.sum()) for k, v in groups.items()},
        "n_variants_reported": len(rows),
        "rows": rows[:top_n],
        "truncated": len(rows) > top_n,
        "method": (
            "Frequencies are internal to this dataset, computed per group over "
            "called samples. They are NOT gnomAD frequencies and must not be "
            "quoted as population reference values — a referral-selected cohort "
            "is enriched for exactly the alleles it was selected on. Groups with "
            "fewer than {} calls at a variant are left out of the spread."
            .format(MIN_CALLS)),
    }


# ------------------------------------------------------- diagnostic yield ----
@register("diagnostic_yield")
def diagnostic_yield_job(spec: Dict[str, Any],
                         context: Dict[str, Any]) -> Dict[str, Any]:
    """Proportion of subjects carrying a qualifying variant, by phenotype group.

    This is the research-mode analogue of the lab's §G03. It is deliberately
    *not* called a diagnostic rate: without curated classifications and a
    phenotype-relevance gate, a qualifying variant is a candidate, not a
    diagnosis. The wording in the result says so.
    """
    from .. import store as research_store

    gm, ph = _load(context)
    if ph is None:
        raise ValueError(
            "Diagnostic yield needs a phenotype file — without a per-subject "
            "indication there is nothing to compute a yield within.")

    annotations = research_store.load_annotations(context["dataset_id"]) or {}
    group_col = spec.get("group") or None
    min_af, max_af = 0.0, float(spec.get("max_af") or 0.05)
    genes_filter = spec.get("genes") or None

    idx = {s: i for i, s in enumerate(ph.sample_ids)}
    cols = np.array([idx.get(s, -1) for s in gm.sample_ids])

    if group_col and group_col in (ph.columns or {}):
        values = ph.get(group_col)
        labels = np.array([str(values[c]) if c >= 0 else "unknown" for c in cols])
    else:
        labels = np.array(["all subjects"] * gm.n_samples)

    # Which variants qualify: rare enough, and in a gene of interest if given.
    qualifying = []
    for i, v in enumerate(gm.variants):
        c = _counts(gm.dosages[i])
        if c["n_called"] < MIN_CALLS:
            continue
        af = _af(c)
        if af is None or af < min_af or af > max_af:
            continue
        gene = _gene_of(_variant_label(v), annotations)
        if genes_filter and (gene not in genes_filter):
            continue
        qualifying.append(i)

    if not qualifying:
        raise ValueError(
            "No variant met the qualifying criteria (allele frequency <= {:.3g}{}). "
            "Raise the frequency ceiling, or check that gene annotations were "
            "supplied with the dataset.".format(
                max_af, " within the selected genes" if genes_filter else ""))

    qual = gm.dosages[np.array(qualifying)]
    has_any = np.any(qual >= 1, axis=0)

    out_groups = []
    for lab in sorted(set(labels.tolist())):
        mask = labels == lab
        n = int(mask.sum())
        k = int(np.sum(has_any & mask))
        out_groups.append({
            "group": lab,
            "subjects": n,
            "with_qualifying_variant": k,
            "yield": (k / n) if n else None,
            "suppressed": 0 < k < MIN_CELL,
        })
    out_groups.sort(key=lambda g: g["yield"] or 0.0, reverse=True)

    total_n = int(len(has_any))
    total_k = int(has_any.sum())
    return {
        "analysis": "diagnostic_yield",
        "title": "Diagnostic yield",
        "grouped_by": group_col,
        "n_qualifying_variants": len(qualifying),
        "max_allele_frequency": max_af,
        "overall": {"subjects": total_n, "with_qualifying_variant": total_k,
                    "yield": (total_k / total_n) if total_n else None},
        "groups": out_groups,
        "method": (
            "A subject counts once if they carry at least one qualifying variant "
            "(allele frequency <= {:.3g}{}). This is a CANDIDATE rate, not a "
            "diagnostic rate: it applies no ACMG classification and no "
            "phenotype-relevance gate, so a variant unrelated to the subject's "
            "indication still counts. Treat it as an upper bound."
            .format(max_af, " in the selected genes" if genes_filter else "")),
    }


# ---------------------------------------------------------- segregation ------
@register("segregation")
def segregation_job(spec: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    """Mendelian-consistency check across declared trios.

    A de-novo call from genotypes alone is unreliable — most apparent de novos
    are genotyping error in one of the three samples — so this reports
    *inconsistencies* and says what they usually mean rather than labelling them
    de novo mutations.
    """
    from .. import store as research_store

    gm, ph = _load(context)
    if ph is None:
        raise ValueError("Segregation needs a phenotype/pedigree file with "
                         "father_id and mother_id columns.")

    cols = ph.columns or {}
    fa_col = next((c for c in ("father_id", "paternal_id", "father") if c in cols), None)
    mo_col = next((c for c in ("mother_id", "maternal_id", "mother") if c in cols), None)
    if not fa_col or not mo_col:
        raise ValueError(
            "No pedigree columns found. Expected father_id and mother_id "
            "(or paternal_id / maternal_id) in the phenotype file; found: {}."
            .format(", ".join(sorted(cols)) or "no columns"))

    annotations = research_store.load_annotations(context["dataset_id"]) or {}
    sidx = gm.sample_index()
    fathers, mothers = ph.get(fa_col), ph.get(mo_col)

    trios = []
    for i, sid in enumerate(ph.sample_ids):
        fa, mo = str(fathers[i]), str(mothers[i])
        if sid in sidx and fa in sidx and mo in sidx:
            trios.append((sid, fa, mo))
    if not trios:
        raise ValueError(
            "No complete trios: a trio needs the child, father and mother all "
            "present in the genotype data. Found {} pedigree rows."
            .format(len(ph.sample_ids)))

    d = gm.dosages
    findings, per_trio = [], []
    for child, fa, mo in trios:
        ci, fi, mi = sidx[child], sidx[fa], sidx[mo]
        cg, fg, mg = d[:, ci], d[:, fi], d[:, mi]
        complete = (cg != MISSING) & (fg != MISSING) & (mg != MISSING)

        # A child carrying ALT when neither parent does is Mendelian-inconsistent.
        inconsistent = complete & (cg >= 1) & (fg == 0) & (mg == 0)
        n_bad = int(inconsistent.sum())
        per_trio.append({
            "child": child, "father": fa, "mother": mo,
            "variants_compared": int(complete.sum()),
            "inconsistent": n_bad,
            "rate": (n_bad / int(complete.sum())) if complete.sum() else None,
        })
        for i in np.flatnonzero(inconsistent)[:50]:
            v = gm.variants[int(i)]
            findings.append({
                "child": child, "variant": _variant_label(v),
                "gene": _gene_of(_variant_label(v), annotations),
                "chrom": v.chrom, "pos": v.pos,
                "child_dosage": int(cg[i]), "father_dosage": int(fg[i]),
                "mother_dosage": int(mg[i]),
            })

    rates = [t["rate"] for t in per_trio if t["rate"] is not None]
    return {
        "analysis": "segregation",
        "title": "Segregation / de novo",
        "n_trios": len(trios),
        "median_inconsistency_rate": float(np.median(rates)) if rates else None,
        "trios": per_trio,
        "findings": findings[:200],
        "n_findings": len(findings),
        "method": (
            "A variant is flagged when the child carries an ALT allele and both "
            "parents are called homozygous reference. This is NOT a de-novo "
            "call: at typical genotyping error rates most such sites are errors "
            "in one of the three samples, and a real de novo needs orthogonal "
            "confirmation. A trio whose rate is far above the others is usually "
            "a sample swap or mislabelled pedigree, not a mutator phenotype."),
    }
