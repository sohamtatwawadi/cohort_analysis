"""Level-1 dataset dashboard — what a registered dataset actually contains.

The gap this fills: registering a dataset took you straight to the capability
matrix, which answers "what can this support?" without ever answering "what is
this?". A reader needs the second question first.

It is modelled on an oncology cohort dashboard — KPI row, top altered genes,
composition breakdowns, distributions — with one substitution that is not
cosmetic. A somatic dashboard leads with tumour mutational burden, stage and
linked therapy. None of those exist for germline data, and inventing them would
mean a dashboard that reads as clinically meaningful while being fabricated.
So the same *slots* are filled from what a germline research dataset really
holds: carriers per gene, ancestry and case/control composition, the allele
frequency spectrum, relatedness.

The one rule carried over verbatim is the important one, and their mockup states
it too: a per-gene rate divides by the samples actually **called** at that
gene's variants, never by the cohort size. A sample with no call at a gene is
not evidence of absence.

Computing top genes needs the genotype matrix, which is ~450 MB for a
genome-wide dataset, so the result is cached. Recomputing it on every page view
would make the dashboard the most expensive screen in the product.
"""
from __future__ import annotations

import json
from datetime import datetime
from typing import Any, Dict, List, Optional

import numpy as np

from .. import db
from .types import MISSING

# Rank this many genes. Enough to fill the panel; not so many that the payload
# becomes a second copy of the annotation table.
TOP_GENES = 14

# Below this, a per-gene rate is being computed from too few called samples to
# mean anything, so the gene is ranked but its rate is withheld.
MIN_CALLED = 20

# Only variants at or below this frequency count toward a gene's carrier rate.
#
# Ranking on ALL variants put every gene at 100%: with common variants included,
# every sample carries something in every gene, and the panel says nothing. An
# oncology dashboard has the same constraint and solves it the same way — it
# counts reportable alterations, not every SNP. Rare is the germline analogue of
# "worth reporting", and the threshold is stated on the panel so the reader
# knows what they are looking at.
GENE_MAX_AF = 0.01

# Small-cell suppression, matching the governance rule on the lab side.
MIN_CELL = 5


# --------------------------------------------------------------- persistence --
def _ensure_table() -> None:
    db.execute("""
        CREATE TABLE IF NOT EXISTS dataset_dashboard (
            dataset_id     VARCHAR PRIMARY KEY,
            computed_at    TIMESTAMP NOT NULL,
            dashboard_json VARCHAR NOT NULL
        )
    """)


def get(dataset_id: str, refresh: bool = False) -> Dict[str, Any]:
    """Cached dashboard, computed on first request."""
    _ensure_table()
    if not refresh:
        row = db.row("SELECT dashboard_json FROM dataset_dashboard WHERE dataset_id = ?",
                     [dataset_id])
        if row:
            return json.loads(row["dashboard_json"])

    payload = compute(dataset_id)
    db.execute("DELETE FROM dataset_dashboard WHERE dataset_id = ?", [dataset_id])
    db.insert_rows("dataset_dashboard",
                   ["dataset_id", "computed_at", "dashboard_json"],
                   [[dataset_id, datetime.utcnow(), json.dumps(payload, default=str)]])
    return payload


def invalidate(dataset_id: str) -> None:
    _ensure_table()
    db.execute("DELETE FROM dataset_dashboard WHERE dataset_id = ?", [dataset_id])


# ------------------------------------------------------------------ helpers --
def _counts(row: np.ndarray) -> Dict[str, int]:
    called = row[row != MISSING]
    return {
        "n_called": int(called.size),
        "het": int(np.sum(called == 1)),
        "hom_alt": int(np.sum(called == 2)),
    }


def _histogram(values: np.ndarray, bins: int = 18) -> List[Dict[str, Any]]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size < 2:
        return []
    lo, hi = float(vals.min()), float(vals.max())
    if not hi > lo:
        return []
    step = (hi - lo) / bins
    out = []
    for i in range(bins):
        a, b = lo + i * step, lo + (i + 1) * step
        in_bin = (vals >= a) & (vals < b) if i < bins - 1 else (vals >= a) & (vals <= b)
        out.append({"lo": a, "hi": b, "count": int(in_bin.sum())})
    return out


def _humanise(label: str, raw: str) -> str:
    """Turn a stored code into something a reader recognises.

    A composition panel reading "1.0 / 2.0" is data the reader has to decode. A
    phenotype file stores numbers, so the display has to do this — but only
    where the meaning is actually known, which is why this is keyed on the
    column name rather than applied to every 1/2 column in the file.
    """
    v = raw[:-2] if raw.endswith(".0") else raw
    key = label.strip().lower()
    if key in ("sex", "gender"):
        # PLINK .fam convention, which is also what the phenotype loader emits.
        return {"1": "Male", "2": "Female", "0": "Unknown", "-9": "Unknown"}.get(v, v)
    if v in ("0", "1"):
        if key in ("affected", "case", "status", "is_case", "phenotype"):
            return {"1": "Case", "0": "Control"}[v]
        return {"1": "Yes", "0": "No"}[v]
    return v


# capitalize() would render these as "Ldl" and "Bmi".
_ACRONYMS = {"ldl", "hdl", "bmi", "af", "maf", "prs", "tmb", "egfr", "bp", "hba1c"}


def _title(label: str) -> str:
    words = label.replace("_", " ").strip().split()
    return " ".join(w.upper() if w.lower() in _ACRONYMS else w.capitalize()
                    for w in words) or label


def _breakdown(values, label: str) -> Optional[Dict[str, Any]]:
    """Category counts, largest first — the composition panels."""
    vals = [_humanise(label, str(v)) for v in values
            if v is not None and str(v) not in ("", "nan", "None")]
    if not vals:
        return None
    counts: Dict[str, int] = {}
    for v in vals:
        counts[v] = counts.get(v, 0) + 1
    total = len(vals)
    # More than this and it is a list, not a composition — a stacked bar with
    # twenty segments communicates nothing.
    if len(counts) > 8:
        ordered = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
        head, tail = ordered[:7], ordered[7:]
        counts = dict(head)
        counts["Other ({} more)".format(len(tail))] = sum(v for _, v in tail)
    return {
        "label": _title(label),
        "column": label,
        "total": total,
        "groups": [{"name": k, "count": v, "pct": v / total}
                   for k, v in sorted(counts.items(), key=lambda kv: kv[1],
                                      reverse=True)],
    }


# ------------------------------------------------------------------ compute --
def compute(dataset_id: str) -> Dict[str, Any]:
    from . import registry, store as research_store

    ds = registry.get_dataset(dataset_id)
    if not ds:
        raise ValueError("unknown dataset: {}".format(dataset_id))
    prof = registry.get_profile(dataset_id) or {}

    anc = prof.get("ancestry") or {}
    rel = prof.get("relatedness") or {}
    phenos = prof.get("phenotypes") or {}

    # ---- KPI row -----------------------------------------------------------
    binary = {k: v for k, v in phenos.items() if v.get("kind") == "binary"}
    lead = max(binary.items(), key=lambda kv: kv[1].get("cases", 0) or 0,
               default=(None, {}))
    lead_name, lead_meta = lead

    kpis = [
        {"key": "Samples", "value": prof.get("n_samples"),
         "detail": "{} unrelated".format(prof.get("n_unrelated") or 0)},
        {"key": "Variants", "value": prof.get("n_variants"),
         "detail": (prof.get("density_class") or "").replace("_", " ")},
        {"key": "Genes annotated", "value": prof.get("n_annotated_genes") or 0,
         "detail": "with a gene symbol" if prof.get("has_gene_annotations")
                   else "no annotations supplied"},
    ]
    if lead_name:
        cases = lead_meta.get("cases") or 0
        controls = lead_meta.get("controls") or 0
        kpis.append({
            "key": "Cases / controls", "value": "{} / {}".format(cases, controls),
            "detail": "on " + lead_name})
    kpis.append({
        "key": "Call rate",
        "value": (round(100 * prof["mean_call_rate"], 1) if prof.get("mean_call_rate")
                  else None),
        "detail": "mean across samples", "unit": "%"})
    kpis.append({
        "key": "Related pairs", "value": rel.get("n_related_pairs", 0),
        "detail": "first-degree or closer" if rel.get("n_related_pairs")
                  else "none detected"})

    # ---- top genes ---------------------------------------------------------
    # The expensive part, and the reason this is cached.
    top_genes: List[Dict[str, Any]] = []
    genes_note = ""
    annotations = research_store.load_annotations(dataset_id) or {}
    if annotations:
        gm = research_store.load_genotypes(dataset_id)

        # Per-variant alternate allele frequency over called samples only.
        d = gm.dosages
        called = d != MISSING
        n_called_v = called.sum(axis=1)
        alt = np.where(called, np.maximum(d, 0), 0).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            af = np.where(n_called_v > 0, alt / (2.0 * np.maximum(n_called_v, 1)), np.nan)
        qualifying = np.isfinite(af) & (af > 0) & (af <= GENE_MAX_AF)

        by_gene: Dict[str, List[int]] = {}
        for i, v in enumerate(gm.variants):
            if not qualifying[i]:
                continue
            gene = (annotations.get(v.key) or {}).get("gene")
            if gene:
                by_gene.setdefault(str(gene), []).append(i)

        rows = []
        for gene, idx in by_gene.items():
            block = gm.dosages[np.array(idx)]
            # A sample counts once for the gene, however many of its variants
            # they carry — a gene-level rate, not a variant-level one.
            called_any = np.any(block != MISSING, axis=0)
            carries = np.any(block >= 1, axis=0) & called_any
            n_called = int(called_any.sum())
            n_carriers = int(carries.sum())
            if n_called < MIN_CALLED:
                continue
            rows.append({
                "gene": gene,
                "n_variants": len(idx),
                "carriers": n_carriers,
                "n_called": n_called,
                "pct": n_carriers / n_called,
                "suppressed": 0 < n_carriers < MIN_CELL,
            })
        rows.sort(key=lambda r: r["pct"], reverse=True)
        top_genes = rows[:TOP_GENES]
        genes_note = (
            "Carriers of a variant at or below {:.2g}% allele frequency, grouped "
            "by gene. The percentage divides by the samples actually CALLED at "
            "that gene, not by cohort size — a sample with no call there is not "
            "evidence of absence. Counting common variants too would put every "
            "gene at 100%. {} of {} annotated genes have a qualifying variant; "
            "the top {} are shown."
            .format(GENE_MAX_AF * 100, len(by_gene),
                    prof.get("n_annotated_genes") or len(by_gene), len(top_genes)))
    else:
        genes_note = ("No gene annotations were supplied with this dataset, so "
                      "variants cannot be grouped by gene. Attach an annotation "
                      "file to populate this panel.")

    # ---- composition -------------------------------------------------------
    composition: List[Dict[str, Any]] = []
    distributions: List[Dict[str, Any]] = []
    ph = None
    try:
        ph = research_store.load_phenotypes(dataset_id)
    except Exception:                                    # noqa: BLE001
        ph = None

    if ph is not None:
        for name, meta in phenos.items():
            kind = meta.get("kind")
            if name in ("father_id", "mother_id"):        # pedigree, not a trait
                continue
            try:
                values = ph.get(name)
            except KeyError:
                continue

            if kind in ("binary", "categorical"):
                b = _breakdown(values.tolist(), name)
                if b and 1 < len(b["groups"]) <= 8:
                    composition.append(b)
            elif kind in ("quantitative", "time_to_event"):
                bins = _histogram(values)
                if bins:
                    distributions.append({
                        "label": _title(name), "column": name, "bins": bins,
                        "mean": meta.get("mean"), "sd": meta.get("sd"),
                        "min": meta.get("min"), "max": meta.get("max"),
                    })

    # ---- variant landscape -------------------------------------------------
    spectrum = prof.get("maf_spectrum") or {}
    per_chrom = prof.get("variants_per_chromosome") or {}

    return {
        "dataset_id": dataset_id,
        "name": ds.get("name"),
        "genome_build": ds.get("genome_build") or prof.get("genome_build"),
        "source_format": ds.get("source_format"),
        "uploaded_at": ds.get("uploaded_at"),
        "computed_at": datetime.utcnow().isoformat(timespec="seconds"),

        "kpis": kpis,
        "top_genes": top_genes,
        "genes_note": genes_note,
        "composition": composition,
        "distributions": distributions,

        "spectrum": [
            {"label": "Monomorphic", "count": spectrum.get("monomorphic", 0)},
            {"label": "Ultra-rare <0.1%", "count": spectrum.get("ultra_rare_lt_0.1pct", 0)},
            {"label": "Rare 0.1–1%", "count": spectrum.get("rare_0.1_1pct", 0)},
            {"label": "Low freq 1–5%", "count": spectrum.get("low_freq_1_5pct", 0)},
            {"label": "Common ≥5%", "count": spectrum.get("common_ge_5pct", 0)},
        ],
        "per_chromosome": [{"chrom": c, "count": n} for c, n in per_chrom.items()],

        "ancestry": {"n_pcs": anc.get("n_pcs", 0),
                     "variance_explained": anc.get("variance_explained")},
        "relatedness": {
            "n_related_pairs": rel.get("n_related_pairs", 0),
            "max_kinship": rel.get("max_kinship"),
            "n_samples_with_relative": rel.get("n_samples_with_relative", 0),
            "n_trios": prof.get("n_trios", 0),
        },
        "ascertained": bool(prof.get("ascertained")),
        "ascertainment_rationale": prof.get("ascertainment_rationale") or "",
        "warnings": prof.get("warnings") or [],

        "caveat": (
            "Descriptive summary of the uploaded dataset. Frequencies are "
            "internal to this cohort and are not population reference values."),
    }
