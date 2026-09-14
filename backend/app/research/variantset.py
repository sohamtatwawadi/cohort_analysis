"""Qualifying-variant builder — Part II §4.4.

    "Qualifying-variant builder: MAF threshold, consequence class, in-silico
     predictor thresholds (CADD, REVEL, SpliceAI, AlphaMissense), ClinVar
     classification, LoF confidence, custom. Saved as a reusable, versioned
     variant set."

Versioned matters. A burden test is only interpretable against the exact
qualifying-variant definition that produced it, and those definitions get
tweaked — a result computed under "REVEL > 0.5" and compared against one
computed under "REVEL > 0.75" is not a comparison. The definition is stored
with a version and referenced by every job that uses it.
"""
from __future__ import annotations

import json
import uuid
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from .. import config, db
from .types import GenotypeMatrix

# Consequence classes, ordered by severity. Used by the presets.
LOF_CONSEQUENCES = ["Nonsense", "Frameshift", "Splice site", "Start loss", "Stop loss",
                    "Exon deletion"]
MISSENSE_CONSEQUENCES = ["Missense"]

PREDICTORS = ["CADD", "REVEL", "SpliceAI", "AlphaMissense"]

# §4.4: "Presets: ultra-rare LoF, rare damaging missense, P/LP only."
PRESETS: Dict[str, Dict[str, Any]] = {
    "ultra_rare_lof": {
        "label": "Ultra-rare loss-of-function",
        "max_maf": 0.001,
        "consequences": LOF_CONSEQUENCES,
        "lof_confidence": "HC",
        "description": "High-confidence LoF at MAF < 0.1%. The most conservative "
                       "and usually best-powered rare-variant definition.",
    },
    "rare_damaging_missense": {
        "label": "Rare damaging missense",
        "max_maf": 0.01,
        "consequences": MISSENSE_CONSEQUENCES,
        "predictors": {"REVEL": {"min": 0.5}},
        "description": "Missense at MAF < 1% predicted damaging by REVEL.",
    },
    "plp_only": {
        "label": "Pathogenic / likely pathogenic only",
        "max_maf": 0.01,
        "clinvar": ["Pathogenic", "Likely pathogenic"],
        "description": "Clinically classified P/LP variants only.",
    },
    "lof_plus_damaging_missense": {
        "label": "LoF + damaging missense",
        "max_maf": 0.01,
        "consequences": LOF_CONSEQUENCES + MISSENSE_CONSEQUENCES,
        "predictors": {"REVEL": {"min": 0.5}},
        "description": "The common composite definition: any LoF, plus missense "
                       "predicted damaging.",
    },
}


def apply_definition(gm: GenotypeMatrix, definition: Dict[str, Any],
                     annotations: Optional[Dict[str, Dict[str, Any]]] = None
                     ) -> Dict[str, Any]:
    """Select qualifying variants and report WHY each filter dropped what it did.

    The per-filter counts are not decoration: a burden test on 2 qualifying
    variants and one on 200 are different analyses, and the researcher needs to
    see which filter did the damage.
    """
    annotations = annotations or {}
    n = gm.n_variants
    mask = np.ones(n, dtype=bool)
    steps: List[Dict[str, Any]] = [{"filter": "all variants", "kept": int(n), "dropped": 0}]

    def apply(name: str, keep: np.ndarray) -> None:
        nonlocal mask
        before = int(mask.sum())
        mask = mask & keep
        after = int(mask.sum())
        steps.append({"filter": name, "kept": after, "dropped": before - after})

    max_maf = definition.get("max_maf")
    if max_maf is not None:
        maf = gm.maf()
        apply("MAF ≤ {:g}".format(max_maf),
              np.isfinite(maf) & (maf <= float(max_maf)) & (maf > 0))

    min_maf = definition.get("min_maf")
    if min_maf is not None:
        maf = gm.maf()
        apply("MAF ≥ {:g}".format(min_maf), np.isfinite(maf) & (maf >= float(min_maf)))

    consequences = definition.get("consequences")
    if consequences:
        want = set(consequences)
        keep = np.array([
            (annotations.get(v.key, {}).get("consequence") in want) for v in gm.variants])
        apply("consequence in {}".format(", ".join(sorted(want))), keep)

    clinvar = definition.get("clinvar")
    if clinvar:
        want = set(clinvar)
        keep = np.array([
            (annotations.get(v.key, {}).get("clinvar_sig") in want) for v in gm.variants])
        apply("ClinVar in {}".format(", ".join(sorted(want))), keep)

    for predictor, bounds in (definition.get("predictors") or {}).items():
        lo = bounds.get("min")
        hi = bounds.get("max")
        vals = np.array([
            _num(annotations.get(v.key, {}).get(predictor)) for v in gm.variants])
        keep = np.isfinite(vals)
        label = predictor
        if lo is not None:
            keep = keep & (vals >= float(lo))
            label += " ≥ {:g}".format(lo)
        if hi is not None:
            keep = keep & (vals <= float(hi))
            label += " ≤ {:g}".format(hi)
        apply(label, keep)

    lof_conf = definition.get("lof_confidence")
    if lof_conf:
        keep = np.array([
            (annotations.get(v.key, {}).get("lof_confidence") == lof_conf)
            for v in gm.variants])
        apply("LoF confidence = {}".format(lof_conf), keep)

    genes = definition.get("genes")
    if genes:
        want = set(genes)
        keep = np.array([
            (annotations.get(v.key, {}).get("gene") in want) for v in gm.variants])
        apply("gene in set ({} genes)".format(len(want)), keep)

    return {
        "mask": mask,
        "n_qualifying": int(mask.sum()),
        "funnel": steps,
        "definition": definition,
    }


def _num(v: Any) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def group_by_gene(gm: GenotypeMatrix, mask: np.ndarray,
                  annotations: Dict[str, Dict[str, Any]]) -> Dict[str, List[int]]:
    """Qualifying variant indices per gene — the unit a burden test operates on."""
    out: Dict[str, List[int]] = {}
    for i, v in enumerate(gm.variants):
        if not mask[i]:
            continue
        gene = (annotations.get(v.key) or {}).get("gene")
        if gene:
            out.setdefault(gene, []).append(i)
    return out


# ------------------------------------------------------------- persistence ---
def save_set(project_id: str, name: str, definition: Dict[str, Any],
             created_by: Optional[str] = None) -> Dict[str, Any]:
    """Store a definition. A re-save under the same name bumps the version
    rather than overwriting, so a past result's definition stays retrievable."""
    prev = db.scalar("""
        SELECT MAX(version) FROM variant_set WHERE project_id = ? AND name = ?
    """, [project_id, name], default=0) or 0
    sid = "vs-" + uuid.uuid4().hex[:10]
    db.insert_rows(
        "variant_set",
        ["set_id", "project_id", "name", "definition_json", "version",
         "created_by", "created_at"],
        [[sid, project_id, name, json.dumps(definition), int(prev) + 1,
          created_by or config.DEFAULT_USER, datetime.utcnow()]])
    return {"set_id": sid, "name": name, "version": int(prev) + 1,
            "definition": definition}


def list_sets(project_id: str) -> List[Dict[str, Any]]:
    rows = db.rows("""
        SELECT set_id, name, definition_json, version, created_by, created_at
        FROM variant_set WHERE project_id = ? ORDER BY created_at DESC
    """, [project_id])
    for r in rows:
        r["definition"] = json.loads(r.pop("definition_json"))
    return rows


def get_set(set_id: str) -> Optional[Dict[str, Any]]:
    r = db.row("SELECT * FROM variant_set WHERE set_id = ?", [set_id])
    if r:
        r["definition"] = json.loads(r.pop("definition_json"))
    return r
