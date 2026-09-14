"""Persistence for uploaded genotype and phenotype data.

Kept out of DuckDB deliberately. A genotype matrix is dense numeric data of
shape (variants x samples); storing it as rows would be tens of millions of
tuples for a dataset that is a few megabytes as an array. It lives on disk as
compressed numpy, with the metadata in the relational store.

Every path is namespaced by dataset_id under the project's directory so that
`registry.delete_dataset` can remove the payload as well as the derived tables
— §8 requires deletion to actually propagate.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from .. import config
from .types import GenotypeMatrix, PhenotypeTable, Variant

RESEARCH_DIR = Path(config.DB_PATH).parent / "research"


def dataset_dir(dataset_id: str) -> Path:
    return RESEARCH_DIR / dataset_id


def save_genotypes(dataset_id: str, gm: GenotypeMatrix) -> str:
    d = dataset_dir(dataset_id)
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / "genotypes.npz", dosages=gm.dosages)
    (d / "variants.json").write_text(json.dumps([
        {"chrom": v.chrom, "pos": v.pos, "ref": v.ref, "alt": v.alt, "vid": v.vid}
        for v in gm.variants]))
    (d / "samples.json").write_text(json.dumps({
        "sample_ids": gm.sample_ids, "build": gm.build,
        "source_format": gm.source_format, "source_files": gm.source_files}))
    return str(d)


def load_genotypes(dataset_id: str) -> GenotypeMatrix:
    d = dataset_dir(dataset_id)
    if not (d / "genotypes.npz").exists():
        raise FileNotFoundError("no stored genotypes for dataset {}".format(dataset_id))
    dosages = np.load(d / "genotypes.npz")["dosages"]
    variants = [Variant(**v) for v in json.loads((d / "variants.json").read_text())]
    meta = json.loads((d / "samples.json").read_text())
    return GenotypeMatrix(
        sample_ids=meta["sample_ids"], variants=variants, dosages=dosages,
        build=meta.get("build"), source_format=meta.get("source_format"),
        source_files=meta.get("source_files", []))


def save_phenotypes(dataset_id: str, ph: PhenotypeTable) -> None:
    d = dataset_dir(dataset_id)
    d.mkdir(parents=True, exist_ok=True)
    # Categorical columns (ancestry, site, batch) must survive as strings.
    # Coercing everything to float turns ancestry labels into NaN, which then
    # silently disables the ancestry-stratified PRS reporting §4.5 makes
    # mandatory — a governance requirement lost to a serialisation shortcut.
    columns: Dict[str, Any] = {}
    for k, v in ph.columns.items():
        arr = np.asarray(v)
        if ph.kinds.get(k) == "categorical" or arr.dtype.kind in "OUS":
            columns[k] = [None if x is None else str(x) for x in arr]
        else:
            columns[k] = np.asarray(arr, dtype=float).tolist()

    payload = {
        "sample_ids": ph.sample_ids,
        "kinds": ph.kinds,
        "labels": ph.labels,
        "columns": columns,
    }
    (d / "phenotypes.json").write_text(json.dumps(payload))


def load_phenotypes(dataset_id: str) -> Optional[PhenotypeTable]:
    p = dataset_dir(dataset_id) / "phenotypes.json"
    if not p.exists():
        return None
    payload = json.loads(p.read_text())
    kinds = payload["kinds"]
    columns: Dict[str, np.ndarray] = {}
    for k, v in payload["columns"].items():
        if kinds.get(k) == "categorical":
            columns[k] = np.asarray(v, dtype=object)
        else:
            columns[k] = np.asarray(v, dtype=float)
    return PhenotypeTable(
        sample_ids=payload["sample_ids"], columns=columns,
        kinds=kinds, labels=payload.get("labels", {}))


def save_annotations(dataset_id: str, annotations: Dict[str, Dict[str, Any]]) -> None:
    """Per-variant annotation (gene, consequence, predictor scores).

    Stored with the dataset rather than passed in each job spec. A gene-based
    test needs a gene for every variant; requiring the caller to post tens of
    thousands of annotations with every request made the analysis unusable from
    the UI, which is why it only ever worked from a script.
    """
    d = dataset_dir(dataset_id)
    d.mkdir(parents=True, exist_ok=True)
    (d / "annotations.json").write_text(json.dumps(annotations))


def load_annotations(dataset_id: str) -> Dict[str, Dict[str, Any]]:
    p = dataset_dir(dataset_id) / "annotations.json"
    return json.loads(p.read_text()) if p.exists() else {}


def annotation_summary(dataset_id: str) -> Dict[str, Any]:
    ann = load_annotations(dataset_id)
    genes = {a.get("gene") for a in ann.values() if a.get("gene")}
    return {"n_annotated_variants": len(ann), "n_genes": len(genes)}


def save_covariates(dataset_id: str, name: str, values: np.ndarray) -> None:
    """Ancestry PCs and similar derived covariates, computed once at profiling
    rather than recomputed per analysis."""
    d = dataset_dir(dataset_id)
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / "cov_{}.npz".format(name), values=values)


def load_covariates(dataset_id: str, name: str) -> Optional[np.ndarray]:
    p = dataset_dir(dataset_id) / "cov_{}.npz".format(name)
    return np.load(p)["values"] if p.exists() else None


def delete_payload(dataset_id: str) -> bool:
    d = dataset_dir(dataset_id)
    if d.exists():
        shutil.rmtree(d)
        return True
    return False
