"""Core data structures shared across Research Mode.

Part II §1: "The platform runs one engine in two modes." These types are the
boundary between ingestion (which differs by mode) and analysis (which does
not) — a VariMAT load, a VCF and a PLINK fileset all land here.

Genotype convention, fixed everywhere:

    dosage = count of ALT alleles = 0 | 1 | 2, with -1 for missing

-1 rather than NaN because the matrix is int8: a biobank-scale genotype matrix
held as float64 is 8x the memory for no gain, and every consumer has to handle
missingness explicitly rather than having it silently propagate through an
arithmetic operation.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

MISSING = -1

# Genome builds we will name. A build is never inferred from coordinates alone
# (Part II §2.3) — GRCh37 and GRCh38 positions overlap in range, so guessing is
# a coin flip that silently corrupts every downstream annotation.
BUILDS = ("GRCh37", "GRCh38", "T2T-CHM13")


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    vid: Optional[str] = None
    qual: Optional[float] = None
    filter_status: Optional[str] = None
    info: Optional[Dict[str, Any]] = None

    @property
    def key(self) -> str:
        return "{}:{}:{}:{}".format(self.chrom, self.pos, self.ref, self.alt)

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1

    def __repr__(self) -> str:  # keep debugging output readable
        return "Variant({})".format(self.key)


@dataclass
class GenotypeMatrix:
    """Variants x samples. The unit every analysis consumes."""

    sample_ids: List[str]
    variants: List[Variant]
    dosages: np.ndarray            # shape (n_variants, n_samples), int8
    build: Optional[str] = None
    source_format: Optional[str] = None
    source_files: List[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.dosages.ndim != 2:
            raise ValueError("dosages must be 2-D (variants x samples)")
        m, n = self.dosages.shape
        if m != len(self.variants):
            raise ValueError(
                "dosage rows ({}) != variants ({})".format(m, len(self.variants)))
        if n != len(self.sample_ids):
            raise ValueError(
                "dosage cols ({}) != samples ({})".format(n, len(self.sample_ids)))

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    @property
    def n_variants(self) -> int:
        return len(self.variants)

    def sample_index(self) -> Dict[str, int]:
        return {s: i for i, s in enumerate(self.sample_ids)}

    def subset_samples(self, keep: Sequence[str]) -> "GenotypeMatrix":
        idx_of = self.sample_index()
        cols = [idx_of[s] for s in keep if s in idx_of]
        return GenotypeMatrix(
            sample_ids=[self.sample_ids[c] for c in cols],
            variants=list(self.variants),
            dosages=self.dosages[:, cols],
            build=self.build, source_format=self.source_format,
            source_files=list(self.source_files))

    def subset_variants(self, mask: np.ndarray) -> "GenotypeMatrix":
        mask = np.asarray(mask, dtype=bool)
        return GenotypeMatrix(
            sample_ids=list(self.sample_ids),
            variants=[v for v, keep in zip(self.variants, mask) if keep],
            dosages=self.dosages[mask, :],
            build=self.build, source_format=self.source_format,
            source_files=list(self.source_files))

    # ---------------------------------------------------------- statistics --
    def call_rate_per_variant(self) -> np.ndarray:
        return (self.dosages != MISSING).mean(axis=1)

    def call_rate_per_sample(self) -> np.ndarray:
        return (self.dosages != MISSING).mean(axis=0)

    def allele_frequency(self) -> np.ndarray:
        """ALT allele frequency per variant, ignoring missing calls."""
        d = self.dosages
        observed = d != MISSING
        n_called = observed.sum(axis=1)
        alt_count = np.where(observed, d, 0).sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            af = np.where(n_called > 0, alt_count / (2.0 * n_called), np.nan)
        return af

    def maf(self) -> np.ndarray:
        af = self.allele_frequency()
        return np.minimum(af, 1.0 - af)


@dataclass
class PhenotypeTable:
    """Sample-indexed phenotypes. Columns are typed because the analysis engine
    picks its model from the type (Part II §4.1 automatic model selection)."""

    sample_ids: List[str]
    columns: Dict[str, np.ndarray]           # name -> values, aligned to sample_ids
    kinds: Dict[str, str]                    # name -> binary | quantitative | categorical
    labels: Dict[str, str] = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    def get(self, name: str) -> np.ndarray:
        if name not in self.columns:
            raise KeyError("no phenotype column '{}'".format(name))
        return self.columns[name]

    def kind(self, name: str) -> str:
        return self.kinds.get(name, "categorical")

    def completeness(self, name: str) -> float:
        v = self.columns.get(name)
        if v is None or len(v) == 0:
            return 0.0
        return float(np.isfinite(_as_float(v)).mean())

    def case_control_counts(self, name: str) -> Dict[str, int]:
        v = _as_float(self.get(name))
        finite = np.isfinite(v)
        return {
            "cases": int(np.sum(finite & (v == 1))),
            "controls": int(np.sum(finite & (v == 0))),
            "missing": int(np.sum(~finite)),
        }


def _as_float(v: np.ndarray) -> np.ndarray:
    """Coerce a column to float with NaN for missing, whatever it arrived as."""
    if v.dtype.kind in "fc":
        return v.astype(float)
    if v.dtype.kind in "iub":
        return v.astype(float)
    out = np.full(len(v), np.nan)
    for i, x in enumerate(v):
        try:
            out[i] = float(x)
        except (TypeError, ValueError):
            out[i] = np.nan
    return out


as_float = _as_float
