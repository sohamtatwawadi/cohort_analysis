"""Dataset profiling — Part II §3.1.

    "On ingestion, the platform profiles the dataset and computes which
     analyses it can support. This replaces feature removal."

The profile is the evidence base for every gating decision in capability.py.
It is computed once at upload and stored, so that what a researcher was shown
can be reconstructed later.

Fields required by §3.1:

    n_samples · n_variants · variant density (genome-wide vs targeted)
    genome build · ancestry composition (PCA vs reference)
    relatedness structure (kinship distribution)
    phenotype completeness · case/control counts per phenotype
    quantitative phenotype distributions · missingness
    coverage/scope confidence · batch structure
    family structure present · longitudinal depth
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

from .types import MISSING, GenotypeMatrix, PhenotypeTable, as_float
from .validate import harmonise_chrom

# A dataset is treated as genome-wide if it has enough variants spread across
# enough chromosomes. Both conditions matter: a deeply sequenced single gene
# can have thousands of variants and is emphatically not genome-wide, and a
# panel covering 22 chromosomes with 40 variants is not either.
GENOME_WIDE_MIN_VARIANTS = 100_000
GENOME_WIDE_MIN_CHROMS = 20
EXOME_MIN_VARIANTS = 20_000

# Above this on one or two chromosomes, the dataset is chromosome-scale rather
# than a single gene. Set at the exome floor deliberately: no assay targeting
# one gene produces twenty thousand sites, so anything above it on a single
# chromosome is a region or a whole chromosome, not a gene.
SINGLE_CHROMOSOME_MIN_VARIANTS = EXOME_MIN_VARIANTS
EXOME_MIN_GENES = 5_000


@dataclass
class DataProfile:
    n_samples: int = 0
    n_variants: int = 0
    genome_build: Optional[str] = None

    # variant density
    n_chromosomes: int = 0
    chromosomes: List[str] = field(default_factory=list)
    variants_per_chromosome: Dict[str, int] = field(default_factory=dict)
    density_class: str = "targeted"        # genome_wide | exome | targeted
                                           # | single_chromosome | single_gene
    density_rationale: str = ""

    # genotype quality
    mean_call_rate: float = 0.0
    sample_call_rate_min: float = 0.0
    variant_call_rate_min: float = 0.0
    missingness: float = 0.0
    n_polymorphic: int = 0
    maf_spectrum: Dict[str, int] = field(default_factory=dict)

    # ancestry & relatedness
    ancestry: Dict[str, Any] = field(default_factory=dict)
    relatedness: Dict[str, Any] = field(default_factory=dict)
    n_unrelated: int = 0
    family_structure_present: bool = False
    n_trios: int = 0

    # phenotypes
    phenotypes: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    n_binary_phenotypes: int = 0
    n_quantitative_phenotypes: int = 0
    n_coded_phenotypes: int = 0
    has_controls: bool = False
    max_cases: int = 0
    has_time_to_event: bool = False
    longitudinal_depth: int = 0

    # provenance
    n_annotated_genes: int = 0
    has_gene_annotations: bool = False
    coverage_confidence: str = "unknown"   # declared | inferred | unknown
    batch_structure: Dict[str, Any] = field(default_factory=dict)
    ascertained: bool = False
    ascertainment_rationale: str = ""

    warnings: List[str] = field(default_factory=list)

    def as_dict(self) -> Dict[str, Any]:
        return asdict(self)


# ----------------------------------------------------------------- density ---
def classify_density(n_variants: int, chromosomes: Sequence[str],
                     n_genes: Optional[int] = None) -> Dict[str, str]:
    n_chrom = len(chromosomes)
    if n_variants >= GENOME_WIDE_MIN_VARIANTS and n_chrom >= GENOME_WIDE_MIN_CHROMS:
        return {"density_class": "genome_wide",
                "density_rationale": "{:,} variants across {} chromosomes".format(
                    n_variants, n_chrom)}
    if n_variants >= EXOME_MIN_VARIANTS and n_chrom >= GENOME_WIDE_MIN_CHROMS:
        return {"density_class": "exome",
                "density_rationale": "{:,} variants across {} chromosomes — consistent "
                                     "with exome capture".format(n_variants, n_chrom)}
    if n_chrom <= 2:
        # Chromosome count alone does not separate a single gene from a whole
        # chromosome, and calling 123,000 variants on chr22 "single_gene" is
        # simply wrong: a gene spans tens to hundreds of kilobases, a
        # chromosome tens of megabases. No assay of one gene yields this many
        # sites, so the variant count is the discriminator available here.
        if n_variants >= SINGLE_CHROMOSOME_MIN_VARIANTS:
            return {"density_class": "single_chromosome",
                    "density_rationale":
                        "{:,} variants on {} chromosome(s) — chromosome-scale, "
                        "but not genome-wide".format(n_variants, n_chrom)}
        return {"density_class": "single_gene",
                "density_rationale": "{:,} variants on {} chromosome(s)".format(
                    n_variants, n_chrom)}
    return {"density_class": "targeted",
            "density_rationale": "{:,} variants across {} chromosomes — consistent with "
                                 "a targeted panel".format(n_variants, n_chrom)}


# ------------------------------------------------------------------ profile --
def profile_dataset(gm: GenotypeMatrix,
                    phenotypes: Optional[PhenotypeTable] = None,
                    annotations: Optional[Dict[str, Dict[str, Any]]] = None,
                    coverage_confidence: str = "unknown",
                    batch: Optional[Dict[str, str]] = None,
                    pedigree: Optional[List[Dict[str, str]]] = None,
                    compute_genetics: bool = True,
                    max_pca_variants: int = 50_000) -> DataProfile:
    """Build the full §3.1 profile. Genetics (PCA, kinship) can be skipped for
    speed during validation and filled in by the job runner afterwards."""
    p = DataProfile()
    p.n_samples = gm.n_samples
    p.n_variants = gm.n_variants
    p.genome_build = gm.build
    p.coverage_confidence = coverage_confidence

    chroms: Dict[str, int] = {}
    for v in gm.variants:
        c = harmonise_chrom(v.chrom)
        chroms[c] = chroms.get(c, 0) + 1
    p.chromosomes = sorted(chroms)
    p.n_chromosomes = len(chroms)
    p.variants_per_chromosome = chroms
    p.__dict__.update(classify_density(p.n_variants, p.chromosomes))

    if gm.n_variants and gm.n_samples:
        observed = gm.dosages != MISSING
        p.missingness = float(1.0 - observed.mean())
        p.mean_call_rate = float(observed.mean())
        p.sample_call_rate_min = float(gm.call_rate_per_sample().min())
        p.variant_call_rate_min = float(gm.call_rate_per_variant().min())

        maf = gm.maf()
        finite = np.isfinite(maf)
        p.n_polymorphic = int(np.sum(finite & (maf > 0)))
        p.maf_spectrum = {
            "monomorphic": int(np.sum(finite & (maf <= 0))),
            "ultra_rare_lt_0.1pct": int(np.sum(finite & (maf > 0) & (maf < 0.001))),
            "rare_0.1_1pct": int(np.sum(finite & (maf >= 0.001) & (maf < 0.01))),
            "low_freq_1_5pct": int(np.sum(finite & (maf >= 0.01) & (maf < 0.05))),
            "common_ge_5pct": int(np.sum(finite & (maf >= 0.05))),
        }

    if annotations:
        p.n_annotated_genes = len({a.get("gene") for a in annotations.values()
                                   if a.get("gene")})
        p.has_gene_annotations = p.n_annotated_genes > 0

    if batch:
        counts: Dict[str, int] = {}
        for b in batch.values():
            counts[str(b)] = counts.get(str(b), 0) + 1
        p.batch_structure = {"n_batches": len(counts), "counts": counts}

    if pedigree:
        p.family_structure_present = True
        p.n_trios = sum(1 for r in pedigree
                        if r.get("pat") and r.get("pat") != "0"
                        and r.get("mat") and r.get("mat") != "0")

    if phenotypes is not None:
        _profile_phenotypes(p, phenotypes)

    if compute_genetics and gm.n_variants > 0 and gm.n_samples > 2:
        _profile_genetics(p, gm, max_pca_variants)

    _assess_ascertainment(p)
    return p


def _profile_phenotypes(p: DataProfile, ph: PhenotypeTable) -> None:
    for name in ph.columns:
        kind = ph.kind(name)
        entry: Dict[str, Any] = {
            "kind": kind,
            "label": ph.labels.get(name, name),
            "completeness": round(ph.completeness(name), 4),
        }
        v = as_float(ph.get(name))
        finite = np.isfinite(v)
        entry["n_present"] = int(finite.sum())

        if kind == "binary":
            counts = ph.case_control_counts(name)
            entry.update(counts)
            p.n_binary_phenotypes += 1
            p.max_cases = max(p.max_cases, counts["cases"])
            if counts["controls"] > 0:
                p.has_controls = True
        elif kind == "quantitative":
            vals = v[finite]
            if len(vals):
                entry.update({
                    "mean": float(np.mean(vals)), "sd": float(np.std(vals, ddof=1))
                    if len(vals) > 1 else 0.0,
                    "min": float(np.min(vals)), "max": float(np.max(vals)),
                    "median": float(np.median(vals)),
                    # Heavy skew changes the recommended model (§4.1 RINT).
                    "skew": float(_skew(vals)),
                })
            p.n_quantitative_phenotypes += 1
        elif kind == "time_to_event":
            p.has_time_to_event = True
        else:
            p.n_coded_phenotypes += 1

        p.phenotypes[name] = entry

    # PheWAS gates on the number of CODED phenotypes available to test.
    p.n_coded_phenotypes += p.n_binary_phenotypes


def _skew(v: np.ndarray) -> float:
    if len(v) < 3:
        return 0.0
    m, s = np.mean(v), np.std(v)
    return 0.0 if s == 0 else float(np.mean(((v - m) / s) ** 3))


def _profile_genetics(p: DataProfile, gm: GenotypeMatrix, max_variants: int) -> None:
    """Ancestry PCs and kinship. Both are prerequisites in the capability
    matrix, so a dataset that cannot support them is not merely missing a
    nice-to-have — it cannot run an association analysis at all."""
    try:
        from .stats import kinship as kin_mod
        from .stats import pca as pca_mod
    except ImportError:
        p.warnings.append("genetics kernels unavailable; ancestry and relatedness "
                          "were not profiled")
        return

    work = gm
    if gm.n_variants > max_variants:
        step = int(np.ceil(gm.n_variants / max_variants))
        mask = np.zeros(gm.n_variants, dtype=bool)
        mask[::step] = True
        work = gm.subset_variants(mask)

    try:
        pca = pca_mod.compute_pca(work, n_components=min(10, max(2, work.n_samples - 1)))
        p.ancestry = {
            "n_pcs": int(pca.components.shape[1]),
            "n_variants_used": int(pca.n_variants_used),
            "explained_variance_ratio": [round(float(x), 5)
                                         for x in pca.explained_variance_ratio[:10]],
            "pc1_spread": float(np.std(pca.components[:, 0])) if pca.components.size else 0.0,
            "warnings": list(pca.warnings),
        }
    except Exception as exc:                      # profiling must never block upload
        p.warnings.append("PCA failed: {}".format(exc))

    try:
        k = kin_mod.king_robust(work)
        p.relatedness = kin_mod.kinship_summary(k)
        p.n_unrelated = len(kin_mod.unrelated_set(k))
        if p.relatedness.get("first_degree", 0) or p.relatedness.get("second_degree", 0):
            p.family_structure_present = True
    except Exception as exc:
        p.warnings.append("kinship failed: {}".format(exc))
        p.n_unrelated = p.n_samples


def _assess_ascertainment(p: DataProfile) -> None:
    """Spec §4.7: "The platform must detect ascertainment from the data profile
    and label the output accordingly rather than presenting a clean curve."

    Penetrance from an ascertained clinical cohort is biased upward, often
    severely — families come to attention BECAUSE they are affected. This is a
    heuristic flag, not a proof, and it is phrased as such wherever it surfaces.
    """
    signals: List[str] = []
    if p.density_class in ("targeted", "single_gene"):
        signals.append("targeted panel rather than genome-wide genotyping")
    if not p.has_controls and p.n_binary_phenotypes:
        signals.append("no unselected control group")
    if p.family_structure_present:
        signals.append("family structure present (relatives recruited together)")
    if p.n_samples < 1000:
        signals.append("cohort size consistent with a clinical referral series")

    p.ascertained = len(signals) >= 2
    if p.ascertained:
        p.ascertainment_rationale = (
            "Profile suggests an ascertained (referral-selected) cohort: {}. "
            "Penetrance and prevalence estimates from such a cohort are biased "
            "upward and must not be quoted as population estimates."
            .format("; ".join(signals)))
