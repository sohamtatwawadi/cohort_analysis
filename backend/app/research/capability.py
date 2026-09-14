"""Capability gating — Part II §3. The core mechanism.

    "Shipping everything unconditionally means a researcher runs a GWAS on 180
     ascertained samples and gets a Manhattan plot full of noise that looks
     publishable. Shipping nothing means the platform is only useful to one
     tenant."

Three rules shape this module, and each one is a deliberate constraint rather
than an implementation convenience:

1. LOCKED IS VISIBLE, NOT HIDDEN (§3.2). An unavailable analysis renders with
   the specific unmet requirement and the observed value, because a researcher
   needs to know what to add. Hiding it just makes the product look incapable.

2. LOCKED CANNOT BE CLICKED THROUGH (§3.3). There is no "run anyway with a
   warning" path. From the spec: "A researcher under deadline will click
   through a warning. A locked control cannot be clicked through — it can only
   be unlocked by supplying the data."

3. OVERRIDE IS SEPARATE, RECORDED, AND STAMPS THE OUTPUT (§3.4). Borderline
   cases get an explicit escape hatch that requires typed justification, is
   audited, and marks every downstream artifact UNDERPOWERED. The stamp is
   derived from the override record, so it cannot be edited off a result.

Thresholds here are the spec's starting positions and are explicitly flagged as
needing calibration (§3.2 note). Where a power calculation is possible it is
better evidence than a raw n, and §5.1 computes one before execution.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

from .profile import DataProfile

# --------------------------------------------------------------- thresholds --
# Spec §3.2: "Thresholds above are starting positions and must be calibrated.
# They should be expressed as *power* rather than raw n wherever possible."
THRESHOLDS = {
    "association_min_n": 500,
    "association_min_cases": 50,
    "phewas_min_phenotypes": 100,
    "phewas_min_n": 2000,
    "gwas_min_variants": 100_000,
    "gwas_min_n": 2000,
    "burden_min_n": 1000,
    "heritability_min_n": 5000,
}


@dataclass
class Requirement:
    """One prerequisite, with the observed value so a locked card can say
    'Requires n >= 2,000 / Your dataset has 486 subjects'."""
    label: str
    met: bool
    observed: str
    remedy: Optional[str] = None

    def as_dict(self) -> Dict[str, Any]:
        return {"label": self.label, "met": self.met,
                "observed": self.observed, "remedy": self.remedy}


@dataclass
class Capability:
    analysis: str
    title: str
    available: bool
    requirements: List[Requirement] = field(default_factory=list)
    remedy: str = ""
    phase: str = ""

    @property
    def unmet(self) -> List[Requirement]:
        return [r for r in self.requirements if not r.met]

    def as_dict(self) -> Dict[str, Any]:
        return {
            "analysis": self.analysis, "title": self.title,
            "available": self.available, "phase": self.phase,
            "requirements": [r.as_dict() for r in self.requirements],
            "unmet": [r.as_dict() for r in self.unmet],
            "remedy": self.remedy,
        }


def _fmt(n: Any) -> str:
    return "{:,}".format(n) if isinstance(n, int) else str(n)


# ------------------------------------------------------------- the matrix ----
# §3.2 capability matrix. Each entry states its prerequisites against the
# profile. Order is the order the UI lists them.
def assess(profile: DataProfile) -> List[Capability]:
    p = profile
    caps: List[Capability] = []

    def cap(analysis: str, title: str, reqs: List[Requirement],
            remedy: str, phase: str = "") -> None:
        caps.append(Capability(
            analysis=analysis, title=title,
            available=all(r.met for r in reqs),
            requirements=reqs, remedy=remedy, phase=phase))

    has_genotypes = p.n_variants > 0 and p.n_samples > 0
    genome_wide = p.density_class in ("genome_wide",)
    exome_or_better = p.density_class in ("genome_wide", "exome")
    n_pcs = int(p.ancestry.get("n_pcs", 0) or 0)

    # ---- always-on descriptive analyses ----------------------------------
    cap("carrier_frequency", "Carrier / allele frequency", [
        Requirement("Genotype data", has_genotypes,
                    "{} variants x {} samples".format(_fmt(p.n_variants), _fmt(p.n_samples))),
        Requirement("Coverage resolvable", p.coverage_confidence != "unknown"
                    or has_genotypes,
                    "coverage {}".format(p.coverage_confidence)),
    ], "Upload genotype data.", "R1")

    cap("diagnostic_yield", "Diagnostic yield", [
        Requirement("Phenotype or indication per subject", bool(p.phenotypes),
                    "{} phenotype column(s)".format(len(p.phenotypes)),
                    "Attach a phenotype file with an indication per subject."),
        Requirement("Gene-disease reference table", True, "curated table present"),
    ], "Attach per-subject phenotypes.", "R1")

    cap("zygosity", "Zygosity / inheritance", [
        Requirement("Genotype data", has_genotypes, _fmt(p.n_variants) + " variants"),
        Requirement("Sex available", bool(p.phenotypes.get("sex")) or True,
                    "sex inferred from X heterozygosity where available"),
    ], "Upload genotype data.", "R1")

    cap("population_frequency", "Population frequency", [
        Requirement("Reference allele frequencies", True, "gnomAD / GA100K available"),
        Requirement("Ancestry labels or PCs", bool(p.ancestry) or bool(p.phenotypes),
                    "{} ancestry PCs computed".format(n_pcs),
                    "Provide self-reported ancestry, or enough variants to compute PCs."),
    ], "Provide ancestry information.", "R1")

    cap("segregation", "Segregation / de novo", [
        Requirement("Trio or family structure", p.family_structure_present,
                    "{} trio(s); family structure {}".format(
                        p.n_trios, "present" if p.family_structure_present else "absent"),
                    "Upload a PED file or a pedigree column."),
    ], "Upload pedigree information.", "R1")

    # ---- association -----------------------------------------------------
    assoc_n_ok = p.n_samples >= THRESHOLDS["association_min_n"]
    assoc_pheno_ok = (p.n_binary_phenotypes + p.n_quantitative_phenotypes) > 0
    assoc_cases_ok = (p.max_cases >= THRESHOLDS["association_min_cases"]
                      or p.n_quantitative_phenotypes > 0)
    cap("association", "Association (single variant)", [
        Requirement("Case/control or quantitative phenotype", assoc_pheno_ok,
                    "{} binary, {} quantitative".format(
                        p.n_binary_phenotypes, p.n_quantitative_phenotypes),
                    "Attach an outcome column."),
        Requirement("Adequate cases", assoc_cases_ok,
                    "{} cases in the largest binary phenotype".format(p.max_cases),
                    "Needs at least {} cases, or a quantitative outcome.".format(
                        THRESHOLDS["association_min_cases"])),
        Requirement("n >= {}".format(_fmt(THRESHOLDS["association_min_n"])), assoc_n_ok,
                    "{} subjects".format(_fmt(p.n_samples))),
        Requirement("Ancestry PCs", n_pcs > 0,
                    "{} PCs computed".format(n_pcs),
                    "Needs enough variants to compute ancestry PCs, or supplied PCs."),
    ], "Attach an outcome and ensure n is adequate.", "R2")

    # ---- PheWAS ----------------------------------------------------------
    cap("phewas", "PheWAS", [
        Requirement("≥ {} coded phenotypes".format(THRESHOLDS["phewas_min_phenotypes"]),
                    p.n_coded_phenotypes >= THRESHOLDS["phewas_min_phenotypes"],
                    "{} coded phenotype(s)".format(_fmt(p.n_coded_phenotypes)),
                    "Upload an EHR/ICD phenotype table."),
        Requirement("Exposure definable", has_genotypes,
                    "genotypes present" if has_genotypes else "no genotypes"),
        Requirement("n >= {}".format(_fmt(THRESHOLDS["phewas_min_n"])),
                    p.n_samples >= THRESHOLDS["phewas_min_n"],
                    "{} subjects".format(_fmt(p.n_samples))),
    ], "Upload a coded phenotype table (ICD/PheCode) for a larger cohort.", "R3")

    # ---- GWAS ------------------------------------------------------------
    cap("gwas", "GWAS", [
        Requirement("Genome-wide genotypes", genome_wide,
                    p.density_rationale or "{} variants".format(_fmt(p.n_variants)),
                    "Upload array, imputed or WGS data."),
        Requirement("≥ {} variants".format(_fmt(THRESHOLDS["gwas_min_variants"])),
                    p.n_variants >= THRESHOLDS["gwas_min_variants"],
                    "{} variants".format(_fmt(p.n_variants))),
        Requirement("n >= {}".format(_fmt(THRESHOLDS["gwas_min_n"])),
                    p.n_samples >= THRESHOLDS["gwas_min_n"],
                    "{} subjects".format(_fmt(p.n_samples))),
        Requirement("Case/control or quantitative phenotype", assoc_pheno_ok,
                    "{} binary, {} quantitative".format(
                        p.n_binary_phenotypes, p.n_quantitative_phenotypes)),
        Requirement("Unrelated control group", p.has_controls,
                    "controls {}".format("identified" if p.has_controls
                                         else "not identified"),
                    "Include unselected controls."),
        Requirement("Ancestry PCs", n_pcs > 0, "{} PCs computed".format(n_pcs)),
    ], "Upload array or WGS data with controls to enable this analysis.", "R4")

    # ---- burden / SKAT ---------------------------------------------------
    cap("burden", "Burden / SKAT / SKAT-O", [
        Requirement("Exome or genome sequencing", exome_or_better,
                    p.density_rationale or p.density_class,
                    "Upload exome or genome data."),
        Requirement("Case-control outcome", p.n_binary_phenotypes > 0,
                    "{} binary phenotype(s)".format(p.n_binary_phenotypes),
                    "Attach a case/control outcome."),
        Requirement("n >= {}".format(_fmt(THRESHOLDS["burden_min_n"])),
                    p.n_samples >= THRESHOLDS["burden_min_n"],
                    "{} subjects".format(_fmt(p.n_samples))),
        Requirement("Gene annotation per variant", p.has_gene_annotations,
                    "{} annotated gene(s)".format(_fmt(p.n_annotated_genes)),
                    "Upload an annotated VCF, or annotate the dataset — a gene-based "
                    "test has nothing to group by without it."),
        Requirement("Qualifying-variant definition", True,
                    "variant-set builder available"),
    ], "Upload exome/genome data with case-control status.", "R5")

    # ---- PRS -------------------------------------------------------------
    cap("prs", "Polygenic risk score", [
        Requirement("Genome-wide genotypes", genome_wide,
                    p.density_rationale or p.density_class,
                    "Upload array, imputed or WGS data."),
        Requirement("External score weights", True,
                    "PGS Catalog import or custom upload"),
        Requirement("Ancestry-matched validation possible", n_pcs > 0 or bool(p.ancestry),
                    "{} PCs computed".format(n_pcs),
                    "Ancestry information is required to report stratified performance."),
    ], "Upload genome-wide genotypes and a score definition.", "R6")

    # ---- survival / penetrance -------------------------------------------
    cap("survival", "Survival / penetrance", [
        Requirement("Time-to-event with follow-up", p.has_time_to_event,
                    "time-to-event column {}".format(
                        "present" if p.has_time_to_event else "absent"),
                    "Attach time and event columns."),
        Requirement("Unascertained carriers for unbiased penetrance", not p.ascertained,
                    "cohort appears {}".format(
                        "ascertained" if p.ascertained else "unselected"),
                    "Penetrance from an ascertained cohort is biased upward; the "
                    "analysis can still run but the output will be labelled."),
    ], "Attach follow-up time and event status.", "R6")

    # ---- summary-statistic methods ---------------------------------------
    cap("fine_mapping", "Fine-mapping", [
        Requirement("A completed GWAS result", False, "no GWAS results in this project",
                    "Run a GWAS first."),
        Requirement("LD reference panel", False, "no LD reference configured",
                    "Configure an ancestry-matched LD reference."),
    ], "Run a GWAS and configure an LD reference.", "R7")

    cap("colocalization", "Colocalization", [
        Requirement("Two association datasets over the same region", False,
                    "fewer than two association results",
                    "Run or upload two association analyses."),
    ], "Provide two association datasets.", "R7")

    cap("heritability", "Heritability / genetic correlation", [
        Requirement("GWAS summary statistics", False, "no summary statistics available",
                    "Run a GWAS or upload summary statistics."),
        Requirement("n >= {}".format(_fmt(THRESHOLDS["heritability_min_n"])),
                    p.n_samples >= THRESHOLDS["heritability_min_n"],
                    "{} subjects".format(_fmt(p.n_samples))),
    ], "Upload GWAS summary statistics from a larger cohort.", "R7")

    return caps


def assess_dict(profile: DataProfile) -> Dict[str, Dict[str, Any]]:
    return {c.analysis: c.as_dict() for c in assess(profile)}


# ---------------------------------------------------------------- override ---
UNDERPOWERED_STAMP = "UNDERPOWERED — threshold overridden"

# Some requirements are not thresholds and cannot be overridden: they are
# statements that the data required for the analysis does not exist. You cannot
# override your way to genome-wide genotypes you did not upload.
NON_OVERRIDABLE = {
    "Genome-wide genotypes",
    "Exome or genome sequencing",
    "Gene annotation per variant",
    "Time-to-event with follow-up",
    "Trio or family structure",
    "A completed GWAS result",
    "LD reference panel",
    "GWAS summary statistics",
    "Two association datasets over the same region",
    "Case/control or quantitative phenotype",
}


def overridable(cap: Capability) -> Dict[str, Any]:
    """Which of an analysis's unmet requirements could be overridden (§3.4).

    Borderline sample size is a judgement call a researcher may reasonably make
    and defend. Missing data is not a judgement call.
    """
    blocking = [r for r in cap.unmet if r.label in NON_OVERRIDABLE]
    soft = [r for r in cap.unmet if r.label not in NON_OVERRIDABLE]
    return {
        "can_override": bool(soft) and not blocking,
        "overridable": [r.as_dict() for r in soft],
        "blocking": [r.as_dict() for r in blocking],
        "reason": ("Missing data cannot be overridden — it has to be supplied."
                   if blocking else
                   "Borderline thresholds may be overridden with a recorded "
                   "justification. Results will be stamped."),
    }
