"""Deterministic synthetic germline cohort generator.

Purpose: make the tool runnable and testable before real VariMAT files land,
and give the hand-count verification test (spec E03.9) a fixed dataset whose
answers can be recomputed independently.

It is seeded, so the same seed always produces the same store, which is what
lets "any saved cohort re-runs to an identical number" (E03.6) be asserted in
CI rather than hoped for.

Three properties are modelled deliberately because the tool exists to handle
them, and code that never meets them is untested code:

  * families           most subjects are singletons, some are trios and
                       extended families -> probands-only has something to do
  * consent variety    including Withdrawn and SF-declined -> the consent
                       gates have something to exclude
  * registry gaps      two test codes with incomplete reportable scope, plus a
                       slice of runs with NO test code at all -> the §3.4
                       fingerprint inference path is actually exercised

Findings are written the way a real file would carry them: `zygosity_raw` holds
Het/Hom ONLY. Hemizygous and compound-het are produced by the derivation pass
in store.py, never written here (spec §3.5).
"""
from __future__ import annotations

import math
import random
from datetime import date, timedelta
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .. import config
from ..reference.genes import GENE_DISEASE, GENE_LIST
from ..reference.tests import (ANCESTRY_WEIGHTS, CONSENT_CLASSES, HPO_BY_INDICATION,
                               PIPELINE_VERSIONS, REFERRAL_SOURCES, SAMPLE_TYPES,
                               TEST_CODES)
from .clinical import age_bucket

TODAY = date(2026, 8, 13)
AAS = "ACDEFGHIKLMNPQRSTVWY"

TEST_WEIGHTS = {
    "MG-HBOC-27": 112, "MG-CARRIER-113": 96, "MG-THAL-1": 74, "MG-CES": 62,
    "MG-NDD-58": 58, "MG-CARDIO-42": 46, "MG-LYNCH-7": 42, "MG-WES": 38,
    "MG-HEARING-96": 32, "MG-WGS": 18,
}
CONSENT_WEIGHTS = {
    "Full research + secondary findings": 46,
    "Clinical only, secondary findings declined": 31,
    "Clinical only, research declined": 18,
    "Withdrawn": 5,
}
FAMILY_SIZE_WEIGHTS = {1: 54, 2: 14, 3: 22, 4: 6, 5: 4}
SAMPLE_TYPE_WEIGHTS = {"Blood (EDTA)": 74, "Saliva": 13, "Buccal swab": 8,
                       "Dried blood spot": 5}
TRUNCATING = {"Frameshift", "Nonsense", "Splice site", "Exon deletion"}

# Fraction of runs whose test code is absent from the LIMS export, forcing the
# fingerprint-clustering path in §3.4 to resolve their coverage.
UNLABELLED_RUN_FRACTION = 0.14

# Relative draw weight for an ACMG SF gene on a broad test. Sets the §G09
# secondary-findings rate; see the comment at the draw site for why it is this
# small. Calibrated to land the SF rate in the low single-digit percents.
SF_GENE_DRAW_WEIGHT = 0.03

# Fraction of each gene's allele catalogue that is pathogenic or likely
# pathogenic. Fixed rather than drawn, so per-test-code yield reflects the
# cohort rather than catalogue luck.
CATALOG_PLP_FRACTION = 0.22

# Substitutions are NOT uniform in real sequence. Transitions (purine<->purine,
# pyrimidine<->pyrimidine) outnumber transversions roughly 2:1 genome-wide and
# ~3:1 in coding regions, largely because deaminated methyl-CpG gives C>T.
# Drawing uniformly produces Ti/Tv ~0.5, which makes the QC metric that exists
# to detect false-positive calls flag every sample in the cohort.
TRANSITIONS = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]
TRANSVERSIONS = [("A", "C"), ("A", "T"), ("C", "A"), ("C", "G"),
                 ("G", "C"), ("G", "T"), ("T", "A"), ("T", "G")]
TARGET_TI_TV = 3.0          # coding-region expectation

# A handful of samples are deliberately poor, so the QC screen has something
# real to catch. A demo where every sample passes does not show that the check
# works — it shows that the check ran.
CONTAMINATED_FRACTION = 0.018
LOW_DEPTH_FRACTION = 0.012


def _simulated_run_qc(rng: "Rng", contaminated: bool, low_depth: bool,
                      sex: str) -> Dict[str, Any]:
    """Per-run QC over the full call set a real pipeline would see."""
    n_called = int(_clamp(rng.normal(2600 if not low_depth else 900, 420), 200, 6000))
    ti_tv = rng.normal(1.55, 0.09) if contaminated else rng.normal(3.02, 0.11)
    het_hom = rng.normal(3.4, 0.35) if contaminated else rng.normal(1.62, 0.14)
    het_vaf = rng.normal(0.345, 0.02) if contaminated else rng.normal(0.502, 0.012)
    # Male X heterozygosity is near zero; a swap is what makes it disagree.
    x_het = rng.normal(0.03, 0.012) if sex == "M" else rng.normal(0.61, 0.06)
    return {
        "qc_n_called": n_called,
        "qc_ti_tv": round(_clamp(ti_tv, 0.4, 4.0), 4),
        "qc_het_hom": round(_clamp(het_hom, 0.3, 8.0), 4),
        "qc_mean_het_vaf": round(_clamp(het_vaf, 0.15, 0.85), 4),
        "qc_x_het_rate": round(_clamp(x_het, 0.0, 1.0), 4),
    }


def _substitution(rng: "Rng"):
    """ref, alt with a realistic transition/transversion balance."""
    if rng.chance(TARGET_TI_TV / (1.0 + TARGET_TI_TV)):
        return rng.pick(TRANSITIONS)
    return rng.pick(TRANSVERSIONS)


class Rng:
    def __init__(self, seed: int):
        self.r = random.Random(seed)

    def pick(self, seq: Sequence):
        return seq[self.r.randrange(len(seq))]

    def weighted(self, weights: Dict[Any, float]):
        total = sum(weights.values())
        x = self.r.random() * total
        for k, w in weights.items():
            x -= w
            if x <= 0:
                return k
        return list(weights)[-1]

    def normal(self, mu: float, sd: float) -> float:
        return self.r.gauss(mu, sd)

    def chance(self, p: float) -> bool:
        return self.r.random() < p


def _clamp(v, lo, hi):
    return max(lo, min(hi, v))


# ------------------------------------------------------------ variant catalog --
def build_catalog(rng: Rng) -> Dict[str, List[Dict[str, Any]]]:
    """A recurrent-variant catalogue per gene.

    Real germline cohorts are dominated by a limited set of recurrent alleles,
    a few of them founder variants. Sampling every observation independently
    from a continuous space would make §G07 (population frequency, founder
    detection) vacuous — nothing would ever recur.
    """
    catalog: Dict[str, List[Dict[str, Any]]] = {}
    high_diversity = {"HBB", "CFTR", "GJB2", "CYP21A2"}

    for gene in GENE_LIST:
        rec = GENE_DISEASE[gene]
        n = rng.r.randrange(12, 20) if gene in high_diversity else rng.r.randrange(3, 10)
        # Moderate/Limited-validity genes are guaranteed a pathogenic allele.
        # This is not a convenience: a gene like CHEK2 is in the reference
        # precisely BECAUSE it has well-established pathogenic variants — what
        # is weak is the gene-disease relationship, not the variant call. If
        # the generator leaves them all benign, §G05's yield-inflation figure
        # is structurally zero and the screen's whole purpose is untestable.
        force_plp = (2 if gene in high_diversity
                     else 1 if rec.validity in ("Moderate", "Limited") else 0)

        # How many of this gene's alleles are pathogenic is FIXED, not drawn.
        # Letting it emerge from per-variant coin flips gives a catalogue of 3-10
        # variants a huge variance in P/LP fraction, and that variance lands
        # directly on the per-test-code yield: a 5-gene panel whose catalogues
        # happened to come up benign reports a 4% yield, one that came up
        # pathogenic reports 70%. Neither is informative about the tool.
        target_plp = max(force_plp, round(n * CATALOG_PLP_FRACTION))
        plp_slots = set(range(target_plp))
        variants: List[Dict[str, Any]] = []

        for i in range(n):
            var_class = rng.weighted({"SNV": 74, "Indel": 15, "CNV": 5, "Splice": 6})
            if var_class == "CNV":
                conseq = rng.pick(["Exon deletion", "Whole-gene duplication"])
            elif var_class == "Indel":
                conseq = rng.pick(["Frameshift", "In-frame indel"])
            elif var_class == "Splice":
                conseq = "Splice site"
            else:
                conseq = rng.pick(["Missense"] * 5 + ["Nonsense"])
            truncating = conseq in TRUNCATING

            plen = rec.protein_len or 800
            aa = 1 + rng.r.randrange(min(plen, 4000))

            if conseq == "Frameshift":
                hgvs_p = "p.{}{}fs*{}".format(rng.pick(AAS), aa, 3 + rng.r.randrange(18))
            elif conseq == "Nonsense":
                hgvs_p = "p.{}{}*".format(rng.pick(AAS), aa)
            elif conseq == "Splice site":
                hgvs_p = None
            elif var_class == "CNV":
                hgvs_p = conseq
            else:
                hgvs_p = "p.{}{}{}".format(rng.pick(AAS), aa, rng.pick(AAS))

            hgvs_c = ("c.{}{}".format(aa * 3, rng.pick(["+1G>A", "-2A>G"]))
                      if conseq == "Splice site"
                      else "c.{}{}>{}".format(aa * 3, rng.pick("ACGT"), rng.pick("ACGT")))

            if i in plp_slots:
                # Truncating changes in a Definitive gene read as Pathogenic;
                # elsewhere Likely pathogenic is the commoner call.
                if truncating and rec.validity == "Definitive":
                    acmg = "Pathogenic" if rng.chance(0.75) else "Likely pathogenic"
                else:
                    acmg = "Pathogenic" if rng.chance(0.4) else "Likely pathogenic"
            else:
                acmg = rng.weighted({"Uncertain significance": 70,
                                     "Likely benign": 22, "Benign": 8})
            is_plp = acmg in ("Pathogenic", "Likely pathogenic")

            # P/LP alleles are rare in population databases; benign ones are not.
            gnomad = round(10 ** -(3.2 + rng.r.random() * 2.6), 8) if is_plp \
                else round(10 ** -(1.6 + rng.r.random() * 2.4), 8)
            founder = is_plp and rng.chance(0.22)

            variants.append({
                "gene": gene,
                "variant_key": "{}:{}:{}:{}".format(
                    rec.chrom, 1_000_000 + rng.r.randrange(200_000_000),
                    *_substitution(rng)),
                "chrom": rec.chrom,
                "pos": 1_000_000 + rng.r.randrange(200_000_000),
                "hgvs_p": hgvs_p, "hgvs_c": hgvs_c, "aa_pos": aa,
                "protein_len": rec.protein_len,
                "var_class": var_class, "consequence": conseq,
                "classification": acmg, "plp": is_plp,
                "gnomad_af": gnomad,
                "gnomad_sas_af": round(gnomad * (0.6 + rng.r.random() * 1.8), 8),
                "ga100k_sas_af": round(gnomad * (0.5 + rng.r.random() * 2.6), 8),
                "founder": founder,
                "founder_pop": rng.pick(["South Asian", "Middle Eastern",
                                         "South Asian", "East Asian"]) if founder else None,
                "clinvar_sig": (rng.pick(["Pathogenic", "Likely pathogenic",
                                          "Pathogenic/Likely pathogenic"]) if is_plp
                                else rng.pick(["Uncertain significance",
                                               "Conflicting interpretations",
                                               "Not reported", "Benign"])),
                "clinvar_id": "VCV{:09d}".format(rng.r.randrange(999_999_999)),
                "acmg_codes": ", ".join(
                    ["PVS1", "PM2", "PP5"] if acmg == "Pathogenic"
                    else ["PM2", "PP3", "PM1"] if acmg == "Likely pathogenic"
                    else ["PM2", "BP4"]),
                # Volume of new public evidence since last interpretation —
                # the reclassification pressure signal for §G08.
                "evidence_delta": (round(rng.r.random() * rng.r.random() * 100)
                                   if acmg == "Uncertain significance" else 0),
            })
        catalog[gene] = variants
    return catalog


# ------------------------------------------------------------------ generator --
def generate(n_families: int = 640, seed: int = 20260813) -> Dict[str, List[Dict[str, Any]]]:
    from ..reference.genes import SF_GENES
    from ..reference.tests import INDICATIONS, relevant_genes

    # "Secondary finding" is defined RELATIVE TO THE INDICATION, not absolutely.
    # MYH7 on a cardiomyopathy referral is the target gene; MYH7 on an
    # undiagnosed-exome referral is an incidental. Only the latter gets
    # downweighted — treating the SF list as globally rare instead starves
    # cardiac and hereditary-cancer panels of findings in their own core genes.
    incidental_genes = {
        ind: set(SF_GENES) - set(relevant_genes(ind)) for ind in INDICATIONS
    }

    def gene_weights(test):
        skip = incidental_genes.get(test.indication, set())
        return {g: SF_GENE_DRAW_WEIGHT if g in skip else 1.0 for g in test.genes}

    rng = Rng(seed)
    catalog = build_catalog(rng)

    bundle: Dict[str, List[Dict[str, Any]]] = {
        "family": [], "subject": [], "sample": [], "run": [],
        "finding": [], "interpretation": [], "run_fingerprint": [],
    }
    sj_n = run_n = find_n = 0

    for f in range(1, n_families + 1):
        family_id = "FAM-{:04d}".format(f)
        ancestry = rng.weighted(ANCESTRY_WEIGHTS)
        test_code = rng.weighted(TEST_WEIGHTS)
        test = TEST_CODES[test_code]
        size = rng.weighted(FAMILY_SIZE_WEIGHTS)
        consent = rng.weighted(CONSENT_WEIGHTS)

        # A segregating family variant, so related members share a finding and
        # probands-only de-duplication has a real inflation to suppress.
        family_variant = None
        if rng.chance(0.42):
            # Same SF downweighting as the per-subject draw below. Missing it
            # here leaks the whole effect: ~40% of families carry a segregating
            # variant and the proband inherits it 94% of the time, so an
            # unweighted draw dominates the secondary-findings rate on its own.
            gene = rng.weighted(gene_weights(test))
            if catalog.get(gene):
                family_variant = rng.pick(catalog[gene])

        bundle["family"].append({
            "family_id": family_id, "tenant_id": config.TENANT_ID,
            "ancestry": ancestry, "consent_class": consent,
        })

        for m in range(size):
            is_proband = m == 0
            if is_proband:
                relation = "Proband"
            elif size == 3:
                relation = "Mother" if m == 1 else "Father"
            else:
                relation = rng.pick(["Mother", "Father", "Sibling", "Sibling", "Child"])

            if is_proband:
                mu, sd = ((6, 5) if test.indication == "Neurodevelopmental delay"
                          else (31, 6) if test.indication == "Reproductive carrier screening"
                          else (43, 15))
            else:
                mu, sd = (56, 12) if relation in ("Mother", "Father") else (28, 14)
            age = int(_clamp(round(rng.normal(mu, sd)), 0, 90))
            sex = ("F" if relation == "Mother" else "M" if relation == "Father"
                   else ("F" if rng.chance(0.5) else "M"))

            days_ago = int(rng.r.random() * rng.r.random() * 700)
            collected = TODAY - timedelta(days=days_ago)

            sj_n += 1
            subject_id = "SJ-{:05d}".format(sj_n)
            bundle["subject"].append({
                "subject_id": subject_id, "tenant_id": config.TENANT_ID,
                "mrn": "MG{}".format(510000 + rng.r.randrange(89999)),
                "sex": sex, "age": age, "age_bucket": age_bucket(age),
                "ancestry": ancestry, "consent_class": consent,
                "family_id": family_id, "relation": relation, "is_proband": is_proband,
                "indication": test.indication,
                "affected_status": ("Affected" if is_proband and rng.chance(0.82)
                                    else "At risk / unaffected" if is_proband
                                    else rng.weighted({"Unaffected": 64, "Affected": 21,
                                                       "At risk / unaffected": 15})),
                "family_history": "Positive" if rng.chance(0.42) else "Negative / unknown",
                "referral_source": rng.pick(REFERRAL_SOURCES),
                "phenotype_hpo": rng.pick(HPO_BY_INDICATION.get(
                    test.indication, ["HP:0000001 All"])),
            })

            run_n += 1
            sample_id = "SP-{:05d}".format(run_n)
            run_id = "GR-{:05d}".format(run_n)
            is_wgs = test_code == "MG-WGS"
            mean_depth = int(_clamp(rng.normal(38 if is_wgs else 180,
                                               8 if is_wgs else 60), 18, 420))
            # A slice of runs has no test code in the LIMS export, so its
            # coverage must be inferred from the fingerprint (§3.4).
            labelled = not rng.chance(UNLABELLED_RUN_FRACTION)

            # A few samples are deliberately poor so the QC screen has real
            # findings. Contamination adds a second individual's alleles, which
            # surface as spurious heterozygotes with off-centre allele balance
            # and a depressed Ti/Tv; low-coverage samples simply lose depth.
            contaminated = rng.chance(CONTAMINATED_FRACTION)
            low_depth = (not contaminated) and rng.chance(LOW_DEPTH_FRACTION)
            if low_depth:
                mean_depth = int(mean_depth * 0.28)

            bundle["sample"].append({
                "sample_id": sample_id, "subject_id": subject_id,
                "sample_type": rng.weighted(SAMPLE_TYPE_WEIGHTS),
                "collection_date": collected,
            })
            bundle["run"].append({
                "run_id": run_id, "sample_id": sample_id, "subject_id": subject_id,
                "test_code": test_code if labelled else None,
                "assay_version": test.assay_version if labelled else None,
                "pipeline_version": rng.pick(PIPELINE_VERSIONS),
                "reference_build": "GRCh38" if rng.chance(0.84) else "GRCh37",
                "caller": "DRAGEN-4.2" if rng.chance(0.8) else "GATK-4.5",
                "implied_panel_id": None, "mean_depth": mean_depth,
                "pct_bases_20x": round(
                    _clamp(rng.normal(70.0, 4.0) if low_depth
                           else rng.normal(98.4, 1.9), 45, 100), 1),
                "qc_status": ("Pass" if rng.chance(0.965)
                              else rng.pick(["Low coverage", "Contamination flag"])),
                "collection_date": collected,
                "source_file": "synthetic",
                # What a pipeline would report over the WHOLE call set for this
                # sample — thousands of variants, of which only the reviewable
                # handful is retained. Simulated here rather than derived from
                # `finding`, because deriving it from a few retained rows is
                # exactly the mistake this column exists to prevent.
                **_simulated_run_qc(rng, contaminated, low_depth, sex),
            })

            # Fingerprint: genes the assay looked at that produced any on-target
            # call. Presence scales with gene length — this is precisely why the
            # §3.4 threshold has to be per-gene.
            for gene in test.genes:
                plen = GENE_DISEASE[gene].protein_len or 500
                p = _clamp(0.18 + 0.42 * math.log10(max(plen, 10)) / 3.0, 0.12, 0.97)
                if rng.chance(p):
                    bundle["run_fingerprint"].append(
                        {"run_id": run_id, "gene_symbol": gene})

            # ------------------------------------------------ observations --
            chosen: List[Dict[str, Any]] = []
            seen = set()

            def add(cv: Dict[str, Any]) -> None:
                if cv["variant_key"] in seen:
                    return
                seen.add(cv["variant_key"])
                chosen.append(cv)

            if family_variant and rng.chance(0.94 if is_proband else 0.52):
                add(family_variant)

            extra = (rng.r.randrange(2, 6) if test_code in ("MG-WGS", "MG-WES", "MG-CES")
                     else rng.r.randrange(1, 4))
            for _ in range(extra):
                if test_code == "MG-HBOC-27" and rng.chance(0.46):
                    gene = rng.pick(["BRCA1", "BRCA2", "BRCA1", "BRCA2",
                                     "PALB2", "ATM", "CHEK2"])
                elif test_code == "MG-LYNCH-7" and rng.chance(0.55):
                    gene = rng.pick(["MLH1", "MSH2", "MSH6", "PMS2"])
                elif test_code == "MG-THAL-1":
                    if not rng.chance(0.62):
                        continue
                    gene = "HBB"
                else:
                    # Broad tests draw widely across their scope — that spread
                    # is what makes exome yield LOWER than a targeted panel once
                    # the phenotype-relevance gate is applied, which is the
                    # calibration the spec expects (23-36% vs 40-60%).
                    #
                    # But SF genes must be downweighted hard. This reference
                    # table holds 43 genes of which 17 are on the ACMG SF list
                    # (40%); a real exome covers ~20,000 genes of which 81 are
                    # SF (0.4%). Drawing uniformly would inflate the incidental
                    # rate by two orders of magnitude and make §G09 report ~70%
                    # secondary findings instead of the few percent a lab sees.
                    gene = rng.weighted(gene_weights(test))
                if catalog.get(gene):
                    add(rng.pick(catalog[gene]))

            # Second P/LP hit in an AR gene, so compound-het derivation has
            # real pairs to find. Note we do NOT write 'Compound heterozygous'
            # here — the derivation pass computes it (§3.5).
            ar_plp_genes = {cv["gene"] for cv in chosen
                            if cv["plp"] and GENE_DISEASE[cv["gene"]].inheritance == "AR"}
            for gene in ar_plp_genes:
                if rng.chance(0.18):
                    alts = [v for v in catalog[gene]
                            if v["plp"] and v["variant_key"] not in seen]
                    if alts:
                        add(rng.pick(alts))

            if contaminated:
                # Extra calls drawn with a UNIFORM substitution spectrum — the
                # point is that they are not real variants, so they carry none
                # of sequence's transition bias and drag Ti/Tv down.
                for _ in range(rng.r.randrange(14, 30)):
                    g = rng.pick(test.genes)
                    if not catalog.get(g):
                        continue
                    noise = dict(rng.pick(catalog[g]))
                    ref, alt = rng.pick("ACGT"), rng.pick("ACGT")
                    if ref == alt:
                        alt = "ACGT".replace(ref, "")[0]
                    noise["variant_key"] = "{}:{}:{}:{}".format(
                        noise["chrom"], 1_000_000 + rng.r.randrange(200_000_000),
                        ref, alt)
                    noise["classification"] = "Uncertain significance"
                    noise["plp"] = False
                    add(noise)

            for cv in chosen:
                inh = GENE_DISEASE[cv["gene"]].inheritance
                sex_linked = cv["chrom"] in ("X", "Y")

                # zygosity_raw: what the FILE says. Het/Hom only (§3.3).
                # X-linked male calls are commonly emitted as homozygous by
                # callers — the derivation pass is what turns them hemizygous.
                if contaminated:
                    # A second individual's alleles read as heterozygous and sit
                    # well off the 0.5 a true het produces.
                    zyg_raw = "Heterozygous" if rng.chance(0.93) else "Homozygous"
                elif sex_linked and sex == "M":
                    zyg_raw = "Homozygous" if rng.chance(0.7) else "Heterozygous"
                elif inh == "AR":
                    zyg_raw = rng.weighted({"Heterozygous": 90, "Homozygous": 10})
                else:
                    zyg_raw = rng.weighted({"Heterozygous": 94, "Homozygous": 6})

                vaf = (round(_clamp(rng.normal(0.98, 0.02), 0.85, 1.0), 3)
                       if zyg_raw == "Homozygous"
                       else round(_clamp(rng.normal(
                           0.34 if contaminated else 0.5, 0.05), 0.15, 0.68), 3))

                find_n += 1
                finding_id = "GO-{:07d}".format(find_n)
                bundle["finding"].append({
                    "finding_id": finding_id, "run_id": run_id,
                    "subject_id": subject_id, "family_id": family_id,
                    "gene_id": GENE_DISEASE[cv["gene"]].gene_id,
                    "gene_symbol": cv["gene"], "variant_key": cv["variant_key"],
                    "chrom": cv["chrom"], "pos": cv["pos"],
                    "ref_allele": cv["variant_key"].split(":")[2],
                    "alt_allele": cv["variant_key"].split(":")[3],
                    "hgvs_c": cv["hgvs_c"], "hgvs_p": cv["hgvs_p"],
                    "aa_pos": cv["aa_pos"], "protein_len": cv["protein_len"],
                    "var_class": cv["var_class"], "consequence": cv["consequence"],
                    "zygosity_raw": zyg_raw, "zygosity": zyg_raw,
                    "vaf": vaf,
                    "depth": int(_clamp(rng.normal(mean_depth, mean_depth * 0.3), 20, 900)),
                    "alt_depth": None, "filter_status": "PASS",
                    "gnomad_af": cv["gnomad_af"],
                    "gnomad_sas_af": cv["gnomad_sas_af"],
                    "ga100k_sas_af": cv["ga100k_sas_af"],
                    "clinvar_sig": cv["clinvar_sig"], "clinvar_id": cv["clinvar_id"],
                    "mane": True, "reviewable": True,
                })
                bundle["interpretation"].append({
                    "finding_id": finding_id, "framework": "ACMG/AMP",
                    "classification": cv["classification"],
                    "acmg_codes": cv["acmg_codes"],
                    "kb_snapshot_id": config.KB_SNAPSHOT_ID,
                    "curated_flag": rng.chance(0.56),
                    "interpreted_at": TODAY - timedelta(days=rng.r.randrange(430)),
                    # Where the registry does not declare full reportable scope,
                    # reportability itself is uncertain (§G02 provisional).
                    "reportable": True if test.scope_complete else rng.chance(0.7),
                    "inherited_from": (rng.pick(["Maternal", "Paternal", "De novo", "Unknown"])
                                       if is_proband and rng.chance(0.5) else "Unknown"),
                    "segregation": (rng.pick(["Segregates with phenotype",
                                              "Does not segregate", "Incomplete data"])
                                    if size > 1 else "Not assessed"),
                    "evidence_delta": cv["evidence_delta"],
                })

    return bundle


def rebuild(n_families: int = 640, seed: int = 20260813) -> Dict[str, Any]:
    from . import store

    store.init(fresh=True)
    store.load_reference()
    bundle = generate(n_families=n_families, seed=seed)
    store.persist(bundle)
    derived = store.rebuild_derivations()

    from .. import db
    db.meta_set("source", "synthetic")
    db.meta_set("seed", str(seed))
    db.meta_set("n_families", str(n_families))
    return {
        "families": len(bundle["family"]),
        "subjects": len(bundle["subject"]),
        "runs": len(bundle["run"]),
        "findings": len(bundle["finding"]),
        "derived": derived,
    }
