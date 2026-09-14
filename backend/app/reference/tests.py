"""Germline test-code registry and the phenotype-relevance map.

The test-code registry is what makes a *gene-specific* denominator possible
(spec §3.7, C01): a cohort has one subject count but as many denominators as it
has genes, because each test code assays and reports a different gene list.

Two registry entries deliberately carry `scope_complete = False`. Real lab
registries have gaps; the tool must count those runs in numerators but flag
them as provisional in denominators (spec §G02 provisional warning) rather than
silently dropping or silently including them.
"""
from __future__ import annotations

from typing import Dict, List, NamedTuple, Tuple

from .genes import GENE_LIST, GENE_SETS


class TestCode(NamedTuple):
    code: str
    name: str
    assay_version: str
    indication: str
    sf_capable: bool       # scope can report ACMG SF v3.3 incidentals
    scope_complete: bool   # registry declares the full reportable gene list
    accredited: bool
    genes: Tuple[str, ...]


def _t(code, name, ver, indication, sf, complete, accredited, genes) -> TestCode:
    return TestCode(code, name, ver, indication, sf, complete, accredited, tuple(sorted(set(genes))))


TEST_CODES: Dict[str, TestCode] = {t.code: t for t in [
    _t("MG-HBOC-27", "Hereditary breast & ovarian cancer panel", "v3.2",
       "Hereditary cancer risk", False, True, True,
       GENE_SETS["HBOC (breast/ovarian)"] + ["TP53", "PTEN", "STK11"]),
    _t("MG-LYNCH-7", "Lynch syndrome panel", "v2.1",
       "Hereditary cancer risk", False, True, True,
       GENE_SETS["Lynch / MMR"] + ["APC"]),
    _t("MG-CARRIER-113", "Expanded carrier screen", "v4.0",
       "Reproductive carrier screening", False, True, True,
       GENE_SETS["Expanded carrier screen"] + ["ALDOB", "ATP7B", "USH2A"]),
    _t("MG-THAL-1", "Beta thalassemia (HBB) targeted", "v1.4",
       "Haemoglobinopathy", False, True, True, ["HBB"]),
    _t("MG-CARDIO-42", "Cardiomyopathy & arrhythmia panel", "v2.0",
       "Cardiac", False, True, True, GENE_SETS["Cardiomyopathy / arrhythmia"]),
    # registry gap: reportable scope not declared
    _t("MG-NDD-58", "Neurodevelopmental disorder panel", "v3.0",
       "Neurodevelopmental delay", False, False, False,
       GENE_SETS["Neurodevelopmental"] + ["NF1", "FMR1"]),
    _t("MG-CES", "Clinical exome (CES)", "v5.1",
       "Undiagnosed / multisystem", True, True, True, GENE_LIST),
    _t("MG-WES", "Whole exome (WES)", "v6.0",
       "Undiagnosed / multisystem", True, True, True, GENE_LIST),
    _t("MG-WGS", "Whole genome (WGS)", "v2.2",
       "Undiagnosed / multisystem", True, True, True, GENE_LIST),
    # registry gap: reportable scope not declared
    _t("MG-HEARING-96", "Hearing loss panel", "v2.3",
       "Sensorineural hearing loss", False, False, True,
       GENE_SETS["Hearing loss"] + ["GJB2", "SLC26A4"]),
]}

TEST_CODE_LIST: List[str] = list(TEST_CODES)
INDICATIONS: List[str] = sorted({t.indication for t in TEST_CODES.values()})

ASSAY_VERSIONS: List[str] = sorted({t.assay_version for t in TEST_CODES.values()})
PIPELINE_VERSIONS = ["LINC-v2.8.3", "LINC-v2.9.0", "LINC-v2.9.1"]
REFERENCE_BUILDS = ["GRCh38", "GRCh37"]
SAMPLE_TYPES = ["Blood (EDTA)", "Saliva", "Buccal swab", "Dried blood spot"]
CALLERS = ["DRAGEN-4.2", "GATK-4.5"]

# ------------------------------------------------------------- demographics --
ANCESTRY_WEIGHTS = {
    "South Asian": 62, "Middle Eastern": 9, "East Asian": 7,
    "European": 11, "African": 4, "Admixed / other": 7,
}
ANCESTRIES: List[str] = list(ANCESTRY_WEIGHTS)

CONSENT_CLASSES = [
    "Full research + secondary findings",
    "Clinical only, secondary findings declined",
    "Clinical only, research declined",
    "Withdrawn",
]
# Only this class permits ACMG SF analysis (spec §G09).
CONSENT_SF_OK = {"Full research + secondary findings"}
# Consent withdrawn is excluded from every cohort, unconditionally (§G01 step 2).
CONSENT_WITHDRAWN = "Withdrawn"

AGE_BUCKETS = ["<1y", "1–11y", "12–17y", "18–39y", "40–59y", "60+"]
AFFECTED_STATUSES = ["Affected", "At risk / unaffected", "Unaffected"]
FAMILY_HISTORY = ["Positive", "Negative / unknown"]
RELATIONS = ["Proband", "Mother", "Father", "Sibling", "Child"]
REFERRAL_SOURCES = [
    "Medical genetics", "Oncology", "Obstetrics",
    "Paediatrics", "Cardiology", "Self-referred",
]

# --------------------------------------------------------- phenotype (HPO) --
HPO_BY_INDICATION: Dict[str, List[str]] = {
    "Hereditary cancer risk": [
        "HP:0100013 Breast neoplasm", "HP:0100615 Ovarian neoplasm",
        "HP:0100273 Colon neoplasm", "HP:0003002 Breast carcinoma"],
    "Reproductive carrier screening": [
        "HP:0000006 Autosomal dominant inheritance",
        "HP:0032322 Preconception screening", "HP:0001939 Metabolic abnormality"],
    "Haemoglobinopathy": [
        "HP:0001903 Anemia", "HP:0011902 Abnormal hemoglobin", "HP:0001744 Splenomegaly"],
    "Cardiac": [
        "HP:0001639 Hypertrophic cardiomyopathy", "HP:0001279 Syncope",
        "HP:0001645 Sudden cardiac death"],
    "Neurodevelopmental delay": [
        "HP:0001263 Global developmental delay", "HP:0001250 Seizure",
        "HP:0000717 Autism", "HP:0001249 Intellectual disability"],
    "Undiagnosed / multisystem": [
        "HP:0001939 Metabolic abnormality", "HP:0001627 Abnormal heart morphology",
        "HP:0000365 Hearing impairment", "HP:0001263 Global developmental delay"],
    "Sensorineural hearing loss": [
        "HP:0000407 Sensorineural hearing impairment",
        "HP:0008527 Congenital sensorineural hearing impairment"],
}

# ------------------------------------------------ phenotype-relevance gate ---
# Spec §G02: "Without relevance gating, a broad exome counts *any* P/LP as
# diagnostic — which pushed yield to 81% in testing. A P/LP finding in an
# unrelated organ system is a secondary finding, not a diagnosis."
INDICATION_RELEVANT_SETS: Dict[str, Tuple[str, ...]] = {
    "Hereditary cancer risk": ("HBOC", "Lynch", "GI", "SF"),
    "Reproductive carrier screening": ("Carrier",),
    "Haemoglobinopathy": ("Heme", "Carrier"),
    "Cardiac": ("Cardio",),
    "Neurodevelopmental delay": ("Neuro",),
    "Sensorineural hearing loss": ("Hearing",),
    "Undiagnosed / multisystem": ("Neuro", "Metabolic", "Endocrine"),
}


def relevant_genes(indication: str) -> List[str]:
    """Genes considered phenotype-relevant for an indication (yield gate)."""
    from .genes import GENE_DISEASE
    keys = set(INDICATION_RELEVANT_SETS.get(indication, ()))
    if not keys:
        return []
    return [g for g, r in GENE_DISEASE.items() if keys & set(r.sets)]
