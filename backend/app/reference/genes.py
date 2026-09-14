"""Gene-disease validity reference (spec §G05).

    "Not available in VariMAT — must be curated, sourced from ClinGen where
     available."

SEED TABLE — REPLACE BEFORE PRODUCTION.
This is a hand-entered starter set covering the 44 genes the prototype models.
Validity and penetrance calls must be replaced with a dated export from the
ClinGen Gene-Disease Validity and Gene Curation Coalition resources, and
`gene_id` with an HGNC export, before any figure leaves the building. The
loader reads whatever is in this table; swapping the source is a one-file
change. `kb_snapshot_id` in the cohort manifest is what pins which version
produced a given number (spec C02).

Fields
    condition/mondo   the disease entity the gene is curated against
    inheritance       AD | AR | AR/AD | XLR | XLD  -- drives §G04 zygosity logic
    validity          Definitive | Strong | Moderate | Limited
                      only Definitive/Strong count toward diagnostic yield (§G02)
    penetrance        High | Moderate | Low
    sets              curated gene-set membership (§G01 "curated gene sets")
    chrom             required to DERIVE hemizygosity (§3.5) -- X/Y + subject sex
    protein_len       for the §G10 lollipop; None hides the plot, keeps the panel
"""
from __future__ import annotations

from typing import Dict, List, NamedTuple, Optional, Tuple


class GeneRecord(NamedTuple):
    symbol: str
    gene_id: str
    condition: str
    mondo: str
    inheritance: str
    validity: str
    penetrance: str
    sets: Tuple[str, ...]
    chrom: str
    protein_len: Optional[int]
    domains: Tuple[Tuple[int, int, str], ...] = ()


def _g(*args, **kw) -> GeneRecord:
    return GeneRecord(*args, **kw)


GENE_DISEASE: Dict[str, GeneRecord] = {r.symbol: r for r in [
    # ---------------------------------------------------- hereditary cancer --
    _g("BRCA1", "HGNC:1100", "Hereditary breast and ovarian cancer", "MONDO:0003582",
       "AD", "Definitive", "High", ("HBOC", "SF"), "17", 1863,
       ((24, 64, "RING"), (1650, 1863, "BRCT"))),
    _g("BRCA2", "HGNC:1101", "Hereditary breast and ovarian cancer", "MONDO:0003582",
       "AD", "Definitive", "High", ("HBOC", "SF"), "13", 3418,
       ((1002, 2085, "BRC repeats"), (2481, 3186, "Helical/OB"))),
    _g("PALB2", "HGNC:26144", "Familial breast cancer", "MONDO:0007254",
       "AD", "Definitive", "Moderate", ("HBOC",), "16", 1186, ((853, 1186, "WD40"),)),
    _g("ATM", "HGNC:795", "Ataxia-telangiectasia / breast cancer susceptibility",
       "MONDO:0008842", "AR/AD", "Definitive", "Moderate", ("HBOC",), "11", 3056,
       ((2712, 2962, "PI3K/PI4K"),)),
    _g("CHEK2", "HGNC:16627", "Familial breast cancer", "MONDO:0007254",
       "AD", "Moderate", "Low", ("HBOC",), "22", 543, ((220, 486, "Kinase"),)),
    _g("RAD51C", "HGNC:9820", "Breast-ovarian cancer susceptibility", "MONDO:0013253",
       "AD", "Strong", "Moderate", ("HBOC",), "17", 376, ((80, 300, "RecA/RAD51"),)),
    _g("BRIP1", "HGNC:20473", "Ovarian cancer susceptibility", "MONDO:0013254",
       "AD", "Moderate", "Low", ("HBOC",), "17", 1249, ((1, 600, "Helicase"),)),
    _g("MLH1", "HGNC:7127", "Lynch syndrome", "MONDO:0005835",
       "AD", "Definitive", "High", ("Lynch", "SF"), "3", 756, ((20, 320, "ATPase"),)),
    _g("MSH2", "HGNC:7325", "Lynch syndrome", "MONDO:0005835",
       "AD", "Definitive", "High", ("Lynch", "SF"), "2", 934, ((1, 124, "MutS I"),)),
    _g("MSH6", "HGNC:7329", "Lynch syndrome", "MONDO:0005835",
       "AD", "Definitive", "Moderate", ("Lynch", "SF"), "2", 1360, ((1080, 1230, "MutS V"),)),
    _g("PMS2", "HGNC:9122", "Lynch syndrome", "MONDO:0005835",
       "AD", "Definitive", "Moderate", ("Lynch", "SF"), "7", 862, ((1, 200, "MutS"),)),
    _g("APC", "HGNC:583", "Familial adenomatous polyposis", "MONDO:0007383",
       "AD", "Definitive", "High", ("GI", "SF"), "5", 2843, ((1020, 1169, "Beta-cat bind"),)),
    _g("TP53", "HGNC:11998", "Li-Fraumeni syndrome", "MONDO:0018875",
       "AD", "Definitive", "High", ("SF",), "17", 393, ((95, 289, "DNA-binding"),)),
    _g("PTEN", "HGNC:9588", "PTEN hamartoma tumour syndrome", "MONDO:0017623",
       "AD", "Definitive", "High", ("SF",), "10", 403, ((14, 185, "Phosphatase"),)),
    _g("STK11", "HGNC:11389", "Peutz-Jeghers syndrome", "MONDO:0008280",
       "AD", "Definitive", "High", ("GI", "SF"), "19", 433, ((49, 309, "Kinase"),)),
    _g("RET", "HGNC:9967", "Multiple endocrine neoplasia type 2", "MONDO:0017411",
       "AD", "Definitive", "High", ("SF",), "10", 1114, ((724, 1016, "Kinase"),)),
    _g("VHL", "HGNC:12687", "von Hippel-Lindau disease", "MONDO:0008667",
       "AD", "Definitive", "High", ("SF",), "3", 213, ((63, 155, "Beta-domain"),)),
    _g("SDHB", "HGNC:10681", "Hereditary paraganglioma-pheochromocytoma", "MONDO:0017366",
       "AD", "Definitive", "Moderate", ("SF",), "1", 280, ((50, 250, "Fe-S"),)),

    # ------------------------------------------------------ carrier / heme --
    _g("HBB", "HGNC:4827", "Beta thalassemia / sickle cell disease", "MONDO:0009709",
       "AR", "Definitive", "High", ("Carrier", "Heme"), "11", 147, ((1, 147, "Globin"),)),
    _g("CFTR", "HGNC:1884", "Cystic fibrosis", "MONDO:0009061",
       "AR", "Definitive", "High", ("Carrier",), "7", 1480,
       ((389, 673, "NBD1"), (1210, 1443, "NBD2"))),
    _g("SMN1", "HGNC:11117", "Spinal muscular atrophy", "MONDO:0001516",
       "AR", "Definitive", "High", ("Carrier",), "5", 294, ((1, 91, "Tudor"),)),
    _g("GJB2", "HGNC:4284", "Nonsyndromic hearing loss", "MONDO:0019497",
       "AR", "Definitive", "High", ("Carrier", "Hearing"), "13", 226, ((1, 226, "Connexin"),)),
    _g("GBA", "HGNC:4177", "Gaucher disease", "MONDO:0018150",
       "AR", "Definitive", "High", ("Carrier", "Metabolic"), "1", 536, ((80, 500, "Glyco hydro 30"),)),
    _g("G6PD", "HGNC:4057", "G6PD deficiency", "MONDO:0009945",
       "XLR", "Definitive", "Moderate", ("Carrier", "Heme"), "X", 515, ((30, 500, "G6PD"),)),
    _g("F8", "HGNC:3546", "Haemophilia A", "MONDO:0010602",
       "XLR", "Definitive", "High", ("Heme",), "X", 2351,
       ((1, 336, "F5/8 A1"), (1690, 2019, "C1"))),
    _g("CYP21A2", "HGNC:2600", "Congenital adrenal hyperplasia", "MONDO:0008724",
       "AR", "Definitive", "High", ("Carrier", "Endocrine"), "6", 495, ((40, 480, "Cytochrome P450"),)),
    _g("PAH", "HGNC:8582", "Phenylketonuria", "MONDO:0009861",
       "AR", "Definitive", "High", ("Carrier", "Metabolic"), "12", 452, ((118, 424, "Biopterin-dep AAH"),)),

    # ----------------------------------------------------- cardiac / SF set --
    _g("MYBPC3", "HGNC:7551", "Hypertrophic cardiomyopathy", "MONDO:0005045",
       "AD", "Definitive", "Moderate", ("Cardio", "SF"), "11", 1274, ((1, 100, "C0 Ig"),)),
    _g("MYH7", "HGNC:7577", "Hypertrophic cardiomyopathy", "MONDO:0005045",
       "AD", "Definitive", "Moderate", ("Cardio", "SF"), "14", 1935, ((181, 937, "Myosin head"),)),
    _g("KCNQ1", "HGNC:6294", "Long QT syndrome", "MONDO:0011185",
       "AD", "Definitive", "Moderate", ("Cardio", "SF"), "11", 676, ((122, 348, "Ion transport"),)),
    _g("LDLR", "HGNC:6547", "Familial hypercholesterolaemia", "MONDO:0007750",
       "AD", "Definitive", "High", ("Cardio", "SF"), "19", 860, ((25, 313, "LDL-R class A"),)),
    _g("TTN", "HGNC:12403", "Dilated cardiomyopathy", "MONDO:0005021",
       "AD", "Moderate", "Low", ("Cardio",), "2", 34350, ()),

    # ----------------------------------------------------- neurodevelopment --
    _g("SCN1A", "HGNC:10585", "Developmental and epileptic encephalopathy", "MONDO:0018614",
       "AD", "Definitive", "High", ("Neuro",), "2", 2009, ((130, 1800, "Ion transport"),)),
    _g("MECP2", "HGNC:6990", "Rett syndrome", "MONDO:0010726",
       "XLD", "Definitive", "High", ("Neuro",), "X", 486, ((90, 162, "MBD"), (302, 306, "TRD"))),
    _g("FMR1", "HGNC:3775", "Fragile X syndrome", "MONDO:0010383",
       "XLD", "Definitive", "High", ("Neuro", "Carrier"), "X", 632, ((216, 280, "KH1"),)),
    _g("NF1", "HGNC:7765", "Neurofibromatosis type 1", "MONDO:0018975",
       "AD", "Definitive", "High", ("Neuro",), "17", 2839, ((1198, 1530, "GAP"),)),
    _g("DMD", "HGNC:2928", "Duchenne/Becker muscular dystrophy", "MONDO:0010679",
       "XLR", "Definitive", "High", ("Neuro",), "X", 3685,
       ((14, 240, "CH domains"), (3080, 3360, "Cys-rich"))),

    # ---------------------------------------------------------- metabolic ---
    _g("ALDOB", "HGNC:417", "Hereditary fructose intolerance", "MONDO:0009535",
       "AR", "Definitive", "High", ("Metabolic",), "9", 364, ((5, 360, "Aldolase"),)),
    _g("ATP7B", "HGNC:870", "Wilson disease", "MONDO:0010200",
       "AR", "Definitive", "High", ("Metabolic",), "13", 1465, ((980, 1230, "E1-E2 ATPase"),)),
    _g("POLG", "HGNC:9179", "Mitochondrial DNA depletion syndrome", "MONDO:0018158",
       "AR", "Definitive", "Moderate", ("Metabolic",), "15", 1239, ((440, 1230, "DNA pol A"),)),
    _g("DUOX2", "HGNC:13273", "Congenital hypothyroidism", "MONDO:0018612",
       "AR", "Strong", "Moderate", ("Endocrine",), "15", 1548, ((1060, 1500, "Ferric reductase"),)),

    # ------------------------------------------------------------ hearing ---
    _g("USH2A", "HGNC:12601", "Usher syndrome type 2", "MONDO:0010168",
       "AR", "Definitive", "High", ("Hearing",), "1", 5202, ((1, 200, "Laminin N"),)),
    _g("SLC26A4", "HGNC:8818", "Pendred syndrome", "MONDO:0010134",
       "AR", "Definitive", "High", ("Hearing", "Carrier"), "7", 780, ((90, 500, "Sulphate transp"),)),
]}

GENE_LIST: List[str] = sorted(GENE_DISEASE)

# Aliases seen in older VariMAT annotations. Spec §3.2: "fall back to GENE_NAME
# with an alias map, and log unresolved symbols".
GENE_ALIASES: Dict[str, str] = {
    "GBA1": "GBA", "BRCC2": "BRCA2", "RNF53": "BRCA1", "FANCS": "BRCA1",
    "FANCD1": "BRCA2", "FANCN": "PALB2", "FANCJ": "BRIP1", "FANCO": "RAD51C",
    "COL4A5": "COL4A5", "HBB1": "HBB", "ABCC7": "CFTR", "PDS": "SLC26A4",
    "CX26": "GJB2", "LKB1": "STK11", "MMAC1": "PTEN", "CYP21": "CYP21A2",
}

# HGNC id -> symbol, so `GENE_ID` (preferred per §3.3) resolves first.
GENE_BY_ID: Dict[str, str] = {r.gene_id: r.symbol for r in GENE_DISEASE.values()}

# ----------------------------------------------------------- curated sets ---
# Spec §G01: HBOC · Lynch/MMR · Expanded carrier screen · Cardiomyopathy/
# arrhythmia · Neurodevelopmental · Hearing loss · Inborn errors of metabolism ·
# ACMG SF v3.3.
_SET_LABELS = [
    ("HBOC (breast/ovarian)", "HBOC"),
    ("Lynch / MMR", "Lynch"),
    ("Expanded carrier screen", "Carrier"),
    ("Cardiomyopathy / arrhythmia", "Cardio"),
    ("Neurodevelopmental", "Neuro"),
    ("Hearing loss", "Hearing"),
    ("Inborn errors of metabolism", "Metabolic"),
    ("ACMG SF v3.3", "SF"),
]

GENE_SETS: Dict[str, List[str]] = {
    label: [g for g in GENE_LIST if key in GENE_DISEASE[g].sets]
    for label, key in _SET_LABELS
}

SF_GENES: List[str] = GENE_SETS["ACMG SF v3.3"]

VALIDITY_ORDER = ["Definitive", "Strong", "Moderate", "Limited"]
PENETRANCE_ORDER = ["High", "Moderate", "Low"]
INHERITANCE_ORDER = ["AD", "AR", "AR/AD", "XLR", "XLD"]

# Only these count toward diagnostic yield (spec §G02).
ESTABLISHED_VALIDITY = {"Definitive", "Strong"}

SEX_CHROMOSOMES = {"X", "Y", "chrX", "chrY"}


def resolve_symbol(gene_id: Optional[str], gene_name: Optional[str]) -> Optional[str]:
    """GENE_ID first (§3.3), then GENE_NAME, then the alias map."""
    if gene_id:
        hit = GENE_BY_ID.get(gene_id.strip())
        if hit:
            return hit
    if gene_name:
        name = gene_name.strip().upper()
        if name in GENE_DISEASE:
            return name
        alias = GENE_ALIASES.get(name)
        if alias in GENE_DISEASE:
            return alias
    return None


def is_sex_linked(symbol: str) -> bool:
    rec = GENE_DISEASE.get(symbol)
    return bool(rec and rec.chrom in ("X", "Y"))
