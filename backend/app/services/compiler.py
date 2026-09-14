"""Query compiler — spec D01 · "Ask a cohort question".

    "Convert a plain-language question into cohort criteria. Criteria only —
     never numbers, never conclusions."

The hard constraints (D01), none of which have exceptions:

    Model emits criteria only, never a number or percentage
        -> a wrong number sounds authoritative and is unverifiable
    Model never sees patient data — only field names and allowed values
        -> PHI containment
    Unmapped terms fail loudly and stop
        -> a wrong-but-plausible cohort is worse than a refusal
    Nothing executes until the user approves the compiled criteria
        -> the user owns the query
    Every compilation logged with its source question
        -> builds the evaluation corpus

This module implements the CONTRACT and the VALIDATOR, plus a deterministic
rule-based compiler that needs no model at all. The validator is the part that
matters: whatever produces the criteria — the rule compiler here, or an LLM
wired in behind `compile_question` — has to pass through `validate`, and
`validate` fails closed.

Wiring an LLM in: have it return the same {profile, criteria, target_module}
shape, pass it to `validate`, and do not touch `execute` unless validation
passed and the user approved. The model must be given `vocabulary()` and
nothing else — that function returns field names and allowed values only, no
patient data.
"""
from __future__ import annotations

import json
import re
import uuid
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from .. import config, db
from ..reference.genes import GENE_LIST, GENE_SETS
from ..reference.tests import (ANCESTRIES, CONSENT_CLASSES, INDICATIONS,
                               TEST_CODE_LIST)
from .cohort import ALLOWED_FIELDS, LIST_FIELDS, SCALAR_FIELDS, normalise
from .governance import MODULE_CLASS, MODULE_TITLES

TARGET_MODULES = [m for m in MODULE_CLASS if m.startswith("g_")]

CLASSIFICATIONS = ["Pathogenic", "Likely pathogenic", "Uncertain significance",
                   "Likely benign", "Benign"]
ZYGOSITIES = ["Heterozygous", "Homozygous", "Compound heterozygous", "Hemizygous"]
INHERITANCES = ["AD", "AR", "AR/AD", "XLR", "XLD"]
VALIDITIES = ["Definitive", "Strong", "Moderate", "Limited"]
CONSEQUENCES = ["Missense", "Nonsense", "Frameshift", "In-frame indel",
                "Splice site", "Exon deletion", "Whole-gene duplication"]
VAR_CLASSES = ["SNV", "Indel", "CNV", "Splice", "SV"]


def vocabulary() -> Dict[str, Any]:
    """The ONLY thing a model is ever shown. Field names and allowed values —
    no subjects, no findings, no counts (spec D01: PHI containment)."""
    return {
        "profile": config.PROFILE,
        "fields": {
            "indication": INDICATIONS,
            "sex": ["F", "M"],
            "age_bucket": ["<1y", "1–11y", "12–17y", "18–39y", "40–59y", "60+"],
            "ancestry": ANCESTRIES,
            "affected_status": ["Affected", "At risk / unaffected", "Unaffected"],
            "family_history": ["Positive", "Negative / unknown"],
            "relation": ["Proband", "Mother", "Father", "Sibling", "Child"],
            "referral_source": ["Medical genetics", "Oncology", "Obstetrics",
                                "Paediatrics", "Cardiology", "Self-referred"],
            "consent_class": CONSENT_CLASSES,
            "test_code": TEST_CODE_LIST,
            "assay_version": [], "pipeline_version": [], "reference_build": [],
            "sample_type": [],
            "gene": GENE_LIST,
            "var_class": VAR_CLASSES,
            "consequence": CONSEQUENCES,
            "classification": CLASSIFICATIONS,
            "zygosity": ZYGOSITIES,
            "inheritance": INHERITANCES,
            "validity": VALIDITIES,
            "penetrance": ["High", "Moderate", "Low"],
            "clinvar_sig": [],
        },
        "scalars": dict(SCALAR_FIELDS),
        "gene_sets": list(GENE_SETS),
        "target_modules": [{"module": m, "title": MODULE_TITLES[m],
                            "output_class": MODULE_CLASS[m]} for m in TARGET_MODULES],
        "rules": [
            "Emit criteria only. Never emit a number, a percentage or a conclusion.",
            "Every field name must appear in `fields` or `scalars`.",
            "Every list value must appear in that field's allowed values, "
            "where the list is non-empty.",
            "If a term in the question does not map, emit it in `unmapped` and stop.",
        ],
    }


# ------------------------------------------------------------------ validator
def validate(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Schema + vocabulary validation. Fails closed.

    "Validate emitted field names against the profile's allowed field list.
     Anything outside it is a hard failure, not a warning."
    """
    errors: List[str] = []
    unmapped: List[str] = []
    vocab = vocabulary()["fields"]

    if payload.get("profile") not in (None, config.PROFILE):
        errors.append("profile must be '{}' — this tool is germline only"
                      .format(config.PROFILE))

    criteria = payload.get("criteria")
    if not isinstance(criteria, dict):
        errors.append("criteria must be an object")
        criteria = {}

    for field, value in criteria.items():
        if field not in ALLOWED_FIELDS:
            unmapped.append(field)
            continue
        if field in LIST_FIELDS:
            values = value if isinstance(value, (list, tuple)) else [value]
            allowed = vocab.get(field) or []
            if allowed:
                for v in values:
                    if v not in allowed:
                        unmapped.append("{}={}".format(field, v))
        elif field == "gene_set":
            if value is not None and value not in GENE_SETS:
                unmapped.append("gene_set={}".format(value))

    target = payload.get("target_module")
    if target is not None and target not in MODULE_CLASS:
        unmapped.append("target_module={}".format(target))

    # A number in the payload means the model tried to answer rather than
    # compile. That is the one failure mode the spec cares most about.
    leaked = _leaked_figures(payload)
    if leaked:
        errors.append("payload contains computed figures ({}) — the compiler emits "
                      "criteria only; the database computes every number"
                      .format(", ".join(leaked)))

    ok = not errors and not unmapped
    return {
        "validation": "PASSED" if ok else "FAILED",
        "ok": ok,
        "errors": errors,
        "unmapped": sorted(set(unmapped)),
        "fields_mapped": len([f for f in criteria if f in ALLOWED_FIELDS]),
        "target_module": target,
        "output_class": MODULE_CLASS.get(target or "g_dashboard", "RESEARCH"),
    }


_FIGURE_KEYS = {"answer", "count", "n", "rate", "percent", "percentage", "pct",
                "result", "value", "yield", "total", "estimate"}


def _leaked_figures(payload: Dict[str, Any]) -> List[str]:
    return sorted(k for k in payload if k.lower() in _FIGURE_KEYS)


# ----------------------------------------------------------- rule compiler ---
# Deterministic mapping from phrase to criteria. No model required, and it is
# the baseline the 100-question acceptance gate (D01) is measured against.
_PHRASES: List[Tuple[str, Dict[str, Any], Optional[str]]] = [
    (r"\bp/?lp\b|\bpathogenic\b",
     {"classification": ["Pathogenic", "Likely pathogenic"]}, None),
    (r"\bvus\b|uncertain significance",
     {"classification": ["Uncertain significance"]}, "g_vus"),
    (r"\bhboc\b|breast and ovarian|breast/ovarian",
     {"gene_set": "HBOC (breast/ovarian)"}, "g_carrier"),
    (r"\blynch\b|\bmmr\b|mismatch repair", {"gene_set": "Lynch / MMR"}, "g_carrier"),
    (r"carrier screen", {"gene_set": "Expanded carrier screen"}, "g_carrier"),
    (r"cardiomyopathy|arrhythmia|\bcardiac\b",
     {"gene_set": "Cardiomyopathy / arrhythmia"}, "g_carrier"),
    (r"neurodevelopment|\bndd\b", {"gene_set": "Neurodevelopmental"}, "g_carrier"),
    (r"hearing loss|deafness", {"gene_set": "Hearing loss"}, "g_carrier"),
    (r"metabolic|inborn error", {"gene_set": "Inborn errors of metabolism"}, "g_carrier"),
    (r"secondary finding|\bacmg sf\b|\bsf v3", {"sf_only": True}, "g_sf"),
    (r"founder", {"founder_only": True}, "g_popfreq"),
    (r"compound het", {"zygosity": ["Compound heterozygous"]}, "g_zygosity"),
    (r"homozygous", {"zygosity": ["Homozygous"]}, "g_zygosity"),
    (r"hemizygous|x-?linked male", {"zygosity": ["Hemizygous"]}, "g_zygosity"),
    (r"biallelic", {"zygosity": ["Homozygous", "Compound heterozygous"]}, "g_zygosity"),
    (r"heterozygous", {"zygosity": ["Heterozygous"]}, "g_zygosity"),
    (r"\bexome\b|\bwes\b", {"test_code": ["MG-CES", "MG-WES"]}, None),
    (r"\bgenome\b|\bwgs\b", {"test_code": ["MG-WGS"]}, None),
    (r"\bmale\b|\bmen\b", {"sex": ["M"]}, None),
    (r"\bfemale\b|\bwomen\b", {"sex": ["F"]}, None),
    (r"\baffected\b", {"affected_status": ["Affected"]}, None),
    (r"family history", {"family_history": ["Positive"]}, None),
    (r"include relatives|all family members|family members",
     {"probands_only": False}, None),
    (r"definitive|established validity", {"established_validity_only": True}, "g_genedisease"),
    (r"diagnostic yield|\byield\b", {}, "g_dashboard"),
    (r"carrier rate|carrier frequency", {}, "g_carrier"),
    (r"gene.?disease|validity", {}, "g_genedisease"),
    (r"phenotype|\bhpo\b", {}, "g_phenotype"),
    (r"denominator|coverage|assayed", {}, "c_denominator"),
    (r"\bqc\b|quality|provenance|pipeline version", {}, "g_runs"),
]

_ANCESTRY_RE = {a: re.compile(re.escape(a.split(" /")[0]), re.I) for a in ANCESTRIES}


def compile_question(question: str) -> Dict[str, Any]:
    """Compile a question to criteria. Deterministic, auditable, no model.

    Returns the D01 payload shape: {profile, criteria, target_module}, plus the
    trace of which phrases fired, so a wrong compilation is debuggable rather
    than mysterious.
    """
    q = question.lower()
    criteria: Dict[str, Any] = {}
    target: Optional[str] = None
    trace: List[str] = []

    for pattern, crit, module in _PHRASES:
        if re.search(pattern, q):
            trace.append(pattern)
            for k, v in crit.items():
                if k in criteria and isinstance(v, list) and isinstance(criteria[k], list):
                    criteria[k] = sorted(set(criteria[k]) | set(v))
                else:
                    criteria[k] = v
            if module and not target:
                target = module

    for gene in GENE_LIST:
        if re.search(r"\b{}\b".format(re.escape(gene)), question, re.I):
            criteria.setdefault("gene", []).append(gene)
            trace.append("gene:" + gene)
    for ind in INDICATIONS:
        if ind.lower().split(" /")[0] in q:
            criteria.setdefault("indication", []).append(ind)
            trace.append("indication:" + ind)
    for anc, rx in _ANCESTRY_RE.items():
        if rx.search(question):
            criteria.setdefault("ancestry", []).append(anc)
            trace.append("ancestry:" + anc)
    for code in TEST_CODE_LIST:
        if code.lower() in q:
            criteria.setdefault("test_code", []).append(code)
            trace.append("test_code:" + code)

    m = re.search(r"(?:last|past|within)\s+(\d+)\s*(month|year)", q)
    if m:
        months = int(m.group(1)) * (12 if m.group(2) == "year" else 1)
        criteria["months"] = months
        trace.append("window:{}mo".format(months))

    return {
        "profile": config.PROFILE,
        "criteria": criteria,
        "target_module": target or "g_dashboard",
        "trace": trace,
        "matched_terms": len(trace),
    }


# --------------------------------------------------------------- compile step
def compile_and_validate(question: str, payload: Optional[Dict[str, Any]] = None,
                         user: Optional[str] = None) -> Dict[str, Any]:
    """The §D01 compile screen payload. NOTHING HAS RUN YET after this call."""
    compiled = payload if payload is not None else compile_question(question)
    result = validate(compiled)

    cid = "cmp-" + uuid.uuid4().hex[:10]
    db.insert_rows(
        "compilation_log",
        ["compilation_id", "asked_at", "question", "profile", "criteria_json",
         "target_module", "validation", "unmapped", "executed"],
        [[cid, datetime.utcnow(), question, config.PROFILE,
          json.dumps(compiled.get("criteria", {})), compiled.get("target_module"),
          result["validation"], json.dumps(result["unmapped"]), False]],
    )

    return {
        "compilation_id": cid,
        "question": question,
        "compiled": compiled,
        "compiled_json": json.dumps(
            {k: compiled[k] for k in ("profile", "criteria", "target_module")
             if k in compiled}, indent=2),
        "validation": result["validation"],
        "ok": result["ok"],
        "errors": result["errors"],
        "unmapped": result["unmapped"],
        "fields_mapped": result["fields_mapped"],
        "target_module": compiled.get("target_module"),
        "target_title": MODULE_TITLES.get(compiled.get("target_module") or "", ""),
        "output_class": result["output_class"],
        "kb_snapshot_id": db.meta_get("kb_snapshot_id", config.KB_SNAPSHOT_ID),
        "status": "Nothing has run yet.",
    }


def mark_executed(compilation_id: str) -> None:
    db.execute("UPDATE compilation_log SET executed = TRUE WHERE compilation_id = ?",
               [compilation_id])


def corpus(limit: int = 200) -> List[Dict[str, Any]]:
    """The evaluation corpus the acceptance gate is measured on (D01)."""
    return db.rows("""
        SELECT compilation_id, asked_at, question, criteria_json, target_module,
               validation, unmapped, executed
        FROM compilation_log ORDER BY asked_at DESC LIMIT ?
    """, [limit])


SUGGESTED_QUESTIONS = [
    "What is the diagnostic yield for neurodevelopmental delay referrals?",
    "How many BRCA1 and BRCA2 P/LP carriers are in the HBOC panel cohort?",
    "Carrier rate for HBB in South Asian subjects",
    "Show me biallelic P/LP findings in recessive genes",
    "Which VUS should we reclassify first?",
    "Founder variant candidates in the cohort",
    "Secondary findings among consented exome subjects",
    "Hemizygous X-linked males with pathogenic findings",
    "Lynch syndrome carriers referred in the last 12 months",
    "Gene-disease validity breakdown for established genes only",
    "Compound heterozygous subjects on the carrier screen",
    "Which genes have provisional coverage in this cohort?",
]
