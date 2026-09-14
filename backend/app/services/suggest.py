"""Analysis recommendation — which questions are worth asking of THIS cohort.

A flat menu of eight questions puts the same options in front of a user whose
cohort is 59 BRCA carriers and a user whose cohort is 600 unselected subjects,
even though most of those options are pointless for one of them. Worse, the
useless ones are indistinguishable from the useful ones until you click.

So we look at what the resolved cohort actually contains and rank the questions
against it, with the reason stated in counts. Three verdicts:

    recommended   the cohort has the data this question needs, and enough of it
    available     it will run, but the cohort was not built for it
    not_useful    it would return nothing, and we say why rather than let the
                  user find out by clicking

This is the same philosophy as Research Mode's capability gating (Part II §3):
never hide an option, always say what is missing. The difference is that here
nothing is *locked* — a lab cohort can always run any descriptive analysis. The
ranking is advice, not a gate, so the wording is "worth asking" not "allowed".
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from .. import db
from .carrier import PLP_SQL, SECONDARY_FINDING_SQL, SF_CONSENT_SQL, VUS_SQL, _obs_where
from .cohort import OBS_JOIN


def cohort_signals(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """Cheap counts that decide what is worth asking. One pass, not eight."""
    obs = _obs_where(criteria)
    base = ("WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject) "
            "AND {obs}".format(obs=obs))

    row = db.row("""
        SELECT
          COUNT(DISTINCT CASE WHEN {plp} THEN f.subject_id END)          AS plp_subjects,
          COUNT(DISTINCT CASE WHEN {plp} THEN f.gene_symbol END)         AS plp_genes,
          COUNT(DISTINCT CASE WHEN {vus} THEN f.variant_key END)         AS vus_variants,
          COUNT(DISTINCT CASE WHEN {plp} AND f.zygosity IN
                ('Homozygous','Compound heterozygous') THEN f.subject_id END) AS biallelic,
          COUNT(DISTINCT CASE WHEN {plp} AND f.zygosity = 'Hemizygous'
                THEN f.subject_id END)                                   AS hemizygous,
          COUNT(DISTINCT CASE WHEN {plp} AND gd.validity NOT IN ('Definitive','Strong')
                THEN f.subject_id END)                                   AS weak_validity,
          COUNT(DISTINCT CASE WHEN {plp} THEN gd.condition END)          AS conditions
        {join} {base}
    """.format(join=OBS_JOIN, base=base, plp=PLP_SQL, vus=VUS_SQL)) or {}

    recurrent = int(db.scalar("""
        SELECT COUNT(*) FROM (
          SELECT f.variant_key {join} {base}
          GROUP BY f.variant_key HAVING COUNT(DISTINCT f.subject_id) >= 2)
    """.format(join=OBS_JOIN, base=base), default=0))

    sf_subjects = int(db.scalar("""
        SELECT COUNT(DISTINCT f.subject_id) {join}
        JOIN subject s ON s.subject_id = f.subject_id
        {base} AND {plp} AND {secondary} AND {consent}
    """.format(join=OBS_JOIN, base=base, plp=PLP_SQL,
               secondary=SECONDARY_FINDING_SQL, consent=SF_CONSENT_SQL), default=0))

    demo = db.row("""
        SELECT COUNT(DISTINCT s.ancestry) AS ancestries,
               COUNT(DISTINCT s.indication) AS indications,
               COUNT(DISTINCT s.phenotype_hpo) AS hpo_terms,
               COUNT(*) AS subjects
        FROM cohort_subject cs JOIN subject s USING (subject_id)
    """) or {}

    multi_family = int(db.scalar("""
        SELECT COUNT(*) FROM (
          SELECT family_id FROM cohort_subject GROUP BY family_id HAVING COUNT(*) > 1)
    """, default=0))

    out = {k: int(v or 0) for k, v in row.items()}
    out.update({k: int(v or 0) for k, v in demo.items()})
    out["recurrent_variants"] = recurrent
    out["sf_subjects"] = sf_subjects
    out["multi_member_families"] = multi_family
    return out


def _n(n: int, one: str, many: Optional[str] = None) -> str:
    return "{:,} {}".format(n, one if n == 1 else (many or one + "s"))


def suggest(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """Rank the question catalogue against this cohort."""
    s = cohort_signals(criteria)
    c = criteria or {}
    picked = lambda k: bool(c.get(k))            # noqa: E731 — reads better inline

    out: List[Dict[str, Any]] = []

    def add(qid: str, verdict: str, reason: str, score: int) -> None:
        out.append({"question": qid, "verdict": verdict, "reason": reason,
                    "score": score})

    # --- diagnostic yield -------------------------------------------------
    if s["plp_subjects"] and s["indications"]:
        add("yield", "recommended",
            "{} carry a pathogenic finding across {} — yield is meaningful here."
            .format(_n(s["plp_subjects"], "subject"),
                    _n(s["indications"], "indication")), 90)
    elif s["plp_subjects"]:
        add("yield", "available",
            "{} carry a pathogenic finding, but every subject shares one indication, "
            "so there is nothing to stratify by.".format(_n(s["plp_subjects"], "subject")), 40)
    else:
        add("yield", "not_useful",
            "No pathogenic findings in this cohort, so the yield is zero by "
            "construction.", 0)

    # --- carrier rate -----------------------------------------------------
    if s["plp_genes"] >= 2:
        add("carrier", "recommended",
            "Pathogenic findings span {} — carrier rates will separate them."
            .format(_n(s["plp_genes"], "gene")),
            95 if picked("gene_set") or picked("gene") else 88)
    elif s["plp_genes"] == 1:
        add("carrier", "available", "Only one gene carries pathogenic findings.", 45)
    else:
        add("carrier", "not_useful", "No pathogenic findings to compute a rate from.", 0)

    # --- zygosity ---------------------------------------------------------
    interesting = s["biallelic"] + s["hemizygous"]
    if interesting:
        add("zygosity", "recommended",
            "{} biallelic and {} hemizygous — these are affected individuals, not "
            "carriers, and the distinction matters clinically."
            .format(s["biallelic"], s["hemizygous"]),
            98 if (picked("zygosity") or picked("inheritance")) else 85)
    elif s["plp_subjects"]:
        add("zygosity", "available",
            "All pathogenic findings here are heterozygous.", 35)
    else:
        add("zygosity", "not_useful", "No pathogenic findings to classify.", 0)

    # --- gene-disease validity -------------------------------------------
    if s["weak_validity"]:
        add("validity", "recommended",
            "{} on a gene with weaker than Strong evidence — this is exactly "
            "the inflation the yield figure hides."
            .format(_n(s["weak_validity"], "finding rests", "findings rest")), 92)
    elif s["conditions"] >= 3:
        add("validity", "available",
            "{} implicated, all on well-established genes."
            .format(_n(s["conditions"], "condition")), 50)
    else:
        add("validity", "not_useful", "Too few conditions implicated to compare.", 5)

    # --- VUS queue --------------------------------------------------------
    if s["vus_variants"] >= 20:
        add("vus", "recommended",
            "{} uncertain variants are waiting on review."
            .format("{:,}".format(s["vus_variants"])),
            97 if picked("classification") and "Uncertain significance" in (
                c.get("classification") or []) else 80)
    elif s["vus_variants"]:
        add("vus", "available",
            "{} to review — a short queue.".format(_n(s["vus_variants"], "variant")), 40)
    else:
        add("vus", "not_useful", "No uncertain variants in this cohort.", 0)

    # --- population frequency --------------------------------------------
    if s["recurrent_variants"] >= 5:
        add("popfreq", "recommended",
            "{} appear in two or more subjects — enough to spot founder alleles "
            "and implausible frequencies."
            .format(_n(s["recurrent_variants"], "variant")),
            96 if picked("founder_only") else 75)
    elif s["recurrent_variants"]:
        add("popfreq", "available",
            "Only {} recur.".format(_n(s["recurrent_variants"], "variant")), 30)
    else:
        add("popfreq", "not_useful",
            "No variant is seen in more than one subject, so there is no "
            "frequency to compare.", 0)

    # --- secondary findings -----------------------------------------------
    if s["sf_subjects"]:
        add("secondary", "recommended",
            "{} have a returnable secondary finding and consented to receive it."
            .format(_n(s["sf_subjects"], "subject")),
            99 if picked("sf_only") else 78)
    else:
        add("secondary", "not_useful",
            "No consented subject on an SF-capable test carries one. Nothing to "
            "return.", 0)

    # --- phenotype --------------------------------------------------------
    if s["hpo_terms"] >= 3 and s["plp_subjects"]:
        add("phenotype", "recommended",
            "{} recorded across the cohort — worth cross-tabulating against genes."
            .format(_n(s["hpo_terms"], "phenotype term")), 70)
    elif s["plp_subjects"]:
        add("phenotype", "available", "Few distinct phenotype terms recorded.", 25)
    else:
        add("phenotype", "not_useful", "No pathogenic findings to cross-tabulate.", 0)

    out.sort(key=lambda x: -x["score"])
    return {"signals": s, "suggestions": out,
            "summary": _summary(out, s)}


def _summary(out: List[Dict[str, Any]], s: Dict[str, Any]) -> str:
    top = [x for x in out if x["verdict"] == "recommended"][:3]
    if not top:
        return ("This cohort has no pathogenic or uncertain findings, so most "
                "analyses would return nothing. Widen the cohort in step 1.")
    from .suggest import QUESTION_TITLES
    names = [QUESTION_TITLES.get(x["question"], x["question"]) for x in top]
    return "Based on what this cohort contains, start with: " + "; ".join(names) + "."


QUESTION_TITLES = {
    "yield": "diagnostic yield",
    "carrier": "carrier rate by gene",
    "zygosity": "zygosity and inheritance",
    "validity": "gene–disease evidence",
    "vus": "the uncertain-variant queue",
    "popfreq": "population frequency",
    "secondary": "secondary findings",
    "phenotype": "phenotype cross-tab",
}
