"""Germline cohort engine — spec §G01.

The gating order is not a style choice. From the spec:

    1. all subjects
    2. CONSENT GATE      -- exclude consent_class = 'Withdrawn' entirely
    3. PROBANDS ONLY     -- if on, keep is_proband = true
    4. clinical criteria
    5. collection window
    6. assay criteria
    7. SF CONSENT GATE   -- if sfOnly, require explicit secondary-findings consent
    8. genomic + KB criteria (subject must have >=1 matching observation)

    "Steps 2, 3 and 7 run BEFORE everything else because they change what the
     cohort *is*, not merely what is displayed. A withdrawn-consent subject
     must not be reachable by any filter combination."

That is enforced structurally here: the consent gates are compiled into the
member SQL as unconditional predicates, so there is no criteria object — valid
or malformed — that can produce a cohort containing a withdrawn subject. It is
not a display filter (spec E04: "Display-filtering consent -> data one bug from
exposure").

Resolving a cohort materialises two temp tables that every analytics module
reads: `cohort_subject` and `cohort_run`. Numbers are computed by SQL over
those, never in the browser and never by a model.
"""
from __future__ import annotations

import hashlib
import json
import threading
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .. import config, db
from ..reference.genes import GENE_SETS
from ..reference.tests import CONSENT_SF_OK, CONSENT_WITHDRAWN

PLP = ("Pathogenic", "Likely pathogenic")
VUS = "Uncertain significance"

# A resolved cohort lives in TEMP TABLES on the DuckDB connection, and the app
# shares one connection across the request threadpool. Two requests resolving
# at once therefore race: the second CREATE lands between the first DROP and
# CREATE and fails with "Table already exists" — and worse, a request that wins
# the race can read another request's cohort mid-analysis and return numbers for
# the wrong cohort, silently.
#
# So the whole resolve-then-read sequence is one critical section. Use
# `cohort_scope()` — it holds the lock across resolve AND the reads that follow.
#
# `resolve()` alone locks only its own DROP/CREATE. That is enough when the
# caller never touches the temp tables afterwards (the CLI, the reproducibility
# check), and NOT enough for anything that then runs an analytic. Relying on a
# comment to tell those two cases apart is how a route gets missed, so the
# context manager is the supported entry point.
COHORT_LOCK = threading.RLock()


@contextmanager
def cohort_scope(raw_criteria: Optional[Dict[str, Any]] = None,
                 name: str = "All germline subjects · probands",
                 as_of: Optional[str] = None):
    """Resolve a cohort and hold the lock for as long as you read it.

        with cohort_scope(criteria, name) as co:
            payload = carrier.dashboard(co.criteria)   # reads cohort_subject

    Re-entrant, so nesting is safe.
    """
    with COHORT_LOCK:
        yield resolve(raw_criteria, name, as_of=as_of)

# ---------------------------------------------------------------- criteria ---
# The allowed field list. The query compiler (§D01) validates emitted field
# names against exactly this — anything outside it is a hard failure.
LIST_FIELDS = [
    # clinical
    "indication", "sex", "age_bucket", "ancestry", "affected_status",
    "family_history", "relation", "referral_source", "consent_class",
    # assay
    "test_code", "assay_version", "pipeline_version", "reference_build", "sample_type",
    # genomic
    "gene", "var_class", "consequence", "classification", "zygosity",
    "inheritance", "validity", "penetrance", "clinvar_sig",
]
SCALAR_FIELDS = {
    "gene_set": None,           # one of GENE_SETS
    "months": 0,                # collection window
    "gnomad_max": 0.0,          # max population AF
    "probands_only": True,      # family de-duplication, ON by default
    "reportable_only": True,
    "sf_only": False,           # ACMG SF v3.3 — consent-gated
    "founder_only": False,
    "curated_only": False,
    "established_validity_only": False,   # §G05 yield-inflation toggle
}
ALLOWED_FIELDS = set(LIST_FIELDS) | set(SCALAR_FIELDS)

CRITERIA_LABELS = {
    "indication": "Indication", "sex": "Sex", "age_bucket": "Age",
    "ancestry": "Ancestry", "affected_status": "Status",
    "family_history": "Family hx", "relation": "Relation",
    "referral_source": "Referral", "consent_class": "Consent",
    "test_code": "Test code", "assay_version": "Assay ver",
    "pipeline_version": "Pipeline", "reference_build": "Reference",
    "sample_type": "Sample", "gene": "Gene", "var_class": "Class",
    "consequence": "Consequence", "classification": "ACMG",
    "zygosity": "Zygosity", "inheritance": "Inheritance",
    "validity": "Gene validity", "penetrance": "Penetrance",
    "clinvar_sig": "ClinVar", "gene_set": "Gene set", "months": "Collected within",
    "gnomad_max": "Max gnomAD AF", "probands_only": "Probands only",
    "reportable_only": "Reportable only", "sf_only": "ACMG SF v3.3 only",
    "founder_only": "Founder candidates only", "curated_only": "Curated only",
    "established_validity_only": "Established validity only",
}

# Subject-level vs observation-level: which table each criterion constrains.
CLINICAL_FIELDS = ["indication", "sex", "age_bucket", "ancestry", "affected_status",
                   "family_history", "relation", "referral_source", "consent_class"]
ASSAY_FIELDS = ["test_code", "assay_version", "pipeline_version",
                "reference_build", "sample_type"]
GENOMIC_FIELDS = ["gene", "var_class", "consequence", "classification", "zygosity",
                  "inheritance", "validity", "penetrance", "clinvar_sig"]

_SUBJECT_COL = {
    "indication": "s.indication", "sex": "s.sex", "age_bucket": "s.age_bucket",
    "ancestry": "s.ancestry", "affected_status": "s.affected_status",
    "family_history": "s.family_history", "relation": "s.relation",
    "referral_source": "s.referral_source", "consent_class": "s.consent_class",
}
_RUN_COL = {
    "test_code": "r.test_code", "assay_version": "r.assay_version",
    "pipeline_version": "r.pipeline_version", "reference_build": "r.reference_build",
    "sample_type": "sm.sample_type",
}
_OBS_COL = {
    "gene": "f.gene_symbol", "var_class": "f.var_class", "consequence": "f.consequence",
    "classification": "i.classification", "zygosity": "f.zygosity",
    "inheritance": "gd.inheritance", "validity": "gd.validity",
    "penetrance": "gd.penetrance", "clinvar_sig": "f.clinvar_sig",
}


def blank() -> Dict[str, Any]:
    c: Dict[str, Any] = {k: [] for k in LIST_FIELDS}
    c.update(SCALAR_FIELDS)
    return c


def normalise(raw: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Merge a partial criteria object onto the blank, dropping unknown keys."""
    c = blank()
    for k, v in (raw or {}).items():
        if k not in ALLOWED_FIELDS:
            continue
        if k in LIST_FIELDS:
            c[k] = list(v) if isinstance(v, (list, tuple)) else ([v] if v else [])
        elif k in ("months",):
            c[k] = int(v or 0)
        elif k in ("gnomad_max",):
            c[k] = float(v or 0)
        elif k == "gene_set":
            c[k] = v if v in GENE_SETS else None
        else:
            c[k] = bool(v)
    return c


def unknown_fields(raw: Dict[str, Any]) -> List[str]:
    return sorted(k for k in (raw or {}) if k not in ALLOWED_FIELDS)


def criteria_hash(c: Dict[str, Any]) -> str:
    canonical = json.dumps(_canonical(c), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()[:8]


def _canonical(c: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k in sorted(ALLOWED_FIELDS):
        v = c.get(k)
        if k in LIST_FIELDS:
            if v:
                out[k] = sorted(v)
        elif v != SCALAR_FIELDS.get(k):
            out[k] = v
    return out


def active_chips(c: Dict[str, Any]) -> List[Dict[str, Any]]:
    """One chip per active filter (spec §2.3).

    "Chips are the only place a user can see the whole cohort definition at a
     glance. Every filter must produce a chip; a filter with no chip is a
     filter users will forget is on."

    That includes probands_only, which defaults to ON — a default that changes
    every rate on the screen is exactly the kind of filter a user forgets.
    """
    chips: List[Dict[str, Any]] = []
    for k in LIST_FIELDS:
        for v in c.get(k) or []:
            chips.append({"field": k, "value": v,
                          "label": "{}: {}".format(CRITERIA_LABELS[k], v)})
    if c.get("gene_set"):
        chips.append({"field": "gene_set", "value": c["gene_set"],
                      "label": "Gene set: {}".format(c["gene_set"])})
    if c.get("months"):
        chips.append({"field": "months", "value": c["months"],
                      "label": "Collected within {} months".format(c["months"])})
    if c.get("gnomad_max"):
        chips.append({"field": "gnomad_max", "value": c["gnomad_max"],
                      "label": "gnomAD AF ≤ {}".format(c["gnomad_max"])})
    for k in ("probands_only", "reportable_only", "sf_only", "founder_only",
              "curated_only", "established_validity_only"):
        if c.get(k):
            chips.append({"field": k, "value": True, "label": CRITERIA_LABELS[k],
                          "sticky": k == "probands_only"})
    return chips


# ------------------------------------------------------------ SQL building ---
def _in_clause(col: str, values: Sequence[Any], params: List[Any]) -> str:
    params.extend(values)
    return "{} IN ({})".format(col, ", ".join(["?"] * len(values)))


def _obs_predicates(c: Dict[str, Any], params: List[Any]) -> List[str]:
    """Predicates on a single observation (finding + interpretation + gene)."""
    preds: List[str] = []
    if c.get("reportable_only"):
        preds.append("i.reportable")
    for k in GENOMIC_FIELDS:
        vals = c.get(k) or []
        if vals:
            preds.append(_in_clause(_OBS_COL[k], vals, params))
    if c.get("gene_set"):
        genes = GENE_SETS.get(c["gene_set"]) or []
        preds.append(_in_clause("f.gene_symbol", genes or ["\x00"], params))
    if c.get("gnomad_max"):
        params.append(c["gnomad_max"])
        preds.append("(f.gnomad_af IS NOT NULL AND f.gnomad_af <= ?)")
    if c.get("established_validity_only"):
        preds.append("gd.validity IN ('Definitive','Strong')")
    if c.get("curated_only"):
        preds.append("i.curated_flag")
    if c.get("sf_only"):
        # SF genes only — the consent side of this gate is applied at subject
        # level in step 7, not here.
        preds.append("gd.gene_sets LIKE '%SF%'")
    if c.get("founder_only"):
        # A founder candidate is a P/LP allele recurrent in the cohort; the
        # recurrence test is the subquery, so it cannot be faked per-row.
        preds.append("""i.classification IN ('Pathogenic','Likely pathogenic')
            AND f.variant_key IN (
                SELECT variant_key FROM finding
                GROUP BY variant_key HAVING COUNT(DISTINCT subject_id) >= 2)""")
    return preds


def observation_filter_sql(c: Dict[str, Any], params: List[Any]) -> str:
    preds = _obs_predicates(c, params)
    return " AND ".join(preds) if preds else "TRUE"


def has_genomic_criteria(c: Dict[str, Any]) -> bool:
    if any(c.get(k) for k in GENOMIC_FIELDS):
        return True
    return bool(c.get("gene_set") or c.get("gnomad_max") or c.get("sf_only")
                or c.get("founder_only") or c.get("curated_only")
                or c.get("established_validity_only"))


OBS_JOIN = """
    FROM finding f
    JOIN interpretation i USING (finding_id)
    JOIN gene_disease gd ON gd.gene_symbol = f.gene_symbol
"""


# ------------------------------------------------------------------ resolve ---
@dataclass
class Cohort:
    criteria: Dict[str, Any]
    name: str
    funnel: List[Dict[str, Any]]
    n_subjects: int
    n_families: int
    n_observations: int
    n_runs: int
    criteria_hash: str
    member_hash: str
    resolved_at: str
    kb_snapshot_id: str
    as_of: Optional[str] = None
    warnings: List[Dict[str, str]] = field(default_factory=list)

    def as_dict(self) -> Dict[str, Any]:
        return {
            "criteria": self.criteria, "name": self.name, "funnel": self.funnel,
            "counts": {"subjects": self.n_subjects, "families": self.n_families,
                       "observations": self.n_observations, "runs": self.n_runs},
            "criteria_hash": self.criteria_hash, "member_hash": self.member_hash,
            "resolved_at": self.resolved_at, "kb_snapshot_id": self.kb_snapshot_id,
            "as_of": self.as_of,
            "chips": active_chips(self.criteria), "warnings": self.warnings,
        }


def collection_anchor() -> Optional[str]:
    """The date a relative collection window is measured back from."""
    v = db.scalar("SELECT MAX(collection_date) FROM sample")
    return str(v)[:10] if v is not None else None


def resolve(raw_criteria: Optional[Dict[str, Any]] = None,
            name: str = "All germline subjects · probands",
            as_of: Optional[str] = None) -> Cohort:
    """Run the §G01 pipeline and materialise the cohort member tables.

    `as_of` pins the collection-window anchor. Left unset it defaults to the
    store's most recent collection, which is stable against wall-clock drift but
    NOT against ingest: loading new samples moves the anchor, so "last 24 months"
    silently covers a different window and the member hash changes. That is real
    membership change rather than a bug, but a reproducibility check has to be
    able to tell the two apart — so the resolved anchor is returned on the cohort
    and carried in the manifest, and re-runs pass it back.
    """
    with COHORT_LOCK:
        return _resolve(raw_criteria, name, as_of)


def _resolve(raw_criteria: Optional[Dict[str, Any]], name: str,
             as_of: Optional[str] = None) -> Cohort:
    c = normalise(raw_criteria)
    anchor = as_of or collection_anchor()
    funnel: List[Dict[str, Any]] = []

    db.execute("DROP TABLE IF EXISTS cohort_subject")
    db.execute("DROP TABLE IF EXISTS cohort_run")

    # ---- step 1: all subjects -------------------------------------------
    total = int(db.scalar("SELECT COUNT(*) FROM subject", default=0))
    funnel.append({"step": "All subjects in ImpactOmics", "group": "base",
                   "count": total, "dropped": 0})
    prev = total

    where: List[str] = []
    params: List[Any] = []

    # ---- step 2: has sequencing data ------------------------------------
    # Every later step counts through a join to run and sample, so a subject
    # registered in the clinical system but not yet sequenced would silently
    # disappear at whatever step came next — and be blamed on it. The funnel
    # promises every exclusion is attributable to exactly one choice, so the
    # join gets its own step rather than hiding inside the consent gate.
    n = _count([], [])
    funnel.append({"step": "Has sequencing data", "group": "base",
                   "count": n, "dropped": prev - n, "mandatory": True})
    prev = n

    # ---- step 3: CONSENT GATE -------------------------------------------
    # Unconditional. Not driven by criteria, so no criteria object can lift it.
    where.append("s.consent_class <> ?")
    params.append(CONSENT_WITHDRAWN)
    n = _count(where, params)
    funnel.append({"step": "Consent gate — withdrawn excluded", "group": "consent",
                   "count": n, "dropped": prev - n, "mandatory": True})
    prev = n

    # ---- step 3: PROBANDS ONLY ------------------------------------------
    if c["probands_only"]:
        where.append("s.is_proband")
        n = _count(where, params)
        funnel.append({"step": "Probands only — family de-duplication",
                       "group": "family", "count": n, "dropped": prev - n})
        prev = n

    # ---- step 4: clinical criteria --------------------------------------
    clinical_preds: List[str] = []
    for k in CLINICAL_FIELDS:
        vals = c.get(k) or []
        if vals:
            clinical_preds.append(_in_clause(_SUBJECT_COL[k], vals, params))
    if clinical_preds:
        where.extend(clinical_preds)
        n = _count(where, params)
        funnel.append({"step": "Clinical criteria", "group": "clinical",
                       "count": n, "dropped": prev - n})
        prev = n

    # ---- step 5: collection window --------------------------------------
    if c["months"]:
        # Anchored on the store's most recent collection, not wall-clock today,
        # so a saved cohort re-runs to an identical number (spec E03.6) instead
        # of silently shrinking as the calendar moves.
        where.append("sm.collection_date >= (CAST(? AS DATE) - INTERVAL (?) DAY)")
        params.extend([anchor, int(round(c["months"] * 30.44))])
        n = _count(where, params)
        funnel.append({"step": "Collected within {} months".format(c["months"]),
                       "group": "clinical", "count": n, "dropped": prev - n})
        prev = n

    # ---- step 6: assay criteria -----------------------------------------
    assay_preds: List[str] = []
    for k in ASSAY_FIELDS:
        vals = c.get(k) or []
        if vals:
            assay_preds.append(_in_clause(_RUN_COL[k], vals, params))
    if assay_preds:
        where.extend(assay_preds)
        n = _count(where, params)
        funnel.append({"step": "Assay criteria", "group": "assay",
                       "count": n, "dropped": prev - n})
        prev = n

    # ---- step 7: SF CONSENT GATE ----------------------------------------
    # Spec §G09: subjects who declined are excluded from the cohort QUERY
    # entirely, not filtered from display.
    if c["sf_only"]:
        where.append(_in_clause("s.consent_class", sorted(CONSENT_SF_OK), params))
        n = _count(where, params)
        funnel.append({"step": "Secondary-findings consent required",
                       "group": "consent", "count": n, "dropped": prev - n,
                       "mandatory": True})
        prev = n

    # ---- step 8: genomic + KB criteria ----------------------------------
    if has_genomic_criteria(c):
        obs_params: List[Any] = []
        obs_sql = observation_filter_sql(c, obs_params)
        where.append("""s.subject_id IN (
            SELECT f.subject_id {join} WHERE {pred})""".format(join=OBS_JOIN, pred=obs_sql))
        params.extend(obs_params)
        n = _count(where, params)
        funnel.append({"step": "Genomic / knowledgebase criteria", "group": "genomic",
                       "count": n, "dropped": prev - n})
        prev = n

    # ---- materialise ----------------------------------------------------
    db.execute("""
        CREATE TEMP TABLE cohort_subject AS
        SELECT DISTINCT s.subject_id, s.family_id FROM subject s
        JOIN run r ON r.subject_id = s.subject_id
        JOIN sample sm ON sm.sample_id = r.sample_id
        WHERE {}
    """.format(" AND ".join(where)), params)

    # Only the runs that MATCH the assay criteria, not every run those subjects
    # ever had. Otherwise filtering to a 1-gene thalassemia panel and finding the
    # subject also had an exome credits that exome's scope to the cohort: genes
    # the selected panel never looks at acquire real denominators, and a gene
    # with no findings reads "0 / 40" — tested and clear — when nothing on the
    # chosen test code ever examined it.
    run_where = ["r.subject_id IN (SELECT subject_id FROM cohort_subject)"]
    run_params: List[Any] = []
    for k in ASSAY_FIELDS:
        vals = c.get(k) or []
        if vals:
            col = _RUN_COL[k]
            run_where.append(_in_clause(col.replace("sm.", "sm."), vals, run_params))
    db.execute("""
        CREATE TEMP TABLE cohort_run AS
        SELECT r.run_id, r.subject_id
        FROM run r JOIN sample sm ON sm.sample_id = r.sample_id
        WHERE {}
    """.format(" AND ".join(run_where)), run_params)

    n_subjects = int(db.scalar("SELECT COUNT(*) FROM cohort_subject", default=0))
    n_families = int(db.scalar(
        "SELECT COUNT(DISTINCT family_id) FROM cohort_subject", default=0))
    n_runs = int(db.scalar("SELECT COUNT(*) FROM cohort_run", default=0))
    funnel.append({"step": "Cohort", "group": "final", "count": n_subjects,
                   "dropped": 0, "final": True})

    obs_params2: List[Any] = []
    obs_sql2 = observation_filter_sql({"reportable_only": c["reportable_only"]}, obs_params2)
    n_obs = int(db.scalar("""
        SELECT COUNT(*) {join}
        WHERE f.subject_id IN (SELECT subject_id FROM cohort_subject) AND {pred}
    """.format(join=OBS_JOIN, pred=obs_sql2), obs_params2, default=0))

    member_hash = _member_hash()
    return Cohort(
        as_of=anchor,
        criteria=c, name=name, funnel=funnel, n_subjects=n_subjects,
        n_families=n_families, n_observations=n_obs, n_runs=n_runs,
        criteria_hash=criteria_hash(c), member_hash=member_hash,
        resolved_at=datetime.utcnow().isoformat(timespec="seconds") + "Z",
        kb_snapshot_id=db.meta_get("kb_snapshot_id", config.KB_SNAPSHOT_ID),
        warnings=collect_warnings(c),
    )


def _count(where: List[str], params: List[Any]) -> int:
    """Subjects surviving `where`, counted through the run/sample join every
    later step also uses. An empty predicate list means "just the join"."""
    return int(db.scalar("""
        SELECT COUNT(DISTINCT s.subject_id) FROM subject s
        JOIN run r ON r.subject_id = s.subject_id
        JOIN sample sm ON sm.sample_id = r.sample_id
        WHERE {}
    """.format(" AND ".join(where) if where else "TRUE"), params, default=0))


def _member_hash() -> str:
    """8-char hash of the sorted member IDs. Detects silent membership drift
    when the same criteria are re-run against a changed store (spec C02)."""
    ids = db.rows("SELECT subject_id FROM cohort_subject ORDER BY subject_id")
    joined = ",".join(r["subject_id"] for r in ids)
    return hashlib.sha256(joined.encode()).hexdigest()[:8]


# ------------------------------------------------------------------ warnings --
def collect_warnings(c: Dict[str, Any]) -> List[Dict[str, str]]:
    """Cohort-level warnings the context bar and dashboards must surface."""
    out: List[Dict[str, str]] = []

    # Probands-only must always state its effect (spec §G01).
    fam = family_stats()
    if c.get("probands_only"):
        out.append({
            "level": "info", "code": "probands_on",
            "text": ("One subject per family. {} related subjects are being excluded "
                     "to keep carrier rates independent.".format(fam["suppressed"])),
        })
    else:
        out.append({
            "level": "warn", "code": "probands_off",
            "text": ("{} related subjects included across {} multi-member families. "
                     "Carrier rates are inflated and must not be quoted."
                     .format(fam["extra_subjects"], fam["multi_member_families"])),
        })

    # Provenance mixing (spec §G13).
    prov = db.row("""
        SELECT COUNT(DISTINCT reference_build) AS builds,
               COUNT(DISTINCT pipeline_version) AS pipelines
        FROM run WHERE run_id IN (SELECT run_id FROM cohort_run)
    """) or {"builds": 0, "pipelines": 0}
    if (prov["builds"] or 0) > 1 or (prov["pipelines"] or 0) > 1:
        out.append({
            "level": "warn", "code": "provenance",
            "text": ("Cohort mixes {} reference builds and {} pipeline versions. "
                     "Acceptable for descriptive counts; not acceptable for frequency "
                     "comparison without stratification."
                     .format(prov["builds"], prov["pipelines"])),
        })

    # Provisional scope (spec §G02).
    prov_scope = provisional_stats()
    if prov_scope["provisional_subjects"]:
        out.append({
            "level": "warn", "code": "provisional",
            "text": ("{} subjects ({}%) ran on a test code with incomplete scope. "
                     "Counted in numerators, flagged as provisional in denominators."
                     .format(prov_scope["provisional_subjects"],
                             prov_scope["provisional_pct"])),
        })
    return out


def family_stats() -> Dict[str, int]:
    """Family clustering for the probands-only messaging (spec §G01)."""
    row = db.row("""
        WITH fam AS (
            SELECT family_id, COUNT(*) AS n FROM cohort_subject GROUP BY family_id)
        SELECT COUNT(*) AS families,
               COALESCE(SUM(CASE WHEN n > 1 THEN 1 ELSE 0 END), 0) AS multi,
               COALESCE(SUM(n - 1), 0) AS extra
        FROM fam
    """) or {}
    # What probands-only is currently suppressing: consented related subjects
    # in the cohort's families who would join if the toggle were off.
    suppressed = int(db.scalar("""
        SELECT COUNT(*) FROM subject s
        WHERE s.family_id IN (SELECT DISTINCT family_id FROM cohort_subject)
          AND NOT s.is_proband
          AND s.consent_class <> ?
          AND s.subject_id NOT IN (SELECT subject_id FROM cohort_subject)
    """, [CONSENT_WITHDRAWN], default=0))
    return {
        "families": int(row.get("families") or 0),
        "multi_member_families": int(row.get("multi") or 0),
        "extra_subjects": int(row.get("extra") or 0),
        "suppressed": suppressed,
    }


def provisional_stats() -> Dict[str, Any]:
    """Subjects whose run has any scope row of unverified confidence."""
    total = int(db.scalar("SELECT COUNT(*) FROM cohort_subject", default=0))
    n = int(db.scalar("""
        SELECT COUNT(DISTINCT cr.subject_id) FROM cohort_run cr
        WHERE cr.run_id IN (
            SELECT run_id FROM run_scope WHERE confidence = 'unknown')
    """, default=0))
    return {"provisional_subjects": n, "cohort_subjects": total,
            "provisional_pct": round(100.0 * n / total, 1) if total else 0.0}


def ensure_resolved() -> None:
    """Guard: every analytics call requires a materialised cohort."""
    exists = db.scalar("""
        SELECT COUNT(*) FROM duckdb_tables()
        WHERE table_name = 'cohort_subject'
    """, default=0)
    if not exists:
        resolve()
