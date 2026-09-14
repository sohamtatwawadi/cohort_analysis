"""Builds the analytics store: reference load, ingest, and the derivation passes.

The store is derived and rebuildable (spec E03.7). Two ingest paths feed it and
both land in exactly the same tables and run exactly the same derivations:

    varimat   — N per-sample VariMAT files + a clinical/LIMS sidecar
    synthetic — a deterministic generator, for development and for the
                hand-count verification test (spec E03.9)

Derivations run AFTER ingest, over the whole store, because they are
cohort-independent facts that must not vary with which cohort is being viewed:

    scope     §3.4  per-run, per-gene assayability -> every denominator
    zygosity  §3.5  hemizygous and compound-het, COMPUTED not read
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .. import config, db
from ..reference.genes import GENE_DISEASE, GENE_LIST
from ..reference.tests import TEST_CODES
from . import fingerprint as fp_mod
from . import varimat as vm

PLP = ("Pathogenic", "Likely pathogenic")


# ------------------------------------------------------------------ lifecycle --
# Lab Mode data is DERIVED and rebuildable from VariMAT. Research Mode data is
# UPLOADED BY SOMEONE ELSE and is not ours to delete (Part II §8). Rebuilding
# the lab store used to drop the whole database file, which took every
# registered project, dataset, job and result with it — and left the genotype
# payloads orphaned on disk. So a rebuild clears the lab tables and nothing else.
LAB_TABLES = [
    "finding", "interpretation", "run_scope", "run_fingerprint", "run",
    "sample", "subject", "family", "gene_disease", "test_code",
    "test_code_gene", "indication_gene", "cohort_snapshot", "analysis_run",
]


def reset_lab() -> None:
    """Clear the derived lab tables, preserving everything a user uploaded."""
    db.connect()
    for table in LAB_TABLES:
        db.execute("DELETE FROM {}".format(table))
    db.meta_set("built_at", datetime.utcnow().isoformat(timespec="seconds"))
    db.meta_set("kb_snapshot_id", config.KB_SNAPSHOT_ID)


def init(fresh: bool = False) -> None:
    """`fresh` rebuilds the LAB store. It does not touch uploaded datasets —
    use `db.connect(fresh=True)` directly if a full wipe is really intended."""
    if fresh:
        reset_lab()
    else:
        db.connect()


def is_empty() -> bool:
    return db.scalar("SELECT COUNT(*) FROM subject", default=0) == 0


def load_reference() -> None:
    """Write the curated reference tables. Idempotent."""
    db.execute("DELETE FROM gene_disease")
    db.insert_rows(
        "gene_disease",
        ["gene_symbol", "gene_id", "condition", "mondo_id", "inheritance",
         "validity", "penetrance", "gene_sets", "protein_len", "chrom"],
        [[r.symbol, r.gene_id, r.condition, r.mondo, r.inheritance, r.validity,
          r.penetrance, ",".join(r.sets), r.protein_len, r.chrom]
         for r in GENE_DISEASE.values()],
    )
    db.execute("DELETE FROM test_code")
    db.insert_rows(
        "test_code",
        ["code", "name", "assay_version", "indication", "sf_capable",
         "scope_complete", "accredited"],
        [[t.code, t.name, t.assay_version, t.indication, t.sf_capable,
          t.scope_complete, t.accredited] for t in TEST_CODES.values()],
    )
    db.execute("DELETE FROM test_code_gene")
    db.insert_rows(
        "test_code_gene",
        ["code", "gene_symbol", "reportable"],
        [[t.code, g, True] for t in TEST_CODES.values() for g in t.genes],
    )

    # The phenotype-relevance gate, materialised so the yield query can join it
    # rather than reimplementing the map (spec §G02).
    from ..reference.tests import INDICATIONS, relevant_genes
    db.execute("DELETE FROM indication_gene")
    db.insert_rows(
        "indication_gene", ["indication", "gene_symbol"],
        [[ind, g] for ind in INDICATIONS for g in relevant_genes(ind)],
    )
    db.meta_set("kb_snapshot_id", config.KB_SNAPSHOT_ID)


# ------------------------------------------------------------------- persist --
FAMILY_COLS = ["family_id", "tenant_id", "ancestry", "consent_class"]
SUBJECT_COLS = ["subject_id", "tenant_id", "mrn", "sex", "age", "age_bucket",
                "ancestry", "consent_class", "family_id", "relation", "is_proband",
                "indication", "affected_status", "family_history", "referral_source",
                "phenotype_hpo"]
SAMPLE_COLS = ["sample_id", "subject_id", "sample_type", "collection_date"]
RUN_COLS = ["run_id", "sample_id", "subject_id", "test_code", "assay_version",
            "pipeline_version", "reference_build", "caller", "implied_panel_id",
            "mean_depth", "pct_bases_20x", "qc_status", "collection_date", "source_file",
            "qc_n_called", "qc_ti_tv", "qc_het_hom", "qc_mean_het_vaf", "qc_x_het_rate"]
FINDING_COLS = ["finding_id", "run_id", "subject_id", "family_id", "gene_id",
                "gene_symbol", "variant_key", "chrom", "pos", "ref_allele", "alt_allele",
                "hgvs_c", "hgvs_p", "aa_pos", "protein_len", "var_class", "consequence",
                "zygosity_raw", "zygosity", "vaf", "depth", "alt_depth", "filter_status",
                "gnomad_af", "gnomad_sas_af", "ga100k_sas_af", "clinvar_sig",
                "clinvar_id", "mane", "reviewable"]
INTERP_COLS = ["finding_id", "framework", "classification", "acmg_codes",
               "kb_snapshot_id", "curated_flag", "interpreted_at", "reportable",
               "inherited_from", "segregation", "evidence_delta"]


def persist(bundle: Dict[str, List[Dict[str, Any]]]) -> None:
    """Insert a full set of records. Keys map to table names."""
    spec = [
        ("family", FAMILY_COLS), ("subject", SUBJECT_COLS), ("sample", SAMPLE_COLS),
        ("run", RUN_COLS), ("finding", FINDING_COLS), ("interpretation", INTERP_COLS),
    ]
    for table, cols in spec:
        records = bundle.get(table) or []
        db.insert_rows(table, cols, [[r.get(c) for c in cols] for r in records])

    prints = bundle.get("run_fingerprint") or []
    db.insert_rows("run_fingerprint", ["run_id", "gene_symbol"],
                   [[r["run_id"], r["gene_symbol"]] for r in prints])


# --------------------------------------------------------------- derivations --
def derive_scope() -> Dict[str, Any]:
    """§3.4. Rebuild run_scope for every run in the store.

    Runs whose test code is known contribute the calibration set; runs without
    one are resolved by fingerprint clustering. A run on a test code whose
    registry entry has `scope_complete = false` is marked 'unknown' — counted
    in numerators, provisional in denominators (§G02).
    """
    prints: Dict[str, Set[str]] = {}
    for r in db.rows("SELECT run_id, gene_symbol FROM run_fingerprint"):
        prints.setdefault(r["run_id"], set()).add(r["gene_symbol"])

    declared: Dict[str, Set[str]] = {}
    complete: Dict[str, bool] = {}
    labels: Dict[str, str] = {}
    for r in db.rows("SELECT run_id, test_code FROM run"):
        prints.setdefault(r["run_id"], set())
        code = r["test_code"]
        if not code:
            continue
        t = TEST_CODES.get(code)
        if not t:
            continue
        declared[r["run_id"]] = set(t.genes)
        complete[r["run_id"]] = t.scope_complete
        labels[r["run_id"]] = code

    calls, clusters, thresholds = fp_mod.infer_scope(
        prints, GENE_LIST, declared_scope=declared,
        scope_complete=complete, labels=labels)

    db.execute("DELETE FROM run_scope")
    db.insert_rows(
        "run_scope",
        ["run_id", "gene_symbol", "variant_classes", "reportable", "confidence"],
        # Inference is valid for SNV/indel only; CNV and fusion absence is
        # uninformative, so the scope row does not claim them (§3.4).
        [[c.run_id, c.gene_symbol, "SNV,Indel,Splice", True, c.confidence]
         for c in calls],
    )

    # Stamp the implied panel back onto the run, for the C01 inspector.
    for c in clusters:
        if not c.run_ids:
            continue
        db.execute(
            "UPDATE run SET implied_panel_id = ? WHERE run_id IN ({})".format(
                ", ".join("?" * len(c.run_ids))),
            [c.implied_panel] + c.run_ids)

    db.meta_set("scope_thresholds", json.dumps(
        {g: round(t, 4) for g, t in sorted(thresholds.items())}))
    return {
        "scope_rows": len(calls),
        "clusters": [{"cluster_id": c.cluster_id, "implied_panel": c.implied_panel,
                      "runs": c.size, "genes": len(c.centroid)} for c in clusters],
    }


def derive_zygosity() -> Dict[str, int]:
    """§3.5. ZYGOSITY carries only Het/Hom. Hemizygous and compound
    heterozygous do not exist in the file and must be computed.

    Getting this wrong is the §E04 trap "reading ZYGOSITY for hemizygous —
    X-linked males misclassified". A hemizygous P/LP in an XLR gene in a male
    is an affected individual; the same call read as heterozygous reads as an
    unaffected carrier.
    """
    db.execute("UPDATE finding SET zygosity = zygosity_raw")

    # Hemizygous: a call on X or Y in a male subject. Single copy — it cannot
    # be heterozygous regardless of what the file says.
    hemi = db.scalar("""
        SELECT COUNT(*) FROM finding f JOIN subject s USING (subject_id)
        WHERE f.chrom IN ('X','Y','chrX','chrY') AND s.sex = 'M'
    """, default=0)
    db.execute("""
        UPDATE finding SET zygosity = 'Hemizygous'
        WHERE chrom IN ('X','Y','chrX','chrY')
          AND subject_id IN (SELECT subject_id FROM subject WHERE sex = 'M')
    """)

    # Compound heterozygous: >= 2 DISTINCT P/LP variants in the same gene in
    # the same subject. Distinct is on variant_key — two annotations of one
    # locus are not two hits. Phase is unconfirmed; §G04 flags it as such.
    db.execute("""
        UPDATE finding SET zygosity = 'Compound heterozygous'
        WHERE zygosity = 'Heterozygous'
          AND (subject_id, gene_symbol) IN (
            SELECT f.subject_id, f.gene_symbol
            FROM finding f JOIN interpretation i USING (finding_id)
            WHERE i.classification IN ('Pathogenic','Likely pathogenic')
              AND f.zygosity = 'Heterozygous'
            GROUP BY f.subject_id, f.gene_symbol
            HAVING COUNT(DISTINCT f.variant_key) >= 2)
          AND finding_id IN (
            SELECT finding_id FROM interpretation
            WHERE classification IN ('Pathogenic','Likely pathogenic'))
    """)
    chet = db.scalar(
        "SELECT COUNT(*) FROM finding WHERE zygosity = 'Compound heterozygous'",
        default=0)
    return {"hemizygous": int(hemi), "compound_het": int(chet)}


def rebuild_derivations() -> Dict[str, Any]:
    scope = derive_scope()
    zyg = derive_zygosity()
    return {"scope": scope, "zygosity": zyg}


# ------------------------------------------------------------ varimat ingest --
def ingest_varimat(
    paths: Sequence[Path],
    clinical: Optional[Dict[str, Dict[str, Any]]] = None,
    kb_snapshot_id: Optional[str] = None,
) -> Dict[str, Any]:
    """Ingest N VariMAT files — one per sample — as a set.

    `clinical` maps sample_id -> the sidecar record (subject, family, consent,
    indication, test code...). VariMAT carries none of that (spec §3.5); a
    sample with no sidecar record gets a stub subject in its own family so the
    genomic data is still loadable and visibly incomplete rather than dropped.
    """
    from .clinical import stub_record

    clinical = clinical or {}
    kb = kb_snapshot_id or config.KB_SNAPSHOT_ID

    # Two files resolving to one sample_id means either a re-sequenced sample or
    # a stale file left in the directory. Either way the operator has to decide
    # which run is authoritative — silently keeping one would put a denominator
    # on the wrong assay, and letting it reach the database surfaces as an
    # opaque primary-key error that says nothing about which files collided.
    by_sample: Dict[str, List[str]] = {}
    for p in paths:
        by_sample.setdefault(vm.parse_filename(Path(p))["sample_id"], []).append(Path(p).name)
    collisions = {s: f for s, f in by_sample.items() if len(f) > 1}
    if collisions:
        detail = "; ".join("{} <- {}".format(s, ", ".join(f))
                           for s, f in sorted(collisions.items()))
        raise ValueError(
            "{} sample ID(s) map to more than one VariMAT file. Remove the "
            "superseded file or pass an explicit file list. {}".format(
                len(collisions), detail))

    bundle: Dict[str, List[Dict[str, Any]]] = {
        "family": [], "subject": [], "sample": [], "run": [],
        "finding": [], "interpretation": [], "run_fingerprint": [],
    }
    seen_families: Set[str] = set()
    seen_subjects: Set[str] = set()
    funnels: List[Dict[str, Any]] = []
    fid = 0

    for sample in vm.parse_many(paths):
        rec = clinical.get(sample.sample_id) or stub_record(sample.sample_id)
        funnels.append(sample.funnel.as_dict())

        if rec["family_id"] not in seen_families:
            seen_families.add(rec["family_id"])
            bundle["family"].append({
                "family_id": rec["family_id"], "tenant_id": config.TENANT_ID,
                "ancestry": rec.get("ancestry"), "consent_class": rec.get("consent_class"),
            })
        if rec["subject_id"] not in seen_subjects:
            seen_subjects.add(rec["subject_id"])
            bundle["subject"].append({
                "subject_id": rec["subject_id"], "tenant_id": config.TENANT_ID,
                "mrn": rec.get("mrn"), "sex": rec.get("sex"), "age": rec.get("age"),
                "age_bucket": rec.get("age_bucket"), "ancestry": rec.get("ancestry"),
                "consent_class": rec["consent_class"], "family_id": rec["family_id"],
                "relation": rec.get("relation"), "is_proband": rec.get("is_proband", True),
                "indication": rec.get("indication"),
                "affected_status": rec.get("affected_status"),
                "family_history": rec.get("family_history"),
                "referral_source": rec.get("referral_source"),
                "phenotype_hpo": rec.get("phenotype_hpo"),
            })

        run_id = "GR-" + sample.sample_id
        bundle["sample"].append({
            "sample_id": sample.sample_id, "subject_id": rec["subject_id"],
            "sample_type": rec.get("sample_type"),
            "collection_date": rec.get("collection_date"),
        })
        test = TEST_CODES.get(rec.get("test_code") or "")
        bundle["run"].append({
            "run_id": run_id, "sample_id": sample.sample_id,
            "subject_id": rec["subject_id"], "test_code": rec.get("test_code"),
            "assay_version": test.assay_version if test else None,
            "pipeline_version": sample.pipeline_version or rec.get("pipeline_version"),
            "reference_build": sample.reference_build or rec.get("reference_build"),
            "caller": sample.caller, "implied_panel_id": None,
            "mean_depth": rec.get("mean_depth"), "pct_bases_20x": rec.get("pct_bases_20x"),
            "qc_status": rec.get("qc_status", "Pass"),
            "collection_date": rec.get("collection_date"),
            "source_file": sample.source_file,
            **sample.qc,
        })
        for g in sorted(sample.fingerprint):
            bundle["run_fingerprint"].append({"run_id": run_id, "gene_symbol": g})

        for v in sample.variants:
            if not v.gene_symbol:
                continue   # unresolved symbol: logged in the funnel, not counted
            fid += 1
            finding_id = "GO-{:07d}".format(fid)
            bundle["finding"].append({
                "finding_id": finding_id, "run_id": run_id,
                "subject_id": rec["subject_id"], "family_id": rec["family_id"],
                "gene_id": v.gene_id, "gene_symbol": v.gene_symbol,
                "variant_key": v.variant_key, "chrom": v.chrom, "pos": v.pos,
                "ref_allele": v.ref, "alt_allele": v.alt, "hgvs_c": v.hgvs_c,
                "hgvs_p": v.hgvs_p, "aa_pos": v.aa_pos, "protein_len": v.protein_len,
                "var_class": v.var_class, "consequence": v.consequence,
                "zygosity_raw": v.zygosity_raw, "zygosity": v.zygosity_raw,
                "vaf": v.vaf, "depth": v.depth, "alt_depth": v.alt_depth,
                "filter_status": v.filter_status, "gnomad_af": v.gnomad_af,
                "gnomad_sas_af": v.gnomad_sas_af, "ga100k_sas_af": v.ga100k_sas_af,
                "clinvar_sig": v.clinvar_sig, "clinvar_id": v.clinvar_id,
                "mane": v.mane, "reviewable": True,
            })
            reportable = bool(test and v.gene_symbol in test.genes) if test else True
            bundle["interpretation"].append({
                "finding_id": finding_id, "framework": "ACMG/AMP",
                "classification": v.acmg, "acmg_codes": v.acmg_codes,
                "kb_snapshot_id": kb, "curated_flag": False,
                "interpreted_at": rec.get("collection_date"),
                "reportable": reportable, "inherited_from": None,
                "segregation": "Not assessed", "evidence_delta": 0.0,
            })

    persist(bundle)
    derived = rebuild_derivations()
    return {
        "files": len(funnels),
        "subjects": len(bundle["subject"]),
        "runs": len(bundle["run"]),
        "findings": len(bundle["finding"]),
        "funnels": funnels,
        "derived": derived,
    }


def rebuild_from_varimat(
    directory: Optional[Path] = None,
    paths: Optional[Sequence[Path]] = None,
    clinical_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    from .clinical import load_sidecar

    files = list(paths) if paths else vm.discover(directory or config.VARIMAT_DIR)
    if not files:
        raise FileNotFoundError(
            "No VariMAT files found in {}".format(directory or config.VARIMAT_DIR))
    init(fresh=True)
    load_reference()
    sidecar = load_sidecar(clinical_dir or config.CLINICAL_DIR)
    result = ingest_varimat(files, clinical=sidecar)
    db.meta_set("source", "varimat")
    db.meta_set("source_files", json.dumps([str(f) for f in files]))
    return result
