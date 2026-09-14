-- ============================================================================
-- Germline Cohort Analytics — derived analytics schema
-- Implements Cohort-Analytics-Developer-Spec.md §3.6, germline subset.
--
-- This store is REBUILDABLE and READ-ONLY to the application (spec E03.7,
-- E04 "Analytics store treated as truth"). It is derived from VariMAT files
-- plus the clinical/LIMS sidecar; it is never the system of record.
-- ============================================================================

-- ---------------------------------------------------------------- reference --
-- Gene-disease validity. NOT available in VariMAT (spec §G05) — curated,
-- sourced from ClinGen where available.
CREATE TABLE IF NOT EXISTS gene_disease (
    gene_symbol   VARCHAR PRIMARY KEY,
    gene_id       VARCHAR,              -- HGNC:nnnnn
    condition     VARCHAR NOT NULL,
    mondo_id      VARCHAR,
    inheritance   VARCHAR NOT NULL,     -- AD | AR | AR/AD | XLR | XLD
    validity      VARCHAR NOT NULL,     -- Definitive | Strong | Moderate | Limited
    penetrance    VARCHAR NOT NULL,     -- High | Moderate | Low
    gene_sets     VARCHAR NOT NULL,     -- comma-joined set membership
    protein_len   INTEGER,
    chrom         VARCHAR               -- for hemizygosity derivation (§3.5)
);

-- Test-code registry (mirrors the lab's All_Test_Codes_Master).
CREATE TABLE IF NOT EXISTS test_code (
    code            VARCHAR PRIMARY KEY,
    name            VARCHAR NOT NULL,
    assay_version   VARCHAR NOT NULL,
    indication      VARCHAR NOT NULL,
    sf_capable      BOOLEAN NOT NULL,   -- can report ACMG SF v3.3
    scope_complete  BOOLEAN NOT NULL,   -- registry declares full reportable scope
    accredited      BOOLEAN NOT NULL
);

-- Phenotype-relevance gate for diagnostic yield (spec §G02). Without it a
-- broad exome counts ANY P/LP as diagnostic and yield goes to 81%. A P/LP
-- finding in an unrelated organ system is a secondary finding, not a diagnosis.
CREATE TABLE IF NOT EXISTS indication_gene (
    indication   VARCHAR NOT NULL,
    gene_symbol  VARCHAR NOT NULL,
    PRIMARY KEY (indication, gene_symbol)
);

CREATE TABLE IF NOT EXISTS test_code_gene (
    code         VARCHAR NOT NULL,
    gene_symbol  VARCHAR NOT NULL,
    reportable   BOOLEAN NOT NULL,
    PRIMARY KEY (code, gene_symbol)
);

-- --------------------------------------------------------------- core facts --
CREATE TABLE IF NOT EXISTS family (
    family_id      VARCHAR PRIMARY KEY,
    tenant_id      INTEGER NOT NULL,
    ancestry       VARCHAR,
    consent_class  VARCHAR
);

CREATE TABLE IF NOT EXISTS subject (
    subject_id      VARCHAR PRIMARY KEY,
    tenant_id       INTEGER NOT NULL,
    mrn             VARCHAR,
    sex             VARCHAR,            -- F | M  (drives hemizygosity, §3.5)
    age             INTEGER,
    age_bucket      VARCHAR,
    ancestry        VARCHAR,
    consent_class   VARCHAR NOT NULL,   -- gates cohort membership (§G01 step 2)
    family_id       VARCHAR NOT NULL,
    relation        VARCHAR,
    is_proband      BOOLEAN NOT NULL,
    indication      VARCHAR,
    affected_status VARCHAR,
    family_history  VARCHAR,
    referral_source VARCHAR,
    phenotype_hpo   VARCHAR             -- comma-joined HPO terms
);

CREATE TABLE IF NOT EXISTS sample (
    sample_id       VARCHAR PRIMARY KEY,
    subject_id      VARCHAR NOT NULL,
    sample_type     VARCHAR,
    collection_date DATE
);

CREATE TABLE IF NOT EXISTS run (
    run_id            VARCHAR PRIMARY KEY,
    sample_id         VARCHAR NOT NULL,
    subject_id        VARCHAR NOT NULL,
    test_code         VARCHAR,
    assay_version     VARCHAR,
    pipeline_version  VARCHAR,
    reference_build   VARCHAR,
    caller            VARCHAR,
    implied_panel_id  VARCHAR,          -- from fingerprint clustering (§3.4)
    mean_depth        DOUBLE,
    pct_bases_20x     DOUBLE,
    qc_status         VARCHAR,
    collection_date   DATE,
    source_file       VARCHAR,          -- provenance back to the VariMAT file

    -- Sample QC, computed at ingest over EVERY PASS call in the source file.
    -- It cannot be derived later: §3.2 has us archive the intronic/intergenic
    -- rows after fingerprinting, and `finding` retains only the reviewable
    -- subset — one to seven rows per sample, far too few for a ratio. These
    -- are the numbers a contaminated sample or a swap shows up in.
    qc_n_called       INTEGER,          -- PASS calls the metrics were computed over
    qc_ti_tv          DOUBLE,
    qc_het_hom        DOUBLE,
    qc_mean_het_vaf   DOUBLE,
    qc_x_het_rate     DOUBLE
);

-- Per-run, per-gene assayability. THE denominator table (§3.4, §3.7, C01).
-- confidence: 'declared'  — the test-code registry declares this gene reportable
--             'inferred'  — derived from fingerprint clustering
--             'unknown'   — in scope but reportable status unverified (provisional)
CREATE TABLE IF NOT EXISTS run_scope (
    run_id         VARCHAR NOT NULL,
    gene_symbol    VARCHAR NOT NULL,
    variant_classes VARCHAR,            -- classes this scope is valid for
    reportable     BOOLEAN NOT NULL,
    confidence     VARCHAR NOT NULL,
    PRIMARY KEY (run_id, gene_symbol)
);

-- Per-run gene-presence fingerprint (§3.2 storage requirement (b)).
CREATE TABLE IF NOT EXISTS run_fingerprint (
    run_id      VARCHAR NOT NULL,
    gene_symbol VARCHAR NOT NULL,
    PRIMARY KEY (run_id, gene_symbol)
);

CREATE TABLE IF NOT EXISTS finding (
    finding_id      VARCHAR PRIMARY KEY,
    run_id          VARCHAR NOT NULL,
    subject_id      VARCHAR NOT NULL,
    family_id       VARCHAR NOT NULL,
    gene_id         VARCHAR,
    gene_symbol     VARCHAR NOT NULL,
    variant_key     VARCHAR NOT NULL,   -- CHROM:POS:REF:ALT — the dedup key (§3.2)
    chrom           VARCHAR,
    pos             BIGINT,
    ref_allele      VARCHAR,
    alt_allele      VARCHAR,
    hgvs_c          VARCHAR,
    hgvs_p          VARCHAR,
    aa_pos          INTEGER,
    protein_len     INTEGER,
    var_class       VARCHAR,            -- SNV | Indel | CNV | Splice | SV
    consequence     VARCHAR,            -- Missense | Nonsense | Frameshift | ...
    zygosity_raw    VARCHAR,            -- as read from VariMAT: Het/Hom ONLY
    zygosity        VARCHAR,            -- DERIVED (§3.5): + Hemizygous / Compound het
    vaf             DOUBLE,             -- 0..1
    depth           INTEGER,
    alt_depth       INTEGER,
    filter_status   VARCHAR,
    gnomad_af       DOUBLE,
    gnomad_sas_af   DOUBLE,
    ga100k_sas_af   DOUBLE,
    clinvar_sig     VARCHAR,
    clinvar_id      VARCHAR,
    mane            BOOLEAN,
    reviewable      BOOLEAN NOT NULL    -- §3.2 reviewable subset
);

CREATE TABLE IF NOT EXISTS interpretation (
    finding_id      VARCHAR PRIMARY KEY,
    framework       VARCHAR NOT NULL,   -- ACMG/AMP
    classification  VARCHAR NOT NULL,   -- Pathogenic | Likely pathogenic | ...
    acmg_codes      VARCHAR,
    kb_snapshot_id  VARCHAR NOT NULL,
    curated_flag    BOOLEAN NOT NULL,
    interpreted_at  DATE,
    reportable      BOOLEAN NOT NULL,   -- reportable under the run's test code
    inherited_from  VARCHAR,
    segregation     VARCHAR,
    evidence_delta  DOUBLE              -- new public evidence since interpretation
);

-- ------------------------------------------------------------- reproducible --
CREATE TABLE IF NOT EXISTS cohort_def (
    cohort_id     VARCHAR PRIMARY KEY,
    tenant_id     INTEGER NOT NULL,
    profile       VARCHAR NOT NULL,
    name          VARCHAR NOT NULL,
    description   VARCHAR,
    criteria_json VARCHAR NOT NULL,
    version       INTEGER NOT NULL,
    created_by    VARCHAR,
    created_at    TIMESTAMP,
    builtin       BOOLEAN DEFAULT FALSE
);

CREATE TABLE IF NOT EXISTS cohort_snapshot (
    snapshot_id    VARCHAR PRIMARY KEY,
    cohort_id      VARCHAR NOT NULL,
    resolved_at    TIMESTAMP NOT NULL,
    member_hash    VARCHAR NOT NULL,
    criteria_hash  VARCHAR NOT NULL,
    kb_snapshot_id VARCHAR NOT NULL,
    member_count   INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS analysis_run (
    analysis_id  VARCHAR PRIMARY KEY,
    snapshot_id  VARCHAR,
    module       VARCHAR NOT NULL,
    output_class VARCHAR NOT NULL,
    executed_at  TIMESTAMP NOT NULL,
    user_id      VARCHAR,
    cohort_size  INTEGER,
    criteria_hash VARCHAR
);

-- Query-compiler corpus (spec D01: "Every compilation logged with its source
-- question" — builds the evaluation set).
CREATE TABLE IF NOT EXISTS compilation_log (
    compilation_id VARCHAR PRIMARY KEY,
    asked_at       TIMESTAMP NOT NULL,
    question       VARCHAR NOT NULL,
    profile        VARCHAR NOT NULL,
    criteria_json  VARCHAR,
    target_module  VARCHAR,
    validation     VARCHAR NOT NULL,   -- PASSED | FAILED
    unmapped       VARCHAR,
    executed       BOOLEAN DEFAULT FALSE
);

-- Store-level metadata (KB snapshot pin, build info).
CREATE TABLE IF NOT EXISTS store_meta (
    key   VARCHAR PRIMARY KEY,
    value VARCHAR
);
