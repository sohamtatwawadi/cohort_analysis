-- ============================================================================
-- Research Mode schema — Part II (R1)
--
-- Lab Mode data is derived from our own pipeline and is rebuildable. Research
-- Mode data is UPLOADED BY SOMEONE ELSE, which changes the obligations
-- entirely (Part II §8): it belongs to the uploader, is scoped to a project,
-- is never pooled, and must be deletable with the deletion propagating to
-- everything derived from it.
--
-- That is why datasets, jobs and results all carry project_id — isolation is a
-- column on every row, not a filter applied in the application layer.
-- ============================================================================

CREATE TABLE IF NOT EXISTS project (
    project_id      VARCHAR PRIMARY KEY,
    tenant_id       INTEGER NOT NULL,
    name            VARCHAR NOT NULL,
    description     VARCHAR,
    owner           VARCHAR NOT NULL,
    created_at      TIMESTAMP NOT NULL,
    retention_days  INTEGER,          -- §8 retention, enforced automatically
    deleted_at      TIMESTAMP
);

CREATE TABLE IF NOT EXISTS dataset (
    dataset_id        VARCHAR PRIMARY KEY,
    project_id        VARCHAR NOT NULL,
    name              VARCHAR NOT NULL,
    mode              VARCHAR NOT NULL,   -- 'research' | 'lab'
    source_format     VARCHAR,            -- VCF | PLINK1 | VariMAT | ...
    source_files      VARCHAR,            -- JSON list
    genome_build      VARCHAR,            -- NEVER inferred from position (§2.3)
    n_samples         INTEGER,
    n_variants        INTEGER,
    uploaded_at       TIMESTAMP NOT NULL,
    uploaded_by       VARCHAR,
    status            VARCHAR NOT NULL,   -- validating|ready|failed|deleted
    failure_reason    VARCHAR,

    -- §2.4 Annotation provenance. Pinned to the dataset; releases are never
    -- silently mixed, so these are recorded even when the uploader supplied
    -- pre-annotated input.
    annotation_source   VARCHAR,
    annotation_version  VARCHAR,
    transcript_set      VARCHAR,

    -- §8 Consent attestation. The uploader asserts they hold consent and
    -- ethics approval; we record who asserted what and when.
    consent_attested    BOOLEAN NOT NULL DEFAULT FALSE,
    consent_statement   VARCHAR,
    consent_attested_by VARCHAR,
    consent_attested_at TIMESTAMP,

    storage_path      VARCHAR,
    checksum          VARCHAR,
    deleted_at        TIMESTAMP
);

-- Sample roster per dataset, with the reconciliation outcome against phenotypes.
CREATE TABLE IF NOT EXISTS dataset_sample (
    dataset_id     VARCHAR NOT NULL,
    sample_id      VARCHAR NOT NULL,
    reported_sex   VARCHAR,
    genetic_sex    VARCHAR,       -- derived; discordance is a §2.3 hard check
    call_rate      DOUBLE,
    in_phenotype   BOOLEAN,
    is_duplicate   BOOLEAN,
    PRIMARY KEY (dataset_id, sample_id)
);

-- §3.1 Data profile. Stored whole as JSON because the capability matrix reads
-- it as a document and the field set will grow.
CREATE TABLE IF NOT EXISTS dataset_profile (
    dataset_id    VARCHAR PRIMARY KEY,
    profiled_at   TIMESTAMP NOT NULL,
    profile_json  VARCHAR NOT NULL
);

-- §3.2/§3.3 Capability state. Persisted rather than recomputed on each request
-- so that what a researcher was shown is auditable after the fact.
CREATE TABLE IF NOT EXISTS dataset_capability (
    dataset_id   VARCHAR NOT NULL,
    analysis     VARCHAR NOT NULL,
    available    BOOLEAN NOT NULL,
    reasons_json VARCHAR,          -- unmet requirements, each with the observed value
    assessed_at  TIMESTAMP NOT NULL,
    PRIMARY KEY (dataset_id, analysis)
);

-- §3.4 Override with justification. The stamp that results carry is derived
-- from the existence of a row here — it cannot be removed by editing a result.
CREATE TABLE IF NOT EXISTS capability_override (
    override_id   VARCHAR PRIMARY KEY,
    dataset_id    VARCHAR NOT NULL,
    analysis      VARCHAR NOT NULL,
    justification VARCHAR NOT NULL,
    requested_by  VARCHAR NOT NULL,
    created_at    TIMESTAMP NOT NULL,
    unmet_json    VARCHAR            -- what was unmet AT THE TIME of override
);

-- §4.4 Variant sets — "Saved as a reusable, versioned variant set."
CREATE TABLE IF NOT EXISTS variant_set (
    set_id          VARCHAR PRIMARY KEY,
    project_id      VARCHAR NOT NULL,
    name            VARCHAR NOT NULL,
    definition_json VARCHAR NOT NULL,
    version         INTEGER NOT NULL,
    created_by      VARCHAR,
    created_at      TIMESTAMP NOT NULL
);

-- §6 Compute. Research analyses cannot run in a web request.
CREATE TABLE IF NOT EXISTS analysis_job (
    job_id         VARCHAR PRIMARY KEY,
    project_id     VARCHAR NOT NULL,
    dataset_id     VARCHAR NOT NULL,
    analysis       VARCHAR NOT NULL,
    spec_json      VARCHAR NOT NULL,
    status         VARCHAR NOT NULL,   -- queued|validating|running|post|complete|failed|cancelled
    submitted_at   TIMESTAMP NOT NULL,
    started_at     TIMESTAMP,
    finished_at    TIMESTAMP,
    submitted_by   VARCHAR,
    error          VARCHAR,
    log            VARCHAR,
    -- §3.4: set when the analysis ran under an override. Every output and
    -- export derived from this job carries the stamp.
    underpowered   BOOLEAN NOT NULL DEFAULT FALSE,
    override_id    VARCHAR,
    -- §6 Reproducibility: software versions, seed, parameters, dataset release.
    reproducibility_json VARCHAR
);

CREATE TABLE IF NOT EXISTS analysis_result (
    job_id       VARCHAR PRIMARY KEY,
    analysis     VARCHAR NOT NULL,
    result_json  VARCHAR NOT NULL,
    diagnostics_json VARCHAR,
    created_at   TIMESTAMP NOT NULL,
    output_class VARCHAR NOT NULL DEFAULT 'RESEARCH'
);

-- Phenotype columns registered per dataset, with the type that drives model
-- selection (§4.1) and the counts the capability matrix gates on (§3.2).
CREATE TABLE IF NOT EXISTS dataset_phenotype (
    dataset_id    VARCHAR NOT NULL,
    name          VARCHAR NOT NULL,
    kind          VARCHAR NOT NULL,     -- binary|quantitative|categorical|time_to_event
    label         VARCHAR,
    n_present     INTEGER,
    n_cases       INTEGER,
    n_controls    INTEGER,
    completeness  DOUBLE,
    PRIMARY KEY (dataset_id, name)
);
