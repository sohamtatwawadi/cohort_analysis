# ImpactOmics Cohort — Platform Documentation

Current state of the build, as of the latest working session.

Built from two internal specifications, which are not distributed with this
repository.

| Spec | What it covers | Status |
|---|---|---|
| Cohort Analytics Developer Spec | Part B (germline), governance, query compiler | Built |
| Germline Platform Part II — Research Mode | R1, R2, R4, R5, R6 | Built (R3, R7, R8 not) |

The somatic profile (S01–S14) is deliberately **not** built.

**Scale:** ~16,200 lines of backend Python, ~2,400 lines of frontend JavaScript,
~4,400 lines of tests. 316 tests passing. 52 API endpoints. 27 database tables.

---

## 1. Running it

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/uvicorn backend.app.main:app --port 8077     # http://localhost:8077
```

The server seeds a synthetic store on first start if the database is empty, so
there is nothing else to do to see it working.

```bash
.venv/bin/python -m backend.cli synthetic      # rebuild lab demo data
.venv/bin/python -m backend.tools.make_research_fixture --samples 2400 --variants 40000
```

The research fixture plants six causal variants at OR 1.9, two populations with
differing allele frequencies, 30 parent-offspring pairs, 40 complete trios and
per-variant gene annotations — so every analysis has a known right answer to
check against.

#### The demo dataset

One command builds a cohort that unlocks **every implemented analysis**:

```bash
.venv/bin/python -m backend.tools.make_research_fixture \
  --samples 3000 --variants 150000 --seed 20260915 \
  --name "Cardiomyopathy case-control · demo cohort"
```

Stop the server first — DuckDB takes a single writer lock, so the generator
cannot run while the app holds the database.

What each part of the fixture is for:

| Ingredient | Unlocks |
|---|---|
| 150,000 variants over 23 chromosomes | `genome_wide` density → GWAS, PRS |
| ~45% rare (MAF 0.0004–0.0012), LoF-annotated, `lof_confidence: HC` | Burden / SKAT / SKAT-O |
| Balanced case/control on `affected` | Association, burden, GWAS |
| `followup_years` + `died` | Survival / penetrance |
| Two populations at Fst 0.05, 10 PCs | Population frequency, ancestry adjustment |
| 40 complete trios with `father_id` / `mother_id` | Segregation |
| Gene annotations in contiguous 60-variant blocks | Gene-based tests, association by gene |

Frequencies are deliberately rare-but-*observed* — 2 to 8 copies in the cohort.
Drawing them any lower produces monomorphic sites that are present in the file
and absent from every sample, which pass no filter and make a gene-based test
look broken when it is behaving correctly.

### 1.1 Getting your own data in

**Research Mode — from the browser.** *Uploaded data → Register a dataset*.

| Field | Accepts |
|---|---|
| Genotype | VCF / VCF.gz (one file) · PLINK1 (select `.bed` `.bim` `.fam` **together**) |
| Phenotype | CSV/TSV keyed by `sample_id`. Optional — without it only descriptive analyses unlock |
| Genome build | Explicit, never inferred. GRCh37 and GRCh38 positions overlap, so a wrong guess silently corrupts every annotation |
| Consent | Attestation, enforced server-side before a byte is parsed |

Files post as multipart to `POST /api/research/projects/{id}/datasets/upload`,
land in `data/uploads/`, and go straight through validation and profiling — a
dataset that fails validation is rejected rather than stored in a broken state.
A collapsed *"the file is already on the server"* panel keeps the older
path-based route (`/datasets/ingest`) for files too large to stream through a
browser.

**Lab Mode — CLI only.** VariMAT ingest has no upload UI (see §12):

```bash
.venv/bin/python -m backend.cli ingest data/varimat --clinical data/clinical
```

It reads `*.tsv|txt|csv` (optionally `.gz`) recursively, requiring `CHROM`,
`START`, `REF`, `ALT`, `VARCLASS` and `VARIANT_FILTER_STATUS`; sample, pipeline,
caller and build come from the filename. The clinical sidecar is a 23-column
CSV/JSON keyed by `sample_id`.

---

## 2. The workflow

The UI is one spine, four steps, used by both data sources. This replaced an
earlier build that put 23 destinations in a left rail and let the user work out
the order.

```
 ① Build cohort  →  ② Review  →  ③ Analyse  →  ④ Share
   Who is included   Can you       Ask your     Export with
                     trust it      question     provenance
```

### ① Build cohort
Opens on real starting points ("Everyone", "HBOC panel · BRCA1/2 P/LP carriers",
"Biallelic P/LP · recessive genes", …) rather than a blank form, because almost
nobody builds a cohort from nothing. Filters are grouped by the question they
answer — *Who is included* / *How they were tested* / *What findings they carry*
— with everything most users never touch behind an **Advanced** disclosure.

The right-hand panel shows the live subject count and a funnel attributing every
excluded subject to exactly one choice.

### ② Review
Four trust checks, each pass/warn and expandable in place:

| Check | Asks |
|---|---|
| **Sample quality** | Is any sample contaminated, swapped or under-covered? |
| **Independence** | Is it one subject per family, or are relatives inflating the rates? |
| **Coverage** | Does every gene have a real denominator, or are some provisional? |
| **Provenance** | Is the cohort mixing reference builds or pipeline versions? |

Sample quality comes first because it invalidates everything under it: a
contaminated sample contributes spurious heterozygous calls to every rate the
other three checks are busy validating. See §5.1.

The coverage check embeds the full per-gene denominator table (assayed /
provisional / not assayed / basis). In the previous build this lived in a modal
behind a rail link, and the warnings were three banners users learned to scroll
past.

### ③ Analyse
Eleven internal modules are presented as **eight plain-English questions**. The
module name never appears in the UI.

The menu is **ranked against the cohort you actually built** (see §6), with the
reason stated in counts, and split into three groups: *worth asking*, *will run
but this cohort was not built for it*, *would return nothing here*.

Drill-down (gene → variants → subject) is reached by clicking a number, so the
evidence is always one click from the figure it supports.

### ④ Share
Export with the manifest shown as part of the export, not hidden behind a
governance menu. Formats: XLSX (data + manifest on a second sheet), CSV
(manifest as header comments), JSON.

**Ask a question (⌘K)** is available throughout: plain English → cohort filters,
never a number. You approve the compiled filters before anything runs.

---

## 3. Lab Mode — the eight questions

| Question shown to the user | Module | What it computes |
|---|---|---|
| How many subjects got a diagnosis? | `g_carrier` | Diagnostic yield with the phenotype-relevance gate |
| Which genes carry pathogenic variants? | `g_carrier` | Carrier rate per gene on gene-specific denominators |
| Are they carriers, or affected? | `g_zygosity` | Het / hom / compound-het / hemizygous, never summed |
| How strong is the gene–disease evidence? | `g_genedisease` | Yield inflation from weak-validity genes |
| Which uncertain variants should we review first? | `g_vus` | Ranked curation queue |
| Is any variant unusually common here? | `g_popfreq` | Internal AF vs gnomAD, founder candidates |
| Are there secondary findings to return? | `g_sf` | ACMG SF v3.3, consented subjects only |
| Which phenotypes go with which genes? | `g_phenotype` | Gene × HPO cross-tab and review queue |

Drill targets: `g_gene` (L4), `g_variants` (L5), `g_subjects` (L2), `g_runs` (L3).

---

## 4. Research Mode — capability gating

Uploaded datasets. The spec calls capability gating *"the most important safety
feature in the product"*, and three properties are enforced structurally:

**Locked is visible, not hidden.** Every unavailable analysis renders with each
unmet requirement *and the observed value* — "Requires n ≥ 2,000 / your dataset
has 486 subjects".

**Locked cannot be clicked through.** A locked card has no run control at all,
not a disabled one. `jobs.submit` calls `registry.can_run` before queuing, so
an analysis cannot start even from a direct API call.

**Override is separate, recorded, and stamps the output.** Borderline thresholds
can be overridden with a typed justification (≥20 chars, audited). Missing data
cannot — you cannot override your way to genotypes you did not upload. The
`UNDERPOWERED — threshold overridden` stamp is applied by the job runner from
the override record, so it cannot be edited off a result.

### Analyses implemented

| Analysis | Phase | Notes |
|---|---|---|
| Carrier / allele frequency | R1 | Called-sample denominator per variant; carrier model is a choice, not a default |
| Zygosity / inheritance | R1 | Per-variant counts and a per-sample het/hom ratio |
| Population frequency | R1 | Internal frequency by group — explicitly not a gnomAD substitute |
| Diagnostic yield | R1 | Candidate rate, stated as an upper bound: no ACMG gate, no relevance gate |
| Segregation / de novo | R1 | Mendelian consistency over complete trios; reports inconsistencies, never "de novo" |
| Association (single variant) | R2 | Linear / logistic / **Firth** / robust, auto model selection; targeted by gene or variant key |
| GWAS | R4 | Mandatory QC, PCA, relatedness pruning, λ_GC guardrail. Two-stage: score scan then exact refit |
| Burden / SKAT / SKAT-O | R5 | Davies exact p-values, min-carrier guardrail |
| PRS | R6 | Ancestry-stratified performance always computed |
| Survival / penetrance | R6 | KM, log-rank, Cox (Efron), PH diagnostics, ascertainment label |

Each analysis has **its own configuration form** — association asks which
variants and an outcome, survival asks for time and event columns, a gene-based
test asks for a qualifying-variant definition, PRS asks for score weights.

The five R1 analyses were advertised by the capability matrix long before they
existed: the cards rendered, turned green on any genotype data, and did nothing
when clicked. `tests/test_analysis_coverage.py` now asserts the three lists —
what the matrix can offer, what the job registry can run, and what the frontend
has a form and a renderer for — agree in both directions.

#### The GWAS scan is two-stage

A scan used to fit a full GLM per variant. The covariate block is identical every
time, so an 82,000-variant scan repeated the same expensive iteration 82,000
times: over half an hour on 3,000 samples.

It now fits the null model **once** and scores every variant against its
residuals (`stats/scoretest.py`), then refits only the top hits — everything
below p < 10⁻⁴ plus the top 500 — with the exact model, because an odds ratio and
a confidence interval must come from a fitted model rather than a score
statistic. Rows that were scanned but not refitted say so in `estimate_source`
instead of carrying an invented effect size.

Score and likelihood-ratio statistics are asymptotically equivalent under the
null, so the ranking is unchanged; `tests/test_scoretest.py` pins that with a
Spearman correlation against the exact fit, a KS test for uniformity of the null
p-values, and a λ_GC check.

### Statistical guardrails

- **Power before execution** (§5.1), with the spec's verdict text and options.
- **Firth where it is required** — model selection routes to it on sparse 2×2
  cells and says why. P-values from the penalised likelihood-ratio test, not
  Wald, because Wald is unreliable in exactly the separation regime Firth exists
  for.
- **λ_GC guardrail with attributed cause** — names unmodelled ancestry, cryptic
  relatedness, batch or differential missingness from the data profile.
- **Mandatory pre-GWAS QC** with no parameter to skip it.
- **Too-few-carrier genes reported as untested**, not given an unstable p-value.
- **Multiple testing always** — raw p, Bonferroni and BH-FDR on every result.

---

## 5. Safety invariants

These are enforced in the code path rather than by convention, and each has a
test that tries to break it.

| Invariant | Where it lives | Why |
|---|---|---|
| Withdrawn consent is unreachable | `cohort.resolve` compiles it as an unconditional SQL predicate | No criteria object, valid or malformed, can surface a withdrawn subject |
| Secondary findings need explicit consent | `cohort.resolve` step 7, plus `SECONDARY_FINDING_SQL` | Display filtering leaves data one bug from exposure |
| Every gene has its own denominator | `denominator.gene_denominators` | Cohort size as denominator makes rates quietly, consistently too low — and quoted |
| Hemizygous / compound-het are computed | `ingest.store.derive_zygosity` | `ZYGOSITY` carries Het/Hom only; reading it directly misclassifies X-linked males as carriers |
| A percentage always carries its counts | `denominator.rate` + `kit.rate` | There is deliberately no helper that formats a bare percentage |
| Yield is phenotype-gated | `carrier.diagnostic_yield` + `indication_gene` table | Without it a broad exome counts any P/LP as diagnostic |
| Deletion propagates | `registry.delete_dataset` | A flag on `dataset` would leave the profile, results and genotype payload in place |
| Project isolation | `project_id` column on every dataset, job and result | Isolation is a column, not a filter a caller can forget |
| Cohort resolution is serialised | `cohort.COHORT_LOCK` | Temp tables are connection-scoped; concurrent requests could return another cohort's numbers |
| One compilation never produces a number | `compiler.validate` rejects payloads containing figures | A wrong number sounds authoritative and is unverifiable |
| QC ratios are computed over the full call set | `varimat.compute_qc`, stored on `run` | `finding` keeps only the reviewable ~2%; a Ti/Tv from that is noise, not a metric |

### 5.1 Sample QC

A contaminated sample or a sample swap invalidates every rate computed over it,
so this runs before analysis rather than after a result looks strange.

| Metric | Catches |
|---|---|
| **Ti/Tv** | Contamination and false-positive inflation. Real calls are transition-biased (~2.0 exome, ~3.0 panel); noise tends toward the random 0.5 |
| **Het/Hom** | Contamination and sample mixture — two genomes make one look heterozygous everywhere |
| **Mean het VAF** | Allele balance. A true heterozygote centres on 0.5; contamination pulls it off-centre |
| **Mean depth / % ≥20×** | Under-coverage, which produces false negatives rather than false positives |
| **X het rate** | Sex concordance — a male sample heterozygous across X is mislabelled or mixed |
| **Pipeline `qc_status`** | The lab's own verdict, which outranks anything derived here |

**Thresholds are cohort-relative, not absolute.** A panel, an exome and a genome
have genuinely different expected values, so a fixed cut-off would fail an entire
assay type for being itself. Outliers are measured in median absolute deviations
from the cohort median — warn at 3, fail at 5 — with MAD rather than standard
deviation because a handful of extreme values inflate an SD enough to hide
themselves. Two absolute floors still apply, for cases a uniformly bad cohort
would otherwise normalise: Ti/Tv below 1.0 and % ≥20× below 80.

Every flag states the observed value and what it was compared against; "this
sample looks odd" is not actionable.

**Where the numbers come from.** `compute_qc` runs in the VariMAT loader over
every PASS call *before* the reviewable filter, and the results are stored as
`qc_*` columns on `run`. This was a loader change, not a query: §3.2 archives
~98% of rows after fingerprinting and `finding` retains only the reviewable
handful, so computing a ratio from stored findings gives a number derived from a
few dozen calls — 597 of 606 samples came back "not assessable" when tried that
way. The service prefers the stored per-run metrics, falls back to findings only
when there are enough of them, and **always reports `metrics_basis`** so a
panel-sized ratio can never be read as a genome-wide one.

---

## 6. Analysis suggestions

`backend/app/services/suggest.py`

Rather than a flat menu, the cohort is profiled and the eight questions ranked
against what it actually contains, with the reason given in counts.

Signals computed in one pass: P/LP subjects and genes, VUS variants, biallelic
and hemizygous counts, weak-validity findings, conditions implicated, recurrent
variants, SF-eligible subjects, ancestries, indications, HPO terms,
multi-member families.

Three verdicts — `recommended`, `available`, `not_useful`. Nothing is *locked*:
a lab cohort can always run any descriptive analysis, so this is advice, not a
gate.

Example — selecting Pathogenic + Likely pathogenic + Homozygous + Compound
heterozygous + AR + AR/AD produces:

> Based on what this cohort contains, start with: zygosity and inheritance;
> gene–disease evidence; diagnostic yield.
>
> **Are they carriers, or affected?** → 26 biallelic and 3 hemizygous — these
> are affected individuals, not carriers, and the distinction matters clinically.

---

## 7. Architecture

```
backend/app/
  config.py              thresholds, KB snapshot pin, paths
  db.py                  DuckDB connection, query helpers
  schema.sql             lab analytics schema (spec §3.6)
  research_schema.sql    research mode: projects, datasets, jobs, capabilities

  reference/             curated tables that are NOT in VariMAT (spec §3.5)
    genes.py             43 gene-disease records — SEED TABLE, replace with ClinGen
    tests.py             10 test codes, phenotype-relevance map, vocabularies

  ingest/                lab data path
    varimat.py           loader: dedup on locus, MANE canonical, reviewable subset
    fingerprint.py       coverage inference by containment clustering (§3.4)
    clinical.py          LIMS sidecar — consent, family, indication
    store.py             persistence + derivation passes (scope, zygosity)
    synthetic.py         deterministic demo generator

  services/              lab analytics — every number is SQL
    cohort.py            the §G01 gating pipeline, COHORT_LOCK
    denominator.py       gene-specific denominators, coverage inspector
    sampleqc.py          per-sample QC, cohort-relative MAD outliers (§5.1)
    carrier.py           diagnostic yield, carrier rates, dashboard
    zygosity.py          §G04
    analysis.py          §G05–G09
    explore.py           §G10–G13 drill levels
    suggest.py           analysis recommendation
    governance.py        output classes, manifest, audit, export
    library.py           saved cohorts
    compiler.py          query compiler (D01) — contract + validator

  research/              uploaded data path
    types.py             GenotypeMatrix / PhenotypeTable — the mode boundary
    formats/             VCF, PLINK1, format detection
    validate.py          §2.3 mandatory validations (fail the upload, never warn)
    profile.py           §3.1 data profile
    capability.py        §3.2–3.4 the capability matrix
    registry.py          dataset registry, project isolation, deletion
    store.py             genotype/phenotype/annotation payload on disk
    jobs.py              job queue, lifecycle, reproducibility record
    variantset.py        qualifying-variant builder, versioned
    stats/               glm (incl. Firth), pca, kinship, qc, skat, survival, power
    analyses/            association, gwas, burden, prs, survival

  api/
    routes.py            lab endpoints
    research_routes.py   research endpoints

frontend/
  index.html             shell
  app.css                one stylesheet
  js/
    flow.js              shell, stepper, router, action dispatcher
    build.js             step 1
    review.js            step 2
    analyse.js           step 3 — questions, answers, drill, research analyses
    share.js             step 4
    kit.js               shared rendering helpers
    api.js               API client

backend/tools/
  make_varimat_fixtures.py    VariMAT test files with realistic dedup pressure
  make_research_fixture.py    research-scale cohort with planted truth
```

### Data flow

```
VariMAT files  ─┐   CLI only
clinical CSV   ─┴─→ ingest ─→ derived store ─→ cohort engine ─→ services ─→ API ─→ UI
                      │            ↑
                   compute_qc ─────┘  (qc_* columns on `run`)

VCF / PLINK ─→ upload ─→ validate ─→ profile ─→ capability matrix ─→ jobs ─→ analyses
  (browser)                  │
                       reject on failure, never store a broken dataset
```

The analytics store is **derived and rebuildable**; it is never the system of
record. `cli ingest` drops and rebuilds it.

---

## 8. Database

27 tables. The ones that matter:

**Lab** — `subject`, `family`, `sample`, `run`, `run_scope` (the denominator
table), `run_fingerprint`, `finding`, `interpretation`, `gene_disease`,
`test_code`, `test_code_gene`, `indication_gene` (the yield relevance gate),
`cohort_def`, `cohort_snapshot`, `analysis_run`, `compilation_log`, `store_meta`.

`run` also carries five per-sample QC columns written by the loader —
`qc_n_called`, `qc_ti_tv`, `qc_het_hom`, `qc_mean_het_vaf`, `qc_x_het_rate` —
alongside `mean_depth`, `pct_bases_20x` and the pipeline's own `qc_status`.
They are computed over the **full PASS call set** at ingest because that
population does not survive into `finding` (§5.1).

**Research** — `project`, `dataset`, `dataset_sample`, `dataset_profile`,
`dataset_capability`, `dataset_phenotype`, `capability_override`, `variant_set`,
`analysis_job`, `analysis_result`.

Genotype matrices live on disk as compressed numpy under
`data/research/<dataset_id>/`, not in DuckDB — a matrix of variants × samples
would be tens of millions of rows for a few megabytes of array.

---

## 9. API

52 endpoints. Grouped:

**Lab cohort** — `POST /api/cohort/resolve`, `POST /api/cohort/suggestions`,
`GET|POST|DELETE /api/cohorts`, `GET /api/cohorts/{id}/verify`

**Lab modules** — `POST /api/module/{g_dashboard|g_carrier|g_zygosity|
g_genedisease|g_phenotype|g_popfreq|g_vus|g_sf|g_gene|g_variants|g_subjects|
g_runs}`

**Drill** — `POST /api/detail/variant/{id}`, `POST /api/detail/subject/{id}`

**Governance** — `POST /api/governance/{denominator|manifest|export}`,
`GET /api/governance/audit`

**Query compiler** — `POST /api/compile`, `GET /api/compile/vocabulary`,
`GET /api/compile/corpus`, `POST /api/compile/{id}/execute`

**Research** — `GET /api/research/meta`, projects, datasets (ingest / upload /
capabilities / override / delete), jobs (submit / poll / result / cancel),
`POST /api/research/power`, `GET /api/research/estimate`, variant sets

`datasets/upload` is multipart and takes a **list** of genotype files, because a
PLINK fileset is only readable as a set; `datasets/ingest` takes server-side
paths instead. Both converge on the same validation and profiling path.

**Admin** — `POST /api/admin/rebuild/{synthetic|varimat}`,
`GET /api/admin/varimat/discover`

Every module endpoint returns `{module, module_title, output_class, cohort, data}`
and writes an audit entry.

---

## 10. CLI

```bash
python -m backend.cli synthetic [--families N] [--seed S]   # rebuild demo store
python -m backend.cli ingest <varimat-dir> [--clinical <dir>]
python -m backend.cli inspect <varimat-file>    # dedup funnel, writes nothing
python -m backend.cli template <path.csv>       # clinical sidecar header
python -m backend.cli status                    # what is loaded
python -m backend.cli verify                    # re-resolve every saved cohort
```

`inspect` prints the §3.2 funnel per file — raw rows → distinct variants (and
the ratio) → PASS → coding → rare → protein-altering — so the ~1.75× dedup
inflation is visible rather than assumed.

`verify` re-resolves every saved cohort and compares member hashes. A saved
cohort stores criteria, never a member list.

---

## 11. Tests

```bash
.venv/bin/python -m pytest tests/ -q      # 316 passing
```

| File | Tests | Covers |
|---|---|---|
| `test_hand_count.py` | 6 | **Independent hand count (spec E03.9)** |
| `test_cohort_engine.py` | 9 | Lock contract under concurrency, funnel attribution, pinned anchor |
| `test_denominators.py` | 12 | Gene denominators, family independence, derived zygosity |
| `test_sample_qc.py` | 12 | Ti/Tv and het/hom maths, MAD outliers, metrics basis |
| `test_governance_gates.py` | 8 | Consent gates, SF enforcement, small-cell suppression |
| `test_upload.py` | 7 | Multipart upload, PLINK filesets, consent, path traversal |
| `test_loader.py` | 26 | VariMAT dedup, reviewable subset, fingerprinting, QC |
| `test_compiler_and_repro.py` | 24 | Query compiler contract, reproducibility, manifest |
| `test_research_gating.py` | 21 | Upload validation, capability matrix, override, deletion |
| `test_research_formats.py` | 39 | VCF and PLINK readers |
| `test_research_glm.py` | 56 | Linear / logistic / Firth / HC3 / BH |
| `test_research_genetics.py` | 35 | PCA, KING kinship, HWE exact, sex check |
| `test_research_skat_survival.py` | 61 | Davies, SKAT calibration, KM, Cox, competing risks |

`test_hand_count.py` is the one that matters. It recomputes cohort membership,
gene denominators, carrier rates and diagnostic yield from the raw tables in
plain Python, importing nothing from the services package — so a shared helper
cannot make both sides wrong in the same way. Mutating the denominator service
to use cohort size fails it with `ALDOB: hand 69 vs tool 189`.

---

## 12. Known gaps

**Reference data**
- `backend/app/reference/genes.py` is a **seed table**. The 43 gene-disease
  records were hand-entered to match the prototype. Validity and penetrance
  calls must be replaced with a dated ClinGen export, and `gene_id` with an HGNC
  export, before any figure leaves the building.

**Research Mode analyses**
- The four permanently-locked analyses — PheWAS, colocalization, fine-mapping,
  heritability — declare a hardcoded unmet requirement, so no dataset can turn
  them green. That is deliberate: a card that can never unlock is honest, one
  that unlocks and then does nothing is not.
- `association` with no gene or variant target scans the whole dataset one model
  at a time and refuses above 5,000 variants, pointing at GWAS instead.

**Ingest surface**
- **Lab Mode has no upload UI.** VariMAT ingest is CLI-only (§1.1). The work is
  not a modal: ingest is multi-file, needs the clinical sidecar alongside it, and
  runs fingerprint clustering over the result, so it wants a progress view rather
  than a spinner.
- **No database connector.** Both modes read files. `ingest.store.persist(bundle)`
  — which takes dicts keyed by table name — is the seam a LIMS or warehouse
  connector would write to; nothing else needs to change to feed the engine from
  a live source.

**Query compiler**
- Rule-based, not an LLM. `compiler.py` implements the D01 contract and the
  validator, which fails closed. Wiring a model in means returning the same
  `{profile, criteria, target_module}` shape through `validate`. The
  100-question acceptance gate has not been run.

**Research Mode**
- **Engine is in-process numpy/scipy**, not REGENIE/SAIGE/PLINK2. Correct and
  tested at this scale (thousands of samples); **not biobank-scale** — an exact
  in-memory SVD and an O(n²) kinship matrix are the binding limits. The job
  contract is the seam to swap at.
- **Not built:** R3 (PheWAS — needs an ICD→PheCode ontology), R7 (fine-mapping,
  colocalization, heritability — need LD reference panels), R8.
- **Formats not supported:** BCF (listed P1 in the spec), PLINK2/PGEN, BGEN,
  Hail MatrixTable.
- **SKAT-O's correction is approximate** — Šidák on Galwey's effective number of
  tests, accurate to roughly a factor of 2, and the result says so. The exact
  Lee et al. integration is not implemented.
- **Left-alignment is trim-based parsimony**, not true left-shifting through
  repeats, which needs a reference FASTA the loader has no handle on.
- **Gene-based tests need gene annotations** supplied with the dataset. Without
  them the capability matrix locks the analysis rather than failing at runtime.

**Discrepancies found in the specs**
- The spec says the prototype models 44 genes; it models 43.
- The §5.1 power example does not reproduce: its stated 11% at n=4,200 / 180
  cases / MAF 0.008 / OR 2.0 / α=5×10⁻⁸ would require OR≈4.0 or α≈0.005. The
  implementation matches the §7 example exactly (94.1% vs 94%) and gives 72.8%
  on a textbook design, so the §5.1 figure appears illustrative rather than
  computed. Worth confirming before anyone calibrates against it.

**Demo data**
- Synthetic calibration is approximate. Overall yield lands near the spec's
  41.5% and the SF rate in the low single digits, but per-test-code yields vary
  more than a real lab would. This affects demo data only, not the computation.
- CNV and fusion show counts only, never a percentage, because absence of a CNV
  call is not evidence of absence and coverage for them cannot be inferred.

---

## 13. Notable bugs found and fixed

Kept because each one is a trap that would recur.

| Bug | Consequence if unfixed |
|---|---|
| Davies' method transcribed `log1(-x)` as `log1p(-x)` | SKAT p-values off by ~15,000× (0.021 vs 1.38e-6) — "nothing here" vs genome-wide significant |
| Genotype matrix passed to SKAT variants × samples | Silently transposed for any gene with as many variants as samples |
| Blocking association on zero cell counts | Refused exactly the analyses Firth exists to handle |
| SF numerator counted subjects whose test was not SF-capable | Rate could exceed 100%; MYBPC3 on a cardiac panel counted as "secondary" |
| `resolve()` racing on shared DuckDB temp tables | Concurrent requests could return another cohort's numbers, silently |
| Step buttons rendered before the cohort resolved | Users could never leave step 1 |
| `import('./flow.js')` without the cache-busting query | Second module instance with its own state — clicks mutated an unrendered state |
| `{...blank_criteria}` shallow copy | "Clear all" and every compiled question inherited stale filters |
| One shared analysis config form | Three of five analyses failed the moment anyone ran them from the UI |
| Burden error said "no annotations" when a filter had removed everything | Sent users to re-annotate when the real fix was a threshold |
| Five analyses advertised by the capability matrix were never implemented | Cards turned green, opened no form, logged nothing — the analysis did not exist |
| GWAS read `QCResult.keep_mask`, which does not exist | A `hasattr` fallback tried subscripting instead of failing, so a typo surfaced as "'QCResult' object is not subscriptable" |
| GWAS fitted a full GLM per variant | Half an hour for a scan whose covariate block is identical every time |
| `association` documented a `gene` key it never implemented | Asking for one gene silently tested every variant in the cohort |
| Upload endpoint took a single `UploadFile` for a PLINK triple | Kept only the last part, discarded the rest, then failed with an error naming the wrong missing file |
| `cli synthetic` deleted the whole database file | Destroyed uploaded research datasets as a side effect of rebuilding lab demo data |
| Ti/Tv computed from the `finding` table | 597 of 606 samples "not assessable" — the reviewable subset is far too small to support a ratio |
| Funnel step 1 left unjoined to `run` | Un-sequenced subjects were blamed on whichever gate came next, usually consent |
| `cohort_run` ignored the assay filter | A gene the chosen panel never looked at counted as "tested and clear" |

---

## 14. Where to start reading

- **Understand the cohort engine:** `backend/app/services/cohort.py` — the gating
  order is the heart of the whole product.
- **Understand the denominators:** `backend/app/services/denominator.py` — the
  single most important module.
- **Understand the safety model:** `backend/app/research/capability.py` and
  `tests/test_governance_gates.py`.
- **Understand the UX:** `frontend/js/flow.js` — four steps, one dispatcher.
- **Understand what makes a rate trustworthy:** `backend/app/services/sampleqc.py`
  and `frontend/js/review.js` — the four checks a cohort has to pass.
