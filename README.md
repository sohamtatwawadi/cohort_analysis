# Germline Cohort Analytics

> **Full platform documentation: [DOCUMENTATION.md](DOCUMENTATION.md)** — current
> state, architecture, API, safety invariants and known gaps. This README is the
> quick start and the mapping back to the source specifications.

Subject-level germline cohort analytics, built from an internal developer specification
(Part B, plus the governance and query-compiler layers). The source specifications are
not distributed with this repository. **Germline only** — the somatic
profile is deliberately not built.

Every number is computed server-side by SQL over a rebuildable analytics store. The
browser renders; it never calculates. A model compiles questions into criteria and is
never shown a patient record or asked for a figure.

```
backend/app/
  schema.sql        derived analytics schema (spec §3.6)
  reference/        curated gene-disease + test-code tables (not in VariMAT, §3.5)
  ingest/           VariMAT loader, clinical sidecar, coverage inference, generator
  services/         cohort engine, denominators, yield, the analysis modules
  api/routes.py     HTTP surface
frontend/           static SPA (no build step)
tests/              including an independent hand count (spec E03.9)
```

## Run it

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m backend.cli synthetic          # build a demo store
.venv/bin/uvicorn backend.app.main:app --port 8077 # http://localhost:8077
```

The server seeds a synthetic store automatically on first start if none exists.

## Loading real VariMAT data

One VariMAT file = one sample. A cohort is many samples, so ingest takes a directory:

```bash
.venv/bin/python -m backend.cli inspect  path/to/sample.tsv   # dedup funnel, writes nothing
.venv/bin/python -m backend.cli ingest   data/varimat --clinical data/clinical
.venv/bin/python -m backend.cli status
.venv/bin/python -m backend.cli verify   # re-resolve every saved cohort, check hashes
```

`ingest` drops and rebuilds the store — it is derived data, never the system of record.

**The clinical sidecar is not optional.** VariMAT carries no subject, family, consent or
indication data (spec §3.5), and consent gates cohort membership while family structure
gates independence. Get the header the loader expects with:

```bash
.venv/bin/python -m backend.cli template data/clinical/subjects.csv
```

A sample with no sidecar row loads with `consent_class = Withdrawn`, so it is inspectable
but cannot enter any cohort until the clinical system supplies a consent class. Failing
closed is the only safe default.

To exercise the loader without real data:

```bash
.venv/bin/python -m backend.tools.make_varimat_fixtures data/varimat --samples 12
```

These fixtures reproduce what makes the real loader hard: ~2x multi-transcript row
duplication, genome-wide intronic VUS noise, non-PASS rows, and gene-symbol aliases.

## What the implementation is careful about

The spec's §E04 lists ten traps. The ones that shaped this code:

**Dedup.** Rows per distinct variant is ~1.75x — the same variant annotated under several
transcripts and gene symbols. Dedup is on `CHROM:START:REF:ALT` with the canonical row
chosen by `MANE`. Skip it and every count is ~75% high. `cli inspect` prints the
per-file funnel so the ratio is visible rather than assumed.

**Gene-specific denominators.** `frequency(gene) = altered / denominator(gene)`, never
over cohort size. A cohort has one subject count and as many denominators as it has
genes, because each test code assays a different gene list. The C01 inspector shows
every gene's assayed / provisional / not-assayed split. Using cohort size instead
produces rates that are quietly, consistently too low — and get quoted.

**Coverage inference (§3.4).** VariMAT has no panel identifier, so coverage is inferred
from per-run gene fingerprints. Two details matter: clustering uses *containment* rather
than Jaccard (a run's fingerprint is a subsample of its panel, so pairwise Jaccard
scatters same-panel runs into singleton clusters, and every singleton produces a
denominator of one); and the presence threshold is **per-gene**, calibrated against runs
whose test code is known, because a 34,350-aa gene surfaces in nearly every fingerprint
while a 213-aa gene does not even when fully assayed.

**Consent.** The withdrawn-consent gate and the secondary-findings gate are compiled into
the member SQL as unconditional predicates, not applied to results. There is no criteria
object — valid or malformed — that produces a cohort containing a withdrawn subject.
`tests/test_governance_gates.py` enumerates the criteria space to check this.

**Derived zygosity (§3.5).** `ZYGOSITY` carries Het/Hom only. Hemizygous is computed from
chromosome plus subject sex; compound-het from two distinct P/LP variants in one gene in
one subject. Reading the file directly misclassifies X-linked males as unaffected
carriers. The G04 screen shows the file value beside the derived one.

**Diagnostic yield.** Gated on phenotype relevance: a P/LP must be in a Definitive or
Strong gene that is relevant to the subject's stated indication. Without that gate a
broad exome counts any P/LP as diagnostic. A P/LP in an unrelated organ system is a
secondary finding, not a diagnosis.

**Reproducibility.** A saved cohort stores criteria, never a member list. `criteria_hash`
plus `kb_snapshot_id` make a figure reproducible; `member_hash` detects silent membership
drift. `cli verify` re-resolves every saved cohort and compares.

## Tests

```bash
.venv/bin/python -m pytest tests/ -q
```

`tests/test_hand_count.py` is the one that matters (spec E03.9: *"Someone counts one
cohort by hand, independently, and the numbers match"*). It recomputes cohort membership,
gene denominators, carrier rates and diagnostic yield from the raw tables in plain Python,
importing nothing from the services package — so a shared helper cannot make both sides
wrong in the same way. Mutating the denominator service to use cohort size fails it with
`ALDOB: hand 69 vs tool 189`.

## Known gaps

- **`backend/app/reference/genes.py` is a seed table.** The 43 gene-disease records were
  hand-entered to match the prototype. Validity and penetrance calls must be replaced with
  a dated ClinGen export, and `gene_id` with an HGNC export, before any figure leaves the
  building. (The spec says the prototype models 44 genes; it models 43.)
- **The query compiler is rule-based, not an LLM.** `compiler.py` implements the D01
  contract and validator; wiring a model in means returning the same
  `{profile, criteria, target_module}` shape and passing it through `validate`, which
  fails closed. The 100-question acceptance gate has not been run.
- **Synthetic calibration is approximate.** Overall yield lands near the spec's 41.5% and
  the SF rate in the low single digits, but per-test-code yields vary more than a real lab
  would. This affects demo data only, not the computation.
- **CNV and fusion show counts only**, with no percentage, because absence of a CNV call
  is not evidence of absence and coverage for them cannot be inferred (spec §3.4, E03.4).
- **Not built:** the somatic profile (S01–S14), and G06/G07/G09 are present but were
  `[P1-opt]` in the spec.

---

# Part II — Research Mode

Built from the Part II specification. One engine, two modes
(§1): Lab Mode analyses our own pipeline output, Research Mode analyses uploaded
datasets. Switch with the Lab/Research toggle in the rail. The cohort engine,
governance, audit and export are shared; mode changes what data flows in and
which analyses unlock.

**Scope built:** R1 (ingestion, profiling, capability gating, registry, project
isolation, and the five descriptive analyses — carrier frequency, zygosity,
population frequency, diagnostic yield, segregation), R2 (association), R4
(GWAS), R5 (burden/SKAT/SKAT-O), R6 (PRS, survival, penetrance). **Not built:** R3 (PheWAS — needs an ICD→PheCode
ontology), R7 (fine-mapping, colocalization, heritability — need LD reference
panels and summary statistics), R8.

```
backend/app/research/
  types.py        GenotypeMatrix / PhenotypeTable — the mode boundary
  formats/        VCF, PLINK1, format detection
  validate.py     §2.3 mandatory validations (fail the upload, never warn)
  profile.py      §3.1 data profile
  capability.py   §3.2–3.4 the capability matrix
  registry.py     §8 dataset registry, project isolation, deletion
  jobs.py         §6 job queue, lifecycle, reproducibility record
  variantset.py   §4.4 qualifying-variant builder, versioned
  stats/          glm (incl. Firth), pca, kinship, qc, skat, survival, power,
                  scoretest (the vectorised GWAS scan)
  analyses/       association, gwas, burden, prs, survival, descriptive
```

## Try it

```bash
.venv/bin/python -m backend.tools.make_research_fixture --samples 2400 --variants 60000
.venv/bin/uvicorn backend.app.main:app --port 8077
```

Then switch to Research in the rail. The fixture plants six causal variants at
OR 1.9, two populations with differing allele frequencies, 30 parent-offspring
pairs, 40 complete trios and a rare-variant burden signal in three genes, so the
analyses have a known right answer. The GWAS scan recovers all six planted
variants inside its top eight hits, at odds ratios of 1.49–1.78 against a planted
1.9, with λ_GC = 1.007 — it finds what was planted and nothing else. The kinship
engine recovers exactly the 30 planted pairs (max κ = 0.252 against a theoretical
0.25).

## Capability gating is the product

§3 calls it "the most important safety feature", and three properties are
enforced structurally rather than by convention:

**Locked is visible, not hidden.** Every unavailable analysis renders with each
unmet requirement *and the observed value* — "Requires n ≥ 2,000 / your dataset
has 486 subjects". Hiding it would just make the product look incapable.

**Locked cannot be clicked through.** A locked card has no run control at all —
not a disabled one. `jobs.submit` calls `registry.can_run` before queuing, so an
analysis cannot start even from a direct API call. From the spec: *"A researcher
under deadline will click through a warning. A locked control cannot be clicked
through."*

**Override is separate, recorded, and stamps the output.** Borderline thresholds
can be overridden with a typed justification (≥20 chars, audited). Missing data
cannot — you cannot override your way to genotypes you did not upload. The
`UNDERPOWERED — threshold overridden` stamp is derived from the override record
by the job runner, so it cannot be edited off a result payload.

## Statistical guardrails (§5)

- **Power before execution.** Every analysis reports estimated power with the
  spec's verdict text and options. Validated against the spec's §7 worked
  example (94.1% vs the spec's 94%).
- **Firth is used where it is required.** §4.1 mandates it because standard
  logistic separates with rare exposures. Model selection routes to Firth on
  sparse 2×2 cells and states why; p-values come from the penalised likelihood
  ratio test, not Wald, because Wald is unreliable in exactly the separation
  regime Firth exists for.
- **λ_GC guardrail with attributed cause.** A GWAS exceeding the inflation
  threshold names the likely cause from the data profile (unmodelled ancestry,
  cryptic relatedness, batch, differential missingness) rather than just
  printing a number.
- **Mandatory pre-GWAS QC.** Call rates, MAF, HWE (exact, in controls), case
  /control differential missingness, relatedness pruning, ancestry outliers,
  sex check. There is no parameter that skips it.
- **Too-few-carrier genes are reported as untested**, not given an unstable
  p-value — an untested gene is not evidence of no effect.
- **Ancestry-stratified PRS performance is always computed.** No flag turns it
  off; if a group is absent, that absence is the finding.
- **Penetrance from an ascertained cohort is labelled**, detected from the data
  profile.
- **Multiple testing always.** Raw p, Bonferroni and BH-FDR on every result.

## Governance for uploaded data (§8)

Isolation is a `project_id` column on every dataset, job and result, and every
read goes through a function that takes one — there is no query path that
returns another project's rows. Deletion propagates to the profile, capability
assessment, phenotype registry, sample roster, overrides, jobs, results *and the
genotype payload on disk*, and reports the counts so it is verifiable rather
than asserted. Consent attestation is enforced in `register_dataset`, not in the
API layer, so no code path can skip it.

## Known gaps and judgement calls

- **Engine is in-process numpy/scipy**, not REGENIE/SAIGE/PLINK2. §4.2 prefers
  delegating to an established tool but also says engine selection "needs
  benchmarking before committing". This is correct and tested at the scale here
  (thousands of samples); it is **not biobank-scale**. An exact in-memory SVD
  and an O(n²) kinship matrix are the binding limits.
- **The spec's §5.1 power example does not reproduce.** Its stated 11% at
  n=4,200 / 180 cases / MAF 0.008 / OR 2.0 / α=5×10⁻⁸ would require OR≈4.0 or
  α≈0.005; the standard formula gives ≈0%. The implementation matches the §7
  example exactly, and gives 72.8% on a textbook design, so the §5.1 figure
  appears illustrative rather than computed. Worth confirming before anyone
  calibrates against it.
- **BCF, PLINK2/PGEN, BGEN and Hail are not implemented** (§2.1 lists BCF as P1).
- **Left-alignment is trim-based parsimony**, not true left-shifting through
  repeats — that needs a reference FASTA the loader has no handle on.
- **SKAT-O's multiple-ρ correction is approximate** and says so in its result.
- **Gene-based tests need gene annotations supplied**; the loader does not
  annotate. Without them the analysis refuses rather than testing one giant
  pseudo-gene.
- **Capability thresholds are the spec's starting positions** and §3.2 says they
  must be calibrated. Where power can be computed it is better evidence than a
  raw n.
