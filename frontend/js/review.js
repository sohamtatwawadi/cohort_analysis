/* Step 2 — Review. "Can you trust this cohort?"
 *
 * Previously this information existed but was scattered: three warning banners
 * stacked above every screen, and a denominator inspector buried in a modal
 * behind a rail link. Users learned to scroll past the banners.
 *
 * Here it is the whole step, expressed as three checks that either pass or
 * warn, each expandable to the evidence. You cannot reach the analysis step
 * without having seen them. */

import { api } from './api.js';
import {
  bars, card, empty, esc, f1, fmt, head, note, rate, stat, table, tag,
} from './kit.js';
import { datasetDashboard } from './datasetdash.js';

export async function render(state) {
  if (state.mode === 'research') return renderResearch(state);
  return renderLab(state);
}

/* ===================================================================== lab == */
async function renderLab(state) {
  const [dash, den, qc] = await Promise.all([
    api.module('g_dashboard', state.criteria, state.name),
    api.denominator(state.criteria, state.name),
    api.sampleQc(state.criteria, state.name),
  ]);
  state.reviewed = true;
  const d = dash.data;
  const k = d.kpis;
  const warnings = dash.cohort.warnings || [];
  const cov = den.data;

  const byCode = (code) => warnings.find((w) => w.code === code);
  const probandsOff = byCode('probands_off');
  const provenance = byCode('provenance');
  const provisional = byCode('provisional');

  const flagged = cov.genes.filter((g) => g.warn);
  const q = qc.data;
  const qs = q.summary;

  return `<div class="wrap">
    ${head('Does this cohort hold up?',
      'Four checks decide whether these numbers can leave the building. Open any one to '
      + 'see the evidence behind it.', 'Step 2 of 4')}

    ${stat([
      { k: 'Subjects', v: fmt(k.subjects), d: `in ${fmt(k.families)} families` },
      { k: 'Findings', v: fmt(k.all_findings), d: `${fmt(k.plp_findings)} pathogenic` },
      { k: 'Genes covered', v: fmt(cov.genes.filter((g) => g.assayed > 0).length),
        d: `of ${fmt(cov.genes.length)} in the reference` },
      { k: 'Samples failing QC', v: fmt(qs.fail),
        d: qs.fail ? 'exclude or investigate' : 'none' },
    ])}

    ${check({
      pass: qs.fail === 0,
      title: 'Sample quality',
      summary: qs.fail
        ? `${qs.fail} of ${qs.total} samples fail QC and ${qs.warn} are borderline. `
          + 'A contaminated sample or a swap invalidates every rate below it.'
        : `All ${qs.total} samples are consistent with the rest of the cohort`
          + (qs.warn ? ` (${qs.warn} borderline)` : '') + '.',
      more: 'See every sample',
      detail: `
        <p class="hint" style="margin-top:0">${esc(q.method)} ${esc(q.caveat)}</p>
        ${q.flagged.length ? table([
          { label: 'Sample', cell: (x) => `<span class="mono">${esc(x.run_id)}</span>` },
          { label: 'Subject', cell: (x) => `<span class="mono">${esc(x.subject_id)}</span>` },
          { label: 'Test', cell: (x) => esc(x.test_code || 'inferred') },
          { label: 'Ti/Tv', n: true, cell: (x) => x.ti_tv === null ? '—'
              : `<span class="mono">${Number(x.ti_tv).toFixed(2)}</span>` },
          { label: 'Het/Hom', n: true, cell: (x) => x.het_hom === null ? '—'
              : `<span class="mono">${Number(x.het_hom).toFixed(2)}</span>` },
          { label: 'Het VAF', n: true, cell: (x) => x.mean_het_vaf === null ? '—'
              : `<span class="mono">${Number(x.mean_het_vaf).toFixed(3)}</span>` },
          { label: 'Depth', n: true, cell: (x) => x.mean_depth === null ? '—'
              : `<span class="mono">${Number(x.mean_depth).toFixed(0)}×</span>` },
          { label: '% ≥20×', n: true, cell: (x) => x.pct_bases_20x === null ? '—'
              : `<span class="mono">${Number(x.pct_bases_20x).toFixed(0)}</span>` },
          { label: 'Status', cell: (x) => tag(x.status,
              x.status === 'fail' ? 'red' : 'amber') },
          { label: 'Why', cell: (x) => x.flags.map((f) =>
              `<div><b>${esc(f.label)}</b> — <span class="hint">${esc(f.detail)}</span></div>`)
              .join('') },
        ], q.flagged) : note('Every sample is consistent with the cohort.', 'good')}
        <p class="hint" style="margin-top:12px">Cohort medians —
          Ti/Tv ${q.cohort.ti_tv.median}, Het/Hom ${q.cohort.het_hom.median},
          het VAF ${q.cohort.mean_het_vaf.median},
          depth ${q.cohort.mean_depth.median}×.
          Metrics basis: ${esc((q.samples[0] || {}).metrics_basis || 'n/a')}.</p>`,
    })}

    ${check({
      pass: !probandsOff,
      title: 'Independence',
      summary: probandsOff
        ? probandsOff.text
        : 'One subject per family. Carrier rates are independent and can be quoted.',
      more: probandsOff ? 'Why this matters' : 'What was excluded',
      detail: probandsOff
        ? note(`Relatives share alleles by descent. Counting them separately makes a family's
            single variant look like several independent carriers, which inflates every
            carrier rate. Turn <b>One subject per family</b> back on in step 1 unless you
            have a specific reason not to.`, 'warn')
        : `<p>${esc((byCode('probands_on') || {}).text || '')}</p>
           <p class="hint">Relatives share alleles by descent, so counting them separately
           would inflate every carrier rate.</p>`,
    })}

    ${check({
      pass: !provisional && flagged.length === 0,
      title: 'Coverage',
      summary: provisional
        ? provisional.text
        : flagged.length
          ? `${flagged.length} gene(s) have a denominator below 30 or largely provisional scope.`
          : 'Every gene has a solid denominator. Rates are fit for external quotation.',
      more: 'See every gene’s denominator',
      detail: `
        <p class="hint" style="margin-top:0">${esc(cov.explain)}</p>
        ${table([
          { label: 'Gene', cell: (g) => `<span class="gene">${esc(g.gene)}</span>` },
          { label: 'Condition', cell: (g) => esc(g.condition) },
          { label: 'Assayed', n: true, cell: (g) => fmt(g.assayed) },
          { label: 'Provisional', n: true, cell: (g) => g.provisional
              ? `<span style="color:var(--warn)">${fmt(g.provisional)}</span>` : '0' },
          { label: 'Not assayed', n: true, cell: (g) => fmt(g.not_assayed) },
          { label: 'Coverage', cell: (g) => `<div class="row2">
              <div class="bar2" style="flex:1"><i style="width:${g.coverage_pct}%"></i></div>
              <span class="mono" style="width:46px;text-align:right">${f1(g.coverage_pct)}%</span>
            </div>` },
          { label: 'Basis', cell: (g) => tag(g.confidence,
              g.confidence === 'declared' ? 'teal' : g.confidence === 'inferred' ? 'blue' : 'amber') },
          { label: '', cell: (g) => g.warn
              ? `<span class="warnic" title="${esc(g.warn_reason)}">⚠</span>` : '' },
        ], cov.genes)}`,
    })}

    ${check({
      pass: !provenance,
      title: 'Provenance',
      summary: provenance
        ? provenance.text
        : 'One reference build and one pipeline version. Safe to compare across subjects.',
      more: 'See the test codes and pipelines',
      detail: table([
        { label: 'Implied panel', cell: (p) => `<span class="mono">${esc(p.implied_panel)}</span>` },
        { label: 'Declared test code', cell: (p) => `<span class="mono">${esc(p.test_code)}</span>` },
        { label: 'Runs', n: true, cell: (p) => fmt(p.runs) },
      ], cov.panels),
    })}

    ${card('What is in this cohort', `
      <div class="cols c2">
        <div>
          <p class="hint" style="margin-top:0">Reason for testing</p>
          ${bars((d.yield_by_indication || []).map((r) => ({
            label: r.bucket, value: r.subjects, display: fmt(r.subjects) })))}
        </div>
        <div>
          <p class="hint" style="margin-top:0">Ancestry</p>
          ${bars((d.ancestry_stack || []).map((r) => ({ label: r.label, value: r.value })))}
        </div>
      </div>`)}

    <div class="actions">
      <button class="btn primary lg" data-act="goto:analyse">Analyse this cohort →</button>
      <button class="btn" data-act="goto:build">← Change who is included</button>
    </div>
  </div>`;
}

function check({ pass, title, summary, more, detail }) {
  return `<details class="check ${pass ? 'pass' : 'warn'}">
    <summary>
      <span class="badge">${pass ? '✓' : '!'}</span>
      <span><span class="ct">${esc(title)}</span>
        <span class="cs">${summary}</span>
        <span class="more">${esc(more)} ▸</span></span>
    </summary>
    <div class="detail">${detail}</div>
  </details>`;
}

/* ================================================================ research == */
async function renderResearch(state) {
  const p = state.profile || {};
  const caps = state.caps || [];
  const available = caps.filter((c) => c.available);
  const locked = caps.filter((c) => !c.available);

  /* The dashboard answers "what IS this data?", which a reader needs before
     "what can it support?". It is fetched rather than derived from the profile
     because the gene panel needs the genotype matrix; the server caches it.
     A failure here must not take the capability matrix down with it — that is
     the part that gates the analyses. */
  let dash = null;
  try {
    if (state.datasetId) {
      dash = (await api.research.dashboard(state.datasetId)).dashboard;
    }
  } catch (err) {
    dash = null;
  }

  return `<div class="wrap">
    ${head(esc((state.dataset || {}).name || 'This dataset'),
      'What the data contains, and which analyses it can actually answer. '
      + 'Locked analyses state exactly what is missing.', 'Step 2 of 4')}

    ${dash ? datasetDashboard(dash) : `${stat([
      { k: 'Samples', v: fmt(p.n_samples), d: `${fmt(p.n_unrelated || 0)} unrelated` },
      { k: 'Variants', v: fmt(p.n_variants), d: esc(p.density_class || '') },
      { k: 'Ready to run', v: String(available.length),
        d: `${locked.length} need more data` },
    ])}
    ${note('The dataset summary could not be built, so only the capability '
         + 'matrix is shown below.', 'warn')}`}

    ${card(`Ready to run — ${available.length}`,
      available.length
        ? `<div class="qgrid">${available.map((c) => `
            <button class="q" data-act="rrun:${esc(c.analysis)}">
              <span class="qt">${esc(c.title)}</span>
              <span class="qd">${c.requirements.map((r) =>
                esc(r.observed)).slice(0, 2).join(' · ')}</span>
              <span class="qa">Configure →</span></button>`).join('')}</div>`
        : empty('Nothing can run yet', 'The requirements below show what to add.'))}

    ${card(`Needs more data — ${locked.length}`, `<div class="qgrid">${locked.map((c) => `
      <div class="q locked">
        <span class="qt">${esc(c.title)}</span>
        ${c.unmet.map((r) => `<span class="why">✗ ${esc(r.label)}
          <span class="obs">— you have: ${esc(r.observed)}</span></span>`).join('')}
        <span class="qd">${esc(c.remedy)}</span>
        ${c.override && c.override.active
          ? `<span class="qa" style="color:var(--warn)">Override active — results stamped</span>`
          : c.override && c.override.can_override
            ? `<button class="btn sm" data-act="rovr:${esc(c.analysis)}">Request override…</button>`
            : `<span class="hint">${esc((c.override || {}).reason || '')}</span>`}
      </div>`).join('')}</div>`)}

    <div class="actions">
      <button class="btn primary lg" data-act="goto:analyse">Go to analyses →</button>
      <button class="btn" data-act="goto:build">← Choose a different dataset</button>
    </div>
  </div>`;
}

export async function handle() { return false; }
