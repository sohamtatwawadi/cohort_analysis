/* Level-1 dataset dashboard — "what is this data?"
 *
 * Registering a dataset used to drop you straight into the capability matrix,
 * which answers "what can this support?" while never answering the question a
 * reader has first.
 *
 * Laid out after an oncology cohort dashboard: a KPI row, top genes with a bar
 * and its denominator, composition breakdowns, then distributions. One
 * substitution matters. A somatic dashboard leads with tumour mutational
 * burden, stage and linked therapy; germline data has none of those, and
 * rendering empty tiles labelled TMB would read as clinically meaningful while
 * being fabricated. The slots are filled from what the dataset really holds.
 *
 * The denominator note under the gene panel is the same rule the lab side
 * enforces, and the same one the reference dashboard states: a per-gene rate
 * divides by the samples actually called at that gene.
 */

import { histogram } from './charts.js';
import { CATEGORICAL, bars, card, esc, f1, fmt, note, stacked, stat } from './kit.js';

const pct1 = (x) => (x === null || x === undefined ? '—' : (x * 100).toFixed(1) + '%');

export function datasetDashboard(d) {
  if (!d) return '';

  /* Which panel earns the wide column depends on the dataset. A VCF with no
     gene annotations — 1000 Genomes, for one — has nothing to put in the gene
     panel, so giving it two thirds of the width to display an explanation
     while 26 real populations are squeezed into the remainder is the wrong way
     round. */
  const hasGenes = (d.top_genes || []).length > 0;
  return `
    ${kpiRow(d)}
    <div class="cols ${hasGenes ? 'c21' : 'c12'}">
      ${hasGenes ? topGenes(d) + composition(d) : composition(d) + topGenes(d)}
    </div>
    <div class="cols c2">
      ${spectrum(d)}
      ${provenance(d)}
    </div>
    ${distributions(d)}
    ${d.ascertained ? note(`<b>This looks like a referral-selected cohort.</b>
      ${esc(d.ascertainment_rationale)}`, 'warn') : ''}
    ${(d.warnings || []).length
      ? note(`<b>From profiling:</b><ul style="margin:6px 0 0">${
          d.warnings.map((w) => `<li>${esc(w)}</li>`).join('')}</ul>`, 'info')
      : ''}`;
}

/* The headline row. Values carry their own unit and a line of context, because
   a bare "3,000" invites the reader to guess what it counts. */
function kpiRow(d) {
  return stat((d.kpis || []).map((k) => ({
    k: k.key,
    v: k.value === null || k.value === undefined
      ? '—'
      : (typeof k.value === 'number' ? fmt(k.value) : esc(String(k.value)))
        + (k.unit || ''),
    d: esc(k.detail || ''),
  })));
}

/* Top genes. The bar is the visual, but the counts beside it are the finding —
   "21.8%" alone is not checkable, "87 / 400" is. */
function topGenes(d) {
  const rows = d.top_genes || [];
  if (!rows.length) {
    return card('Top genes by carrier rate',
      note(esc(d.genes_note || 'No gene-level summary available.'), 'info'));
  }
  const max = Math.max(...rows.map((r) => r.pct)) || 1;
  return card('Top genes by carrier rate', `
    <div class="genebars">
      ${rows.map((r) => `
        <div class="gb">
          <div class="gb-name"><span class="gene">${esc(r.gene)}</span></div>
          <div class="gb-track"><i style="width:${Math.max(1.5, 100 * r.pct / max)}%"></i></div>
          <div class="gb-pct mono">${r.suppressed ? '—' : pct1(r.pct)}</div>
          <div class="gb-n mono">${r.suppressed
            ? '<span class="hint">&lt;5</span>' : fmt(r.carriers)}
            <span class="hint">/ ${fmt(r.n_called)}</span></div>
          <div class="gb-v hint">${fmt(r.n_variants)} var</div>
        </div>`).join('')}
    </div>`, { sub: 'carriers / samples called', foot: esc(d.genes_note) });
}

function composition(d) {
  const panels = d.composition || [];
  if (!panels.length) {
    return card('Cohort composition',
      note('No categorical phenotype columns were supplied, so there is nothing '
         + 'to break the cohort down by. Attach a phenotype file with columns '
         + 'such as ancestry, sex or case status.', 'info'));
  }
  return card('Cohort composition', panels.map((p) => `
    <div class="comp">
      <p class="comp-h">${esc(p.label)}</p>
      ${stacked(p.groups.map((g) => ({ label: g.name, value: g.count })),
                { colors: CATEGORICAL })}
    </div>`).join(''));
}

/* The allele-frequency spectrum is the shape of the call set, and it is the
   single most diagnostic thing about an uploaded dataset: a file with no rare
   variants is a genotyping array, whatever its header claims. */
function spectrum(d) {
  const rows = (d.spectrum || []).filter((s) => s.count > 0);
  if (!rows.length) return '';
  const total = rows.reduce((a, r) => a + r.count, 0);
  return card('Variant frequency spectrum', bars(rows.map((r) => ({
    label: r.label, value: r.count,
    display: `${fmt(r.count)}  ${f1(100 * r.count / total)}%`,
  }))), { sub: `${fmt(total)} variants`,
    foot: 'A call set with no rare variants is an array, not sequencing — '
        + 'whatever the file header says.' });
}

function provenance(d) {
  const r = d.relatedness || {};
  const a = d.ancestry || {};
  return card('Structure and provenance', `
    <dl class="kv">
      <dt>Genome build</dt><dd class="mono">${esc(d.genome_build || 'not declared')}</dd>
      <dt>Source format</dt><dd class="mono">${esc(d.source_format || '—')}</dd>
      <dt>Ancestry PCs</dt><dd>${fmt(a.n_pcs || 0)} computed</dd>
      <dt>Related pairs</dt><dd>${fmt(r.n_related_pairs || 0)}${
        r.max_kinship ? ` · max κ ${f1(r.max_kinship * 100) / 100}` : ''}</dd>
      <dt>Samples with a relative</dt><dd>${fmt(r.n_samples_with_relative || 0)}</dd>
      <dt>Complete trios</dt><dd>${fmt(r.n_trios || 0)}</dd>
      <dt>Summary computed</dt><dd class="mono">${esc(
        String(d.computed_at || '').replace('T', ' '))}</dd>
    </dl>`, {
    foot: 'Relatedness matters before any analysis assuming independent '
        + 'samples — relatives share alleles by descent and inflate every rate.' });
}

function distributions(d) {
  const ds = d.distributions || [];
  if (!ds.length) return '';
  return card('Distributions', `<div class="cols c2">
    ${ds.slice(0, 4).map((x) => histogram(x.bins, {
      title: x.label,
      sub: x.mean !== null && x.mean !== undefined
        ? `mean ${f1(x.mean)}${x.sd ? ` · SD ${f1(x.sd)}` : ''}` : '',
      xLabel: x.label,
      yLabel: 'samples',
    })).join('')}
  </div>`);
}
