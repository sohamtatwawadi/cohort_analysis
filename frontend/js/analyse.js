/* Step 3 — Analyse.
 *
 * The old rail listed eleven module names: "G05 Gene–disease association",
 * "G07 Population frequency". Those are our filing system, not the user's
 * question. Here the same eleven modules are eight questions in the words an
 * analyst would actually use, and the module name never appears.
 *
 * Drill-down (gene → variant → subject) is reached by clicking a number, not
 * by a separate navigation tree, so the evidence is always one click from the
 * figure it supports. */

import { api } from './api.js';
import {
  dotplot, dumbbell, histogram, kaplanMeier, manhattan, mountCharts, qq,
  scatter, volcano,
} from './charts.js';
import {
  acmg, bars, card, closeModal, crumb, empty, esc, f1, fmt, head, heatmap, kv,
  modal, note, rate, stacked, stat, table, tag, toast, ZYG_TONE,
} from './kit.js';

/* The catalogue. `module` is an implementation detail the user never sees. */
export const QUESTIONS = [
  { id: 'yield', module: 'g_carrier', title: 'How many subjects got a diagnosis?',
    desc: 'Diagnostic yield, counting only pathogenic findings in well-established genes '
        + 'relevant to why the subject was tested.' },
  { id: 'carrier', module: 'g_carrier', title: 'Which genes carry pathogenic variants?',
    desc: 'Carrier rate per gene, each on its own denominator — the subjects whose test '
        + 'actually covers that gene.' },
  { id: 'zygosity', module: 'g_zygosity', title: 'Are they carriers, or affected?',
    desc: 'A heterozygous and a biallelic finding in the same recessive gene are different '
        + 'clinical entities. Includes derived hemizygous and compound-het calls.' },
  { id: 'validity', module: 'g_genedisease', title: 'How strong is the gene–disease evidence?',
    desc: 'What the findings implicate, and how much the yield would inflate if weak-evidence '
        + 'genes were counted.' },
  { id: 'vus', module: 'g_vus', title: 'Which uncertain variants should we review first?',
    desc: 'A ranked curation queue, scored on new evidence, how many subjects are affected, '
        + 'and how stale the interpretation is.' },
  { id: 'popfreq', module: 'g_popfreq', title: 'Is any variant unusually common here?',
    desc: 'Internal frequency against gnomAD. Finds founder candidates — and implausible '
        + 'frequencies, which are a data-quality signal.' },
  { id: 'secondary', module: 'g_sf', title: 'Are there secondary findings to return?',
    desc: 'Medically actionable incidental findings, among consented subjects only.' },
  { id: 'phenotype', module: 'g_phenotype', title: 'Which phenotypes go with which genes?',
    desc: 'Hypothesis generation across recorded HPO terms, plus a phenotype-matched '
        + 'review queue.' },
];

/* Part II §10: the statistical suite is "Restored. Locked in Lab Mode, active
   in Research Mode." Hiding it in Lab Mode would make the product look
   incapable and leave a user wondering where GWAS went; showing it greyed with
   the real reason tells them exactly what to do — upload a study with controls.

   The reason is the one Part II §0 gives for why v2 removed these: our own
   diagnostic data is ascertained, panel-based and has no unselected controls.
   That judgement holds for Lab Mode and does not hold for an upload. */
const LAB_LOCKED = [
  { title: 'Association (single variant)',
    why: 'Needs an unselected control group. Lab cohorts are referral-selected — '
       + 'everyone was tested because something prompted it.' },
  { title: 'Genome-wide association (GWAS)',
    why: 'Needs genome-wide genotypes. Lab data is panel and exome, typically '
       + 'a few thousand variants over selected genes.' },
  { title: 'PheWAS',
    why: 'Needs a coded phenome (ICD or PheCode) per subject. Lab records carry '
       + 'a requisition indication and HPO terms, not an EHR phenome.' },
  { title: 'Burden / SKAT / SKAT-O',
    why: 'Needs case-control status with unselected controls.' },
  { title: 'Polygenic risk score',
    why: 'Needs genome-wide genotypes and an ancestry-matched validation set.' },
  { title: 'Survival / penetrance',
    why: 'Needs time-to-event follow-up. Lab records are a point-in-time result, '
       + 'not a followed cohort.' },
];

export async function render(state) {
  if (state.mode === 'research') return renderResearch(state);
  if (state.gene) return renderGene(state);
  if (!state.question) return renderMenu(state);
  return renderAnswer(state);
}

/* ------------------------------------------------------------------ menu -- */
async function renderMenu(state) {
  const co = state.cohort;
  // The catalogue is ranked against what this cohort actually contains, so a
  // question that would return nothing is not sitting next to one that answers
  // the user's actual problem looking identical.
  const sug = await api.suggestions(state.criteria, state.name);
  const by = {};
  (sug.suggestions || []).forEach((x) => { by[x.question] = x; });
  const ordered = (sug.suggestions || []).map((x) =>
    QUESTIONS.find((q) => q.id === x.question)).filter(Boolean);
  const group = (v) => ordered.filter((q) => (by[q.id] || {}).verdict === v);

  const qcard = (q, muted) => {
    const s = by[q.id] || {};
    return `<button class="q ${muted ? 'locked' : ''}" data-act="ask:${q.id}"
      ${muted ? 'aria-disabled="true"' : ''}>
      <span class="qt">${esc(q.title)}</span>
      <span class="qd">${esc(q.desc)}</span>
      ${s.reason ? `<span class="${muted ? 'why' : 'qd'}"
        style="${muted ? '' : 'color:var(--accent);font-weight:600'}">${
        muted ? '' : '→ '}${esc(s.reason)}</span>` : ''}
      ${muted ? '' : '<span class="qa">Show me →</span>'}
    </button>`;
  };

  const section = (title, list, muted) => list.length ? `
    ${card(`${title} — ${list.length}`,
      `<div class="qgrid">${list.map((q) => qcard(q, muted)).join('')}</div>`)}` : '';

  return `<div class="wrap">
    ${head('What do you want to know?',
      `Every answer is computed on your cohort of ${fmt(co.counts.subjects)} subjects. `
      + 'Click any number in a result to see the subjects behind it.', 'Step 3 of 4')}
    ${note(`<b>${esc(sug.summary)}</b>`, 'info')}
    ${section('Worth asking of this cohort', group('recommended'), false)}
    ${section('Will run, but this cohort was not built for it', group('available'), false)}
    ${section('Would return nothing here', group('not_useful'), true)}
    ${card(`Needs uploaded data — ${LAB_LOCKED.length}`,
      `<p class="hint" style="margin:0 0 12px">These run on uploaded study data, not on
        our own laboratory records. Switch to <b>Uploaded data</b> and register a cohort
        with controls to use them.</p>
      <div class="qgrid">${LAB_LOCKED.map((x) => `
        <div class="q locked">
          <span class="qt">${esc(x.title)}</span>
          <span class="why">✗ ${esc(x.why)}</span>
        </div>`).join('')}</div>`)}

    <div class="actions">
      <button class="btn" data-act="goto:review">← Back to review</button>
      <span class="spacer"></span>
      <button class="btn" data-act="goto:share">Skip to export →</button>
    </div>
  </div>`;
}

/* ---------------------------------------------------------------- answers -- */
async function renderAnswer(state) {
  const q = QUESTIONS.find((x) => x.id === state.question);
  const res = await api.module(q.module, state.criteria, state.name);
  const d = res.data;
  const body = ANSWERS[q.id](d, state);
  return `<div class="wrap">
    ${crumb([{ label: 'All questions', act: 'ask:' }, { label: q.title }])}
    ${head(q.title, q.desc)}
    ${body}
    <div class="actions">
      <button class="btn" data-act="ask:">← Ask something else</button>
      <span class="spacer"></span>
      <button class="btn primary" data-act="goto:share">Export this →</button>
    </div>
  </div>`;
}

const ANSWERS = {
  yield: (d) => {
    const y = d.yield;
    return `
      ${headline(`${f1(y.rate.pct)}%`,
        `${fmt(y.solved)} of ${fmt(y.total)} subjects have a diagnostic finding`)}
      ${note(esc(y.caveat), 'info')}
      <div class="cols c2">
        ${card('Outcome for every subject', stacked([
          { label: 'Diagnostic (P/LP)', value: y.solved, tone: 'plp' },
          { label: 'Uncertain only', value: y.vus_only, tone: 'vus' },
          { label: 'No reportable finding', value: y.negative, tone: 'neg' },
        ]), { foot: 'A subject is diagnostic only if the finding is pathogenic, in a '
              + 'Definitive or Strong gene, reportable, and relevant to why they were tested.' })}
        ${card('Yield by test code', table([
          { label: 'Test code', cell: (r) => `<span class="mono">${esc(r.bucket)}</span>` },
          { label: 'Subjects', n: true, cell: (r) => fmt(r.subjects) },
          { label: 'Yield', n: true, cell: (r) => rate(r.rate) },
        ], d.by_test_code), { flush: true })}
      </div>
      ${card('Yield by ancestry', table([
        { label: 'Ancestry', cell: (r) => esc(r.bucket) },
        { label: 'Subjects', n: true, cell: (r) => fmt(r.subjects) },
        { label: 'Solved', n: true, cell: (r) => fmt(r.solved) },
        { label: 'Yield', n: true, cell: (r) => rate(r.rate) },
      ], d.by_ancestry), { flush: true, foot: esc(d.ancestry_caveat) })}`;
  },

  carrier: (d) => `
    ${headline(fmt(d.table.length), 'genes with at least one pathogenic carrier')}
    ${card('Carrier rate by gene', table([
      { label: 'Gene', cell: (r) => `<a class="gene" data-act="gene:${esc(r.gene)}">${esc(r.gene)}</a>` },
      { label: 'Condition', cell: (r) => esc(r.condition) },
      { label: 'Inh', cell: (r) => esc(r.inheritance) },
      { label: 'Carriers', n: true, cell: (r) => fmt(r.plp_subjects) },
      { label: 'Families', n: true, cell: (r) => `${fmt(r.plp_families)}${
          r.family_clustered ? ' <span class="warnic" title="fewer families than carriers — '
          + 'the same allele is being counted more than once from one pedigree">⚠</span>' : ''}` },
      { label: 'Tested for it', n: true, cell: (r) =>
          `<span class="mono">${esc(r.denominator_label)}</span>` },
      { label: 'Carrier rate', n: true,
        cell: (r) => rate(r.carrier_rate, { warn: r.warn, warnReason: r.warn_reason }) },
    ], d.table, { click: true, attrs: (r) => `data-act="gene:${esc(r.gene)}"` }), {
      flush: true,
      foot: esc(d.warn_legend) + ' Click a gene for its variants and per-ancestry breakdown.',
    })}`,

  zygosity: (d) => {
    const k = d.kpis;
    return `
      ${stat([
        { k: 'Biallelic', v: fmt(k.biallelic_subjects), d: 'consistent with affected' },
        { k: 'AR carriers', v: fmt(k.ar_carriers), d: 'heterozygous — not affected' },
        { k: 'Compound het', v: fmt(k.compound_het_subjects), d: 'phase unconfirmed' },
        { k: 'Hemizygous', v: fmt(k.hemizygous_subjects), d: 'derived, not read from file' },
        { k: 'Dominant het', v: fmt(k.ad_het_subjects), d: 'one copy is sufficient' },
      ])}
      ${note(`These five numbers are deliberately never added together. A heterozygous
        carrier of a recessive condition and a biallelic case are different people with
        different clinical meaning.`, 'info')}
      <div class="cols c2">
        ${card('Inheritance × zygosity', heatmap(d.heatmap.rows, d.heatmap.cols,
          d.heatmap.matrix, d.heatmap.max, { shortCol: (c) => c.split(' ')[0] }),
          { sub: 'subjects with a pathogenic finding' })}
        ${card('Compound heterozygous pairs', table([
          { label: 'Subject', cell: (r) => `<a class="mono" data-act="subject:${esc(r.subject_id)}">${esc(r.subject_id)}</a>` },
          { label: 'Gene', cell: (r) => `<a class="gene" data-act="gene:${esc(r.gene)}">${esc(r.gene)}</a>` },
          { label: 'Variants', cell: (r) => `<span class="mono">${esc(r.variants)}</span>` },
          { label: 'Phase', cell: (r) => tag(r.phase, 'amber') },
        ], d.compound_het, { emptyTitle: 'No compound-heterozygous pairs' }),
          { flush: true, foot: esc(d.compound_het_footer) })}
      </div>
      ${card('By gene', table([
        { label: 'Gene', cell: (r) => `<a class="gene" data-act="gene:${esc(r.gene)}">${esc(r.gene)}</a>` },
        { label: 'Condition', cell: (r) => esc(r.condition) },
        { label: 'Inh', cell: (r) => esc(r.inheritance) },
        { label: 'Het', n: true, cell: (r) => fmt(r.het) },
        { label: 'Hom', n: true, cell: (r) => fmt(r.hom) },
        { label: 'Comp het', n: true, cell: (r) => fmt(r.compound_het) },
        { label: 'Hemi', n: true, cell: (r) => fmt(r.hemizygous) },
        { label: 'What it means', cell: (r) => r.flag ? tag(r.flag, r.flag_tone) : '—' },
      ], d.by_gene), { flush: true })}
      ${card('Hemizygous calls — derived, not read', table([
        { label: 'Subject', cell: (r) => `<a class="mono" data-act="subject:${esc(r.subject_id)}">${esc(r.subject_id)}</a>` },
        { label: 'Gene', cell: (r) => `<span class="gene">${esc(r.gene)}</span>` },
        { label: 'Variant', cell: (r) => `<span class="mono">${esc(r.hgvs_p || r.hgvs_c)}</span>` },
        { label: 'File said', cell: (r) => tag(r.file_says, 'grey') },
        { label: 'Actually', cell: (r) => tag(r.derived, 'amber') },
        { label: 'Status', cell: (r) => esc(r.affected_status) },
      ], d.hemizygous, { emptyTitle: 'No hemizygous calls in this cohort' }),
        { flush: true, foot: esc(d.hemizygous_footer) })}`;
  },

  validity: (d) => {
    const k = d.kpis;
    return `
      ${headline(`${f1(k.inflation_pct)}%`,
        `yield would be inflated by counting weak-evidence genes — ${fmt(k.any_plp)} `
        + `subjects become ${fmt(k.established)} when restricted to established genes`)}
      ${note(esc(d.guidance), 'info')}
      <div class="cols c21">
        ${card('Conditions implicated', table([
          { label: 'Condition', cell: (r) => esc(r.condition) },
          { label: 'Genes', cell: (r) => `<span class="mono">${esc(r.genes)}</span>` },
          { label: 'Evidence', cell: (r) => tag(r.validity,
              r.validity === 'Definitive' ? 'teal' : r.validity === 'Strong' ? 'blue' : 'amber') },
          { label: 'Penetrance', cell: (r) => esc(r.penetrance) },
          { label: 'Subjects', n: true, cell: (r) => fmt(r.subjects) },
        ], d.conditions), { flush: true })}
        ${card('Evidence strength', stacked(d.validity_dist.map((r) =>
          ({ label: r.label, value: r.value }))), {
          foot: 'Only Definitive and Strong count toward the diagnostic yield.' })}
      </div>`;
  },

  vus: (d) => {
    const k = d.kpis;
    return `
      ${headline(fmt(k.unique_vus), `uncertain variants to review, affecting `
        + `${fmt(k.subjects_affected)} subjects`)}
      ${note(esc(d.scoping_note), 'info')}
      ${stat([
        { k: 'High priority', v: fmt(k.high_priority), d: 'score ≥ 40' },
        { k: 'New evidence', v: fmt(k.new_evidence), d: 'public data has moved' },
        { k: 'Stale', v: fmt(k.stale), d: 'over 180 days old' },
        { k: 'Already curated', v: fmt(k.curated), d: 'reviewed by a person' },
      ])}
      ${card('Review queue', table([
        { label: 'Priority', n: true, cell: (r) => `<b class="mono">${f1(r.score)}</b>` },
        { label: 'Gene', cell: (r) => `<a class="gene" data-act="gene:${esc(r.gene)}">${esc(r.gene)}</a>` },
        { label: 'Variant', cell: (r) => `<span class="mono">${esc(r.variant)}</span>` },
        { label: 'Condition', cell: (r) => esc(r.condition) },
        { label: 'Subjects', n: true, cell: (r) => fmt(r.subjects) },
        { label: 'Families', n: true, cell: (r) => fmt(r.families) },
        { label: 'Age', n: true, cell: (r) => `<span class="mono">${fmt(r.age_days)}d</span>` },
        { label: 'Flags', cell: (r) => [
            r.new_evidence ? tag('new evidence', 'teal') : '',
            r.stale ? tag('stale', 'amber') : '',
          ].filter(Boolean).join(' ') || '—' },
      ], d.rows), { flush: true, sub: 'ranked by review priority', foot: esc(d.footer) })}`;
  },

  popfreq: (d) => {
    const k = d.kpis;
    return `
      ${note(`<b>${esc(d.banner)}</b>`, 'warn')}
      ${stat([
        { k: 'Recurrent variants', v: fmt(k.recurrent_variants), d: 'seen in ≥2 subjects' },
        { k: 'Founder candidates', v: fmt(k.founder_candidates), d: 'enriched, unrelated' },
        { k: 'Implausible for P/LP', v: fmt(k.af_too_high), d: 'review the classification' },
        { k: 'Family clustered', v: fmt(k.family_clustered), d: 'not independent' },
      ])}
      ${card('Most enriched over gnomAD', table([
        { label: 'Gene', cell: (r) => `<a class="gene" data-act="gene:${esc(r.gene)}">${esc(r.gene)}</a>` },
        { label: 'Variant', cell: (r) => `<span class="mono">${esc(r.variant)}</span>` },
        { label: 'Class', cell: (r) => acmg(r.classification) },
        { label: 'Carriers', n: true, cell: (r) => fmt(r.carriers) },
        { label: 'Families', n: true, cell: (r) => fmt(r.families) },
        { label: 'Here', n: true, cell: (r) => r.internal_af
            ? `<span class="mono">${r.internal_af.toExponential(1)}</span>` : '—' },
        { label: 'gnomAD', n: true, cell: (r) => r.gnomad_af
            ? `<span class="mono">${r.gnomad_af.toExponential(1)}</span>` : '—' },
        { label: 'Ratio', n: true, cell: (r) => r.ratio
            ? `<span class="mono">${fmt(r.ratio)}×</span>` : '—' },
        { label: 'Flags', cell: (r) => r.flags.map((f) =>
            `<span class="tag t-${f.tone}" title="${esc(f.why)}">${esc(f.label)}</span>`)
            .join(' ') || '—' },
      ], d.rows), { flush: true, foot: esc(d.footer) + ' ' + esc(d.sources) })}`;
  },

  secondary: (d) => {
    const k = d.kpis;
    return `
      ${note(`<b>Consent is enforced when the cohort is built, not when it is displayed.</b>
        ${esc(d.banner)}`, 'bad')}
      ${headline(`${f1(k.rate.pct)}%`,
        `${fmt(k.sf_subjects)} of ${fmt(k.sf_capable)} eligible subjects have a secondary finding`)}
      ${card('By gene', table([
        { label: 'Gene', cell: (r) => `<a class="gene" data-act="gene:${esc(r.gene)}">${esc(r.gene)}</a>` },
        { label: 'Condition', cell: (r) => esc(r.condition) },
        { label: 'Inh', cell: (r) => esc(r.inheritance) },
        { label: 'Penetrance', cell: (r) => esc(r.penetrance) },
        { label: 'Subjects', n: true, cell: (r) => `<span class="mono">${esc(r.subjects_display)}</span>` },
        { label: 'Families', n: true, cell: (r) => `<span class="mono">${esc(r.families_display)}</span>` },
      ], d.by_gene, { emptyTitle: 'No secondary findings among consented subjects' }),
        { flush: true, foot: esc(d.suppression_note) })}`;
  },

  phenotype: (d) => `
    ${card('Gene × phenotype', d.genes.length
      ? heatmap(d.genes, d.terms.map((t) => t.term), d.matrix, d.max,
          { shortCol: (c) => c.replace(/^HP:\d+\s*/, '').split(' ')[0] })
      : empty('No pathogenic findings to cross-tabulate'), {
      sub: 'subjects with both', foot: esc(d.heatmap_footer) })}
    ${card('Phenotype-matched uncertain variants', table([
      { label: 'Subject', cell: (r) => `<a class="mono" data-act="subject:${esc(r.subject_id)}">${esc(r.subject_id)}</a>` },
      { label: 'Gene', cell: (r) => `<a class="gene" data-act="gene:${esc(r.gene)}">${esc(r.gene)}</a>` },
      { label: 'Variant', cell: (r) => `<span class="mono">${esc(r.variant)}</span>` },
      { label: 'Gene condition', cell: (r) => esc(r.condition) },
      { label: 'Subject phenotype', cell: (r) => esc(r.term) },
    ], d.queue, { emptyTitle: 'No phenotype-matched uncertain variants' }),
      { flush: true, sub: 'candidate supporting evidence', foot: esc(d.queue_footer) })}`,
};

const headline = (big, sub) => `
  <div class="card"><div class="body" style="text-align:center;padding:26px 16px">
    <div style="font-family:'IBM Plex Mono',monospace;font-size:52px;font-weight:500;
      letter-spacing:-.035em;line-height:1;color:var(--accent)">${big}</div>
    <div style="color:var(--muted);margin-top:8px;font-size:14.5px">${sub}</div>
  </div></div>`;

/* -------------------------------------------------------------- gene drill */
async function renderGene(state) {
  const res = await api.module('g_gene', state.criteria, state.name, { gene: state.gene });
  const g = res.data.detail;
  if (!g || g.error) return `<div class="wrap">${empty('Gene not found')}</div>`;
  const gd = g.gene_disease;
  const den = g.denominator || {};
  const q = QUESTIONS.find((x) => x.id === state.question);

  return `<div class="wrap">
    ${crumb([
      { label: 'All questions', act: 'ask:' },
      ...(q ? [{ label: q.title, act: 'ask:' + q.id }] : []),
      { label: g.gene },
    ])}
    ${head(g.gene, gd.condition)}
    ${den.warn ? note(`<b>Treat this rate with care.</b> ${esc(den.warn_reason)}. It is
      shown, but it is not fit for external quotation.`, 'warn') : ''}
    ${stat([
      { k: 'Carrier rate', v: `${f1(g.kpis.carrier_rate.pct)}%`,
        d: `<span class="mono">${g.kpis.carrier_rate.label}</span> tested for this gene` },
      { k: 'Carriers', v: fmt(g.kpis.plp_subjects), d: `in ${fmt(g.kpis.plp_families)} families` },
      { k: 'Uncertain', v: fmt(g.kpis.vus_subjects), d: 'subjects with a VUS' },
      { k: 'Distinct variants', v: fmt(g.kpis.unique_variants), d: 'seen in this cohort' },
      { k: 'Evidence', v: `<span style="font-size:17px">${esc(gd.validity)}</span>`,
        d: `${esc(gd.penetrance)} penetrance · ${esc(gd.inheritance)}` },
    ])}
    <div class="cols c21">
      ${card('Variants found', table([
        { label: 'Variant', cell: (r) => `<span class="mono">${esc(r.variant)}</span>` },
        { label: 'Consequence', cell: (r) => esc(r.consequence) },
        { label: 'Class', cell: (r) => acmg(r.classification) },
        { label: 'Zygosity', cell: (r) => esc(r.zygosities || '—') },
        { label: 'Subjects', n: true, cell: (r) => fmt(r.subjects) },
        { label: 'Families', n: true, cell: (r) => fmt(r.families) },
        { label: 'gnomAD', n: true, cell: (r) => r.gnomad_af
            ? `<span class="mono">${r.gnomad_af.toExponential(1)}</span>` : '—' },
      ], g.findings), { flush: true,
        actions: `<button class="btn sm" data-act="variants:${esc(g.gene)}">All observations →</button>` })}
      <div>
        ${card('Carrier rate by ancestry', table([
          { label: 'Ancestry', cell: (r) => esc(r.ancestry) },
          { label: 'Rate', n: true, cell: (r) => rate(r.rate,
              { warn: r.warn, warnReason: 'fewer than 30 subjects tested' }) },
        ], g.ancestry), { flush: true, foot: esc(g.ancestry_caveat) })}
        ${card('Who was tested for it', kv([
          ['Declared by test code', fmt(den.declared)],
          ['Inferred from coverage', fmt(den.inferred)],
          ['Provisional', fmt(den.provisional)],
          ['Not tested', fmt(den.not_assayed)],
        ]), { foot: 'This is the denominator behind the rate above.' })}
      </div>
    </div>
    <div class="actions"><button class="btn" data-act="ask:${q ? q.id : ''}">← Back</button></div>
  </div>`;
}

/* ================================================================ research == */
async function renderResearch(state) {
  const jobs = (await api.research.jobs(state.datasetId)).jobs || [];
  const available = (state.caps || []).filter((c) => c.available);

  if (state.jobId) {
    const job = await api.research.job(state.jobId);
    if (job.status === 'complete') {
      const payload = await api.research.result(state.jobId);
      return `<div class="wrap">
        ${crumb([{ label: 'Analyses', act: 'rback:' },
                 { label: analysisTitle(job.analysis) }])}
        ${head(analysisTitle(job.analysis),
               `Run on ${esc(String(job.submitted_at).slice(0, 10))}.`)}
        ${renderResult(payload)}
        <div class="actions"><button class="btn" data-act="rback:">← Back to analyses</button>
          <span class="spacer"></span>
          <button class="btn primary" data-act="goto:share">Export →</button></div>
      </div>`;
    }
    return `<div class="wrap">
      ${head('Running…', 'Analyses run in the background. This page updates itself.')}
      ${card(job.analysis, `
        ${note(`<b>${esc(job.status)}</b> — submitted ${esc(String(job.submitted_at).slice(0, 19))}`,
          job.status === 'failed' ? 'bad' : 'info')}
        ${job.error ? `<div class="codeblock">${esc(job.error)}</div>` : ''}
        <pre class="codeblock">${esc(job.log || 'waiting…')}</pre>`)}
      <div class="actions"><button class="btn" data-act="rback:">← Back to analyses</button></div>
    </div>`;
  }

  return `<div class="wrap">
    ${head('What do you want to run?',
      'Only analyses this dataset can support are shown. Each one estimates its power '
      + 'before it runs.', 'Step 3 of 4')}
    <div class="qgrid">${available.map((c) => `
      <button class="q" data-act="rrun:${esc(c.analysis)}">
        <span class="qt">${esc(c.title)}</span>
        <span class="qd">${c.requirements.map((r) => esc(r.observed)).slice(0, 2).join(' · ')}</span>
        <span class="qa">Configure →</span></button>`).join('')}</div>
    ${jobs.length ? card('Previous runs', table([
      { label: 'Analysis', cell: (j) => esc(j.analysis) },
      { label: 'Status', cell: (j) => tag(j.status, j.status === 'complete' ? 'green'
          : j.status === 'failed' ? 'red' : 'blue') },
      { label: '', cell: (j) => j.underpowered ? tag('UNDERPOWERED', 'amber') : '' },
      { label: 'Submitted', cell: (j) => `<span class="mono hint">${esc(String(j.submitted_at).slice(0, 19))}</span>` },
      { label: '', cell: (j) => j.status === 'complete'
          ? `<button class="btn sm" data-act="rres:${esc(j.job_id)}">View</button>` : '' },
    ], jobs), { flush: true }) : ''}
    <div class="actions"><button class="btn" data-act="goto:review">← Back to capabilities</button></div>
  </div>`;
}

/* The registry key is an identifier, not a label. Showing it twice — once in
   the breadcrumb and again as the heading — told the reader nothing and looked
   like a bug, which it effectively was. */
function analysisTitle(key) {
  const form = FORMS[key];
  if (form && form.title) return form.title;
  return String(key || '').replace(/_/g, ' ').replace(/^./, (c) => c.toUpperCase());
}

function renderResult(payload) {
  const r = payload.result || {};
  const stamp = r.underpowered_stamp;
  const head2 = stamp
    ? note(`<b>${esc(stamp)}</b><br>${esc(r.override_justification || '')}`, 'warn') : '';
  const body = RESULTS[r.analysis];
  return head2 + (body ? body(r, payload) : genericResult(r));
}


/* ---- the R1 descriptive results -------------------------------------------
 * Each one leads with the sentence a reader would say out loud, then the
 * evidence, then the caveat. A number without its denominator never appears. */
function pct(x, dp) { return x === null || x === undefined ? '—'
  : (Number(x) * 100).toFixed(dp === undefined ? 1 : dp) + '%'; }

/* kit.rate takes {n, d, pct} — passing two bare numbers rendered
   "undefined / 0 · no denominator" on every row. */
function rateOf(n, d) {
  if (!d) return '—';
  return rate({ n: fmt(n), d: fmt(d), pct: (Math.round(1000 * n / d) / 10) });
}

const CARRIER_LABELS = {
  dominant: 'het or hom',
  het_only: 'heterozygotes',
  recessive: 'homozygotes',
};

/* The frequency spectrum is the shape of the data: rare variants dominate any
   real call set, and seeing that is worth more than reading the top of a table. */
function afSpectrum(rows) {
  if (!rows || rows.length < 8) return '';
  const vals = rows.map((r) => r.allele_freq).filter((v) => isFinite(v));
  if (!vals.length) return '';
  const lo = Math.min(...vals);
  const hi = Math.max(...vals);
  if (!(hi > lo)) return '';
  const n = 24;
  const step = (hi - lo) / n;
  const bins = Array.from({ length: n }, (_, i) => ({
    lo: lo + i * step, hi: lo + (i + 1) * step, count: 0,
  }));
  for (const v of vals) {
    bins[Math.min(n - 1, Math.floor((v - lo) / step))].count += 1;
  }
  return histogram(bins, {
    title: 'Allele frequency spectrum',
    sub: `${fmt(vals.length)} variants shown`,
    xLabel: 'alternate allele frequency',
    yLabel: 'variants',
    foot: 'The shape of the reported set, not of the whole dataset — the frequency '
      + 'filter you chose decides which part of the spectrum appears here.',
  });
}

function carrierFreqResult(d) {
  const top = d.rows[0];
  return `
    ${stat([
      { k: 'Variants reported', v: fmt(d.n_variants_reported),
        d: `of ${fmt(d.n_variants_scanned)} scanned` },
      { k: 'Samples', v: fmt(d.n_samples), d: 'in the dataset' },
      { k: 'Excluded, low call rate', v: fmt(d.skipped_low_call_rate),
        d: `under ${d.min_calls_required} calls` },
      { k: 'Carrier definition', v: esc(CARRIER_LABELS[d.carrier_model] || d.carrier_model),
        d: 'as configured' },
    ])}
    ${afSpectrum(d.rows)}
    ${top ? `<p>The most frequent allele is <b>${esc(top.variant)}</b>${
      top.gene ? ` in <span class="gene">${esc(top.gene)}</span>` : ''} at
      <b>${pct(top.allele_freq, 2)}</b>, carried by ${fmt(top.carriers)} of
      ${fmt(top.n_called)} genotyped samples.</p>` : ''}
    ${table([
      { label: 'Variant', cell: (r) => `<span class="mono">${esc(r.variant)}</span>` },
      { label: 'Gene', cell: (r) => r.gene ? `<span class="gene">${esc(r.gene)}</span>` : '—' },
      { label: 'Carriers', n: true, cell: (r) => r.suppressed
          ? `<span class="hint" title="suppressed: fewer than 5">&lt;5</span>`
          : fmt(r.carriers) },
      { label: 'Called', n: true, cell: (r) => fmt(r.n_called) },
      { label: 'Carrier rate', n: true, cell: (r) => r.suppressed ? '—'
          : rateOf(r.carriers, r.n_called) },
      { label: 'Allele freq', n: true, cell: (r) => `<span class="mono">${pct(r.allele_freq, 3)}</span>` },
      { label: 'Het / HomAlt', cell: (r) => `<span class="mono">${fmt(r.het)} / ${fmt(r.hom_alt)}</span>` },
      { label: 'Call rate', n: true, cell: (r) => pct(r.call_rate, 0) },
    ], d.rows)}
    ${d.truncated ? note(`Showing the ${d.rows.length} most frequent. Narrow the
      frequency range to see the rest.`, 'info') : ''}
    ${note(esc(d.method), 'info')}`;
}

/* Heterozygous against homozygous count, per sample. A contaminated sample sits
   visibly off the cloud in a way no column of numbers conveys — which is the
   entire argument for plotting it. */
function zygScatter(d) {
  const ss = (d.samples || []).filter((s) => isFinite(s.het) && isFinite(s.hom_alt));
  if (ss.length < 5) return '';
  const med = d.median_het_hom_ratio;
  return scatter(ss.map((s) => ({
    x: s.hom_alt, y: s.het,
    flag: med && s.het_hom_ratio && (s.het_hom_ratio > med * 1.5
                                     || s.het_hom_ratio < med * 0.6),
    label: `${s.sample_id} · het ${fmt(s.het)} · hom ${fmt(s.hom_alt)}`
      + (s.het_hom_ratio ? ` · ratio ${Number(s.het_hom_ratio).toFixed(2)}` : ''),
  })), {
    title: 'Heterozygous against homozygous calls',
    sub: med ? `cohort median ratio ${Number(med).toFixed(2)}` : '',
    xLabel: 'homozygous ALT', yLabel: 'heterozygous',
    flagLabel: 'far from the cohort ratio',
    baseLabel: 'consistent with the cohort',
    foot: 'Samples sit on a line whose slope is the het/hom ratio. One sitting off '
      + 'that line is the finding.',
  });
}

function zygosityResult(d) {
  return `
    ${stat([
      { k: 'Samples', v: fmt(d.n_samples), d: 'genotyped' },
      { k: 'Variants', v: fmt(d.n_variants), d: 'in the dataset' },
      { k: 'Median het/hom', v: d.median_het_hom_ratio === null ? '—'
          : Number(d.median_het_hom_ratio).toFixed(2), d: 'across samples' },
      { k: 'Sex column', v: d.sex_available ? 'present' : 'absent',
        d: d.sex_available ? 'X calls interpretable' : 'X reads as autosomal' },
    ])}
    <p>A sample far from the median het/hom ratio is worth investigating before
    its genotypes are used: high suggests contamination or sample mixture, low
    suggests a homozygosity-rich background or a failed call set.</p>
    ${zygScatter(d)}
    <div class="cols c2">
      <div><p class="hint" style="margin-top:0">Per sample</p>
      ${table([
        { label: 'Sample', cell: (r) => `<span class="mono">${esc(r.sample_id)}</span>` },
        { label: 'Het', n: true, cell: (r) => fmt(r.het) },
        { label: 'Hom ALT', n: true, cell: (r) => fmt(r.hom_alt) },
        { label: 'Het/Hom', n: true, cell: (r) => r.het_hom_ratio === null ? '—'
            : `<span class="mono">${Number(r.het_hom_ratio).toFixed(2)}</span>` },
        { label: 'Call rate', n: true, cell: (r) => pct(r.call_rate, 0) },
      ], d.samples)}</div>
      <div><p class="hint" style="margin-top:0">Per variant</p>
      ${table([
        { label: 'Variant', cell: (r) => `<span class="mono">${esc(r.variant)}</span>` },
        { label: 'Gene', cell: (r) => r.gene ? `<span class="gene">${esc(r.gene)}</span>` : '—' },
        { label: 'HomRef', n: true, cell: (r) => fmt(r.hom_ref) },
        { label: 'Het', n: true, cell: (r) => fmt(r.het) },
        { label: 'HomAlt', n: true, cell: (r) => fmt(r.hom_alt) },
        { label: '% carriers hom', n: true,
          cell: (r) => pct(r.pct_carriers_homozygous, 0) },
      ], d.variants)}</div>
    </div>
    ${note(esc(d.method), 'info')}`;
}

/* Before-and-after per item is a dumbbell: one hue in two shades, joined. Two
   hues would spend the identity channel on what position already shows. */
function groupDumbbell(d) {
  const labels = Object.keys(d.groups || {});
  if (labels.length !== 2) return '';          // a dumbbell compares exactly two
  const [a, b] = labels;
  const rows = (d.rows || [])
    .filter((r) => r.by_group && r.by_group[a] && r.by_group[b]
                   && isFinite(r.by_group[a].allele_freq)
                   && isFinite(r.by_group[b].allele_freq))
    .slice(0, 12);
  if (rows.length < 3) return '';
  return dumbbell(rows.map((r) => ({
    label: r.gene ? `${r.gene}` : r.variant,
    a: r.by_group[a].allele_freq, b: r.by_group[b].allele_freq,
    aLabel: a, bLabel: b,
    note: r.max_group_difference !== null && r.max_group_difference !== undefined
      ? `Δ ${(r.max_group_difference * 100).toFixed(2)}pp` : '',
  })), {
    title: 'Where the groups disagree most',
    sub: `top ${rows.length} by frequency difference`,
    fmt: (v) => (v * 100).toFixed(1) + '%',
    foot: 'These are internal frequencies in this cohort, not population reference '
      + 'values. A referral-selected cohort is enriched for exactly the alleles it '
      + 'was selected on.',
  });
}

function popFreqResult(d) {
  const groups = Object.entries(d.groups || {});
  const labels = groups.map(([k]) => k);
  return `
    ${stat([
      { k: 'Grouped by', v: esc(d.grouped_by || 'nothing'),
        d: d.grouped_by ? `${groups.length} groups` : 'whole cohort' },
      { k: 'Variants', v: fmt(d.n_variants_reported), d: 'with enough calls' },
      ...groups.slice(0, 3).map(([k, v]) => ({ k: esc(k), v: fmt(v), d: 'samples' })),
    ])}
    ${groupDumbbell(d)}
    ${d.grouped_by ? `<p>Sorted by the largest frequency difference between groups —
      the variants at the top are the ones whose rarity depends on which population
      you ask about.</p>` : ''}
    ${table([
      { label: 'Variant', cell: (r) => `<span class="mono">${esc(r.variant)}</span>` },
      { label: 'Gene', cell: (r) => r.gene ? `<span class="gene">${esc(r.gene)}</span>` : '—' },
      { label: 'Overall AF', n: true,
        cell: (r) => `<span class="mono">${pct(r.overall_allele_freq, 3)}</span>` },
      ...labels.map((lab) => ({
        label: lab, n: true,
        cell: (r) => {
          const g = r.by_group[lab];
          if (!g || g.allele_freq === null) return '—';
          return `<span class="mono">${pct(g.allele_freq, 3)}</span>
            <span class="hint">n=${fmt(g.n_called)}</span>`;
        },
      })),
      { label: 'Spread', n: true, cell: (r) => r.max_group_difference === null ? '—'
          : `<span class="mono">${pct(r.max_group_difference, 2)}</span>` },
    ], d.rows)}
    ${note(esc(d.method), 'warn')}`;
}

function yieldResult(d) {
  const o = d.overall;
  return `
    ${stat([
      { k: 'Subjects', v: fmt(o.subjects), d: 'with genotypes' },
      { k: 'Carrying a qualifying variant', v: fmt(o.with_qualifying_variant),
        d: 'at least one' },
      { k: 'Candidate rate', v: pct(o.yield), d: 'upper bound' },
      { k: 'Qualifying variants', v: fmt(d.n_qualifying_variants),
        d: `AF ≤ ${pct(d.max_allele_frequency, 2)}` },
    ])}
    <p><b>${fmt(o.with_qualifying_variant)} of ${fmt(o.subjects)} subjects</b>
    (${pct(o.yield)}) carry at least one variant below
    ${pct(d.max_allele_frequency, 2)} allele frequency.</p>
    ${d.grouped_by ? `${bars(d.groups.filter((g) => !g.suppressed).map((g) => ({
      label: `${g.group} — ${g.with_qualifying_variant}/${g.subjects}`,
      value: (g.yield || 0) * 100,
      display: pct(g.yield),
    })))}
    ${table([
      { label: esc(d.grouped_by), cell: (g) => esc(g.group) },
      { label: 'Subjects', n: true, cell: (g) => fmt(g.subjects) },
      { label: 'With a variant', n: true, cell: (g) => g.suppressed ? '&lt;5'
          : fmt(g.with_qualifying_variant) },
      { label: 'Rate', n: true, cell: (g) => g.suppressed ? '—'
          : rateOf(g.with_qualifying_variant, g.subjects) },
    ], d.groups)}` : ''}
    ${note(esc(d.method), 'warn')}`;
}

/* One dot per trio against the cohort median. The question this answers is
   "does one trio behave unlike the others", which is a shape question. */
function trioPlot(d) {
  const ts = (d.trios || []).filter((t) => isFinite(t.rate));
  if (ts.length < 3) return '';
  const med = d.median_inconsistency_rate;
  return dotplot(ts.map((t) => ({
    value: t.rate,
    flag: med && t.rate > med * 3,
    label: `${t.child} · ${fmt(t.inconsistent)} of ${fmt(t.variants_compared)} `
      + `(${(t.rate * 100).toFixed(3)}%)`,
  })), {
    title: 'Mendelian inconsistency by trio',
    sub: `${fmt(ts.length)} complete trios`,
    yLabel: 'inconsistent share',
    xLabel: 'one dot per trio',
    fmtY: (v) => (v * 100).toFixed(2) + '%',
    reference: med,
    referenceLabel: 'cohort median',
    foot: 'A trio well above the median is usually a sample swap or a mislabelled '
      + 'pedigree, not a mutator phenotype.',
  });
}

function segregationResult(d) {
  return `
    ${stat([
      { k: 'Complete trios', v: fmt(d.n_trios), d: 'child + both parents' },
      { k: 'Inconsistent sites', v: fmt(d.n_findings), d: 'across all trios' },
      { k: 'Median rate', v: d.median_inconsistency_rate === null ? '—'
          : pct(d.median_inconsistency_rate, 3), d: 'per trio' },
    ])}
    <p>A trio whose rate sits well above the median is the finding here — that
    pattern is a sample swap or a mislabelled pedigree far more often than it is
    biology.</p>
    ${trioPlot(d)}
    ${table([
      { label: 'Child', cell: (t) => `<span class="mono">${esc(t.child)}</span>` },
      { label: 'Father', cell: (t) => `<span class="mono">${esc(t.father)}</span>` },
      { label: 'Mother', cell: (t) => `<span class="mono">${esc(t.mother)}</span>` },
      { label: 'Compared', n: true, cell: (t) => fmt(t.variants_compared) },
      { label: 'Inconsistent', n: true, cell: (t) => fmt(t.inconsistent) },
      { label: 'Rate', n: true, cell: (t) => t.rate === null ? '—'
          : `<span class="mono">${pct(t.rate, 3)}</span>` },
    ], d.trios)}
    ${d.findings.length ? `<p class="hint" style="margin-top:16px">Individual sites</p>
    ${table([
      { label: 'Child', cell: (f) => `<span class="mono">${esc(f.child)}</span>` },
      { label: 'Variant', cell: (f) => `<span class="mono">${esc(f.variant)}</span>` },
      { label: 'Gene', cell: (f) => f.gene ? `<span class="gene">${esc(f.gene)}</span>` : '—' },
      { label: 'Child', n: true, cell: (f) => f.child_dosage },
      { label: 'Father', n: true, cell: (f) => f.father_dosage },
      { label: 'Mother', n: true, cell: (f) => f.mother_dosage },
    ], d.findings)}` : ''}
    ${note(esc(d.method), 'warn')}`;
}

const RESULTS = {
  association: assocResult,
  gwas: assocResult,
  burden: burdenResult,
  survival: survivalResult,
  prs: prsResult,
  carrier_frequency: carrierFreqResult,
  zygosity: zygosityResult,
  population_frequency: popFreqResult,
  diagnostic_yield: yieldResult,
  segregation: segregationResult,
};

function assocResult(r) {
  const rows = r.results || r.top_hits || [];
  const infl = r.inflation_guardrail;
  const hits = rows.filter((x) => (x.p_bonferroni ?? 1) < 0.05).length;

  /* A genome-wide scan gets the Manhattan it was computed for, paired with the
     Q–Q that says whether the scan can be believed at all. The pair belongs
     together: a tall peak means nothing if the whole distribution is lifted. */
  const lam = (r.inflation || {}).lambda_gc;
  const charts = (r.manhattan && r.manhattan.length)
    ? `<div class="chart-grid c-manhattan">
        ${manhattan(r.manhattan, {
          sub: `${fmt(r.n_variants_tested)} variants · ${fmt(r.n_samples)} samples`,
          foot: 'Each point is one variant, ordered along the genome. Points above '
            + 'the dashed red line cleared genome-wide significance. Hover any point '
            + 'for its identity and p-value.',
        })}
        ${qq((r.qq || {}).observed || [], (r.qq || {}).expected || [], {
          sub: lam ? `λ_GC = ${Number(lam).toFixed(3)}` : '',
          foot: 'Points should follow the diagonal until the real associations '
            + 'lift off at the top. A curve that departs early is inflation, not '
            + 'discovery.',
        })}
      </div>`
    : volcano(rows.map((x) => ({
        x: x.beta, y: -Math.log10(Math.max(x.pvalue || 1, 1e-300)),
        name: x.variant, label: `${x.variant} · p = ${(x.pvalue ?? 1).toExponential(2)}`
          + (x.effect ? ` · ${x.effect_label || 'effect'} ${Number(x.effect).toFixed(2)}` : ''),
      })), {
        title: 'Effect against significance',
        sub: `${fmt(rows.length)} variants tested`,
        xLabel: 'log odds (β)',
        threshold: -Math.log10(r.alpha && r.alpha.alpha ? r.alpha.alpha : 0.05),
        foot: 'Points above the dashed line cleared the significance threshold for '
          + 'this many tests. Distance from the centre line is effect size.',
      });

  return charts + (infl && infl.status !== 'ok'
      ? note(`<b>Possible population stratification.</b> ${esc(infl.text)}`,
             infl.status === 'severe' ? 'bad' : 'warn') : '')
    + (r.power ? powerNote(r.power) : '')
    + headline(fmt(hits), `variant${hits === 1 ? '' : 's'} significant after correcting for `
        + `${fmt(r.n_variants_tested)} tests`)
    + card('Strongest associations', table([
      { label: 'Variant', cell: (x) => `<span class="mono">${esc(x.variant)}</span>` },
      { label: 'Model', cell: (x) => esc(x.model) },
      { label: 'N', n: true, cell: (x) => fmt(x.n) },
      { label: 'Effect', n: true, cell: (x) =>
          `<span class="mono">${esc(x.effect_label)} ${x.effect.toFixed(3)}</span>` },
      { label: '95% CI', n: true, cell: (x) =>
          `<span class="mono">${x.effect_ci_low.toFixed(2)}–${x.effect_ci_high.toFixed(2)}</span>` },
      { label: 'p', n: true, cell: (x) => `<span class="mono">${x.pvalue.toExponential(2)}</span>` },
      { label: 'p corrected', n: true, cell: (x) =>
          `<span class="mono">${(x.p_bonferroni ?? 1).toExponential(2)}</span>` },
      { label: '', cell: (x) => x.converged ? '' : tag('did not converge', 'red') },
    ], rows.slice(0, 200)), { flush: true, foot: esc(r.multiple_testing_note || '') });
}

function burdenResult(r) {
  const sig = (r.results || []).filter((x) => (x.p_fdr_bh ?? 1) < 0.05).length;
  const gv = (r.results || []).filter((x) => x.p !== null && x.p !== undefined);
  return volcano(gv.map((x) => ({
      // Burden collapses a gene to one direction, so the x-axis is the excess of
      // carriers in cases — the quantity the test is actually about.
      x: (x.n_carriers_cases || 0) - (x.n_carriers_controls || 0),
      y: -Math.log10(Math.max(x.p, 1e-300)),
      name: x.gene,
      label: `${x.gene} · ${fmt(x.n_carriers)} carriers `
        + `(${fmt(x.n_carriers_cases)} case / ${fmt(x.n_carriers_controls)} control) `
        + `· p = ${x.p.toExponential(2)}`,
    })), {
      title: 'Genes by carrier excess and significance',
      sub: `${fmt(r.n_genes_tested)} genes · ${fmt(r.n_qualifying_variants)} qualifying variants`,
      xLabel: 'carriers in cases − controls',
      threshold: gv.length ? -Math.log10(0.05 / gv.length) : undefined,
      foot: 'Right of centre means the gene carries more variant carriers among '
        + 'cases. Height is significance; the dashed line is the Bonferroni '
        + 'threshold for the number of genes tested.',
    })
    + headline(fmt(sig), `gene${sig === 1 ? '' : 's'} significant at 5% FDR, of `
      + `${fmt(r.n_genes_tested)} tested`)
    + note(esc(r.guardrail_note), 'info')
    + card('Genes', table([
      { label: 'Gene', cell: (x) => `<span class="gene">${esc(x.gene)}</span>` },
      { label: 'Variants', n: true, cell: (x) => fmt(x.n_variants) },
      { label: 'Carriers', n: true, cell: (x) => fmt(x.carriers) },
      { label: 'In cases', n: true, cell: (x) => fmt(x.carriers_cases) },
      { label: 'In controls', n: true, cell: (x) => fmt(x.carriers_controls) },
      { label: 'p', n: true, cell: (x) => `<span class="mono">${(x.p ?? 1).toExponential(2)}</span>` },
      { label: 'FDR', n: true, cell: (x) => `<span class="mono">${(x.p_fdr_bh ?? 1).toExponential(2)}</span>` },
    ], r.results || []), { flush: true,
      sub: `${esc(r.test)} · ${esc(r.variant_definition)}`,
      foot: esc(r.multiple_testing_note || '') })
    + (r.skipped && r.skipped.length ? card('Not tested — too few carriers', table([
        { label: 'Gene', cell: (x) => `<span class="gene">${esc(x.gene)}</span>` },
        { label: 'Carriers', n: true, cell: (x) => fmt(x.carriers) },
      ], r.skipped), { flush: true,
        foot: 'An untested gene is not evidence of no effect.' }) : '');
}

/* Survival used to render as a raw JSON dump of the Kaplan-Meier arrays —
   thousands of float literals, which is not a result. The numbers a reader
   actually needs are the medians, the log-rank p, the hazard ratios and
   whether the proportional-hazards assumption held. */
/* Axis labels are read by people; a raw column name like followup_years is a
   database identifier that happened to leak onto the page. */
function humanCol(name) {
  if (!name) return '';
  return String(name).replace(/[_-]+/g, ' ').trim();
}

function survivalResult(r) {
  const km = r.kaplan_meier || {};
  const groups = km.groups || {};
  const names = km.group_names || Object.keys(groups);
  const lr = r.logrank || {};
  const cox = r.cox || {};
  const ph = r.ph_assumption || {};
  const pen = r.penetrance || {};

  const medians = names.map((n) => {
    const g = groups[n] || {};
    return { group: n, n: g.n, events: g.n_events,
             median: g.median, ci: g.median_ci };
  });

  const coxRows = (cox.names || []).map((name, i) => ({
    name,
    hr: (cox.hazard_ratio || [])[i],
    beta: (cox.beta || [])[i],
    se: (cox.se || [])[i],
    p: (cox.p || [])[i],
    phP: ((ph.covariates || {})[name] || {}).p,
  }));

  return kaplanMeier(groups, {
      title: names.length > 1 ? 'Survival by group' : 'Survival',
      sub: lr && lr.p !== undefined && lr.p !== null
        ? `log-rank p = ${Number(lr.p).toExponential(2)}` : '',
      xLabel: humanCol(r.time_column) || 'time',
      foot: 'Each drop is an event. The curve is drawn as a step because survival '
        + 'changes only when an event occurs — a smooth line would imply deaths '
        + 'between events that did not happen.',
    })
    + (pen.ascertained
      ? note(`<b>${esc(pen.stamp)}</b><br>${esc(pen.caveat)}`, 'warn')
      : note(esc(pen.caveat), 'info'))
    + headline(fmt(r.n_events), `events among ${fmt(r.n)} subjects followed`)
    + card('Median survival by group', table([
        { label: 'Group', cell: (x) => esc(x.group) },
        { label: 'Subjects', n: true, cell: (x) => fmt(x.n) },
        { label: 'Events', n: true, cell: (x) => fmt(x.events) },
        { label: 'Median', n: true, cell: (x) => x.median === null || x.median === undefined
            ? '<span class="hint">not reached</span>'
            : `<span class="mono">${Number(x.median).toFixed(2)}</span>` },
        { label: '95% CI', n: true, cell: (x) => Array.isArray(x.ci) && x.ci[0] != null
            ? `<span class="mono">${Number(x.ci[0]).toFixed(2)}–${
                x.ci[1] == null ? '∞' : Number(x.ci[1]).toFixed(2)}</span>`
            : '—' },
      ], medians, { emptyTitle: 'No groups' }), {
      flush: true,
      foot: lr.p !== undefined
        ? `Log-rank test: χ² = ${Number(lr.chi2).toFixed(3)} on ${lr.df} df, `
          + `p = ${Number(lr.p).toExponential(3)}. `
          + (lr.p < 0.05 ? 'The groups differ.'
             : 'No detectable difference between the groups.')
        : 'No grouping was applied, so no log-rank test was run.' })
    + (coxRows.length ? card('Cox proportional hazards', table([
        { label: 'Covariate', cell: (x) => esc(x.name) },
        { label: 'Hazard ratio', n: true, cell: (x) =>
            `<span class="mono">${Number(x.hr).toFixed(3)}</span>` },
        { label: 'log HR', n: true, cell: (x) =>
            `<span class="mono">${Number(x.beta).toFixed(4)}</span>` },
        { label: 'SE', n: true, cell: (x) =>
            `<span class="mono">${Number(x.se).toFixed(4)}</span>` },
        { label: 'p', n: true, cell: (x) =>
            `<span class="mono">${Number(x.p).toExponential(2)}</span>` },
        { label: 'PH holds?', cell: (x) => x.phP === undefined ? '—'
            : x.phP < 0.05 ? tag('violated (p=' + Number(x.phP).toFixed(3) + ')', 'amber')
                           : tag('ok', 'green') },
      ], coxRows), {
      flush: true,
      sub: `${esc(cox.ties || 'efron')} ties · ${cox.converged ? 'converged' : 'did not converge'}`
           + ` in ${fmt(cox.n_iter)} iterations`,
      foot: `Model likelihood-ratio test: χ² = ${Number(cox.lr_chi2).toFixed(3)}, `
          + `p = ${Number(cox.lr_p).toExponential(3)}. `
          + (ph.p_global !== undefined
             ? `Global proportional-hazards test p = ${Number(ph.p_global).toFixed(3)}`
               + (ph.p_global < 0.05
                  ? ' — the assumption is violated, so these hazard ratios are an average over follow-up.'
                  : ' — the assumption holds.')
             : '') }) : '');
}

function prsResult(r) {
  const cal = r.calibration || {};
  const strat = r.ancestry_stratified || {};
  const d = r.distribution || {};
  return histogram(d.histogram || [], {
      title: 'Score distribution',
      sub: `mean ${Number(d.mean).toFixed(3)} · SD ${Number(d.sd).toFixed(3)}`,
      xLabel: 'polygenic score',
      foot: 'A polygenic score is only interpretable against the distribution it '
        + 'sits in — an absolute value carries no meaning on its own.',
    })
    + note(esc(r.coverage_note), 'info')
    + headline(fmt(r.n_variants_used), `of ${fmt(r.n_variants_requested)} score variants `
        + `found in this dataset`)
    + card('Score distribution', kv([
        ['Mean', Number(d.mean).toFixed(4)],
        ['Standard deviation', Number(d.sd).toFixed(4)],
        ['Range', `${Number(d.min).toFixed(3)} to ${Number(d.max).toFixed(3)}`],
        ['Variants allele-flipped', fmt(r.n_flipped)],
        ['Variants not matched', fmt(r.n_unmatched)],
      ]))
    + (Object.keys(cal).length ? card('Performance', kv(
        Object.entries(cal).filter(([, v]) => typeof v !== 'object')
          .map(([k, v]) => [k.replace(/_/g, ' '),
                            typeof v === 'number' ? v.toFixed(4) : esc(String(v))]))) : '')
    + card('Performance by ancestry', strat.available
        ? table([
            { label: 'Group', cell: (x) => esc(x.group) },
            { label: 'N', n: true, cell: (x) => fmt(x.n) },
            { label: 'AUC', n: true, cell: (x) => x.auc
                ? `<span class="mono">${x.auc.toFixed(3)}</span>` : '—' },
            { label: 'OR per SD', n: true, cell: (x) => x.or_per_sd
                ? `<span class="mono">${x.or_per_sd.toFixed(3)}</span>` : '—' },
            { label: 'Mean score', n: true, cell: (x) => x.mean_z !== undefined
                ? `<span class="mono">${x.mean_z.toFixed(3)}</span>` : '—' },
            { label: '', cell: (x) => x.note ? `<span class="hint">${esc(x.note)}</span>` : '' },
          ], strat.groups || [])
        : note(esc(strat.reason || ''), 'warn'),
      { flush: strat.available, foot: esc(r.ancestry_note || '') });
}

/* Last resort for an analysis with no dedicated renderer. Still a table of
   scalars rather than a JSON dump — long numeric arrays are summarised, never
   printed. */
function genericResult(r) {
  const scalars = Object.entries(r)
    .filter(([, v]) => v === null || ['string', 'number', 'boolean'].includes(typeof v))
    .map(([k, v]) => [k.replace(/_/g, ' '), esc(String(v))]);
  const arrays = Object.entries(r)
    .filter(([, v]) => Array.isArray(v))
    .map(([k, v]) => [k.replace(/_/g, ' '), `${fmt(v.length)} rows`]);
  return card('Summary', kv(scalars.concat(arrays)), {
    foot: 'This analysis has no dedicated view yet; the values above are its scalar outputs.',
  });
}

function powerNote(p) {
  if (p.power === null || p.power === undefined) return '';
  const tone = p.verdict.level === 'adequate' ? 'good'
    : p.verdict.level === 'marginal' ? 'warn' : 'bad';
  return note(`<b>Estimated power: ${p.power_pct}%</b><br>
    <span class="mono hint">${esc(p.summary)}</span><br>${esc(p.verdict.text)}
    ${p.verdict.options.length
      ? `<br><span class="hint">Options: ${p.verdict.options.map(esc).join(' · ')}</span>` : ''}`,
    tone);
}

/* ---------------------------------------------------------------- handlers */
export async function handle(verb, arg, rest, state, ctx) {
  const { render: rerender, go } = ctx;
  switch (verb) {
    case 'ask':
      state.question = arg || null;
      state.gene = null;
      await rerender();
      return true;
    case 'gene':
      state.gene = arg;
      await rerender();
      return true;
    case 'variants':
      await variantsModal(state, arg);
      return true;
    case 'subject':
      await subjectModal(state, arg);
      return true;
    case 'rback':
      state.jobId = null;
      await rerender();
      return true;
    case 'rres':
      state.jobId = arg;
      await rerender();
      return true;
    case 'rrun':
      await configureModal(state, arg);
      return true;
    case 'rovr':
      await overrideModal(state, arg);
      return true;
    default: return false;
  }
}

async function variantsModal(state, gene) {
  modal(`${gene} — every observation`, '<div class="loading">Loading…</div>', { wide: true });
  const res = await api.module('g_variants', state.criteria, state.name, { gene });
  const d = res.data;
  modal(`${gene} — every observation`, table([
    { label: 'Variant', cell: (r) => `<span class="mono">${esc(r.variant)}</span>` },
    { label: 'Consequence', cell: (r) => esc(r.consequence) },
    { label: 'Class', cell: (r) => acmg(r.classification) },
    { label: 'ACMG codes', cell: (r) => `<span class="mono hint">${esc(r.acmg_codes || '—')}</span>` },
    { label: 'Zygosity', cell: (r) => tag(r.zygosity, ZYG_TONE[r.zygosity] || 'grey') },
    { label: 'Subject', cell: (r) => `<a class="mono" data-act="subject:${esc(r.subject_id)}">${esc(r.subject_id)}</a>` },
    { label: 'Reportable', cell: (r) => r.reportable ? tag('yes', 'green') : tag('no', 'amber') },
  ], d.rows), { wide: true, footer: `<span class="hint">${esc(d.footer)} ${esc(d.handoff_note)}</span>` });
}

async function subjectModal(state, id) {
  modal(`Subject ${id}`, '<div class="loading">Loading…</div>', { wide: true });
  const r = await api.subjectDetail(id, state.criteria, state.name);
  modal(`Subject ${id}`, `
    <div class="cols c2">
      <div>${kv([
        ['Family', `${esc(r.family_id)} · ${r.family.length} member(s)`],
        ['Relation', r.is_proband ? `<b>${esc(r.relation)}</b>` : esc(r.relation)],
        ['Sex / age', `${esc(r.sex)} · ${esc(r.age)}`],
        ['Ancestry', esc(r.ancestry)],
        ['Reason for testing', esc(r.indication)],
        ['Affected status', esc(r.affected_status)],
        ['Consent', esc(r.consent_class)],
        ['Result', tag(r.result, r.result === 'P/LP' ? 'red'
          : r.result === 'VUS-only' ? 'grey' : 'green')],
      ])}</div>
      <div>${kv([
        ['Test code', `${esc(r.test_code || '— inferred —')} ${esc(r.test_name || '')}`],
        ['Scope', tag(r.scope_status, r.scope_status === 'Complete' ? 'teal' : 'amber')],
        ['Sample type', esc(r.sample_type)],
        ['Pipeline', esc(r.pipeline_version)],
        ['Reference', esc(r.reference_build)],
        ['Mean depth', `${fmt(Math.round(r.mean_depth || 0))}×`],
        ['QC', esc(r.qc_status)],
        ['Collected', esc(String(r.collection_date).slice(0, 10))],
      ])}</div>
    </div>
    <h4 style="margin:18px 0 8px">Findings</h4>
    ${table([
      { label: 'Gene', cell: (f) => `<span class="gene">${esc(f.gene)}</span>` },
      { label: 'Variant', cell: (f) => `<span class="mono">${esc(f.variant)}</span>` },
      { label: 'Class', cell: (f) => acmg(f.classification) },
      { label: 'Zygosity', cell: (f) => esc(f.zygosity) },
      { label: 'Condition', cell: (f) => esc(f.condition) },
      { label: 'Reportable', cell: (f) => f.reportable ? tag('yes', 'green')
          : tag('outside scope', 'amber') },
    ], r.findings, { emptyTitle: 'No findings' })}`, { wide: true });
}

/* Per-analysis configuration.
 *
 * There was one shared form here, sending {outcome, covariates, n_pcs, min_maf}
 * to every analysis. Association wanted exactly that; survival needs time and
 * event columns; a gene-based test needs a qualifying-variant definition; a PRS
 * needs score weights. So three of the five analyses failed the moment anyone
 * ran them from the UI — the form offered choices the analysis ignored and
 * omitted the ones it required.
 *
 * Each analysis now declares its own fields and builds its own spec. */

/* A variant key carries colons (chrom:pos:ref:alt); a gene symbol does not.
   That is enough to tell the two apart without asking the user which they
   typed. */
function parseTarget(raw) {
  const parts = String(raw || '').split(/[\s,]+/).filter(Boolean);
  const variants = parts.filter((x) => x.includes(':'));
  const genes = parts.filter((x) => !x.includes(':'));
  return { targetVariants: variants, targetGene: genes[0] || '' };
}

const FORMS = {
  association: {
    title: 'Test specific variants against an outcome',
    fields: ['target', 'outcome', 'covariates', 'pcs', 'maf'],
    spec: (f) => ({ outcome: f.outcome, covariates: f.covariates, n_pcs: f.n_pcs,
                    min_maf: f.maf, genetic_model: 'additive',
                    gene: f.targetGene || null,
                    exposure_variants: f.targetVariants.length ? f.targetVariants : null,
                    // Fitting one model per variant across a genome-wide set is
                    // a GWAS run the slow way; GWAS is the analysis for that.
                    max_variants: 5000 }),
  },
  gwas: {
    title: 'Genome-wide scan',
    fields: ['outcome', 'covariates', 'pcs', 'maf'],
    spec: (f) => ({ outcome: f.outcome, covariates: f.covariates, n_pcs: f.n_pcs,
                    min_maf: f.maf, prune_related: true }),
  },
  burden: {
    title: 'Collapse rare variants into genes',
    fields: ['outcome', 'covariates', 'pcs', 'test', 'preset'],
    spec: (f) => ({ outcome: f.outcome, covariates: f.covariates, n_pcs: f.n_pcs,
                    test: f.test, preset: f.preset, min_carriers: 2 }),
  },
  survival: {
    title: 'Time to event',
    fields: ['time', 'event', 'group', 'covariates'],
    spec: (f) => ({ time: f.time, event: f.event, group_by: f.group || null,
                    covariates: f.covariates }),
  },
  prs: {
    title: 'Polygenic risk score',
    fields: ['weights', 'outcome_optional'],
    spec: (f) => ({ weights: f.weights, outcome: f.outcome || null,
                    score_name: f.score_name || 'custom score',
                    ancestry_column: 'ancestry' }),
  },

  /* The R1 descriptive analyses. These take little configuration — that is the
     point of them — but they still need a form, because "run it with whatever
     defaults" is how a carrier rate ends up quoted against the wrong
     denominator. */
  carrier_frequency: {
    title: 'Carrier and allele frequency',
    fields: ['carrier_model', 'af_window'],
    spec: (f) => ({ carrier_model: f.carrier_model, min_af: f.af_min,
                    max_af: f.af_max, top_n: 200 }),
  },
  zygosity: {
    title: 'Zygosity and inheritance pattern',
    fields: ['zyg_note'],
    spec: () => ({ top_n: 200 }),
  },
  population_frequency: {
    title: 'Allele frequency between groups',
    fields: ['group_any'],
    spec: (f) => ({ group: f.groupAny || null, top_n: 200 }),
  },
  diagnostic_yield: {
    title: 'Subjects carrying a qualifying variant',
    fields: ['group_any', 'af_ceiling'],
    spec: (f) => ({ group: f.groupAny || null, max_af: f.afCeiling }),
  },
  segregation: {
    title: 'Mendelian consistency across trios',
    fields: ['seg_note'],
    spec: () => ({}),
  },
};

async function configureModal(state, analysis) {
  const form = FORMS[analysis];
  if (!form) { toast('No configuration form for ' + analysis); return; }

  const p = state.profile || {};
  const phenos = Object.entries(p.phenotypes || {});
  const pick = (kinds) => phenos.filter(([, v]) => kinds.includes(v.kind));
  const opts = (list) => list.map(([k, v]) =>
    `<option value="${esc(k)}">${esc(v.label || k)}</option>`).join('');

  const outcomes = pick(['binary', 'quantitative']);
  const binaries = pick(['binary']);
  const times = pick(['time_to_event', 'quantitative']);
  const quant = pick(['quantitative']).map(([k]) => k);
  const groups = pick(['binary', 'categorical']);

  const est = await api.research.estimate(analysis, state.dataset.n_samples,
                                          state.dataset.n_variants);

  const F = {
    outcome: `<div class="field"><label>Outcome</label>
      <select id="cfOut">${opts(analysis === 'burden' ? binaries : outcomes)}</select></div>`,
    outcome_optional: `<div class="field"><label>Outcome to validate against (optional)</label>
      <select id="cfOut"><option value="">none — just compute scores</option>
        ${opts(outcomes)}</select></div>`,
    time: `<div class="field"><label>Follow-up time</label>
      <select id="cfTime">${opts(times)}</select></div>`,
    event: `<div class="field"><label>Event occurred (1 = event, 0 = censored)</label>
      <select id="cfEvent">${opts(binaries)}</select></div>`,
    group: `<div class="field"><label>Compare groups by (optional)</label>
      <select id="cfGroup"><option value="">no grouping</option>${opts(groups)}</select></div>`,
    covariates: `<div class="field"><label>Adjust for</label>
      <div class="pills" id="cfCovs">${quant.map((c) =>
        `<button class="pill" data-cov="${esc(c)}">${esc(c)}</button>`).join('')}</div></div>`,
    pcs: `<div class="field"><label>Ancestry adjustment</label>
      <select id="cfPcs"><option value="10">first 10 principal components</option>
        <option value="5">first 5</option><option value="0">none</option></select></div>`,
    maf: `<div class="field"><label>Minimum allele frequency</label>
      <select id="cfMaf"><option value="0.05">5% — common variants</option>
        <option value="0.01">1%</option><option value="0.001">0.1%</option></select></div>`,
    test: `<div class="field"><label>Test</label>
      <select id="cfTest"><option value="skat_o">SKAT-O — combines both below</option>
        <option value="burden">Burden — when effects point the same way</option>
        <option value="skat">SKAT — when effects point both ways</option></select></div>`,
    preset: `<div class="field"><label>Which variants qualify</label>
      <select id="cfPreset">${Object.entries(state.rmeta.variant_set_presets || {})
        .map(([k, v]) => `<option value="${esc(k)}">${esc(v.label)} — ${esc(v.description || '')}</option>`)
        .join('')}</select></div>`,
    weights: `<div class="field"><label>Score weights</label>
      <textarea id="cfWeights" rows="5"
        placeholder="One per line:  1:1000:A:G  0.12"></textarea>
      <p class="hint" style="margin:6px 0 0">Variant key and weight, whitespace separated.
        Paste from a PGS Catalog file.</p></div>`,

    target: `<div class="field"><label>Which variants</label>
      <input type="text" id="cfTarget" placeholder="a gene symbol, or variant keys">
      <p class="hint" style="margin:6px 0 0">A gene symbol (<span class="mono">GENE0042</span>)
        tests every common variant in it. Or paste variant keys
        (<span class="mono">1:10000:A:G</span>), comma or space separated.
        Leave empty to scan everything — which is a GWAS run one variant at a
        time, so prefer the <b>Genome-wide association</b> analysis for that.</p></div>`,

    carrier_model: `<div class="field"><label>Who counts as a carrier</label>
      <select id="cfModel">
        <option value="dominant">Anyone with an ALT allele — het or hom</option>
        <option value="het_only">Heterozygotes only — recessive carrier screening</option>
        <option value="recessive">Homozygotes only — affected under a recessive model</option>
      </select>
      <p class="hint" style="margin:6px 0 0">This is a choice, not a fact. For carrier
        screening the carrier is the heterozygote; for a dominant condition it is
        anyone carrying the allele.</p></div>`,

    af_window: `<div class="field"><label>Allele frequency range</label>
      <div class="row2">
        <select id="cfAfMin" style="flex:1">
          <option value="0">no lower bound</option>
          <option value="0.0001">0.01%</option><option value="0.001">0.1%</option>
        </select>
        <select id="cfAfMax" style="flex:1">
          <option value="1">no upper bound</option>
          <option value="0.05" selected>5% and below — rare</option>
          <option value="0.01">1% and below</option>
          <option value="0.001">0.1% and below — ultra-rare</option>
        </select>
      </div></div>`,

    af_ceiling: `<div class="field"><label>A variant qualifies below</label>
      <select id="cfAfCeil">
        <option value="0.05">5% allele frequency</option>
        <option value="0.01" selected>1% allele frequency</option>
        <option value="0.001">0.1% — ultra-rare only</option>
      </select>
      <p class="hint" style="margin:6px 0 0">No ACMG classification is applied, so this
        is a candidate rate and an upper bound, not a diagnostic rate.</p></div>`,

    group_any: `<div class="field"><label>Split by (optional)</label>
      <select id="cfGroupAny"><option value="">no grouping — whole cohort</option>
        ${opts(phenos)}</select>
      <p class="hint" style="margin:6px 0 0">Any phenotype column: ancestry, indication,
        case/control, site.</p></div>`,

    zyg_note: note(`Counts het, homozygous-ALT and call rate for every variant, and a
      het/hom ratio per sample. No configuration needed — but note that hemizygous
      calls on X and Y read as homozygous in the genotype encoding.`, 'info'),

    seg_note: note(`Compares every complete trio and reports variants where the child
      carries an ALT allele neither parent does. These are <b>not</b> de-novo calls:
      at ordinary error rates most are genotyping error, and a trio far above the
      others is usually a sample swap.`, 'info'),
  };

  modal(form.title, `
    ${form.fields.map((f) => F[f] || '').join('')}
    <div id="cfPower"></div>
    ${note(`<b>Estimated compute:</b> ${esc(est.estimated_label)}. It runs in the
      background — you can leave this page.`, 'info')}
    <div id="cfOutMsg"></div>`, {
    footer: `${form.fields.includes('maf')
        ? '<button class="btn" id="cfPow">Check power first</button>' : ''}
      <button class="btn primary" id="cfGo">Run analysis</button>`,
  });

  const chosen = new Set();
  document.querySelectorAll('#cfCovs .pill').forEach((b) => {
    b.onclick = () => {
      const c = b.dataset.cov;
      chosen.has(c) ? chosen.delete(c) : chosen.add(c);
      b.classList.toggle('on');
    };
  });

  const val = (id, dflt) => {
    const el = document.getElementById(id);
    return el ? el.value : dflt;
  };

  const powBtn = document.getElementById('cfPow');
  if (powBtn) powBtn.onclick = async () => {
    const outcome = val('cfOut', '');
    const meta = p.phenotypes[outcome] || {};
    const maf = parseFloat(val('cfMaf', '0.05'));
    const body = meta.kind === 'binary'
      ? { kind: 'binary', n_cases: meta.cases || 0, n_controls: meta.controls || 0,
          maf, odds_ratio: 1.5, alpha_scope: 'genome' }
      : { kind: 'quantitative', n: meta.n_present || 0, maf, beta_sd: 0.2,
          alpha_scope: 'genome' };
    document.getElementById('cfPower').innerHTML = powerNote(await api.research.power(body));
  };

  document.getElementById('cfGo').onclick = async () => {
    const out = document.getElementById('cfOutMsg');
    const fields = {
      outcome: val('cfOut', ''),
      time: val('cfTime', ''),
      event: val('cfEvent', ''),
      group: val('cfGroup', ''),
      covariates: [...chosen],
      n_pcs: parseInt(val('cfPcs', '0'), 10),
      maf: parseFloat(val('cfMaf', '0.01')),
      test: val('cfTest', 'skat_o'),
      preset: val('cfPreset', 'ultra_rare_lof'),
      weights: parseWeights(val('cfWeights', '')),
      carrier_model: val('cfModel', 'dominant'),
      af_min: parseFloat(val('cfAfMin', '0')),
      af_max: parseFloat(val('cfAfMax', '1')),
      afCeiling: parseFloat(val('cfAfCeil', '0.01')),
      groupAny: val('cfGroupAny', ''),
      ...parseTarget(val('cfTarget', '')),
    };
    if (analysis === 'prs' && !Object.keys(fields.weights).length) {
      out.innerHTML = note('Paste at least one variant and weight.', 'bad');
      return;
    }
    try {
      const job = await api.research.submit(state.datasetId,
                                            { analysis, spec: form.spec(fields) });
      closeModal();
      state.jobId = job.job_id;
      toast('Running in the background');
      await ctxRender();
      pollJob(state, job.job_id);
    } catch (e) {
      out.innerHTML = note(esc(e.message), 'bad');
    }
  };
}

/* "1:1000:A:G  0.12" per line — the shape a PGS Catalog file pastes as. */
function parseWeights(text) {
  const out = {};
  (text || '').split(/\n/).forEach((line) => {
    const m = line.trim().split(/[\s,\t]+/);
    if (m.length >= 2 && m[0]) {
      const w = parseFloat(m[m.length - 1]);
      if (Number.isFinite(w)) out[m[0]] = w;
    }
  });
  return out;
}

function pollJob(state, jobId) {
  const timer = setInterval(async () => {
    const job = await api.research.job(jobId);
    if (['complete', 'failed', 'cancelled'].includes(job.status)) {
      clearInterval(timer);
      toast(`Analysis ${job.status}`);
      if (state.jobId === jobId) {
        await ctxRender();
      }
    } else if (state.jobId === jobId && state.step === 'analyse') {
      await ctxRender();
    }
  }, 2500);
}

async function overrideModal(state, analysis) {
  const cap = (state.caps || []).find((c) => c.analysis === analysis);
  modal(`Override · ${cap.title}`, `
    ${note(`<b>This analysis does not meet its requirements.</b> An override lets it run,
      but every result and export will permanently carry
      <span class="mono">${esc(state.rmeta.underpowered_stamp)}</span>. The stamp cannot
      be removed.`, 'warn')}
    <p class="hint">Unmet:</p>
    <ul class="hint">${cap.unmet.map((u) =>
      `<li><b>${esc(u.label)}</b> — you have: ${esc(u.observed)}</li>`).join('')}</ul>
    <div class="field"><label>Why is this defensible? (recorded in the audit log)</label>
      <textarea id="ovT" rows="4"></textarea></div>
    <div id="ovOut"></div>`, {
    footer: '<button class="btn primary" id="ovGo">Record and unlock</button>',
  });
  document.getElementById('ovGo').onclick = async () => {
    try {
      await api.research.override(state.datasetId, {
        analysis, justification: document.getElementById('ovT').value.trim(),
      });
      closeModal();
      toast('Override recorded — results will be stamped');
      const d = await api.research.dataset(state.datasetId);
      state.caps = d.capabilities;
      await ctxRender();
    } catch (e) {
      document.getElementById('ovOut').innerHTML = note(esc(e.message), 'bad');
    }
  };
}
