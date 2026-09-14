/* The shell: one workflow spine, four steps.
 *
 *     ① Build cohort → ② Review → ③ Analyse → ④ Share
 *
 * Both data sources use the same four steps. That is the whole point of the
 * redesign — previously Lab Mode and Research Mode each had their own
 * navigation, so learning one taught you nothing about the other.
 *
 * Every screen answers: where am I, what do I do here, what is next. */

import { api } from './api.js';
import { closeModal, esc, fmt, modal, note, toast } from './kit.js';
import * as build from './build.js';
import * as review from './review.js';
import * as analyse from './analyse.js';
import * as share from './share.js';

export const state = {
  mode: 'lab',                 // 'lab' | 'research'
  step: 'build',
  meta: null,
  rmeta: null,

  // lab
  criteria: {},
  name: 'Everyone in the database',
  cohort: null,
  presets: [],

  // research
  projectId: null,
  datasetId: null,
  dataset: null,
  profile: null,
  caps: [],

  // step 3 drill state
  question: null,
  gene: null,
  finding: null,
  subject: null,
  jobId: null,

  reviewed: false,
};

const STEPS = [
  { id: 'build', n: 1, label: 'Build cohort', sub: 'Who is included' },
  { id: 'review', n: 2, label: 'Review', sub: 'Can you trust it' },
  { id: 'analyse', n: 3, label: 'Analyse', sub: 'Ask your question' },
  { id: 'share', n: 4, label: 'Share', sub: 'Export with provenance' },
];

const RENDERERS = { build, review, analyse, share };

export const blankCriteria = () => JSON.parse(JSON.stringify(state.meta.blank_criteria));

/* Everything a step module is allowed to call on the shell. */
export const ctx = { render, resolveCohort, blankCriteria, go, state };

/* --------------------------------------------------------------------- boot */
async function boot() {
  bindChrome();
  state.meta = await api.meta();
  state.criteria = blankCriteria();
  state.presets = (await api.cohorts()).cohorts;
  await go('build');
}

function bindChrome() {
  document.querySelectorAll('#modes button').forEach((b) => {
    b.onclick = () => setMode(b.dataset.mode);
  });
  document.getElementById('askBtn').onclick = openAsk;
  document.getElementById('helpBtn').onclick = showHelp;

  window.addEventListener('keydown', (e) => {
    if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === 'k') {
      e.preventDefault(); openAsk();
    }
    if (e.key === 'Escape') closeModal();
  });
}

/* ------------------------------------------------------------------- steps */
export function stepAvailable(id) {
  if (state.mode === 'research') {
    // Nothing downstream of "pick a dataset" makes sense without one.
    if (id !== 'build' && !state.datasetId) return false;
    return true;
  }
  if (id === 'build') return true;
  return !!state.cohort;
}

function renderSteps() {
  const order = STEPS.map((s) => s.id);
  const here = order.indexOf(state.step);
  document.getElementById('steps').innerHTML = STEPS.map((s, i) => {
    const on = s.id === state.step;
    const done = i < here;
    const can = stepAvailable(s.id);
    return `<button class="step ${on ? 'on' : ''} ${done ? 'done' : ''}"
      data-step="${s.id}" ${can ? '' : 'disabled'}
      aria-current="${on ? 'step' : 'false'}">
      <span class="dot">${done ? '✓' : s.n}</span>
      <span class="txt"><span class="lbl">${esc(labelFor(s))}</span>
        <span class="sub">${esc(subFor(s))}</span></span>
    </button>`;
  }).join('');
  document.querySelectorAll('.step').forEach((b) => {
    b.onclick = () => go(b.dataset.step);
  });
}

/* The steps keep their meaning across modes but not their wording — "Build
   cohort" is wrong when the job is "choose a dataset". */
function labelFor(s) {
  if (state.mode !== 'research') return s.label;
  return { build: 'Choose data', review: 'Review', analyse: 'Analyse',
           share: 'Share' }[s.id];
}
function subFor(s) {
  if (state.mode !== 'research') return s.sub;
  return { build: 'Upload or pick a dataset', review: 'What it supports',
           analyse: 'Run an analysis', share: 'Results and provenance' }[s.id];
}

export async function go(step) {
  if (!stepAvailable(step)) return;
  state.step = step;
  if (step !== 'analyse') { state.question = null; state.gene = null; }
  await render();
}

export async function render() {
  const stage = document.getElementById('stage');
  stage.innerHTML = '<div class="loading">Working…</div>';
  renderSteps();
  try {
    const html = await RENDERERS[state.step].render(state);
    stage.innerHTML = html;
    renderCohortBar();
  } catch (e) {
    // Clear the context bar too. Leaving the previous render's cohort line
    // above an error message shows lab counts while the user is looking at
    // Research Mode, which reads as "this is the data that failed".
    document.getElementById('cohortbar').innerHTML = '';
    stage.innerHTML = `<div class="wrap">${note(
      `<b>Something went wrong.</b> ${esc(e.message)}`, 'bad')}</div>`;
  }
  // Again, now that the step has run. On first load the cohort does not exist
  // until build.render() resolves it, so the first renderSteps() above sees
  // state.cohort === null and disables every later step — permanently, because
  // nothing re-rendered them. The user could never leave step 1.
  renderSteps();
  stage.scrollTop = 0;
  stage.focus({ preventScroll: true });
}

/* --------------------------------------------------------------- cohort bar */
function renderCohortBar() {
  const bar = document.getElementById('cohortbar');

  if (state.mode === 'research') {
    bar.innerHTML = state.dataset ? `
      <span class="name">${esc(state.dataset.name)}</span>
      <span class="counts">${fmt(state.dataset.n_samples)} samples ·
        ${fmt(state.dataset.n_variants)} variants ·
        ${esc(state.dataset.genome_build || 'build not declared')}</span>
      <span class="chip mut">${esc(state.dataset.source_format)}</span>` : '';
    return;
  }

  if (!state.cohort || state.step === 'build') { bar.innerHTML = ''; return; }
  const c = state.cohort;
  bar.innerHTML = `
    <span class="name">${esc(state.name)}</span>
    <span class="counts">${fmt(c.counts.subjects)} subjects ·
      ${fmt(c.counts.families)} families · ${fmt(c.counts.observations)} findings</span>
    ${c.chips.map((x) => `<span class="chip ${x.sticky ? 'flag' : ''}">${esc(x.label)}
      <button data-act="rmchip:${esc(x.field)}:${esc(String(x.value))}"
        aria-label="Remove ${esc(x.label)}">×</button></span>`).join('')}
    <span class="grow" style="flex:1"></span>
    <button class="btn sm" data-act="goto:build">Edit cohort</button>`;
}

/* --------------------------------------------------------------------- mode */
async function setMode(mode) {
  if (mode === state.mode) return;
  state.mode = mode;
  document.querySelectorAll('#modes button').forEach((b) => {
    const on = b.dataset.mode === mode;
    b.classList.toggle('on', on);
    b.setAttribute('aria-selected', String(on));
  });
  state.step = 'build';
  state.question = null;
  if (mode === 'research') {
    // Re-fetch every time. A cached project id can outlive the project it
    // names — rebuilding the store, or switching tenant, leaves the old id
    // pointing at nothing and every request 404s with "unknown project".
    if (!state.rmeta) state.rmeta = await api.research.meta();
    const { projects } = await api.research.projects();
    const known = projects.some((p) => p.project_id === state.projectId);
    if (!known) {
      state.projectId = projects[0] && projects[0].project_id;
      state.datasetId = null;
      state.dataset = null;
      state.profile = null;
      state.caps = [];
    }
  }
  await render();
}

/* ------------------------------------------------------------- cohort calls */
export async function resolveCohort() {
  state.cohort = await api.resolve(state.criteria, state.name);
  return state.cohort;
}

/* -------------------------------------------------------------- ask (⌘K) -- */
function openAsk() {
  if (state.mode === 'research') {
    toast('Ask a question works on lab cohorts. Switch to "Our lab data".');
    return;
  }
  modal('Ask a question', `
    <p class="hint" style="margin-top:0">Type a question in plain English. It is turned
    into cohort filters — never into a number. You see the filters and approve them
    before anything runs.</p>
    <input type="text" id="askIn" placeholder="e.g. carrier rate for HBB in South Asian subjects"
      autocomplete="off">
    <div id="askOut"></div>
    <p class="hint" style="margin:16px 0 0;font-weight:650">Or start from one of these</p>
    <ul class="asklist">${state.meta.suggested_questions.map((q) =>
      `<li><button data-q="${esc(q)}">${esc(q)}</button></li>`).join('')}</ul>`, {
    footer: `<button class="btn primary" id="askGo">Compile</button>
      <span class="hint">Nothing runs until you approve it.</span>`,
  });

  const input = document.getElementById('askIn');
  input.focus();
  const compile = async () => {
    const q = input.value.trim();
    if (!q) return;
    const out = document.getElementById('askOut');
    out.innerHTML = '<div class="loading">Compiling…</div>';
    const r = await api.compile(q);
    out.innerHTML = `
      ${r.validation === 'PASSED'
        ? note(`<b>Understood.</b> ${r.fields_mapped} filter(s) recognised.`, 'good')
        : note(`<b>I could not map that.</b> ${
            r.unmapped.length ? 'Unrecognised: <span class="mono">'
              + esc(r.unmapped.join(', ')) + '</span>. ' : ''}${esc(r.errors.join(' '))}
            Rephrase, or pick a starting point below.`, 'bad')}
      <div class="codeblock">${esc(r.compiled_json)}</div>
      ${r.validation === 'PASSED'
        ? `<button class="btn primary" id="askRun" style="margin-top:12px">
             Apply these filters →</button>` : ''}`;
    const run = document.getElementById('askRun');
    if (run) run.onclick = async () => {
      state.criteria = { ...blankCriteria(), ...r.compiled.criteria };
      state.name = r.question;
      closeModal();
      await resolveCohort();
      state.question = QUESTION_FOR_MODULE[r.target_module] || null;
      await go(state.question ? 'analyse' : 'review');
      toast('Filters applied — every figure below was computed by the database');
    };
  };
  document.getElementById('askGo').onclick = compile;
  input.onkeydown = (e) => { if (e.key === 'Enter') compile(); };
  document.querySelectorAll('.asklist button').forEach((b) => {
    b.onclick = () => { input.value = b.dataset.q; compile(); };
  });
}

const QUESTION_FOR_MODULE = {
  g_carrier: 'carrier', g_dashboard: null, g_zygosity: 'zygosity',
  g_genedisease: 'validity', g_vus: 'vus', g_popfreq: 'popfreq',
  g_sf: 'secondary', g_phenotype: 'phenotype',
};

function showHelp() {
  modal('How this works', `
    <p style="margin-top:0">Four steps, in order. You can go back at any time.</p>
    <ol style="line-height:1.9;padding-left:20px">
      <li><b>Build cohort</b> — choose who is included. Consent and family
        de-duplication are applied automatically; you will see what they removed.</li>
      <li><b>Review</b> — three checks tell you whether the cohort can be quoted:
        independence, coverage, and provenance.</li>
      <li><b>Analyse</b> — pick a question. Click any number to see the subjects,
        genes or variants behind it.</li>
      <li><b>Share</b> — export with a manifest so someone else can reproduce it.</li>
    </ol>
    ${note(`Every percentage in this tool shows its counts, and every gene uses its own
      denominator — the subjects whose test actually covers that gene. That is why one
      cohort has many denominators.`, 'info')}`);
}

/* ------------------------------------------------------- action dispatcher */
document.addEventListener('click', async (e) => {
  const el = e.target.closest('[data-act]');
  if (!el) return;
  const [verb, ...rest] = el.dataset.act.split(':');
  const arg = rest.join(':');
  const handled = await dispatch(verb, arg, rest, el);
  if (handled === false) return;
});

async function dispatch(verb, arg, rest, el) {
  switch (verb) {
    case 'goto': return go(arg);
    case 'rmchip': {
      const field = rest[0];
      const value = rest.slice(1).join(':');
      const C = state.criteria;
      if (Array.isArray(C[field])) C[field] = C[field].filter((x) => String(x) !== value);
      else if (field === 'gene_set') C.gene_set = null;
      else if (field === 'months' || field === 'gnomad_max') C[field] = 0;
      else C[field] = false;
      await resolveCohort();
      return render();
    }
    default:
      // Steps own their own verbs; ask each in turn.
      //
      // `ctx` carries the shell functions a step needs. The steps used to
      // reach back with `await import('./flow.js')`, which looked harmless but
      // resolved WITHOUT the ?v= cache-busting query the page loaded flow.js
      // with — a different URL is a different module instance, so the step was
      // mutating a second, unrendered `state`. Passing the functions in makes
      // that impossible and removes the circular import.
      for (const mod of [build, review, analyse, share]) {
        if (mod.handle && await mod.handle(verb, arg, rest, state, ctx, el)) return true;
      }
      return false;
  }
}

boot();
