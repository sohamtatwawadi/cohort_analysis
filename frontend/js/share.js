/* Step 4 — Share.
 *
 * Previously export, manifest and audit were three separate rail entries and a
 * top-bar button, so the thing that makes a result defensible was as far from
 * the result as everything else. Here it is the last step of the workflow: you
 * arrive after you have an answer, and the manifest is shown as part of the
 * export rather than hidden behind a governance menu. */

import { api } from './api.js';
import { card, closeModal, esc, fmt, head, kv, modal, note, table, tag, toast } from './kit.js';

let shellState = null;

const DATASETS = [
  ['g_carrier', 'Carrier rate by gene', 'One row per gene with its own denominator.'],
  ['g_vus', 'Uncertain-variant review queue', 'Ranked curation worklist.'],
  ['g_variants', 'Every observation', 'One row per finding, with ACMG codes.'],
  ['g_subjects', 'Subject list', 'One row per subject with their result.'],
  ['g_runs', 'Assay and QC list', 'Provenance per run.'],
  ['g_popfreq', 'Recurrent variants', 'Internal frequency against gnomAD.'],
  ['g_zygosity', 'Zygosity by gene', 'Het / hom / compound-het / hemizygous.'],
  ['g_genedisease', 'Gene–disease reference', 'The curated evidence table.'],
  ['c_denominator', 'Coverage by gene', 'Who was tested for what.'],
];

export async function render(state) {
  shellState = state;
  if (state.mode === 'research') return renderResearch(state);

  const mf = await api.manifest(state.criteria, state.name, 'g_carrier');
  return `<div class="wrap narrow">
    ${head('Take it with you',
      'Every export carries a manifest, so whoever receives it can reproduce the numbers '
      + 'without asking you.', 'Step 4 of 4')}

    ${card('What to export', `
      <div class="field"><label>Dataset</label>
        <select id="exMod">${DATASETS.map(([id, label, desc]) =>
          `<option value="${id}">${esc(label)} — ${esc(desc)}</option>`).join('')}</select></div>
      <div class="field"><label>Format</label>
        <select id="exFmt">
          <option value="xlsx">Excel — data plus the manifest on a second sheet</option>
          <option value="csv">CSV — manifest as header comments</option>
          <option value="json">JSON — manifest and rows together</option>
        </select></div>
      <div id="exOut"></div>
      <div class="actions"><button class="btn primary lg" id="exGo">Generate export</button></div>
    `, { foot: 'Counts below 5 are written as “&lt;5”. A single rare-condition family in a '
          + 'named ancestry is re-identifying.' })}

    ${card('What travels with it', kv([
      ['Cohort', esc(mf.cohort_name)],
      ['Subjects', fmt(mf.member_count) + ` in ${fmt(mf.family_count)} families`],
      ['Criteria fingerprint', `<span class="mono">${esc(mf.criteria_hash)}</span>`],
      ['Membership fingerprint', `<span class="mono">${esc(mf.member_hash)}</span>`],
      ['Knowledgebase version', esc(mf.kb_snapshot_id)],
      ['Reference builds', esc((mf.reference_builds || []).join(', '))],
      ['Pipelines', esc((mf.pipeline_versions || []).join(', '))],
      ['Small-cell suppression', `on, below ${mf.suppression_threshold}`],
    ]), {
      sub: 'the manifest',
      actions: '<button class="btn sm" data-act="manifest">View full manifest</button>',
      foot: 'The two fingerprints are what make a figure reproducible. Re-running the same '
          + 'criteria against the same knowledgebase must produce the same membership '
          + 'fingerprint — if it does not, the underlying data changed.',
    })}

    <div class="actions">
      <button class="btn" data-act="goto:analyse">← Back to analyses</button>
      <span class="spacer"></span>
      <button class="btn" data-act="save">Save this cohort</button>
      <button class="btn" data-act="audit">View activity log</button>
    </div>
  </div>`;
}

function renderResearch(state) {
  return `<div class="wrap narrow">
    ${head('Results and provenance',
      'Every analysis records the software versions, parameters, seed and dataset release '
      + 'that produced it.', 'Step 4 of 4')}
    ${note(`Research results stay inside the project. Individual-level export is governed by
      policy; aggregate exports pass small-cell suppression.`, 'info')}
    ${card('Reproducibility', `<p class="hint" style="margin-top:0">Open any completed
      analysis from the previous step to see its full reproducibility record — container
      image, numpy and scipy versions, reference build, annotation release, parameters and
      random seed.</p>`)}
    <div class="actions">
      <button class="btn" data-act="goto:analyse">← Back to analyses</button>
    </div>
  </div>`;
}

export async function handle(verb, arg, rest, state, ctx) {
  shellState = state;
  switch (verb) {
    case 'manifest': {
      const mf = await api.manifest(state.criteria, state.name, 'g_carrier');
      modal('Cohort manifest', `
        <p class="hint" style="margin-top:0">This is embedded in every export.</p>
        <div class="codeblock">${esc(JSON.stringify(mf, null, 2))}</div>`, {
        footer: '<button class="btn" id="cpMf">Copy</button>',
      });
      document.getElementById('cpMf').onclick = () => {
        navigator.clipboard.writeText(JSON.stringify(mf, null, 2));
        toast('Manifest copied');
      };
      return true;
    }
    case 'audit': {
      modal('Activity log', '<div class="loading">Loading…</div>', { wide: true });
      const d = await api.audit();
      modal('Activity log', table([
        { label: 'Time', cell: (r) => `<span class="mono hint">${
            esc(String(r.executed_at).replace('T', ' ').slice(0, 19))}</span>` },
        { label: 'What was viewed', cell: (r) => esc(r.module_title) },
        { label: 'Class', cell: (r) => tag(r.output_class,
            r.output_class === 'OPERATIONAL' ? 'blue' : 'amber') },
        { label: 'Cohort size', n: true, cell: (r) => fmt(r.cohort_size) },
        { label: 'Criteria', cell: (r) => `<span class="mono">${esc(r.criteria_hash)}</span>` },
      ], d.rows), { wide: true, footer: `<span class="hint">${esc(d.footer)}</span>` });
      return true;
    }
    case 'save': {
      modal('Save this cohort', `
        <p class="hint" style="margin-top:0">A saved cohort stores the filters, not a list of
        people — so re-running it later reflects the data as it is then, and tells you if
        membership has changed.</p>
        <div class="field"><label>Name</label>
          <input type="text" id="svName" value="${esc(state.name)}"></div>
        <div class="field"><label>What question does it answer?</label>
          <input type="text" id="svDesc" placeholder="optional"></div>`, {
        footer: '<button class="btn primary" id="svGo">Save</button>',
      });
      document.getElementById('svGo').onclick = async () => {
        await api.saveCohort({
          name: document.getElementById('svName').value.trim() || 'Untitled cohort',
          criteria: state.criteria,
          description: document.getElementById('svDesc').value.trim(),
        });
        state.presets = (await api.cohorts()).cohorts;
        closeModal();
        toast('Cohort saved — it will appear as a starting point in step 1');
      };
      return true;
    }
    default: return false;
  }
}

/* Export is wired by delegation because the controls live in the page, not a
   modal. `shellState` is captured from the last handle() call rather than
   re-imported — importing './flow.js' here resolves to an unversioned URL and
   would hand back a second, unrendered module instance. */
document.addEventListener('click', async (e) => {
  if (!e.target.closest('#exGo')) return;
  const state = shellState;
  if (!state) return;
  const out = document.getElementById('exOut');
  out.innerHTML = '<div class="loading">Writing…</div>';
  try {
    const r = await api.export({
      criteria: state.criteria, name: state.name,
      module: document.getElementById('exMod').value,
      format: document.getElementById('exFmt').value,
    });
    out.innerHTML = note(`<b>Written.</b> ${esc(r.filename)} — ${fmt(r.rows)} rows.
      ${r.watermarked ? 'Marked research-use-only.' : ''}
      <div class="mono hint" style="margin-top:5px">${esc(r.path)}</div>`, 'good');
  } catch (err) {
    out.innerHTML = note('Export failed: ' + esc(err.message), 'bad');
  }
});
