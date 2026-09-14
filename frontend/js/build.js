/* Step 1 — Build cohort (lab) / Choose data (research).
 *
 * The old builder opened on four accordions containing ~40 controls and no
 * indication of where to start. This opens on a short list of real starting
 * points, because almost nobody builds a cohort from nothing — they start from
 * "everyone" or from a known group and narrow it.
 *
 * Filters are grouped by the question they answer (Who / How tested / What
 * findings) rather than by our internal schema, and anything most users never
 * touch is behind a disclosure. */

import { api } from './api.js';
import {
  card, closeModal, empty, esc, f1, fmt, head, modal, note, stat, table, toast,
} from './kit.js';

const GROUPS = [
  {
    id: 'who', title: 'Who is included',
    fields: [
      ['indication', 'Reason for testing'],
      ['affected_status', 'Affected status'],
      ['family_history', 'Family history'],
      ['sex', 'Sex'],
      ['age_bucket', 'Age'],
      ['ancestry', 'Ancestry'],
    ],
  },
  {
    id: 'test', title: 'How they were tested',
    fields: [
      ['test_code', 'Test code'],
      ['sample_type', 'Sample type'],
      ['reference_build', 'Reference build'],
      ['pipeline_version', 'Pipeline version'],
    ],
  },
  {
    id: 'findings', title: 'What findings they carry',
    fields: [
      ['classification', 'Classification'],
      ['zygosity', 'Zygosity'],
      ['inheritance', 'Inheritance'],
      ['validity', 'Gene–disease validity'],
    ],
  },
];

const ADVANCED = [
  ['gene', 'Specific genes'],
  ['consequence', 'Consequence'],
  ['var_class', 'Variant class'],
  ['clinvar_sig', 'ClinVar significance'],
  ['relation', 'Relation to proband'],
  ['referral_source', 'Referral source'],
  ['penetrance', 'Penetrance'],
  ['consent_class', 'Consent class'],
  ['assay_version', 'Assay version'],
];

export async function render(state) {
  if (state.mode === 'research') return renderResearch(state);
  if (!state.cohort) await resolve(state);
  return renderLab(state);
}

async function resolve(state) {
  state.cohort = await api.resolve(state.criteria, state.name);
  return state.cohort;
}

/* ===================================================================== lab == */
function renderLab(state) {
  const c = state.criteria;
  const co = state.cohort;
  const v = state.meta.vocabulary;
  const fam = co.family_stats || {};

  const activeIn = (fields) => fields.reduce((n, [f]) => {
    const val = c[f];
    return n + (Array.isArray(val) ? val.length : (val ? 1 : 0));
  }, 0);

  const pillRow = (field, values) => `<div class="pills">${(values || []).map((val) =>
    `<button class="pill ${(c[field] || []).includes(val) ? 'on' : ''}"
      data-act="tog:${esc(field)}:${esc(val)}"
      aria-pressed="${(c[field] || []).includes(val)}">${esc(val)}</button>`).join('')}</div>`;

  const fieldBlock = ([field, label]) => `
    <div class="field"><label>${esc(label)}</label>${pillRow(field, v[field])}</div>`;

  const group = (g) => {
    const n = activeIn(g.fields);
    return `<details class="fgroup" ${n ? 'open' : ''}>
      <summary>${esc(g.title)}<span class="grow"></span>
        ${n ? `<span class="count">${n} selected</span>` : '<span class="hint">any</span>'}
      </summary>
      <div class="fbody">${g.fields.map(fieldBlock).join('')}</div>
    </details>`;
  };

  const geneSets = Object.keys(state.meta.gene_sets);

  return `<div class="wrap">
    ${head('Build your cohort',
      'Start from a common group, then narrow it. The panel on the right shows what each '
      + 'choice costs you.', 'Step 1 of 4')}

    <div class="cols c21">
      <div>
        ${card('Start from', `
          <div class="pills">
            ${state.presets.slice(0, 8).map((p) => `
              <button class="pill ${state.name === p.name ? 'on' : ''}"
                data-act="preset:${esc(p.cohort_id)}" title="${esc(p.description || '')}">
                ${esc(p.name)}</button>`).join('')}
          </div>
          <p class="hint" style="margin:12px 0 0">Or leave it as everyone and add filters below.</p>
        `)}

        ${card('Narrow it down', `
          <div class="field">
            <label>Curated gene set</label>
            <div class="pills">${geneSets.map((s) => `
              <button class="pill ${c.gene_set === s ? 'on' : ''}"
                data-act="geneset:${esc(s)}">${esc(s)}</button>`).join('')}</div>
          </div>
          ${GROUPS.map(group).join('')}
          <details class="fgroup">
            <summary>Advanced<span class="grow"></span>
              <span class="count">${activeIn(ADVANCED) || ''}</span></summary>
            <div class="fbody">
              ${ADVANCED.map(fieldBlock).join('')}
              <div class="field"><label>Collected within</label>
                <div class="pills">${[0, 6, 12, 24].map((m) => `
                  <button class="pill ${c.months === m ? 'on' : ''}"
                    data-act="months:${m}">${m ? m + ' months' : 'Any time'}</button>`).join('')}
                </div></div>
              <label class="switch"><input type="checkbox" ${c.reportable_only ? 'checked' : ''}
                data-act="flag:reportable_only">
                <span><span class="st">Reportable findings only</span>
                <span class="sd">Exclude findings outside what the test code reports.</span></span>
              </label>
              <label class="switch"><input type="checkbox"
                ${c.established_validity_only ? 'checked' : ''}
                data-act="flag:established_validity_only">
                <span><span class="st">Established gene–disease validity only</span>
                <span class="sd">Definitive and Strong genes only. Excludes weaker
                  evidence that would inflate a yield figure.</span></span></label>
            </div>
          </details>
        `)}
      </div>

      <div>
        ${card('Your cohort', `
          <div style="text-align:center;padding:6px 0 14px">
            <div style="font-family:'IBM Plex Mono',monospace;font-size:42px;
              font-weight:500;letter-spacing:-.03em;line-height:1">${fmt(co.counts.subjects)}</div>
            <div class="hint" style="font-size:13px">subjects in ${fmt(co.counts.families)} families</div>
          </div>
          ${co.counts.subjects === 0
            ? note('<b>No one matches.</b> The list below shows which choice emptied it — '
                 + 'remove that filter to get subjects back.', 'bad')
            : ''}
          <label class="switch"><input type="checkbox" ${c.probands_only ? 'checked' : ''}
            data-act="flag:probands_only">
            <span><span class="st">One subject per family</span>
            <span class="sd">${c.probands_only
              ? `On. ${fmt(fam.suppressed)} relatives excluded so carrier rates stay independent.`
              : `<b>Off. ${fmt(fam.extra_subjects)} relatives included across `
                + `${fmt(fam.multi_member_families)} families — rates are inflated and `
                + `must not be quoted.</b>`}</span></span></label>
        `, { foot: 'Consent is always applied. Withdrawn subjects cannot appear in any cohort.' })}

        ${card('What each choice removed', funnel(co.funnel), {
          foot: 'Every excluded subject is attributed to exactly one choice.',
        })}

        <div class="actions">
          <button class="btn primary lg" data-act="goto:review"
            ${co.counts.subjects === 0 ? 'disabled' : ''}>Review this cohort →</button>
          <button class="btn" data-act="clear">Clear all</button>
        </div>
      </div>
    </div>
  </div>`;
}

function funnel(steps) {
  const max = Math.max(1, ...steps.map((s) => s.count));
  return `<div class="funnel">${steps.map((s) => `
    <div class="fr ${s.final ? 'final' : ''}">
      <div>${esc(s.step)}${s.mandatory ? ' <span class="lk">always applied</span>' : ''}</div>
      <div class="bar"><i style="width:${Math.max(2, 100 * s.count / max)}%"></i></div>
      <div class="n">${fmt(s.count)}</div>
      <div class="dp">${s.dropped ? '−' + fmt(s.dropped) : ''}</div>
    </div>`).join('')}</div>`;
}

/* ================================================================ research == */
async function renderResearch(state) {
  const d = await api.research.datasets(state.projectId);
  const rows = d.datasets || [];
  return `<div class="wrap">
    ${head('Choose a dataset',
      'Upload your own genotype data, or pick a dataset you have already registered. '
      + 'What you can analyse depends on what the data supports.', 'Step 1 of 4')}
    ${note(esc(state.rmeta.ownership_notice), 'info')}
    ${rows.length ? card('Your datasets', table([
      { label: 'Dataset', cell: (r) => `<a data-act="pick:${esc(r.dataset_id)}">${esc(r.name)}</a>` },
      { label: 'Format', cell: (r) => `<span class="mono">${esc(r.source_format)}</span>` },
      { label: 'Build', cell: (r) => `<span class="mono">${esc(r.genome_build || '—')}</span>` },
      { label: 'Samples', n: true, cell: (r) => fmt(r.n_samples) },
      { label: 'Variants', n: true, cell: (r) => fmt(r.n_variants) },
      { label: '', cell: (r) => `<button class="btn sm" data-act="rmds:${esc(r.dataset_id)}">delete</button>` },
    ], rows, { click: true, attrs: (r) => `data-act="pick:${esc(r.dataset_id)}"` }), { flush: true })
      : card('Your datasets', empty('No datasets yet',
          'Register a VCF or PLINK fileset to begin.'))}
    <div class="actions">
      <button class="btn primary lg" data-act="upload">Register a dataset</button>
    </div>
  </div>`;
}

/* ---------------------------------------------------------------- handlers */
export async function handle(verb, arg, rest, state, ctx) {
  const { render: rerender, resolveCohort, blankCriteria, go } = ctx;
  ctxRender = rerender;
  const C = state.criteria;

  switch (verb) {
    case 'tog': {
      const field = rest[0];
      const value = rest.slice(1).join(':');
      const list = C[field] || (C[field] = []);
      const i = list.indexOf(value);
      if (i < 0) list.push(value); else list.splice(i, 1);
      break;
    }
    case 'flag': C[arg] = !C[arg]; break;
    case 'months': C.months = Number(arg); break;
    case 'geneset': C.gene_set = C.gene_set === arg ? null : arg; break;
    case 'clear':
      state.criteria = blankCriteria();
      state.name = 'Everyone in the database';
      break;
    case 'preset': {
      const p = state.presets.find((x) => x.cohort_id === arg);
      if (p) { state.criteria = { ...blankCriteria(), ...p.criteria }; state.name = p.name; }
      break;
    }
    case 'pick':
      state.datasetId = arg;
      { const d = await api.research.dataset(arg);
        state.dataset = d.dataset; state.profile = d.profile; state.caps = d.capabilities; }
      return go('review');
    case 'rmds':
      if (!window.confirm('Delete this dataset? Deletion removes the profile, every '
        + 'stored result and the genotype data itself. This cannot be undone.')) return true;
      await api.research.del(arg, state.projectId);
      toast('Dataset and all derived data deleted');
      if (state.datasetId === arg) { state.datasetId = null; state.dataset = null; }
      break;
    case 'upload': uploadModal(state); return true;
    default: return false;
  }

  if (state.mode === 'lab') await resolveCohort();
  await rerender();
  return true;
}

let ctxRender = async () => {};

function uploadModal(state) {
  modal('Register a dataset', `
    ${note(esc(state.rmeta.ownership_notice), 'info')}
    <div class="field"><label>Name</label>
      <input type="text" id="upName" placeholder="Cardiomyopathy case-control 2026"></div>

    <div class="field"><label>Genotype file</label>
      <input type="file" id="upGeno" multiple
        accept=".vcf,.gz,.bed,.bim,.fam,.tsv,.txt,.csv">
      <p class="hint" style="margin:6px 0 0">VCF or VCF.gz — one file. PLINK — select
        all three of <span class="mono">.bed .bim .fam</span> together. VariMAT — one
        file per sample.</p></div>

    <div class="field"><label>Phenotype file (optional)</label>
      <input type="file" id="upPheno" accept=".csv,.tsv,.txt">
      <p class="hint" style="margin:6px 0 0">CSV or TSV with a
        <span class="mono">sample_id</span> column. Without it only descriptive
        analyses unlock.</p></div>

    <div class="field"><label>Genome build</label>
      <select id="upBuild"><option value="">Select the build…</option>
        <option>GRCh38</option><option>GRCh37</option><option>T2T-CHM13</option></select>
      <p class="hint" style="margin:6px 0 0">We never guess this. GRCh37 and GRCh38
        positions overlap, so a wrong guess silently corrupts every annotation.</p></div>

    <label class="switch"><input type="checkbox" id="upConsent">
      <span><span class="st">I hold the necessary consent and ethics approval</span>
      <span class="sd">${esc(state.rmeta.consent_statement)}</span></span></label>

    <details class="fgroup"><summary>The file is already on the server<span class="grow"></span>
      <span class="hint">for very large datasets</span></summary>
      <div class="fbody">
        <div class="field"><label>Path on the server</label>
          <input type="text" id="upPath" placeholder="data/uploads/cohort.vcf.gz"></div>
        <div class="field"><label>Phenotype path</label>
          <input type="text" id="upPhenoPath" placeholder="data/uploads/phenotypes.csv"></div>
        <p class="hint">Browsers stream a multi-gigabyte upload poorly. For a biobank
          extract, put it on the server and give the path instead.</p>
      </div>
    </details>

    <div id="upOut"></div>`, {
    footer: '<button class="btn primary" id="upGo">Validate and register</button>',
  });

  document.getElementById('upGo').onclick = async () => {
    const out = document.getElementById('upOut');
    const name = document.getElementById('upName').value.trim() || 'Untitled dataset';
    const build = document.getElementById('upBuild').value || null;
    const consent = document.getElementById('upConsent').checked;
    const files = document.getElementById('upGeno').files;
    const phenoFile = document.getElementById('upPheno').files[0];
    const serverPath = document.getElementById('upPath').value.trim();

    if (!consent) {
      out.innerHTML = note('Consent attestation is required before a dataset can be '
        + 'registered.', 'bad');
      return;
    }
    if (!files.length && !serverPath) {
      out.innerHTML = note('Choose a genotype file, or give a path on the server.', 'bad');
      return;
    }
    out.innerHTML = '<div class="loading">Uploading and validating…</div>';

    try {
      if (files.length) {
        const form = new FormData();
        form.append('name', name);
        if (build) form.append('genome_build', build);
        form.append('consent_attested', 'true');
        for (const f of files) form.append('genotypes', f);
        if (phenoFile) form.append('phenotypes', phenoFile);
        await api.research.upload(state.projectId, form);
      } else {
        await api.research.ingest(state.projectId, {
          name,
          genotype_path: serverPath,
          phenotype_path: document.getElementById('upPhenoPath').value.trim() || null,
          genome_build: build,
          consent_attested: true,
        });
      }
      closeModal();
      toast('Dataset registered');
      await ctxRender();
    } catch (err) {
      const rep = err.detail && err.detail.report;
      out.innerHTML = note(`<b>${esc(err.message)}</b>${rep
        ? '<ul style="margin:8px 0 0">' + rep.issues.filter((i) => i.severity === 'fail')
            .map((i) => `<li>${esc(i.message)}</li>`).join('') + '</ul>' : ''}`, 'bad');
    }
  };
}
