/* Shared rendering helpers.
 *
 * One rule is enforced structurally here: `rate()` is the only way to render a
 * percentage, and it always prints the counts beside it. There is deliberately
 * no helper that formats a bare percentage, because a bare percentage is the
 * thing that gets quoted. */

export const esc = (s) => String(s === null || s === undefined ? '' : s)
  .replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

export const fmt = (n) => (n === null || n === undefined || n === '')
  ? '—' : Number(n).toLocaleString('en-US');

export const f1 = (n) => (n === null || n === undefined) ? '—' : Number(n).toFixed(1);

export function rate(r, opts = {}) {
  if (!r) return '—';
  if (!r.d) return `<span class="rate"><span class="of">${r.n} / 0 · no denominator</span></span>`;
  const warn = opts.warn
    ? ` <span class="warnic" title="${esc(opts.warnReason || '')}">⚠</span>` : '';
  return `<span class="rate"><span class="pc">${r.pct}%</span>
    <span class="of">${r.n}/${r.d}</span></span>${warn}`;
}

export const pct = (n, d) => (d ? Math.round(1000 * n / d) / 10 : 0);

/* ------------------------------------------------------------------ blocks */
export const card = (title, body, opts = {}) => `
  <section class="card">
    ${title ? `<h3>${esc(title)}<span class="grow"></span>
      ${opts.sub ? `<span class="sub">${esc(opts.sub)}</span>` : ''}
      ${opts.actions || ''}</h3>` : ''}
    ${opts.flush ? body : `<div class="body">${body}</div>`}
    ${opts.foot ? `<div class="foot">${opts.foot}</div>` : ''}
  </section>`;

export const note = (text, level = 'info') => {
  const ic = { info: 'i', good: '✓', warn: '!', bad: '!' }[level] || 'i';
  return `<div class="note ${level}"><span class="ic">${ic}</span><div>${text}</div></div>`;
};

export const head = (title, desc, eyebrow) => `
  <div class="pagehead">
    ${eyebrow ? `<div class="eyebrow">${esc(eyebrow)}</div>` : ''}
    <h1>${esc(title)}</h1>
    ${desc ? `<p>${esc(desc)}</p>` : ''}
  </div>`;

export const stat = (cells) => `<div class="stat">${cells.map((c) => `
  <div><div class="k">${esc(c.k)}</div><div class="v">${c.v}</div>
    ${c.d ? `<div class="d">${c.d}</div>` : ''}</div>`).join('')}</div>`;

export const tag = (label, tone = 'grey') =>
  `<span class="tag t-${tone}">${esc(label)}</span>`;

const ACMG_TONE = {
  'Pathogenic': 'p', 'Likely pathogenic': 'lp', 'Uncertain significance': 'vus',
  'Likely benign': 'lb', 'Benign': 'b',
};
export const acmg = (c) => `<span class="tag t-${ACMG_TONE[c] || 'grey'}">${
  esc({ 'Pathogenic': 'P', 'Likely pathogenic': 'LP', 'Uncertain significance': 'VUS',
        'Likely benign': 'LB', 'Benign': 'B' }[c] || c || '—')}</span>`;

export const ZYG_TONE = {
  'Heterozygous': 'grey', 'Homozygous': 'red',
  'Compound heterozygous': 'amber', 'Hemizygous': 'amber',
};

export function table(cols, rows, opts = {}) {
  if (!rows || !rows.length) {
    return `<div class="empty"><b>${esc(opts.emptyTitle || 'Nothing to show')}</b>
      ${esc(opts.emptyBody || '')}</div>`;
  }
  const head = cols.map((c) => `<th class="${c.n ? 'n' : ''}">${esc(c.label)}</th>`).join('');
  const body = rows.map((r, i) =>
    `<tr class="${opts.click ? 'click' : ''}" ${opts.attrs ? opts.attrs(r, i) : ''}>` +
    cols.map((c) => `<td class="${c.n ? 'n' : ''}">${c.cell(r, i)}</td>`).join('') +
    '</tr>').join('');
  return `<div class="scroll"><table><thead><tr>${head}</tr></thead>
    <tbody>${body}</tbody></table></div>`;
}

export const empty = (title, body = '') =>
  `<div class="empty"><b>${esc(title)}</b>${esc(body)}</div>`;

/* ------------------------------------------------------------------ charts */
/* Two palettes, because stacked() serves two different jobs.
 *
 * ORDINAL (the default) — a one-hue teal ramp, light to dark. Correct when the
 * order carries meaning: Definitive → Strong → Moderate → Limited reads as a
 * ramp because it IS one, and a reader sees the ordering in the colour.
 *
 * CATEGORICAL — for identity, where swapping the order would change nothing:
 * ancestry, sex, case/control. The teal ramp fails here and it is measurable,
 * not a matter of taste: adjacent steps sit at ΔE 9.0 under normal vision,
 * below the 15 floor, so two segments of a composition bar are genuinely hard
 * to tell apart. These eight clear it (worst adjacent pair ΔE 19.6 normal,
 * 9.1 protanopia) and carry visible labels and counts as the second channel.
 */
const PALETTE = ['#0E5C63', '#12787F', '#4E9AA0', '#8FBEC2', '#C2A15A', '#B4651A',
                 '#6B3FA0', '#2B6CB0', '#7A7F87', '#2E7D4F'];

export const CATEGORICAL = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100',
                            '#e87ba4', '#008300', '#4a3aa7', '#e34948'];
const TONES = { plp: '#A32B22', vus: '#8A949E', neg: '#B4D9C4' };

export function stacked(parts, opts = {}) {
  // opts.colors lets the caller choose the palette its data's job calls for;
  // the ordinal ramp stays the default.
  const ramp = opts.colors || PALETTE;
  const total = parts.reduce((a, p) => a + (p.value || 0), 0) || 1;
  const seg = parts.map((p, i) => `<i style="width:${100 * (p.value || 0) / total}%;
    background:${TONES[p.tone] || ramp[i % ramp.length]}"
    title="${esc(p.label)}: ${p.value}"></i>`).join('');
  const legend = parts.map((p, i) => `<span><i style="background:${
    TONES[p.tone] || ramp[i % ramp.length]}"></i>${esc(p.label)}
    <span class="hint mono">${fmt(p.value)} · ${Math.round(100 * (p.value || 0) / total)}%</span>
  </span>`).join('');
  return `<div class="stack">${seg}</div>${opts.noLegend ? '' : `<div class="legend">${legend}</div>`}`;
}

export function bars(rows, opts = {}) {
  const max = opts.max || Math.max(1, ...rows.map((r) => r.value));
  return rows.map((r) => `
    <div class="hbar">
      <div class="l ${r.act ? 'link' : ''}" ${r.act ? `data-act="${esc(r.act)}"` : ''}
        title="${esc(r.title || r.label)}">${esc(r.label)}</div>
      <div class="bar2"><i style="width:${Math.max(2, 100 * r.value / max)}%"></i></div>
      <div class="v">${r.display !== undefined ? r.display : fmt(r.value)}</div>
    </div>`).join('');
}

export function heatmap(rows, cols, matrix, max, opts = {}) {
  const head = '<tr><th></th>' + cols.map((c) =>
    `<th title="${esc(c)}">${esc(opts.shortCol ? opts.shortCol(c) : c)}</th>`).join('') + '</tr>';
  const body = rows.map((r) => '<tr><th style="background:none;position:static">' +
    `<span class="gene">${esc(r)}</span></th>` +
    cols.map((c) => {
      const v = (matrix[r] || {})[c] || 0;
      const o = v ? (0.10 + 0.9 * (v / max)) : 0;
      return `<td class="n" title="${esc(r)} × ${esc(c)}: ${v}"
        style="background:${v ? `rgba(14,92,99,${o.toFixed(3)})` : '#F7F9FA'};
        color:${o > 0.55 ? '#fff' : 'var(--text)'};font-family:'IBM Plex Mono',monospace">
        ${v || '·'}</td>`;
    }).join('') + '</tr>').join('');
  return `<div style="overflow-x:auto"><table>${head}${body}</table></div>`;
}

/* ------------------------------------------------------------------- chrome */
export function toast(msg) {
  const d = document.createElement('div');
  d.className = 'toast';
  d.textContent = msg;
  document.getElementById('toasts').appendChild(d);
  setTimeout(() => d.remove(), 3000);
}

export function modal(title, body, opts = {}) {
  const root = document.getElementById('modalRoot');
  root.innerHTML = `
    <div class="modal-wrap" data-close="1">
      <div class="modal ${opts.wide ? 'wide' : ''}" role="dialog" aria-modal="true"
        aria-label="${esc(title)}">
        <header><h3>${esc(title)}</h3><div class="grow"></div>${opts.badge || ''}
          <button class="btn sm" data-close="1">Close</button></header>
        <div class="body">${body}</div>
        ${opts.footer ? `<footer>${opts.footer}</footer>` : ''}
      </div>
    </div>`;
  root.querySelectorAll('[data-close]').forEach((el) => {
    el.addEventListener('click', (e) => { if (e.target === el) closeModal(); });
  });
}
export const closeModal = () => { document.getElementById('modalRoot').innerHTML = ''; };

export const kv = (pairs) => `<dl class="kv">${pairs
  .filter(([, v]) => v !== undefined && v !== null && v !== '')
  .map(([k, v]) => `<dt>${esc(k)}</dt><dd>${v}</dd>`).join('')}</dl>`;

export const crumb = (parts) => `<div class="crumb">${parts.map((p, i) =>
  (i ? '<span>›</span>' : '') + (p.act
    ? `<a data-act="${esc(p.act)}">${esc(p.label)}</a>`
    : `<span>${esc(p.label)}</span>`)).join('')}</div>`;
