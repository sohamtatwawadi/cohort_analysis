/* Inline-SVG charts for the research results.
 *
 * No chart library: the app has no build step, and a CDN dependency in a tool
 * that runs inside a lab network is a liability. Everything here is SVG built
 * from a string, which also means a chart survives being copied into a report.
 *
 * Conventions this file holds to, so the charts read as one system:
 *
 *   - Marks carry colour; text never does. Axis labels, values and legends use
 *     the ink tokens. Identity comes from a swatch beside the text.
 *   - Two hues only, and they are measured, not chosen: #2a78d6 / #eb6834 clear
 *     ΔE 24.7 under protanopia and 33.6 under normal vision on white. The app's
 *     own teal and status colours FAIL as a categorical palette — they sit below
 *     the chroma floor and red↔amber collapse to ΔE 3.2 under deuteranopia — so
 *     they are used for status and emphasis only, never to tell series apart.
 *   - Most of these charts are single-series or emphasis (one hue + grey), which
 *     is the honest form when the story is "this one point matters", not "these
 *     categories differ".
 *   - Gridlines are hairline and recessive. The data is the only loud thing.
 *   - Every chart has a hover layer; a chart that cannot be interrogated is a
 *     picture. Each also sits beside the table it was drawn from, so nothing is
 *     gated behind colour.
 */

import { esc, fmt } from './kit.js';

/* Measured categorical slots — see validate_palette.js in the dataviz skill. */
export const C1 = '#2a78d6';
export const C2 = '#eb6834';
const INK = 'var(--text)';
const MUTED = 'var(--muted)';
const GRID = 'var(--rule)';
const SURFACE = 'var(--panel)';
const EMPH = 'var(--bad)';          // status: the thing that crossed a threshold
const DIM = '#B6BEC7';              // de-emphasis grey for context marks

let uid = 0;
const nextId = () => `c${++uid}`;

/* ------------------------------------------------------------------ scales */
function scale(dLo, dHi, rLo, rHi) {
  const span = (dHi - dLo) || 1;
  const f = (v) => rLo + ((v - dLo) / span) * (rHi - rLo);
  f.invert = (p) => dLo + ((p - rLo) / (rHi - rLo)) * span;
  f.domain = [dLo, dHi];
  return f;
}

/* Axis ticks on round numbers. An axis labelled 0 / 2.5 / 5 is read instantly;
   one labelled 0 / 2.37 / 4.74 makes the reader do arithmetic. */
function ticks(lo, hi, count = 5) {
  const span = (hi - lo) || 1;
  const raw = span / count;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = (norm >= 7.5 ? 10 : norm >= 3.5 ? 5 : norm >= 1.5 ? 2 : 1) * mag;
  const out = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) {
    out.push(Math.abs(t) < step * 1e-9 ? 0 : t);
  }
  return out;
}

const nice = (v) => {
  if (v === null || v === undefined || !isFinite(v)) return '—';
  const a = Math.abs(v);
  if (a !== 0 && (a < 0.001 || a >= 1e6)) return Number(v).toExponential(1);
  if (a >= 1000) return fmt(Math.round(v));
  if (a >= 10) return String(Math.round(v * 10) / 10);
  if (a >= 1) return String(Math.round(v * 100) / 100);
  return String(Math.round(v * 1000) / 1000);
};

/* ------------------------------------------------------------------ frame */
function frame(opts) {
  const { w = 720, h = 300, pad = { t: 14, r: 16, b: 34, l: 52 } } = opts;
  return {
    w, h, pad,
    iw: w - pad.l - pad.r,
    ih: h - pad.t - pad.b,
    x0: pad.l, y0: h - pad.b, x1: w - pad.r, y1: pad.t,
  };
}

function axes(f, x, y, opts = {}) {
  const xt = opts.xTicks || ticks(x.domain[0], x.domain[1], opts.nx || 6);
  const yt = opts.yTicks || ticks(y.domain[0], y.domain[1], opts.ny || 5);
  const fx = opts.fmtX || nice;
  const fy = opts.fmtY || nice;
  return `
    <g class="grid">
      ${yt.map((t) => `<line x1="${f.x0}" x2="${f.x1}" y1="${y(t)}" y2="${y(t)}"/>`).join('')}
    </g>
    <g class="axis">
      <line x1="${f.x0}" x2="${f.x1}" y1="${f.y0}" y2="${f.y0}"/>
      ${opts.hideXTicks ? '' : xt.map((t) => `
        <text x="${x(t)}" y="${f.y0 + 16}" text-anchor="middle">${esc(fx(t))}</text>`).join('')}
      ${yt.map((t) => `
        <text x="${f.x0 - 8}" y="${y(t) + 3.5}" text-anchor="end">${esc(fy(t))}</text>`).join('')}
      ${opts.xLabel ? `<text class="ttl" x="${f.x0 + f.iw / 2}" y="${f.h - 2}"
        text-anchor="middle">${esc(opts.xLabel)}</text>` : ''}
      ${opts.yLabel ? `<text class="ttl" transform="rotate(-90 12 ${f.y1 + f.ih / 2})"
        x="12" y="${f.y1 + f.ih / 2}" text-anchor="middle">${esc(opts.yLabel)}</text>` : ''}
    </g>`;
}

/* The wrapper carries the stylesheet so every chart inherits the same recessive
   grid, the same type scale, and the same hover affordance.
 *
 * The svg is capped at its design width inline. Without that cap, a CSS
 * max-width of 100% lets a small chart scale UP to fill whatever column it
 * lands in — a 330px Q-Q became a full-width poster with 40px axis labels. */
function svg(f, body, opts = {}) {
  const id = nextId();
  return `<figure class="chart-fig" data-cid="${id}">
    ${opts.title ? `<figcaption>${esc(opts.title)}
      ${opts.sub ? `<span>${esc(opts.sub)}</span>` : ''}</figcaption>` : ''}
    ${opts.legend || ''}
    <div class="chart-scroll">
      <svg class="chart" viewBox="0 0 ${f.w} ${f.h}" preserveAspectRatio="xMidYMid meet"
           style="max-width:${f.w}px" role="img"
           aria-label="${esc(opts.alt || opts.title || 'chart')}">
        ${body}
      </svg>
      <div class="chart-tip" hidden></div>
    </div>
    ${opts.foot ? `<p class="chart-foot">${opts.foot}</p>` : ''}
  </figure>`;
}

export function legend(items) {
  return `<div class="chart-legend">${items.map((i) => `
    <span class="lg"><i style="background:${i.color}"></i>${esc(i.label)}</span>`).join('')}</div>`;
}

/* A mark declares its own tooltip text; one delegated listener serves them all
   (see mountCharts). Attaching a listener per point would be thousands of them
   on a Manhattan plot. */
const tip = (text) => `data-tip="${esc(text)}"`;

/* ------------------------------------------------------------- Manhattan -- */
/* Position × significance. The alternating tone is a positional separator for
   chromosome boundaries, NOT series identity — so it is two steps of one grey
   rather than two hues, leaving colour free to mean "this crossed the line". */
export function manhattan(points, opts = {}) {
  const f = frame({ w: 860, h: 290, pad: { t: 14, r: 16, b: 40, l: 54 } });
  if (!points.length) return '';

  const order = [];
  const seen = new Set();
  for (const p of points) {
    if (!seen.has(p.chrom)) { seen.add(p.chrom); order.push(p.chrom); }
  }
  const byChrom = new Map(order.map((c) => [c, []]));
  for (const p of points) byChrom.get(p.chrom).push(p);

  // Lay chromosomes end to end, each scaled to its own span.
  const gap = 4;
  const totalGap = gap * (order.length - 1);
  const spans = order.map((c) => {
    const ps = byChrom.get(c);
    const lo = Math.min(...ps.map((p) => p.pos));
    const hi = Math.max(...ps.map((p) => p.pos));
    return { c, lo, hi, n: ps.length };
  });
  const totalN = spans.reduce((a, s) => a + s.n, 0);
  let cursor = f.x0;
  const placed = [];
  for (const s of spans) {
    const width = ((f.iw - totalGap) * s.n) / totalN;
    placed.push({ ...s, x0: cursor, x1: cursor + width });
    cursor += width + gap;
  }

  const maxY = Math.max(5, Math.ceil(Math.max(...points.map((p) => p.neglog10p)) + 0.5));
  const y = scale(0, maxY, f.y0, f.y1);
  const gw = opts.threshold || 7.301;      // -log10(5e-8)
  const sug = 5;

  let marks = '';
  let bands = '';
  placed.forEach((s, i) => {
    const px = scale(s.lo, s.hi, s.x0, s.x1);
    if (i % 2 === 1) {
      bands += `<rect x="${s.x0 - gap / 2}" y="${f.y1}" width="${s.x1 - s.x0 + gap}"
        height="${f.ih}" fill="var(--rule2)" opacity=".55"/>`;
    }
    for (const p of byChrom.get(s.c)) {
      const hit = p.neglog10p >= gw;
      marks += `<circle cx="${px(p.pos).toFixed(1)}" cy="${y(p.neglog10p).toFixed(1)}"
        r="${hit ? 3.6 : 1.7}" fill="${hit ? EMPH : (i % 2 ? '#7E8A96' : '#5E6975')}"
        ${hit ? `stroke="${SURFACE}" stroke-width="1.5"` : ''}
        ${tip(`${p.variant} · chr${p.chrom}:${fmt(p.pos)} · p = ${Number(Math.pow(10, -p.neglog10p)).toExponential(1)}`)}/>`;
    }
  });

  const labels = placed.map((s) => `<text x="${(s.x0 + s.x1) / 2}" y="${f.y0 + 15}"
    text-anchor="middle">${esc(s.c)}</text>`).join('');

  return svg(f, `
    ${bands}
    ${axes(f, scale(0, 1, f.x0, f.x1), y, { hideXTicks: true, yLabel: '−log₁₀ p' })}
    <g class="axis">${labels}</g>
    <line class="thresh sug" x1="${f.x0}" x2="${f.x1}" y1="${y(sug)}" y2="${y(sug)}"/>
    <line class="thresh gw" x1="${f.x0}" x2="${f.x1}" y1="${y(gw)}" y2="${y(gw)}"/>
    <text class="thresh-lbl" x="${f.x1 - 2}" y="${y(gw) - 5}" text-anchor="end">genome-wide 5×10⁻⁸</text>
    ${marks}`, {
    title: opts.title || 'Association across the genome',
    sub: opts.sub,
    alt: `Manhattan plot of ${fmt(points.length)} variants`,
    foot: opts.foot,
  });
}

/* --------------------------------------------------------------- QQ plot -- */
/* The calibration check. Points on the diagonal mean the null is behaving;
   a whole curve lifting off it is inflation, not discovery. */
export function qq(observed, expected, opts = {}) {
  const f = frame({ w: 330, h: 290, pad: { t: 14, r: 16, b: 40, l: 48 } });
  const n = Math.min(observed.length, expected.length);
  if (!n) return '';
  const hi = Math.max(Math.max(...observed), Math.max(...expected));
  const x = scale(0, hi, f.x0, f.x1);
  const y = scale(0, hi, f.y0, f.y1);

  let pts = '';
  for (let i = 0; i < n; i++) {
    pts += `<circle cx="${x(expected[i]).toFixed(1)}" cy="${y(observed[i]).toFixed(1)}"
      r="1.8" fill="${C1}" ${tip(`expected ${nice(expected[i])} · observed ${nice(observed[i])}`)}/>`;
  }
  return svg(f, `
    ${axes(f, x, y, { nx: 4, ny: 4, xLabel: 'expected −log₁₀ p', yLabel: 'observed' })}
    <line class="thresh diag" x1="${x(0)}" y1="${y(0)}" x2="${x(hi)}" y2="${y(hi)}"/>
    ${pts}`, {
    title: opts.title || 'Calibration (Q–Q)',
    sub: opts.sub,
    alt: 'Quantile-quantile plot of observed against expected p-values',
    foot: opts.foot,
  });
}

/* --------------------------------------------------------------- volcano -- */
/* Effect × significance. Emphasis form: one hue for the marks that cleared the
   threshold, grey for the rest — the story is "these few", not "these groups". */
export function volcano(rows, opts = {}) {
  const f = frame({ w: 560, h: 300, pad: { t: 14, r: 18, b: 40, l: 56 } });
  const pts = rows.filter((r) => r.x !== null && r.x !== undefined
                                 && isFinite(r.x) && isFinite(r.y));
  if (!pts.length) return '';
  const xs = pts.map((p) => p.x);
  const lim = Math.max(Math.abs(Math.min(...xs)), Math.abs(Math.max(...xs))) || 1;
  const x = scale(-lim * 1.05, lim * 1.05, f.x0, f.x1);
  const maxY = Math.max(...pts.map((p) => p.y)) * 1.08 || 1;
  const y = scale(0, maxY, f.y0, f.y1);
  const cut = opts.threshold;

  const marks = pts.map((p) => {
    const hit = cut !== undefined && p.y >= cut;
    return `<circle cx="${x(p.x).toFixed(1)}" cy="${y(p.y).toFixed(1)}"
      r="${hit ? 4 : 2.4}" fill="${hit ? C2 : DIM}"
      ${hit ? `stroke="${SURFACE}" stroke-width="1.5"` : ''}
      ${tip(p.label || '')}/>`;
  }).join('');

  // Direct-label the extremes only. A label on every point is unreadable.
  const top = [...pts].sort((a, b) => b.y - a.y).slice(0, opts.labels === 0 ? 0 : 3);
  const labels = top.map((p) => {
    const anchor = x(p.x) > f.x0 + f.iw * 0.7 ? 'end' : 'start';
    const dx = anchor === 'end' ? -7 : 7;
    return `<text class="pt-lbl" x="${x(p.x) + dx}" y="${y(p.y) + 3.5}"
      text-anchor="${anchor}">${esc(p.name || '')}</text>`;
  }).join('');

  return svg(f, `
    ${axes(f, x, y, { nx: 5, ny: 5, xLabel: opts.xLabel || 'effect',
                      yLabel: opts.yLabel || '−log₁₀ p' })}
    <line class="thresh zero" x1="${x(0)}" x2="${x(0)}" y1="${f.y0}" y2="${f.y1}"/>
    ${cut !== undefined ? `<line class="thresh sug" x1="${f.x0}" x2="${f.x1}"
      y1="${y(cut)}" y2="${y(cut)}"/>` : ''}
    ${marks}${labels}`, {
    title: opts.title, sub: opts.sub, foot: opts.foot,
    alt: opts.title || 'volcano plot',
  });
}

/* ----------------------------------------------------- Kaplan–Meier curves */
/* Change over time, and the series ARE the subject — the one genuinely
   categorical form here, so it takes the two measured hues and a legend. */
export function kaplanMeier(groups, opts = {}) {
  const f = frame({ w: 560, h: 300, pad: { t: 14, r: 18, b: 40, l: 56 } });
  const names = Object.keys(groups);
  if (!names.length) return '';
  const colors = [C1, C2, '#1baf7a'];

  let tMax = 0;
  let sMin = 1;
  for (const n of names) {
    const g = groups[n];
    const times = g.time || g.times || [];
    const surv = g.survival || g.s || [];
    if (times.length) tMax = Math.max(tMax, times[times.length - 1]);
    for (const v of surv) sMin = Math.min(sMin, v);
  }
  const x = scale(0, tMax || 1, f.x0, f.x1);
  const y = scale(Math.max(0, Math.floor(sMin * 10) / 10 - 0.05), 1, f.y0, f.y1);

  const paths = names.map((n, i) => {
    const g = groups[n];
    const times = g.time || g.times || [];
    const surv = g.survival || g.s || [];
    if (!times.length) return '';
    // Survival is a STEP function — it changes only at an event. Drawing it as
    // a smooth line implies deaths between events that never happened.
    let d = `M ${x(0)} ${y(1)}`;
    for (let k = 0; k < times.length; k++) {
      d += ` L ${x(times[k]).toFixed(1)} ${y(k ? surv[k - 1] : 1).toFixed(1)}`;
      d += ` L ${x(times[k]).toFixed(1)} ${y(surv[k]).toFixed(1)}`;
    }
    d += ` L ${x(tMax).toFixed(1)} ${y(surv[surv.length - 1]).toFixed(1)}`;
    return `<path class="km" d="${d}" stroke="${colors[i % colors.length]}"/>`;
  }).join('');

  return svg(f, `
    ${axes(f, x, y, { nx: 5, ny: 5, xLabel: opts.xLabel || 'time',
                      yLabel: 'survival probability',
                      fmtY: (v) => (v * 100).toFixed(0) + '%' })}
    ${paths}`, {
    title: opts.title || 'Survival by group',
    sub: opts.sub, foot: opts.foot,
    alt: 'Kaplan-Meier survival curves',
    legend: names.length > 1
      ? legend(names.map((n, i) => ({ label: n, color: colors[i % colors.length] })))
      : '',
  });
}

/* ------------------------------------------------------------- histogram -- */
export function histogram(bins, opts = {}) {
  const f = frame({ w: 520, h: 250, pad: { t: 14, r: 16, b: 40, l: 54 } });
  if (!bins.length) return '';
  const maxC = Math.max(...bins.map((b) => b.count)) || 1;
  const x = scale(bins[0].lo, bins[bins.length - 1].hi, f.x0, f.x1);
  const y = scale(0, maxC, f.y0, f.y1);
  const gapPx = 2;   // the surface gap; adjacent bars separate by air, not stroke

  const bars = bins.map((b) => {
    const bx = x(b.lo);
    const bw = Math.max(1, x(b.hi) - x(b.lo) - gapPx);
    const bh = f.y0 - y(b.count);
    if (bh <= 0) return '';
    const r = Math.min(4, bw / 2, bh);
    return `<path d="M${bx} ${f.y0} V${y(b.count) + r}
      q0 ${-r} ${r} ${-r} h${bw - 2 * r} q${r} 0 ${r} ${r} V${f.y0} Z"
      fill="${opts.color || C1}"
      ${tip(`${nice(b.lo)} – ${nice(b.hi)} · ${fmt(b.count)}`)}/>`;
  }).join('');

  return svg(f, `
    ${axes(f, x, y, { nx: 5, ny: 4, xLabel: opts.xLabel, yLabel: opts.yLabel || 'samples' })}
    ${bars}`, {
    title: opts.title, sub: opts.sub, foot: opts.foot, alt: opts.title || 'histogram',
  });
}

/* --------------------------------------------------------------- scatter -- */
/* Emphasis form: outliers in the status hue, the cohort in grey. */
export function scatter(points, opts = {}) {
  const f = frame({ w: 520, h: 290, pad: { t: 14, r: 18, b: 40, l: 56 } });
  const pts = points.filter((p) => isFinite(p.x) && isFinite(p.y));
  if (!pts.length) return '';
  const xs = pts.map((p) => p.x);
  const ys = pts.map((p) => p.y);
  const padX = (Math.max(...xs) - Math.min(...xs)) * 0.06 || 1;
  const padY = (Math.max(...ys) - Math.min(...ys)) * 0.08 || 1;
  const x = scale(Math.min(...xs) - padX, Math.max(...xs) + padX, f.x0, f.x1);
  const y = scale(Math.min(...ys) - padY, Math.max(...ys) + padY, f.y0, f.y1);

  const marks = pts.map((p) => `<circle cx="${x(p.x).toFixed(1)}" cy="${y(p.y).toFixed(1)}"
    r="${p.flag ? 4.5 : 3}" fill="${p.flag ? EMPH : DIM}"
    ${p.flag ? `stroke="${SURFACE}" stroke-width="2"` : ''}
    ${tip(p.label || '')}/>`).join('');

  return svg(f, `
    ${axes(f, x, y, { nx: 5, ny: 5, xLabel: opts.xLabel, yLabel: opts.yLabel })}
    ${opts.diagonal ? `<line class="thresh diag"
      x1="${f.x0}" y1="${f.y0}" x2="${f.x1}" y2="${f.y1}"/>` : ''}
    ${marks}`, {
    title: opts.title, sub: opts.sub, foot: opts.foot, alt: opts.title || 'scatter plot',
    legend: opts.flagLabel
      ? legend([{ label: opts.flagLabel, color: EMPH },
                { label: opts.baseLabel || 'within range', color: DIM }])
      : '',
  });
}

/* -------------------------------------------------------------- dumbbell -- */
/* Before → after per item: one hue, two shades, joined by a rule. Two hues here
   would spend the identity channel on what position already shows. */
export function dumbbell(rows, opts = {}) {
  const rowH = 22;
  const f = frame({ w: 560, h: 34 + rows.length * rowH + 30,
                    pad: { t: 14, r: 92, b: 34, l: 128 } });
  if (!rows.length) return '';
  const all = rows.flatMap((r) => [r.a, r.b]).filter(isFinite);
  const x = scale(Math.min(...all), Math.max(...all), f.x0, f.x1);

  // The tooltip uses the axis formatter. A tooltip reading 0.085 beside an axis
  // reading 60.0% makes the reader do a unit conversion to check the chart.
  const fv = opts.fmt || nice;
  const body = rows.map((r, i) => {
    const cy = f.y1 + 12 + i * rowH;
    return `
      <text class="row-lbl" x="${f.x0 - 10}" y="${cy + 3.5}" text-anchor="end">${esc(r.label)}</text>
      <line class="dbl" x1="${x(r.a)}" x2="${x(r.b)}" y1="${cy}" y2="${cy}"/>
      <circle cx="${x(r.a)}" cy="${cy}" r="4.5" fill="#9EC5E8" stroke="${SURFACE}"
        stroke-width="2" ${tip(`${r.label} · ${r.aLabel}: ${fv(r.a)}`)}/>
      <circle cx="${x(r.b)}" cy="${cy}" r="4.5" fill="${C1}" stroke="${SURFACE}"
        stroke-width="2" ${tip(`${r.label} · ${r.bLabel}: ${fv(r.b)}`)}/>
      <text class="row-val" x="${f.x1 + 10}" y="${cy + 3.5}">${esc(r.note || '')}</text>`;
  }).join('');

  const xt = ticks(x.domain[0], x.domain[1], 4);
  return svg(f, `
    <g class="grid">${xt.map((t) => `<line x1="${x(t)}" x2="${x(t)}"
      y1="${f.y1}" y2="${f.y1 + rows.length * rowH + 6}"/>`).join('')}</g>
    <g class="axis">${xt.map((t) => `<text x="${x(t)}"
      y="${f.y1 + rows.length * rowH + 22}" text-anchor="middle">${esc(opts.fmt
        ? opts.fmt(t) : nice(t))}</text>`).join('')}</g>
    ${body}`, {
    title: opts.title, sub: opts.sub, foot: opts.foot, alt: opts.title || 'dumbbell chart',
    legend: legend([{ label: rows[0].aLabel, color: '#9EC5E8' },
                    { label: rows[0].bLabel, color: C1 }]),
  });
}

/* -------------------------------------------------------------- dot plot -- */
/* One value per item against a reference line — built for spotting the one that
   sits away from the rest. */
export function dotplot(rows, opts = {}) {
  const f = frame({ w: 560, h: 210, pad: { t: 14, r: 18, b: 42, l: 56 } });
  const pts = rows.filter((r) => isFinite(r.value));
  if (!pts.length) return '';
  const vals = pts.map((p) => p.value);
  const y = scale(0, Math.max(...vals) * 1.12 || 1, f.y0, f.y1);
  const x = scale(0, pts.length - 1 || 1, f.x0, f.x1);

  const marks = pts.map((p, i) => `<circle cx="${x(i).toFixed(1)}" cy="${y(p.value).toFixed(1)}"
    r="${p.flag ? 5 : 3.4}" fill="${p.flag ? EMPH : C1}" stroke="${SURFACE}" stroke-width="1.6"
    ${tip(p.label || '')}/>`).join('');

  return svg(f, `
    ${axes(f, x, y, { hideXTicks: true, ny: 4, yLabel: opts.yLabel,
                      fmtY: opts.fmtY, xLabel: opts.xLabel })}
    ${opts.reference !== undefined ? `
      <line class="thresh med" x1="${f.x0}" x2="${f.x1}"
        y1="${y(opts.reference)}" y2="${y(opts.reference)}"/>
      <text class="thresh-lbl" x="${f.x1 - 2}" y="${y(opts.reference) - 5}"
        text-anchor="end">${esc(opts.referenceLabel || 'median')}</text>` : ''}
    ${marks}`, {
    title: opts.title, sub: opts.sub, foot: opts.foot, alt: opts.title || 'dot plot',
  });
}

/* ----------------------------------------------------------- interaction -- */
/* One delegated listener for every chart on the page. A Manhattan plot has tens
   of thousands of marks; per-mark listeners would be a memory leak with a
   rendering cost. */
export function mountCharts(root = document) {
  if (root.__chartsMounted) return;
  root.__chartsMounted = true;

  const show = (e) => {
    const mark = e.target.closest('[data-tip]');
    const fig = e.target.closest('.chart-scroll');
    if (!fig) return;
    const tipEl = fig.querySelector('.chart-tip');
    if (!tipEl) return;
    if (!mark || !mark.dataset.tip) { tipEl.hidden = true; return; }
    tipEl.textContent = mark.dataset.tip;
    tipEl.hidden = false;
    const box = fig.getBoundingClientRect();
    const m = mark.getBoundingClientRect();
    const left = m.left - box.left + m.width / 2;
    tipEl.style.left = `${Math.max(4, Math.min(left, box.width - 4))}px`;
    tipEl.style.top = `${m.top - box.top - 10}px`;
  };

  root.addEventListener('mousemove', show, { passive: true });
  root.addEventListener('mouseleave', (e) => {
    const fig = e.target.closest && e.target.closest('.chart-scroll');
    if (fig) { const t = fig.querySelector('.chart-tip'); if (t) t.hidden = true; }
  }, true);
}
