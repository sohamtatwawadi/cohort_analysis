/* API client. The browser NEVER computes a figure — it asks the server. */

async function req(path, opts = {}) {
  const res = await fetch(path, {
    headers: { 'Content-Type': 'application/json' }, ...opts,
  });
  if (!res.ok) {
    let detail = res.statusText;
    try {
      const d = await res.json();
      detail = typeof d.detail === 'string' ? d.detail
        : (d.detail && d.detail.message) || detail;
      const err = new Error(detail);
      err.detail = d.detail;
      throw err;
    } catch (e) {
      if (e instanceof Error && e.message !== res.statusText) throw e;
      throw new Error(detail);
    }
  }
  return res.json();
}

const post = (p, b) => req(p, { method: 'POST', body: JSON.stringify(b || {}) });
const get = (p) => req(p);

export const api = {
  meta: () => get('/api/meta'),
  resolve: (criteria, name) => post('/api/cohort/resolve', { criteria, name }),
  suggestions: (criteria, name) => post('/api/cohort/suggestions', { criteria, name }),
  module: (module, criteria, name, query) =>
    post(`/api/module/${module}${query ? '?' + new URLSearchParams(query) : ''}`,
         { criteria, name }),
  variantDetail: (id, criteria, name) => post(`/api/detail/variant/${id}`, { criteria, name }),
  subjectDetail: (id, criteria, name) => post(`/api/detail/subject/${id}`, { criteria, name }),

  sampleQc: (criteria, name) => post('/api/qc/samples', { criteria, name }),
  denominator: (criteria, name) => post('/api/governance/denominator', { criteria, name }),
  manifest: (criteria, name, module) =>
    post(`/api/governance/manifest${module ? '?module=' + module : ''}`, { criteria, name }),
  audit: () => get('/api/governance/audit'),
  export: (payload) => post('/api/governance/export', payload),

  cohorts: () => get('/api/cohorts'),
  saveCohort: (payload) => post('/api/cohorts', payload),

  compile: (question) => post('/api/compile', { question }),

  research: {
    meta: () => get('/api/research/meta'),
    projects: () => get('/api/research/projects'),
    datasets: (pid) => get(`/api/research/projects/${pid}/datasets`),
    dataset: (did) => get(`/api/research/datasets/${did}`),
    ingest: (pid, body) => post(`/api/research/projects/${pid}/datasets/ingest`, body),
    upload: async (pid, form) => {
      // No Content-Type header — the browser sets the multipart boundary, and
      // overriding it produces a body the server cannot parse.
      const res = await fetch(`/api/research/projects/${pid}/datasets/upload`,
                              { method: 'POST', body: form });
      if (!res.ok) {
        const d = await res.json().catch(() => ({}));
        const err = new Error(typeof d.detail === 'string' ? d.detail
          : (d.detail && d.detail.message) || res.statusText);
        err.detail = d.detail;
        throw err;
      }
      return res.json();
    },
    del: (did, pid) => req(`/api/research/datasets/${did}?project_id=${pid}`,
                           { method: 'DELETE' }),
    override: (did, body) => post(`/api/research/datasets/${did}/override`, body),
    jobs: (did) => get(`/api/research/datasets/${did}/jobs`),
    submit: (did, body) => post(`/api/research/datasets/${did}/jobs`, body),
    job: (jid) => get(`/api/research/jobs/${jid}`),
    result: (jid) => get(`/api/research/jobs/${jid}/result`),
    estimate: (a, n, v) =>
      get(`/api/research/estimate?analysis=${a}&n_samples=${n}&n_variants=${v}`),
    power: (body) => post('/api/research/power', body),
  },
};
