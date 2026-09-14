"""Sample QC — does the data support any analysis at all?

A contaminated sample or a sample swap invalidates every carrier rate in the
cohort no matter how careful the denominators are. Coverage answers "was this
gene looked at"; this answers "should we believe what came back".

The metrics come from `run`, computed AT INGEST over every PASS call in the
source file. That turned out to be necessary rather than optional: §3.2 has the
loader archive the intronic/intergenic rows after fingerprinting, and `finding`
retains only the reviewable subset — one to seven rows per sample in practice.
A Ti/Tv computed from three reportable findings is noise wearing a number's
clothes. Where the stored metrics are absent the service falls back to the
retained findings and reports which basis it used, so nobody reads a
panel-sized ratio as a genome-wide one.

    Ti/Tv        transitions / transversions. Real coding sequence sits near
                 3.0 (genome-wide ~2.0) because deaminated methyl-CpG gives
                 C>T. A depressed ratio means the extra calls are random, i.e.
                 false positives — the classic signature of a noisy or
                 contaminated sample.
    Het/Hom      heterozygous / homozygous-alt. Contamination adds a second
                 individual's alleles, which read as spurious heterozygotes and
                 push this up. A very low value suggests consanguinity or a
                 partly-uncalled sample.
    Het VAF      allele balance at heterozygous sites. A true het sits at 0.5;
                 contamination drags it off-centre.
    Depth        mean depth and % bases >=20x, straight from the run record.
    Call count   an outlier count either way is a processing problem.
    Sex          X heterozygosity AS THE FILE REPORTED IT against recorded sex.

COHORT-RELATIVE, NOT ABSOLUTE. A panel, an exome and a genome have completely
different expected values for most of these, so a fixed threshold either passes
everything or fails a whole assay type. Samples are compared against the median
of their own cohort using median absolute deviation, which — unlike mean and SD
— is not dragged around by the very outliers we are looking for.

Metrics computed from too few calls are reported as unavailable rather than
given a number nobody should act on, the same discipline as the min-carriers
guardrail in the burden test.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from .. import db

# Below this many usable SNV calls a per-sample ratio is noise. A 1-gene
# thalassemia panel will never clear it, and that is the correct outcome.
MIN_CALLS_FOR_RATIO = 20

# How many scaled MADs from the cohort median before a sample is flagged.
MAD_WARN = 3.0
MAD_FAIL = 5.0

# Floors that hold regardless of assay type. A Ti/Tv near 0.5 is what you get
# from calling noise uniformly across the genome; no real assay produces it.
TI_TV_FLOOR = 1.0
PCT20X_FLOOR = 80.0

TRANSITIONS = ("AG", "GA", "CT", "TC")


def _metrics_sql() -> str:
    """Per-sample metrics over the resolved cohort.

    `zygosity_raw` is used for the sex check, not the derived `zygosity`:
    derive_zygosity already rewrote X calls in males to Hemizygous using the
    recorded sex, so checking the derived value against recorded sex would be
    circular and always agree.
    """
    ti = " OR ".join("(f.ref_allele || f.alt_allele) = '{}'".format(t)
                     for t in TRANSITIONS)
    return """
        SELECT
          cr.run_id, cr.subject_id, r.test_code, r.mean_depth, r.pct_bases_20x,
          r.qc_status, r.reference_build, r.pipeline_version, s.sex AS reported_sex,
          r.qc_n_called, r.qc_ti_tv, r.qc_het_hom, r.qc_mean_het_vaf, r.qc_x_het_rate,
          COUNT(f.finding_id) AS n_calls,
          COUNT(*) FILTER (WHERE length(f.ref_allele) = 1
                             AND length(f.alt_allele) = 1
                             AND f.ref_allele <> f.alt_allele)          AS n_snv,
          COUNT(*) FILTER (WHERE ({ti}))                                AS n_ti,
          COUNT(*) FILTER (WHERE length(f.ref_allele) = 1
                             AND length(f.alt_allele) = 1
                             AND f.ref_allele <> f.alt_allele
                             AND NOT ({ti}))                            AS n_tv,
          COUNT(*) FILTER (WHERE f.zygosity_raw = 'Heterozygous')       AS n_het,
          COUNT(*) FILTER (WHERE f.zygosity_raw = 'Homozygous')         AS n_hom,
          AVG(f.vaf) FILTER (WHERE f.zygosity_raw = 'Heterozygous'
                               AND f.vaf IS NOT NULL)                   AS mean_het_vaf,
          AVG(f.depth) FILTER (WHERE f.depth IS NOT NULL)               AS mean_call_depth,
          COUNT(*) FILTER (WHERE f.chrom IN ('X','chrX'))               AS n_x,
          COUNT(*) FILTER (WHERE f.chrom IN ('X','chrX')
                             AND f.zygosity_raw = 'Heterozygous')       AS n_x_het
        FROM cohort_run cr
        JOIN run r USING (run_id)
        JOIN subject s ON s.subject_id = cr.subject_id
        LEFT JOIN finding f ON f.run_id = cr.run_id
        GROUP BY 1,2,3,4,5,6,7,8,9,10,11,12,13,14
    """.format(ti=ti)


def _median(vals: List[float]) -> Optional[float]:
    v = sorted(x for x in vals if x is not None)
    if not v:
        return None
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0


def _mad(vals: List[float], med: Optional[float]) -> Optional[float]:
    if med is None:
        return None
    dev = [abs(x - med) for x in vals if x is not None]
    m = _median(dev)
    # 1.4826 makes MAD comparable to a standard deviation for normal data, so
    # "3 MADs" reads like "3 sigma" without inheriting SD's sensitivity to the
    # outliers being hunted.
    return (m * 1.4826) if m else None


def _deviation(value: Optional[float], med: Optional[float],
               mad: Optional[float]) -> Optional[float]:
    if value is None or med is None or not mad:
        return None
    return (value - med) / mad


def sample_qc(criteria: Dict[str, Any]) -> Dict[str, Any]:
    """Per-sample QC for the resolved cohort, with cohort-relative outliers."""
    rows = db.rows(_metrics_sql())

    samples: List[Dict[str, Any]] = []
    for r in rows:
        n_ti, n_tv = int(r["n_ti"] or 0), int(r["n_tv"] or 0)
        n_het, n_hom = int(r["n_het"] or 0), int(r["n_hom"] or 0)
        n_snv = int(r["n_snv"] or 0)

        # Prefer the metrics computed at ingest over the whole call set. Falling
        # back to the retained findings is only defensible when there are enough
        # of them, which for a reportable subset there usually are not — so the
        # source is reported alongside the number.
        stored = r["qc_ti_tv"] is not None or r["qc_het_hom"] is not None
        if stored:
            ti_tv = r["qc_ti_tv"]
            het_hom = r["qc_het_hom"]
            het_vaf = r["qc_mean_het_vaf"]
            x_het = r["qc_x_het_rate"]
            basis = "full call set ({:,} calls)".format(int(r["qc_n_called"] or 0))
            enough = True
        else:
            enough = n_snv >= MIN_CALLS_FOR_RATIO
            ti_tv = round(n_ti / n_tv, 3) if (enough and n_tv) else None
            het_hom = round(n_het / n_hom, 3) if (enough and n_hom) else None
            het_vaf = (round(float(r["mean_het_vaf"]), 4)
                       if r["mean_het_vaf"] is not None else None)
            x_het = (round(int(r["n_x_het"]) / int(r["n_x"]), 3)
                     if int(r["n_x"] or 0) >= 5 else None)
            basis = ("retained findings ({} calls)".format(n_snv) if enough
                     else "too few retained calls")

        samples.append({
            "run_id": r["run_id"],
            "subject_id": r["subject_id"],
            "test_code": r["test_code"],
            "reported_sex": r["reported_sex"],
            "qc_status": r["qc_status"],
            "reference_build": r["reference_build"],
            "pipeline_version": r["pipeline_version"],
            "n_calls": int(r["n_calls"] or 0),
            "n_snv": n_snv,
            "ti_tv": ti_tv,
            "het_hom": het_hom,
            "mean_het_vaf": het_vaf,
            "metrics_basis": basis,
            "mean_depth": round(float(r["mean_depth"]), 1)
                          if r["mean_depth"] is not None else None,
            "pct_bases_20x": round(float(r["pct_bases_20x"]), 1)
                             if r["pct_bases_20x"] is not None else None,
            "mean_call_depth": round(float(r["mean_call_depth"]), 1)
                               if r["mean_call_depth"] is not None else None,
            "x_het_rate": x_het,
            "metrics_available": enough,
        })

    cohort = _cohort_stats(samples)
    for s in samples:
        s.update(_flag(s, cohort))

    flagged = [s for s in samples if s["status"] != "pass"]
    flagged.sort(key=lambda s: (s["status"] != "fail", -len(s["flags"])))

    return {
        "samples": samples,
        "cohort": cohort,
        "summary": {
            "total": len(samples),
            "pass": sum(1 for s in samples if s["status"] == "pass"),
            "warn": sum(1 for s in samples if s["status"] == "warn"),
            "fail": sum(1 for s in samples if s["status"] == "fail"),
            "not_assessable": sum(1 for s in samples if not s["metrics_available"]),
        },
        "flagged": flagged,
        "method": (
            "Samples are compared against the median of this cohort using median "
            "absolute deviation, not against fixed thresholds — a panel, an exome "
            "and a genome have different expected values for most of these. "
            "Flagged at {:g} MADs, failed at {:g}.".format(MAD_WARN, MAD_FAIL)),
        "caveat": (
            "A sample flagged here has not been shown to be wrong. It is unlike "
            "the rest of this cohort, which is a reason to look before quoting "
            "any rate that includes it."),
        "not_assessable_note": (
            "Ratio metrics need at least {} single-nucleotide calls to be stable. "
            "Small panels will not reach that, and are reported as not assessable "
            "rather than given a number.".format(MIN_CALLS_FOR_RATIO)),
    }


def _cohort_stats(samples: List[Dict[str, Any]]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for key in ("ti_tv", "het_hom", "mean_het_vaf", "mean_depth",
                "pct_bases_20x", "n_calls", "mean_call_depth"):
        vals = [s[key] for s in samples if s.get(key) is not None]
        med = _median(vals)
        out[key] = {"median": round(med, 3) if med is not None else None,
                    "mad": round(_mad(vals, med) or 0, 3),
                    "n": len(vals)}
    return out


def _flag(s: Dict[str, Any], cohort: Dict[str, Any]) -> Dict[str, Any]:
    """Decide a sample's status and say, in numbers, why."""
    flags: List[Dict[str, str]] = []

    def dev(key):
        return _deviation(s.get(key), cohort[key]["median"], cohort[key]["mad"])

    def add(level, label, detail):
        flags.append({"level": level, "label": label, "detail": detail})

    # --- the lab's own verdict comes first --------------------------------
    if s["qc_status"] and s["qc_status"] != "Pass":
        add("fail", s["qc_status"], "The sequencing pipeline already flagged this run.")

    # --- Ti/Tv ------------------------------------------------------------
    d = dev("ti_tv")
    if s["ti_tv"] is not None:
        if s["ti_tv"] < TI_TV_FLOOR:
            add("fail", "Ti/Tv {:.2f}".format(s["ti_tv"]),
                "Below {:.1f} in any assay means the extra calls are close to random — "
                "these are false positives, not biology.".format(TI_TV_FLOOR))
        elif d is not None and d <= -MAD_FAIL:
            add("fail", "Ti/Tv {:.2f}".format(s["ti_tv"]),
                "{:.1f} MADs below the cohort median of {:.2f}.".format(
                    abs(d), cohort["ti_tv"]["median"]))
        elif d is not None and d <= -MAD_WARN:
            add("warn", "Ti/Tv {:.2f}".format(s["ti_tv"]),
                "{:.1f} MADs below the cohort median of {:.2f} — possible "
                "false-positive excess.".format(abs(d), cohort["ti_tv"]["median"]))

    # --- Het/Hom ----------------------------------------------------------
    d = dev("het_hom")
    if s["het_hom"] is not None and d is not None:
        if d >= MAD_FAIL:
            add("fail", "Het/Hom {:.2f}".format(s["het_hom"]),
                "{:.1f} MADs above the cohort median of {:.2f}. Excess heterozygosity "
                "is the signature of a contaminated or mixed sample.".format(
                    d, cohort["het_hom"]["median"]))
        elif d >= MAD_WARN:
            add("warn", "Het/Hom {:.2f}".format(s["het_hom"]),
                "{:.1f} MADs above the cohort median of {:.2f}.".format(
                    d, cohort["het_hom"]["median"]))
        elif d <= -MAD_WARN:
            add("warn", "Het/Hom {:.2f}".format(s["het_hom"]),
                "{:.1f} MADs below the cohort median — consanguinity, or an "
                "incompletely called sample.".format(abs(d)))

    # --- allele balance ---------------------------------------------------
    if s["mean_het_vaf"] is not None:
        off = abs(s["mean_het_vaf"] - 0.5)
        if off >= 0.10:
            add("fail", "Het allele balance {:.2f}".format(s["mean_het_vaf"]),
                "A true heterozygote sits near 0.50. This far off-centre indicates "
                "contamination or systematic allelic dropout.")
        elif off >= 0.06:
            add("warn", "Het allele balance {:.2f}".format(s["mean_het_vaf"]),
                "Drifting from the expected 0.50.")

    # --- coverage ---------------------------------------------------------
    if s["pct_bases_20x"] is not None and s["pct_bases_20x"] < PCT20X_FLOOR:
        add("fail", "{:.0f}% bases at 20x".format(s["pct_bases_20x"]),
            "Below {:.0f}% coverage, absence of a call is uninformative — this "
            "sample cannot contribute to a denominator.".format(PCT20X_FLOOR))
    d = dev("mean_depth")
    if d is not None and d <= -MAD_FAIL:
        add("fail", "Mean depth {:.0f}x".format(s["mean_depth"]),
            "{:.1f} MADs below the cohort median of {:.0f}x.".format(
                abs(d), cohort["mean_depth"]["median"]))
    elif d is not None and d <= -MAD_WARN:
        add("warn", "Mean depth {:.0f}x".format(s["mean_depth"]),
            "{:.1f} MADs below the cohort median of {:.0f}x.".format(
                abs(d), cohort["mean_depth"]["median"]))

    # --- call count -------------------------------------------------------
    d = dev("n_calls")
    if d is not None and abs(d) >= MAD_FAIL and s["n_calls"] > 0:
        add("warn", "{} calls".format(s["n_calls"]),
            "{:.1f} MADs {} the cohort median of {:.0f}. Compare the test code "
            "before reading anything into it — a different panel legitimately "
            "yields a different count.".format(
                abs(d), "above" if d > 0 else "below", cohort["n_calls"]["median"]))

    # --- sex concordance --------------------------------------------------
    if s["x_het_rate"] is not None and s["reported_sex"]:
        inferred = "M" if s["x_het_rate"] < 0.10 else "F"
        if inferred != s["reported_sex"]:
            add("fail", "Sex mismatch",
                "X heterozygosity of {:.2f} reads as {}, but the record says {}. "
                "Usually a sample swap or a data-entry error — resolve it before "
                "this subject contributes to anything.".format(
                    s["x_het_rate"], inferred, s["reported_sex"]))

    status = ("fail" if any(f["level"] == "fail" for f in flags)
              else "warn" if flags else "pass")
    return {"flags": flags, "status": status}
