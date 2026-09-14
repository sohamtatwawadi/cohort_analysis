"""Survival, penetrance and age-of-onset — Part II §4.7.

    "Penetrance guardrail: penetrance estimated from an ascertained clinical
     cohort is biased upward, often severely. The platform must detect
     ascertainment from the data profile and label the output accordingly
     rather than presenting a clean curve."

The detection lives in profile._assess_ascertainment; this module refuses to
present a penetrance estimate without it. That is the whole point: a
cumulative-incidence curve from a cohort of families who came to attention
BECAUSE they were affected is not a penetrance estimate, and it is
indistinguishable from one by eye.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from ..jobs import register
from ..stats import survival as surv
from ..types import MISSING, as_float


ASCERTAINMENT_LABEL = "ASCERTAINED COHORT — penetrance biased upward"


@register("survival")
def survival_job(spec: Dict[str, Any], context: Dict[str, Any]) -> Dict[str, Any]:
    """Kaplan-Meier / Cox / cumulative incidence, with the ascertainment label.

    Spec keys: time, event, group_by (optional carrier definition or phenotype
    column), covariates, event_type (for competing risks), analysis_kind.
    """
    from .. import registry as reg
    from .. import store as research_store

    log = context.get("log", lambda m: None)
    dataset_id = context["dataset_id"]
    ph = research_store.load_phenotypes(dataset_id)
    if ph is None:
        raise ValueError("Survival analysis requires phenotypes with time and event.")

    time_col = spec["time"]
    event_col = spec["event"]
    for col in (time_col, event_col):
        if col not in ph.columns:
            raise ValueError("Phenotype column '{}' not found.".format(col))

    time = as_float(ph.get(time_col))
    event = as_float(ph.get(event_col))

    groups = None
    group_label = None
    if spec.get("group_by"):
        group_label = spec["group_by"]
        if group_label in ph.columns:
            groups = np.asarray(ph.get(group_label))
        else:
            groups = _carrier_groups(dataset_id, ph, spec)

    finite = np.isfinite(time) & np.isfinite(event)
    if finite.sum() < 10:
        raise ValueError("Fewer than 10 subjects have both time and event recorded.")

    log("fitting Kaplan-Meier on {} subjects".format(int(finite.sum())))
    km = surv.kaplan_meier(time[finite], event[finite],
                           groups[finite] if groups is not None else None)

    result: Dict[str, Any] = {
        "analysis": "survival",
        "time_column": time_col,
        "event_column": event_col,
        "group_by": group_label,
        "n": int(finite.sum()),
        "n_events": int(np.nansum(event[finite])),
        "kaplan_meier": km,
    }

    if groups is not None and len(set(str(g) for g in groups[finite])) > 1:
        result["logrank"] = surv.logrank_test(
            time[finite], event[finite], groups[finite])

    # ---- Cox -------------------------------------------------------------
    cov_names: List[str] = []
    cols: List[np.ndarray] = []
    if groups is not None:
        levels = sorted(set(str(g) for g in groups[finite]))
        if len(levels) == 2:
            cols.append((np.array([str(g) for g in groups[finite]]) == levels[1]).astype(float))
            cov_names.append("{}={}".format(group_label, levels[1]))
    for cov in spec.get("covariates", []):
        if cov in ph.columns:
            cols.append(as_float(ph.get(cov))[finite])
            cov_names.append(cov)

    if cols:
        X = np.column_stack(cols)
        log("fitting Cox model with {} covariate(s)".format(len(cov_names)))
        cox = surv.cox_ph(time[finite], event[finite], X, cov_names)
        result["cox"] = cox
        beta = cox.get("beta") or cox.get("betas")
        if beta is not None:
            result["ph_assumption"] = surv.ph_assumption_test(
                time[finite], event[finite], X, beta, cov_names)

    # ---- competing risks --------------------------------------------------
    if spec.get("event_type") and spec["event_type"] in ph.columns:
        et = np.asarray(ph.get(spec["event_type"]))[finite]
        result["cumulative_incidence"] = surv.cumulative_incidence(
            time[finite], event[finite], et,
            groups[finite] if groups is not None else None,
            cause=spec.get("cause"))

    # ---- the penetrance guardrail (§4.7) ----------------------------------
    profile = reg.get_profile(dataset_id) or {}
    ascertained = bool(profile.get("ascertained"))
    result["penetrance"] = {
        "reported_as_penetrance": not ascertained,
        "ascertained": ascertained,
        "stamp": ASCERTAINMENT_LABEL if ascertained else None,
        "rationale": profile.get("ascertainment_rationale", ""),
        "caveat": (
            "This cohort appears to be ascertained (referral-selected). Families "
            "enter such a cohort BECAUSE they are affected, so the cumulative "
            "incidence shown is an upper bound on penetrance, not an estimate of "
            "it — often severely so. An unbiased penetrance estimate needs "
            "carriers identified independently of their disease status."
            if ascertained else
            "No ascertainment signal was detected in the data profile. Confirm "
            "that carriers were identified independently of disease status before "
            "quoting these as population penetrance estimates."),
    }

    result["diagnostics"] = {
        "n_analysed": int(finite.sum()),
        "n_dropped_missing": int(len(time) - finite.sum()),
        "ascertained": ascertained,
        "ph_assumption": result.get("ph_assumption"),
    }
    return result


def _carrier_groups(dataset_id: str, ph, spec: Dict[str, Any]) -> Optional[np.ndarray]:
    """Carrier vs non-carrier for a named variant or gene, aligned to phenotypes."""
    from .. import store as research_store

    variants = spec.get("carrier_variants")
    if not variants:
        return None
    gm = research_store.load_genotypes(dataset_id)
    index = {v.key: i for i, v in enumerate(gm.variants)}
    rows = [index[k] for k in variants if k in index]
    if not rows:
        raise ValueError("None of the carrier-defining variants are in this dataset.")

    carrier = (gm.dosages[rows, :] > 0).any(axis=0)
    by_sample = {s: bool(c) for s, c in zip(gm.sample_ids, carrier)}
    return np.array(["Carrier" if by_sample.get(s) else "Non-carrier"
                     for s in ph.sample_ids])
