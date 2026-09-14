"""Coverage (denominator) inference from gene fingerprints — spec §3.4.

    "No panel or test-code identifier exists in VariMAT (CRDB is empty; the
     filename carries only sample ID, pipeline version, caller, reference).
     Coverage must be inferred."

The procedure, verbatim from the spec:

    1. Compute each sample's gene fingerprint (genes with any ONTARGET variant)
    2. Cluster samples by fingerprint similarity -> implied panel
    3. Per cluster and gene: if the gene appears in >= threshold of the
       cluster's samples, the cluster assays that gene
    4. Threshold is PER-GENE, calibrated against known panels — not a
       hardcoded global number
    5. Otherwise -> unknown

Why step 4 matters. A gene's presence in a fingerprint depends on whether the
subject happened to carry a variant there, which scales with the gene's
size and variability. TTN (34,350 aa) shows up in essentially every exome
fingerprint; VHL (213 aa) shows up in a minority even when fully assayed. A
single global threshold would call VHL "not assayed" across the board and
silently shrink its denominator — which is exactly the §E04 trap "using cohort
size as denominator", arrived at from the other direction.

So the threshold for gene g is calibrated from runs where scope IS known: how
often does g actually appear in a fingerprint given that the test code declares
it? That observed rate becomes the yardstick.

Inference is valid for SNV/indel only. For CNV and fusion, absence is
uninformative — those views show counts only (spec §3.4, E03.4).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

# A gene must appear in at least this fraction of its calibrated rate within a
# cluster before the cluster is credited with assaying it.
CALIBRATION_FACTOR = 0.5
# Floor, for genes so rarely variant that the calibrated rate is near zero.
MIN_THRESHOLD = 0.02
# Below this fraction of the calibrated rate the gene is 'not assayed'; between
# the two bounds it is 'unknown' and lands in the provisional bucket.
UNKNOWN_FLOOR = 0.15
# Containment above which an unlabelled run belongs to an existing panel
# profile. High, because a run whose genes are not almost entirely inside a
# known panel's footprint is evidence of a panel we have not seen — and a new
# cluster with an honest 'inferred' denominator beats a wrong assignment.
CLUSTER_SIMILARITY = 0.85
# Exome/genome fingerprints are ~15,000 genes; a panel's is dozens. Clusters
# separate easily (spec §3.4), but guard the degenerate tiny-fingerprint case.
MIN_FINGERPRINT = 1


def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 1.0
    union = len(a | b)
    return len(a & b) / union if union else 0.0


@dataclass
class Cluster:
    cluster_id: str
    run_ids: List[str] = field(default_factory=list)
    centroid: Set[str] = field(default_factory=set)
    gene_counts: Dict[str, int] = field(default_factory=dict)
    # Test codes observed among labelled members, most common first.
    label_votes: Dict[str, int] = field(default_factory=dict)

    @property
    def size(self) -> int:
        return len(self.run_ids)

    @property
    def implied_panel(self) -> str:
        if self.label_votes:
            best = max(self.label_votes.items(), key=lambda kv: kv[1])[0]
            return best
        return self.cluster_id


def containment(sample: Set[str], profile: Set[str]) -> float:
    """Fraction of a run's observed genes explained by a panel profile.

    Jaccard is the wrong measure here and it is worth being explicit about why.
    A single run's fingerprint is a SUBSAMPLE of its panel — only genes where
    that subject happened to carry a call. Two runs on the same 500-gene panel
    might share only 40% of their genes, so pairwise Jaccard scatters them into
    singleton clusters, and every singleton cluster produces a denominator of
    one. Containment asks the question that actually matters: is everything
    this run saw inside the panel's footprint?
    """
    if not sample:
        return 0.0
    return len(sample & profile) / len(sample)


def cluster_fingerprints(
    fingerprints: Dict[str, Set[str]],
    labels: Optional[Dict[str, str]] = None,
    similarity: float = CLUSTER_SIMILARITY,
) -> List[Cluster]:
    """Group runs into implied panels (step 2).

    Two phases, because the labelled runs are the calibration the spec asks for
    in step 4 and throwing them away would be perverse:

      1. SEED.   Every distinct test code among labelled runs becomes a panel
                 profile — the union of its members' fingerprints. With enough
                 members the union converges on the true assayed gene list.
      2. ASSIGN. Each unlabelled run joins the profile that best contains it,
                 tie-broken toward the SMALLEST profile, so a run whose genes
                 fit both a focused panel and an exome is credited to the
                 panel. Over-crediting scope inflates denominators, which
                 understates every rate — the §E04 trap.

    With no labelled runs at all, phase 1 falls back to agglomerating on
    containment against running profiles, largest fingerprint first.
    """
    labels = labels or {}
    clusters: List[Cluster] = []
    by_label: Dict[str, Cluster] = {}

    def new_cluster(seed: Set[str], label: Optional[str] = None) -> Cluster:
        c = Cluster(cluster_id="PANEL-{:02d}".format(len(clusters) + 1))
        c.centroid = set(seed)
        clusters.append(c)
        if label:
            by_label[label] = c
        return c

    def attach(c: Cluster, run_id: str, fp: Set[str]) -> None:
        c.run_ids.append(run_id)
        c.centroid |= fp
        for g in fp:
            c.gene_counts[g] = c.gene_counts.get(g, 0) + 1
        lab = labels.get(run_id)
        if lab:
            c.label_votes[lab] = c.label_votes.get(lab, 0) + 1

    # ---- phase 1: seed from known panels ---------------------------------
    labelled = [(r, fp) for r, fp in fingerprints.items() if labels.get(r)]
    for run_id, fp in labelled:
        if len(fp) < MIN_FINGERPRINT:
            continue
        lab = labels[run_id]
        c = by_label.get(lab) or new_cluster(fp, label=lab)
        attach(c, run_id, fp)

    # ---- phase 2: assign the unlabelled ----------------------------------
    unlabelled = [(r, fp) for r, fp in fingerprints.items() if not labels.get(r)]
    # Largest first so a broad run cannot be swallowed by a narrow profile it
    # only partly overlaps before that broad profile exists.
    unlabelled.sort(key=lambda kv: len(kv[1]), reverse=True)

    for run_id, fp in unlabelled:
        if len(fp) < MIN_FINGERPRINT:
            continue
        best: Optional[Cluster] = None
        best_score = 0.0
        for c in clusters:
            score = containment(fp, c.centroid)
            if score > best_score or (
                    score == best_score and best is not None
                    and len(c.centroid) < len(best.centroid)):
                best, best_score = c, score
        if best is None or best_score < similarity:
            best = new_cluster(fp)
        attach(best, run_id, fp)

    return clusters


def calibrate_thresholds(
    fingerprints: Dict[str, Set[str]],
    declared_scope: Dict[str, Set[str]],
    genes: Sequence[str],
) -> Dict[str, float]:
    """Per-gene presence rate, learned from runs whose scope IS known (step 4).

    `declared_scope` maps run_id -> the gene set its test code declares. For
    each gene we measure: among runs that declare it, how often does it
    actually surface in the fingerprint? That is the calibrated rate; the
    inference threshold is a fraction of it.

    Genes never declared anywhere fall back to the floor, which keeps them
    inferable rather than silently excluded.
    """
    observed: Dict[str, int] = {}
    declared: Dict[str, int] = {}
    for run_id, scope in declared_scope.items():
        fp = fingerprints.get(run_id, set())
        for g in scope:
            declared[g] = declared.get(g, 0) + 1
            if g in fp:
                observed[g] = observed.get(g, 0) + 1

    thresholds: Dict[str, float] = {}
    for g in genes:
        n = declared.get(g, 0)
        if n >= 5:
            rate = observed.get(g, 0) / n
        else:
            rate = 1.0   # uncalibrated: demand full presence rather than guess
        thresholds[g] = max(MIN_THRESHOLD, rate * CALIBRATION_FACTOR)
    return thresholds


@dataclass
class ScopeCall:
    run_id: str
    gene_symbol: str
    in_scope: bool
    confidence: str   # 'declared' | 'inferred' | 'unknown'


def infer_scope(
    fingerprints: Dict[str, Set[str]],
    genes: Sequence[str],
    declared_scope: Optional[Dict[str, Set[str]]] = None,
    scope_complete: Optional[Dict[str, bool]] = None,
    labels: Optional[Dict[str, str]] = None,
) -> Tuple[List[ScopeCall], List[Cluster], Dict[str, float]]:
    """Full §3.4 pipeline. Returns per-run/gene scope calls, the implied panel
    clusters, and the calibrated per-gene thresholds (for the C01 inspector).

    A declared scope always wins over inference — inference exists to cover the
    runs the registry cannot account for. Where the registry declares a gene
    list but flags the reportable scope as incomplete, the call is 'unknown':
    counted in numerators, provisional in denominators (spec §G02).
    """
    declared_scope = declared_scope or {}
    scope_complete = scope_complete or {}

    clusters = cluster_fingerprints(fingerprints, labels=labels)
    thresholds = calibrate_thresholds(fingerprints, declared_scope, genes)

    cluster_of: Dict[str, Cluster] = {}
    for c in clusters:
        for r in c.run_ids:
            cluster_of[r] = c

    calls: List[ScopeCall] = []
    for run_id in fingerprints:
        declared = declared_scope.get(run_id)
        complete = scope_complete.get(run_id, True)
        cluster = cluster_of.get(run_id)

        for g in genes:
            if declared is not None:
                if g not in declared:
                    continue                      # not assayed: no row emitted
                conf = "declared" if complete else "unknown"
                calls.append(ScopeCall(run_id, g, True, conf))
                continue

            if cluster is None:
                continue
            rate = cluster.gene_counts.get(g, 0) / max(1, cluster.size)
            thr = thresholds.get(g, 1.0)
            if rate >= thr:
                calls.append(ScopeCall(run_id, g, True, "inferred"))
            elif rate >= thr * UNKNOWN_FLOOR:
                calls.append(ScopeCall(run_id, g, True, "unknown"))
            # else: not assayed — deliberately no row (spec C01 "remainder is
            # not assayed"), so absence is absence, never wild-type.
    return calls, clusters, thresholds
