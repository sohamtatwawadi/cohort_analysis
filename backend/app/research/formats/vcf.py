"""VCF / VCF.gz reader (Part II §2.1, §2.3).

Produces a `GenotypeMatrix` in the dosage convention fixed by `types.py`:
dosage = count of ALT alleles, MISSING (-1) for no-call.

Two §2.3 mandatory validations are implemented here rather than downstream,
because doing them later means every consumer has to guess:

*   **Multi-allelic split.** A record with `ALT=A,T` becomes two rows, one per
    ALT allele. Unsplit multi-allelics have no well-defined ALT dosage — a 1/2
    call is neither 1 nor 2 copies of "the" ALT — so allele frequency, burden
    and association would all silently use a number that means nothing.
*   **Chromosome naming.** `chr1` and `1` are the same contig; a leading `chr`
    is stripped so a GRCh38 VCF and a PLINK fileset join on the same key.

Build is read from the header or left as None. It is never inferred from the
coordinate range: GRCh37 and GRCh38 positions occupy the same intervals, so a
guess is a coin flip that mislabels every annotation downstream.
"""
from __future__ import annotations

import gzip
import io
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from ..types import MISSING, GenotypeMatrix, Variant

# Build aliases as they appear in ##reference paths and ##contig assembly
# fields. Matched with alphanumeric boundaries rather than plain substring, so
# ".../GRCh38_full_analysis_set.fa" resolves but "b37" inside an accession does
# not. Only ##reference and ##contig are ever consulted.
_BUILD_ALIASES = (
    ("GRCh38", ("grch38", "grch38p13", "hg38", "b38", "hs38dh", "hs38d1")),
    ("GRCh37", ("grch37", "grch37p13", "hg19", "b37", "hs37d5", "human_g1k_v37")),
    ("T2T-CHM13", ("chm13", "t2t", "hs1")),
)

_GT_SPLIT = re.compile(r"[/|]")


# Rows converted to int8 per block. Small enough that the transient Python
# representation stays well under a hundred megabytes at 2,500 samples.
VCF_BLOCK_ROWS = 4096


def read_vcf(path, max_variants: Optional[int] = None) -> GenotypeMatrix:
    """Read a VCF (plain or gzip/bgzip) into a GenotypeMatrix.

    `max_variants` caps the number of *records read from the file*, applied
    before multi-allelic splitting — so a file whose records are bi-allelic
    yields exactly that many rows, and one with multi-allelic sites may yield
    more. Capping the input rather than the output keeps the truncation point
    reproducible and never leaves a site half-split.
    """
    path = Path(path)
    sample_ids: List[str] = []
    variants: List[Variant] = []
    build: Optional[str] = None
    seen_chrom_line = False
    n_records = 0

    # Accumulate in int8 blocks rather than one growing list of lists.
    #
    # A list of N lists of M Python ints costs ~8 bytes per pointer plus the
    # list headers, so 150,000 x 2,504 — a thinned single chromosome of 1000
    # Genomes — is about 3 GB before it becomes a 358 MB int8 array. That is
    # the whole budget of a small instance spent on a representation that gets
    # thrown away. Converting every few thousand rows keeps the transient cost
    # to one block.
    rows: List[List[int]] = []
    blocks: List[np.ndarray] = []

    def flush() -> None:
        if rows:
            blocks.append(np.asarray(rows, dtype=np.int8))
            del rows[:]

    with _open_text(path) as fh:
        for line in fh:
            if not line or line == "\n":
                continue
            if line.startswith("##"):
                build = build or _build_from_meta(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").rstrip("\r").split("\t")
                sample_ids = [s.strip() for s in fields[9:]]
                seen_chrom_line = True
                continue
            if line.startswith("#"):
                continue
            if not seen_chrom_line:
                raise ValueError(
                    "{}: data line before the #CHROM header — not a valid VCF".format(path))
            if max_variants is not None and n_records >= max_variants:
                break
            n_records += 1
            _parse_record(line, len(sample_ids), variants, rows)
            if len(rows) >= VCF_BLOCK_ROWS:
                flush()
    flush()

    if not seen_chrom_line:
        raise ValueError(
            "{}: no #CHROM header line found (file empty or not a VCF)".format(path))

    if blocks:
        dosages = blocks[0] if len(blocks) == 1 else np.concatenate(blocks, axis=0)
        del blocks[:]
    else:
        dosages = np.empty((0, len(sample_ids)), dtype=np.int8)

    return GenotypeMatrix(
        sample_ids=sample_ids,
        variants=variants,
        dosages=dosages,
        build=build,
        source_format="VCF",
        source_files=[str(path)],
    )


# ------------------------------------------------------------------ records --
def _parse_record(line: str, n_samples: int,
                  variants: List[Variant], rows: List[List[int]]) -> None:
    f = line.rstrip("\n").rstrip("\r").split("\t")
    if len(f) < 8:
        return  # truncated line; nothing usable

    chrom = normalise_chrom(f[0])
    pos = int(f[1])
    vid = None if f[2] in (".", "") else f[2]
    ref = f[3]
    alts = [a for a in f[4].split(",") if a != ""]
    qual = _maybe_float(f[5])
    filt = None if f[6] in (".", "") else f[6]
    info = _parse_info(f[7])

    if not alts or alts == ["."]:
        return  # non-variant / reference block record has no ALT to count

    # Per-sample allele calls, parsed once and reused for every ALT allele of
    # this site: alleles[s] is a list of ints, or None if any allele is a
    # no-call. Re-parsing per split would be O(n_alt) passes for no benefit.
    if n_samples:
        gt_idx = _gt_index(f[8]) if len(f) > 8 else None
        calls = [_parse_gt(f[9 + s] if 9 + s < len(f) else ".", gt_idx)
                 for s in range(n_samples)]
    else:
        calls = []

    multiallelic = len(alts) > 1
    for k, alt in enumerate(alts, start=1):
        variants.append(Variant(
            chrom=chrom, pos=pos, ref=ref, alt=alt, vid=vid, qual=qual,
            filter_status=filt,
            # Split rows get their own INFO dict so that a later per-allele
            # edit (AC, AF) on one row cannot reach through to its siblings.
            info=dict(info) if (multiallelic and info) else info))
        # Dosage for this split counts allele k only. GT 1/2 is one copy of
        # allele 1 and one of allele 2, so it is dosage 1 on *both* rows;
        # GT 2/2 is dosage 0 on the allele-1 row and 2 on the allele-2 row.
        rows.append([MISSING if c is None else c.count(k) for c in calls])


def _gt_index(fmt: str) -> Optional[int]:
    """Position of GT within a FORMAT string.

    The VCF spec requires GT first *when present*, but files written by
    filtering and merging tools routinely violate it, so the index is looked up
    rather than assumed to be 0.
    """
    keys = fmt.split(":")
    try:
        return keys.index("GT")
    except ValueError:
        return None


def _parse_gt(field: str, gt_idx: Optional[int]) -> Optional[List[int]]:
    """Return allele indices, or None if the call is missing.

    Handles `0/1`, `0|1`, `./.`, `.`, and haploid `0` / `1`. A partial no-call
    (`./1`) is treated as missing outright: half a genotype cannot be turned
    into a dosage without inventing the other allele.
    """
    if gt_idx is None:
        return None
    if field in (".", "", "./.", ".|."):
        return None
    parts = field.split(":")
    if gt_idx >= len(parts):
        return None
    gt = parts[gt_idx]
    if gt in (".", ""):
        return None
    out: List[int] = []
    for a in _GT_SPLIT.split(gt):
        if a == "." or a == "":
            return None
        try:
            out.append(int(a))
        except ValueError:
            return None
    return out or None


# ------------------------------------------------------------------- header --
def _build_from_meta(line: str) -> Optional[str]:
    """Extract a genome build from `##reference=` or `##contig=<...assembly=>`.

    Returns None when the header does not say. A None build is a correct,
    actionable answer (§2.3 fails the upload and asks the user); an inferred
    one is not.
    """
    if line.startswith("##reference="):
        return _match_build(line.split("=", 1)[1])
    if line.startswith("##contig=") and "assembly" in line:
        m = re.search(r"assembly=([^,>]+)", line)
        if m:
            return _match_build(m.group(1))
    return None


def _match_build(text: str) -> Optional[str]:
    low = text.lower()
    for build, aliases in _BUILD_ALIASES:
        for alias in aliases:
            if re.search(r"(?<![0-9a-z])" + re.escape(alias) + r"(?![0-9a-z])", low):
                return build
    return None


def normalise_chrom(chrom: str) -> str:
    """`chr1` / `CHR1` / `1` all become `1` (§2.3 naming harmonisation)."""
    c = chrom.strip()
    if c[:3].lower() == "chr":
        c = c[3:]
    return c


def _parse_info(raw: str) -> Optional[Dict[str, Any]]:
    if raw in (".", ""):
        return None
    out: Dict[str, Any] = {}
    for item in raw.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = True  # flag
    return out or None


def _maybe_float(raw: str) -> Optional[float]:
    if raw in (".", ""):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


def _open_text(path: Path) -> io.TextIOBase:
    # bgzip is a gzip member stream, so the stdlib reads it fine sequentially.
    # Random access by tabix index would need the BGZF block structure, which
    # we do not use — ingestion is a single forward pass.
    if path.suffix == ".gz" or path.suffix == ".bgz":
        return io.TextIOWrapper(gzip.open(str(path), "rb"), encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")
