"""PLINK 1 binary fileset reader — .bed/.bim/.fam (Part II §2.1).

Array data is the most common research upload, and the binary fileset is its
lingua franca. Two details in here are worth more care than their line count
suggests, because both fail *silently*:

1.  **Allele orientation.** In PLINK 1 the .bim columns are A1 then A2, and A1
    is conventionally the minor / effect allele — it is the allele the .bed
    encoding counts. Our dosage convention counts ALT, so the counted allele
    must land in `Variant.alt`: **ref = A2, alt = A1**. Getting this backwards
    produces a matrix where every dosage is 2-d, which passes every shape and
    range check and flips the sign of every effect estimate.

2.  **Bit order within a byte.** Samples are packed two bits each starting from
    the *low* bits of each byte, so the first sample is bits 0-1, not 7-6.
    Reading it high-first shuffles samples within each group of four while
    keeping allele frequencies identical, which means QC will not catch it.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..types import MISSING, GenotypeMatrix, Variant

_BED_MAGIC = b"\x6c\x1b"

# PLINK 1 two-bit genotype codes -> our ALT (= A1) dosage.
#   00 -> homozygous A1      -> 2
#   01 -> missing            -> -1
#   10 -> heterozygous       -> 1
#   11 -> homozygous A2      -> 0
# Note the deliberate irregularity: 01 is the missing code, not 1 copy. The
# codes are not a count and must not be treated as one.
_CODE_TO_DOSAGE = (2, MISSING, 1, 0)

# PLINK stores the sex chromosomes as integer codes. Harmonising them here
# (§2.3 chromosome naming) keeps a PLINK fileset joinable with a VCF, which
# spells them X / Y / MT.
_PLINK_CHROM_CODES = {"23": "X", "24": "Y", "25": "XY", "26": "MT"}


def read_plink1(prefix) -> GenotypeMatrix:
    """Read a PLINK 1 binary fileset.

    `prefix` may be the shared stem (`/data/cohort`) or the .bed path itself;
    users supply both and neither is wrong.
    """
    bed, bim, fam = resolve_prefix(prefix)

    samples = read_fam(fam)
    sample_ids = [s["iid"] for s in samples]
    variants = _read_bim(bim)

    dosages = _read_bed(bed, n_samples=len(sample_ids), n_variants=len(variants))

    return GenotypeMatrix(
        sample_ids=sample_ids,
        variants=variants,
        dosages=dosages,
        build=None,  # a PLINK fileset carries no build; §2.3 forbids guessing
        source_format="PLINK1",
        source_files=[str(bed), str(bim), str(fam)],
    )


def resolve_prefix(prefix) -> Tuple[Path, Path, Path]:
    p = Path(prefix)
    if p.suffix.lower() == ".bed":
        p = p.with_suffix("")
    bed, bim, fam = (Path(str(p) + ext) for ext in (".bed", ".bim", ".fam"))
    missing = [str(f) for f in (bed, bim, fam) if not f.exists()]
    if missing:
        raise ValueError(
            "incomplete PLINK 1 fileset for prefix '{}': missing {}".format(
                p, ", ".join(missing)))
    return bed, bim, fam


# ---------------------------------------------------------------------- fam --
def read_fam(path) -> List[Dict[str, Any]]:
    """Parse a .fam into dicts: fid, iid, pat, mat, sex, pheno.

    Sex and phenotype are carried out of ingestion because §2.3 requires a
    reported-vs-genetic sex concordance check, and the reported value only
    exists here.
    """
    out: List[Dict[str, Any]] = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 6:
                raise ValueError(
                    "{}:{}: expected 6 whitespace-separated columns "
                    "(FID IID PAT MAT SEX PHENO), got {}".format(path, lineno, len(fields)))
            out.append({
                "fid": fields[0],
                "iid": fields[1],
                "pat": fields[2],
                "mat": fields[3],
                "sex": _fam_sex(fields[4]),
                "pheno": _fam_pheno(fields[5]),
            })
    if not out:
        raise ValueError("{}: empty .fam — no samples".format(path))
    return out


def _fam_sex(raw: str) -> Optional[int]:
    """1 = male, 2 = female; 0 / -9 / anything else = unknown (None)."""
    try:
        code = int(raw)
    except ValueError:
        return None
    return code if code in (1, 2) else None


def _fam_pheno(raw: str) -> Optional[float]:
    # -9 and NA are PLINK's missing markers. Case/control coding (1/2) is left
    # as-is rather than recoded to 0/1: the .fam alone cannot distinguish a
    # 1/2 case-control column from a quantitative trait that happens to take
    # those values, so the caller decides.
    if raw in ("NA", "na", "-9", "."):
        return None
    try:
        return float(raw)
    except ValueError:
        return None


# ---------------------------------------------------------------------- bim --
def _read_bim(path) -> List[Variant]:
    variants: List[Variant] = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, start=1):
            f = line.split()
            if not f:
                continue
            if len(f) < 6:
                raise ValueError(
                    "{}:{}: expected 6 columns (CHROM VID CM POS A1 A2), "
                    "got {}".format(path, lineno, len(f)))
            chrom, vid, _cm, pos, a1, a2 = f[0], f[1], f[2], f[3], f[4], f[5]
            variants.append(Variant(
                chrom=_normalise_plink_chrom(chrom),
                pos=int(pos),
                ref=a2,   # A2 = non-counted allele
                alt=a1,   # A1 = counted allele; see module docstring
                vid=None if vid in (".", "") else vid,
            ))
    if not variants:
        raise ValueError("{}: empty .bim — no variants".format(path))
    return variants


def _normalise_plink_chrom(chrom: str) -> str:
    c = chrom.strip()
    if c[:3].lower() == "chr":
        c = c[3:]
    return _PLINK_CHROM_CODES.get(c, c)


# ---------------------------------------------------------------------- bed --
def _dosage_lookup() -> np.ndarray:
    """256 x 4 table: byte value -> dosages of the four samples it packs.

    A table beats bit arithmetic per genotype because decoding then becomes a
    single fancy-index over the whole file; a biobank .bed has 10^10 genotypes
    and a Python-level loop over them does not finish.
    """
    table = np.empty((256, 4), dtype=np.int8)
    for byte in range(256):
        for slot in range(4):
            table[byte, slot] = _CODE_TO_DOSAGE[(byte >> (2 * slot)) & 0b11]
    return table


_LUT = _dosage_lookup()


def _read_bed(path, n_samples: int, n_variants: int) -> np.ndarray:
    path = Path(path)
    with path.open("rb") as fh:
        header = fh.read(3)
        if len(header) < 3 or header[:2] != _BED_MAGIC:
            raise ValueError(
                "{}: missing PLINK .bed magic bytes 0x6c 0x1b — not a PLINK 1 "
                "binary genotype file".format(path))
        mode = header[2]
        if mode == 0x00:
            raise ValueError(
                "{}: individual-major .bed (mode byte 0x00). Convert it first: "
                "plink --bfile <prefix> --make-bed --out <new_prefix>".format(path))
        if mode != 0x01:
            raise ValueError(
                "{}: unrecognised .bed mode byte 0x{:02x} (expected 0x01 "
                "SNP-major)".format(path, mode))
        raw = np.frombuffer(fh.read(), dtype=np.uint8)

    # SNP-major: one contiguous run of ceil(n_samples/4) bytes per variant, the
    # last byte of each run zero-padded when n_samples is not a multiple of 4.
    # The padding is real data to numpy, so it is sliced off after decoding.
    bytes_per_variant = (n_samples + 3) // 4
    expected = bytes_per_variant * n_variants
    if raw.size != expected:
        raise ValueError(
            "{}: genotype block is {} bytes, expected {} for {} variants x {} "
            "samples — .bed does not match its .bim/.fam".format(
                path, raw.size, expected, n_variants, n_samples))
    if n_variants == 0:
        return np.empty((0, n_samples), dtype=np.int8)

    block = raw.reshape(n_variants, bytes_per_variant)
    return _LUT[block].reshape(n_variants, bytes_per_variant * 4)[:, :n_samples]
