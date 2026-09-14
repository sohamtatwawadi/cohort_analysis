"""Format detection and dispatch — the first step of the §2.3 pipeline.

Detection sniffs *content*, not the extension. Research uploads arrive with
`.txt`, `.dat`, no suffix at all, or a suffix that lies; deciding what a file
is from its name is how a PLINK fileset gets parsed as TSV and silently yields
zero variants.
"""
from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import List, Optional

from ..types import GenotypeMatrix
from . import plink, vcf

VCF = "VCF"
PLINK1 = "PLINK1"
VARIMAT = "VariMAT"
UNKNOWN = "unknown"

# A VariMAT is identified by its required analytical columns rather than by any
# single marker column, because the real files carry ~230 columns and vary.
# VARCLASS is what separates it from a generic variant TSV.
_VARIMAT_REQUIRED = ({"CHROM"}, {"START", "POS"}, {"REF"}, {"ALT"}, {"VARCLASS"})


def detect_format(path) -> str:
    """Return one of VCF / PLINK1 / VariMAT / unknown."""
    p = Path(path)

    # A bare fileset prefix ("/data/cohort") is a legitimate thing to hand us,
    # and it names no file at all.
    if not p.exists():
        if Path(str(p) + ".bed").exists() and Path(str(p) + ".bim").exists():
            return PLINK1
        return UNKNOWN

    if p.is_dir():
        return UNKNOWN

    if p.suffix.lower() in (".bed", ".bim", ".fam"):
        stem = p.with_suffix("")
        if Path(str(stem) + ".bim").exists() or Path(str(stem) + ".fam").exists():
            return PLINK1

    head = _read_head_bytes(p, 2)
    if head == b"\x6c\x1b":
        return PLINK1

    first = _first_nonempty_line(p)
    if first is None:
        return UNKNOWN
    if first.startswith("##fileformat=VCF"):
        return VCF

    columns = {c.strip().upper() for c in _split_header(first)}
    if all(alternatives & columns for alternatives in _VARIMAT_REQUIRED):
        return VARIMAT

    return UNKNOWN


def load_any(path) -> GenotypeMatrix:
    """Detect and load. VariMAT is deliberately not handled here."""
    fmt = detect_format(path)
    if fmt == VCF:
        return vcf.read_vcf(path)
    if fmt == PLINK1:
        return plink.read_plink1(path)
    if fmt == VARIMAT:
        # VariMAT is one file per sample and carries annotation, classification
        # and coverage fingerprints that a GenotypeMatrix has nowhere to put.
        # backend.app.ingest.varimat owns it.
        raise NotImplementedError(
            "{} is VariMAT — load it with backend.app.ingest.varimat, which "
            "handles the per-sample fileset and annotation columns".format(path))
    raise ValueError("{}: unrecognised genomic format".format(path))


# ------------------------------------------------------------------ sniffing --
def _read_head_bytes(path: Path, n: int) -> bytes:
    try:
        with path.open("rb") as fh:
            return fh.read(n)
    except OSError:
        return b""


def _first_nonempty_line(path: Path) -> Optional[str]:
    try:
        with _open_text(path) as fh:
            for line in fh:
                if line.strip():
                    return line.rstrip("\n").rstrip("\r")
    except (OSError, EOFError, UnicodeDecodeError, gzip.BadGzipFile):
        return None
    return None


def _split_header(line: str) -> List[str]:
    delim = "," if line.count(",") > line.count("\t") else "\t"
    return line.split(delim)


def _open_text(path: Path) -> io.TextIOBase:
    if _read_head_bytes(path, 2) == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(str(path), "rb"), encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")
