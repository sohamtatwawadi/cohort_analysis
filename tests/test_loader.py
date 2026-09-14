"""VariMAT loader — spec §3.2. Dedup is the highest-risk item in the build."""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

from backend.app.ingest import varimat
from backend.app.ingest.varimat import (ColumnMap, normalise_acmg,
                                        normalise_consequence, normalise_var_class,
                                        parse_file, parse_filename)

HEADER = ["CHROM", "START", "REF", "ALT", "GENE_ID", "GENE_NAME", "VARCLASS",
          "VARTYPE", "AA_CHG", "AA_POS", "ZYGOSITY", "ALT_ALLELE_PERCENTAGE",
          "OVERALL_READ_DEPTH", "VARIANT_FILTER_STATUS", "VARIANT_LOCATION",
          "autoACMGPrediction", "gnomAD_AF", "MANE"]


def _write(path: Path, rows):
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({**{c: "NA" for c in HEADER}, **r})


def _variant(pos, gene="BRCA1", **kw):
    base = {"CHROM": "17", "START": pos, "REF": "A", "ALT": "G",
            "GENE_ID": "HGNC:1100", "GENE_NAME": gene, "VARCLASS": "MISSENSE",
            "VARTYPE": "SNV", "AA_POS": 100, "ZYGOSITY": "Heterozygous",
            "ALT_ALLELE_PERCENTAGE": "48.20", "OVERALL_READ_DEPTH": "120",
            "VARIANT_FILTER_STATUS": "PASS", "VARIANT_LOCATION": "ONTARGET",
            "autoACMGPrediction": "Pathogenic", "gnomAD_AF": "0.00001",
            "MANE": "MANE_Select"}
    base.update(kw)
    return base


# ------------------------------------------------------------------- dedup ---
def test_multi_transcript_rows_collapse_to_one_variant(tmp_path):
    """The 1.75x trap. Same locus annotated under three transcripts and two
    gene symbols is ONE variant, not three."""
    f = tmp_path / "S1.varimat.tsv"
    _write(f, [
        _variant(43000000, "BRCA1", MANE="MANE_Select"),
        _variant(43000000, "BRCA1", MANE="NA", GENE_ID="NA"),
        _variant(43000000, "RNF53", MANE="NA", GENE_ID="NA"),   # alias symbol
        _variant(43000500, "BRCA1"),
    ])
    s = parse_file(f)
    assert s.funnel.raw_rows == 4
    assert s.funnel.distinct_variants == 2
    assert len(s.variants) == 2


def test_canonical_row_selected_by_mane(tmp_path):
    """When a locus has several annotations, the MANE row's interpretation is
    the one that survives — not whichever happened to be read first."""
    f = tmp_path / "S2.varimat.tsv"
    _write(f, [
        _variant(43000000, VARCLASS="INTRONIC", MANE="NA", GENE_ID="NA",
                 autoACMGPrediction="Benign"),
        _variant(43000000, VARCLASS="NONSENSE", MANE="MANE_Select",
                 autoACMGPrediction="Pathogenic"),
    ])
    s = parse_file(f)
    assert len(s.variants) == 1
    assert s.variants[0].consequence == "Nonsense"
    assert s.variants[0].acmg == "Pathogenic"
    assert s.variants[0].mane is True


def test_rows_per_variant_ratio_is_reported(tmp_path):
    f = tmp_path / "S3.varimat.tsv"
    _write(f, [_variant(43000000)] * 2 + [_variant(43001000)] * 2)
    s = parse_file(f)
    assert s.funnel.rows_per_variant == pytest.approx(2.0)


# -------------------------------------------------------- reviewable subset --
def test_genome_wide_vus_noise_is_excluded(tmp_path):
    """Spec §3.2: autoACMGPrediction marks 309,331 variants VUS in a single
    sample because it runs genome-wide. Intronic VUS must not reach the store."""
    f = tmp_path / "S4.varimat.tsv"
    rows = [_variant(1_000_000 + i, VARCLASS="INTRONIC",
                     autoACMGPrediction="Uncertain significance",
                     gnomAD_AF="0.3") for i in range(500)]
    rows.append(_variant(43000000, VARCLASS="MISSENSE"))
    _write(f, rows)
    s = parse_file(f)
    assert s.funnel.distinct_variants == 501
    assert len(s.variants) == 1, "intronic genome-wide VUS leaked into the store"


def test_non_pass_and_common_variants_are_dropped(tmp_path):
    f = tmp_path / "S5.varimat.tsv"
    _write(f, [
        _variant(43000000, VARIANT_FILTER_STATUS="LowQual"),
        _variant(43001000, gnomAD_AF="0.25"),          # common
        _variant(43002000),                            # keeper
    ])
    s = parse_file(f)
    assert s.funnel.pass_filter == 2
    assert s.funnel.rare == 1
    assert len(s.variants) == 1


def test_missing_gnomad_af_counts_as_rare(tmp_path):
    """Spec: (gnomAD_AF = NA OR gnomAD_AF < 0.01). Absent is not common."""
    f = tmp_path / "S6.varimat.tsv"
    _write(f, [_variant(43000000, gnomAD_AF="NA")])
    s = parse_file(f)
    assert len(s.variants) == 1


# --------------------------------------------------------------- fingerprint --
def test_fingerprint_includes_genes_with_only_non_coding_calls(tmp_path):
    """An ONTARGET call of ANY kind is evidence the assay looked at the gene.
    Fingerprinting must happen before consequence filtering, or the denominator
    shrinks to only genes that happened to carry a coding variant."""
    f = tmp_path / "S7.varimat.tsv"
    _write(f, [
        _variant(43000000, "BRCA1", VARCLASS="INTRONIC", gnomAD_AF="0.4"),
        _variant(32000000, "BRCA2", GENE_ID="HGNC:1101", VARCLASS="MISSENSE"),
    ])
    s = parse_file(f)
    assert s.fingerprint == {"BRCA1", "BRCA2"}
    assert len(s.variants) == 1


def test_off_target_calls_do_not_enter_the_fingerprint(tmp_path):
    f = tmp_path / "S8.varimat.tsv"
    _write(f, [_variant(43000000, "BRCA1", VARIANT_LOCATION="OFFTARGET")])
    s = parse_file(f)
    assert s.fingerprint == set()


# ------------------------------------------------------------ field mapping --
def test_vaf_is_converted_from_percent_to_fraction(tmp_path):
    """VariMAT stores ALT_ALLELE_PERCENTAGE as a percent, e.g. 47.73."""
    f = tmp_path / "S9.varimat.tsv"
    _write(f, [_variant(43000000, ALT_ALLELE_PERCENTAGE="47.73")])
    s = parse_file(f)
    assert s.variants[0].vaf == pytest.approx(0.4773)


def test_gene_id_preferred_over_symbol_and_aliases_resolve():
    from backend.app.reference.genes import resolve_symbol
    assert resolve_symbol("HGNC:1100", "WRONG") == "BRCA1"     # GENE_ID wins
    assert resolve_symbol(None, "BRCA1") == "BRCA1"
    assert resolve_symbol(None, "GBA1") == "GBA"               # alias map
    assert resolve_symbol(None, "NOT_A_GENE") is None


@pytest.mark.parametrize("varclass,expected", [
    ("MISSENSE", "Missense"), ("NONSENSE", "Nonsense"),
    ("FRAMESHIFT-DEL", "Frameshift"), ("FRAMESHIFT-INS", "Frameshift"),
    ("INFRAME-INS", "In-frame indel"), ("ACCEPTOR-SS-VARIANT", "Splice site"),
    ("DONOR-SS-LOSS", "Splice site"), ("SILENT", "Synonymous"),
    ("5UTR", "5' UTR"), ("EXONIC-NC", "Non-coding exonic"),
    ("INTRONIC", None), ("INTERGENIC", None),
])
def test_consequence_vocabulary(varclass, expected):
    assert normalise_consequence(varclass) == expected


def test_filename_metadata_parsing():
    m = parse_filename(Path("MG01001.LINC-v2.9.0.DRAGEN-4.2.GRCh38.varimat.tsv"))
    assert m["sample_id"] == "MG01001"
    assert m["reference_build"] == "GRCh38"
    assert m["caller"] == "DRAGEN-4.2"


def test_non_varimat_file_is_rejected_clearly(tmp_path):
    f = tmp_path / "random.tsv"
    f.write_text("a\tb\tc\n1\t2\t3\n")
    with pytest.raises(ValueError, match="missing required columns"):
        parse_file(f)


def test_empty_file_is_rejected(tmp_path):
    f = tmp_path / "empty.tsv"
    f.write_text("")
    with pytest.raises(ValueError, match="empty"):
        parse_file(f)


# ------------------------------------------------------------ multi-sample ---
def test_duplicate_sample_ids_fail_loudly(tmp_path):
    """A cohort is many files. Two files for one sample is an operator decision,
    not something to resolve silently."""
    from backend.app.ingest import store

    a = tmp_path / "MG1.LINC-v2.9.0.DRAGEN-4.2.GRCh38.varimat.tsv"
    b = tmp_path / "MG1.LINC-v2.9.1.GATK-4.5.GRCh37.varimat.tsv"
    _write(a, [_variant(43000000)])
    _write(b, [_variant(43000000)])
    with pytest.raises(ValueError, match="map to more than one VariMAT file"):
        store.ingest_varimat([a, b])
