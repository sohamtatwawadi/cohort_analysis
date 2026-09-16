"""Research Mode ingestion: VCF, PLINK 1 and format detection (Part II §2).

Every fixture is written inline from a genotype matrix the test already knows,
so a failure points at the reader and not at a checked-in file drifting.
"""
from __future__ import annotations

import gzip
from pathlib import Path
import numpy as np
import pytest

from backend.app.research.formats import detect, plink, vcf
from backend.app.research.types import MISSING, GenotypeMatrix

# ============================================================== VCF fixtures ==

VCF_HEADER = """##fileformat=VCFv4.2
##reference=file:///ref/GRCh38_full_analysis_set.fa
##contig=<ID=chr1,length=248956422,assembly=GRCh38>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
"""

# Rows exercise, in order: chr-prefixed contig + unphased/phased/no-call;
# bare contig + a bare "." genotype; a multi-allelic site whose FORMAT puts GT
# second; a haploid call on X.
VCF_BODY = """chr1\t100\trs1\tA\tG\t50\tPASS\tDP=10\tGT:DP\t0/1:10\t1|1:12\t./.:0
1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/0\t0|1\t.
chr2\t300\trs3\tG\tA,T\t99\tPASS\tAC=2,1\tDP:GT\t8:1/2\t9:2/2\t7:0/1
chrX\t400\t.\tA\tC\t.\tPASS\t.\tGT\t1\t0\t.
"""

# What the four records must become after the multi-allelic split.
EXPECTED_KEYS = ["1:100:A:G", "1:200:C:T", "2:300:G:A", "2:300:G:T", "X:400:A:C"]
EXPECTED_DOSAGES = np.array([
    [1, 2, MISSING],        # 0/1, 1|1, ./.
    [0, 1, MISSING],        # 0/0, 0|1, .
    [1, 0, 1],              # allele A (index 1) of 1/2, 2/2, 0/1
    [1, 2, 0],              # allele T (index 2) of 1/2, 2/2, 0/1
    [1, 0, MISSING],        # haploid 1, 0, .
], dtype=np.int8)


def write_vcf(tmp_path: Path, name: str = "cohort.vcf",
              header: str = VCF_HEADER, body: str = VCF_BODY) -> Path:
    p = tmp_path / name
    p.write_text(header + body)
    return p


# ================================================================= VCF tests ==

def test_vcf_samples_and_shape(tmp_path):
    gm = vcf.read_vcf(write_vcf(tmp_path))
    assert gm.sample_ids == ["S1", "S2", "S3"]
    assert gm.source_format == "VCF"
    assert gm.dosages.dtype == np.int8
    assert_invariants(gm)


def test_vcf_multiallelic_split_and_dosages(tmp_path):
    gm = vcf.read_vcf(write_vcf(tmp_path))
    assert [v.key for v in gm.variants] == EXPECTED_KEYS
    np.testing.assert_array_equal(gm.dosages, EXPECTED_DOSAGES)


def test_vcf_multiallelic_rows_share_site_but_not_dosage(tmp_path):
    """1/2 is one copy of each allele; 2/2 is zero copies of allele 1."""
    gm = vcf.read_vcf(write_vcf(tmp_path))
    a_row = gm.variants.index(next(v for v in gm.variants if v.key == "2:300:G:A"))
    t_row = gm.variants.index(next(v for v in gm.variants if v.key == "2:300:G:T"))
    assert gm.dosages[a_row][0] == 1 and gm.dosages[t_row][0] == 1   # 1/2
    assert gm.dosages[a_row][1] == 0 and gm.dosages[t_row][1] == 2   # 2/2


def test_vcf_gt_not_first_in_format(tmp_path):
    """FORMAT is DP:GT on the multi-allelic record — index 0 would read depths."""
    gm = vcf.read_vcf(write_vcf(tmp_path))
    assert set(np.unique(gm.dosages)) <= {MISSING, 0, 1, 2}


def test_vcf_chrom_harmonised(tmp_path):
    gm = vcf.read_vcf(write_vcf(tmp_path))
    assert [v.chrom for v in gm.variants] == ["1", "1", "2", "2", "X"]


def test_vcf_fields_parsed(tmp_path):
    gm = vcf.read_vcf(write_vcf(tmp_path))
    v0 = gm.variants[0]
    assert v0.vid == "rs1" and v0.qual == 50.0 and v0.filter_status == "PASS"
    assert v0.info == {"DP": "10"}
    v1 = gm.variants[1]
    assert v1.vid is None and v1.qual is None and v1.filter_status is None
    assert v1.info is None


def test_vcf_build_from_reference(tmp_path):
    gm = vcf.read_vcf(write_vcf(tmp_path))
    assert gm.build == "GRCh38"


def test_vcf_build_from_contig_assembly(tmp_path):
    header = ('##fileformat=VCFv4.2\n'
              '##contig=<ID=1,length=249250621,assembly=b37>\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')
    p = write_vcf(tmp_path, "b37.vcf", header, "1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\n")
    assert vcf.read_vcf(p).build == "GRCh37"


def test_vcf_build_none_when_header_silent(tmp_path):
    """§2.3: build is never inferred from coordinates — unknown stays unknown."""
    header = ('##fileformat=VCFv4.2\n'
              '##reference=file:///ref/mystery.fa\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')
    p = write_vcf(tmp_path, "nobuild.vcf", header, "1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\n")
    assert vcf.read_vcf(p).build is None


def test_vcf_gzip_round_trip(tmp_path):
    p = tmp_path / "cohort.vcf.gz"
    with gzip.open(str(p), "wt") as fh:
        fh.write(VCF_HEADER + VCF_BODY)
    gm = vcf.read_vcf(p)
    np.testing.assert_array_equal(gm.dosages, EXPECTED_DOSAGES)
    assert gm.sample_ids == ["S1", "S2", "S3"]


def test_vcf_max_variants_caps_records(tmp_path):
    gm = vcf.read_vcf(write_vcf(tmp_path), max_variants=2)
    assert [v.key for v in gm.variants] == EXPECTED_KEYS[:2]
    assert_invariants(gm)


def test_vcf_empty_file_raises(tmp_path):
    p = tmp_path / "empty.vcf"
    p.write_text("")
    with pytest.raises(ValueError, match="#CHROM"):
        vcf.read_vcf(p)


def test_vcf_without_chrom_line_raises(tmp_path):
    p = tmp_path / "headeronly.vcf"
    p.write_text("##fileformat=VCFv4.2\n##reference=hg19\n")
    with pytest.raises(ValueError, match="#CHROM"):
        vcf.read_vcf(p)


def test_vcf_sites_only(tmp_path):
    header = ('##fileformat=VCFv4.2\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    p = write_vcf(tmp_path, "sites.vcf", header, "chr7\t50\t.\tA\tG\t.\t.\t.\n")
    gm = vcf.read_vcf(p)
    assert gm.sample_ids == [] and gm.n_variants == 1
    assert gm.dosages.shape == (1, 0)


# ============================================================ PLINK fixtures ==

# Six samples (not a multiple of four, so each variant pads its last byte) and
# three variants, including every code: 2, 1, 0 and missing.
PLINK_DOSAGES = np.array([
    [2, 1, 0, MISSING, 2, 1],
    [0, 0, 1, 2, MISSING, 0],
    [MISSING, 2, 2, 1, 0, 1],
], dtype=np.int8)

PLINK_IIDS = ["S1", "S2", "S3", "S4", "S5", "S6"]

# Independent forward encoder: dosage of A1 -> PLINK two-bit code.
_DOSAGE_TO_CODE = {2: 0b00, MISSING: 0b01, 1: 0b10, 0: 0b11}


def pack_bed(dosages: np.ndarray, snp_major: bool = True) -> bytes:
    out = bytearray(b"\x6c\x1b" + (b"\x01" if snp_major else b"\x00"))
    n_variants, n_samples = dosages.shape
    for v in range(n_variants):
        byte = 0
        for s in range(n_samples):
            slot = s % 4
            byte |= _DOSAGE_TO_CODE[int(dosages[v, s])] << (2 * slot)
            if slot == 3:
                out.append(byte)
                byte = 0
        if n_samples % 4:
            out.append(byte)      # trailing partial byte, high slots left zero
    return bytes(out)


def write_plink(tmp_path: Path, prefix: str = "array",
                dosages: np.ndarray = PLINK_DOSAGES,
                snp_major: bool = True) -> Path:
    stem = tmp_path / prefix
    # A1 first then A2; the reader must map A1 -> alt because A1 is counted.
    stem.with_suffix(".bim").write_text(
        "1\trs1\t0\t1000\tA\tG\n"
        "chr2\trs2\t0\t2000\tT\tC\n"
        "23\trs3\t0\t3000\tG\tA\n")
    stem.with_suffix(".fam").write_text(
        "FAM1 S1 0 0 1 2\n"
        "FAM1 S2 0 0 2 1\n"
        "FAM2 S3 S1 S2 1 -9\n"
        "FAM3 S4 0 0 0 1\n"
        "FAM4 S5 0 0 2 2\n"
        "FAM5 S6 0 0 1 1\n")
    stem.with_suffix(".bed").write_bytes(pack_bed(dosages, snp_major=snp_major))
    return stem


# =============================================================== PLINK tests ==

def test_bed_bit_order_is_low_bits_first():
    """Pins the packing the reader assumes, independent of the reader itself.

    Samples 0..3 of variant 0 are dosages 2,1,0,missing -> codes 00,10,11,01,
    which occupy bits 0-1, 2-3, 4-5, 6-7 in that order -> 0x78.
    """
    assert pack_bed(PLINK_DOSAGES)[3] == 0x78


def test_plink_round_trip_exact_dosages(tmp_path):
    gm = plink.read_plink1(write_plink(tmp_path))
    np.testing.assert_array_equal(gm.dosages, PLINK_DOSAGES)
    assert gm.dosages.dtype == np.int8
    assert gm.source_format == "PLINK1"
    assert_invariants(gm)


def test_plink_padding_does_not_leak_extra_samples(tmp_path):
    """6 samples is 2 bytes per variant; the 2 padding slots decode to 0b00
    (= dosage 2) and would appear as phantom homozygotes if not sliced off."""
    gm = plink.read_plink1(write_plink(tmp_path))
    assert gm.dosages.shape == (3, 6)


def test_plink_accepts_bed_path_directly(tmp_path):
    stem = write_plink(tmp_path)
    gm = plink.read_plink1(str(stem) + ".bed")
    np.testing.assert_array_equal(gm.dosages, PLINK_DOSAGES)


def test_plink_ref_is_a2_and_alt_is_a1(tmp_path):
    gm = plink.read_plink1(write_plink(tmp_path))
    assert [(v.ref, v.alt) for v in gm.variants] == [("G", "A"), ("C", "T"), ("A", "G")]


def test_plink_counted_allele_matches_dosage_direction(tmp_path):
    """Sample S1 is homozygous A1 at rs1, so its dosage of alt must be 2."""
    gm = plink.read_plink1(write_plink(tmp_path))
    assert gm.variants[0].alt == "A"
    assert gm.dosages[0, gm.sample_index()["S1"]] == 2


def test_plink_chrom_harmonised(tmp_path):
    gm = plink.read_plink1(write_plink(tmp_path))
    assert [v.chrom for v in gm.variants] == ["1", "2", "X"]


def test_plink_build_is_none(tmp_path):
    """A PLINK fileset records no build, and §2.3 forbids inventing one."""
    assert plink.read_plink1(write_plink(tmp_path)).build is None


def test_read_fam_fields(tmp_path):
    rows = plink.read_fam(str(write_plink(tmp_path)) + ".fam")
    assert [r["iid"] for r in rows] == PLINK_IIDS
    assert rows[0] == {"fid": "FAM1", "iid": "S1", "pat": "0", "mat": "0",
                       "sex": 1, "pheno": 2.0}
    assert rows[2]["pat"] == "S1" and rows[2]["mat"] == "S2"
    assert rows[2]["pheno"] is None        # -9 is PLINK's missing phenotype
    assert rows[3]["sex"] is None          # 0 is unknown sex


def test_plink_individual_major_rejected(tmp_path):
    stem = write_plink(tmp_path, snp_major=False)
    with pytest.raises(ValueError, match="make-bed"):
        plink.read_plink1(stem)


def test_plink_bad_magic_rejected(tmp_path):
    stem = write_plink(tmp_path)
    stem.with_suffix(".bed").write_bytes(b"\x00\x00\x01" + b"\x00" * 6)
    with pytest.raises(ValueError, match="magic"):
        plink.read_plink1(stem)


def test_plink_size_mismatch_rejected(tmp_path):
    stem = write_plink(tmp_path)
    truncated = pack_bed(PLINK_DOSAGES)[:-1]
    stem.with_suffix(".bed").write_bytes(truncated)
    with pytest.raises(ValueError, match="does not match"):
        plink.read_plink1(stem)


def test_plink_missing_sibling_file(tmp_path):
    stem = write_plink(tmp_path)
    stem.with_suffix(".fam").unlink()
    with pytest.raises(ValueError, match="incomplete PLINK 1 fileset"):
        plink.read_plink1(stem)


def test_plink_allele_frequency_uses_a1(tmp_path):
    """AF must be the A1 frequency; an A2-oriented read would give 1 - af."""
    gm = plink.read_plink1(write_plink(tmp_path))
    af = gm.allele_frequency()
    # Variant 0 called dosages: 2,1,0,2,1 over 5 samples -> 6 / 10.
    assert af[0] == pytest.approx(0.6)


# ============================================================== detect tests ==

def test_detect_vcf(tmp_path):
    assert detect.detect_format(write_vcf(tmp_path)) == "VCF"


def test_detect_vcf_gz(tmp_path):
    p = tmp_path / "c.vcf.gz"
    with gzip.open(str(p), "wt") as fh:
        fh.write(VCF_HEADER + VCF_BODY)
    assert detect.detect_format(p) == "VCF"


def test_detect_vcf_ignores_extension(tmp_path):
    p = write_vcf(tmp_path, "calls.txt")
    assert detect.detect_format(p) == "VCF"


def test_detect_plink_by_prefix(tmp_path):
    stem = write_plink(tmp_path)
    assert detect.detect_format(stem) == "PLINK1"


def test_detect_plink_by_bed_path(tmp_path):
    stem = write_plink(tmp_path)
    assert detect.detect_format(str(stem) + ".bed") == "PLINK1"


def test_detect_varimat(tmp_path):
    p = tmp_path / "sample.varimat.tsv"
    p.write_text("CHROM\tSTART\tREF\tALT\tVARCLASS\tGENE_NAME\n"
                 "1\t100\tA\tG\tMISSENSE\tBRCA1\n")
    assert detect.detect_format(p) == "VariMAT"


def test_detect_unknown(tmp_path):
    p = tmp_path / "notes.txt"
    p.write_text("this is not a genomic file\n")
    assert detect.detect_format(p) == "unknown"


def test_detect_missing_path(tmp_path):
    assert detect.detect_format(tmp_path / "nope") == "unknown"


def test_load_any_dispatch(tmp_path):
    vcf_gm = detect.load_any(write_vcf(tmp_path))
    assert vcf_gm.source_format == "VCF"
    plink_gm = detect.load_any(write_plink(tmp_path))
    assert plink_gm.source_format == "PLINK1"
    assert_invariants(vcf_gm)
    assert_invariants(plink_gm)


def test_load_any_varimat_defers(tmp_path):
    p = tmp_path / "s.tsv"
    p.write_text("CHROM\tSTART\tREF\tALT\tVARCLASS\n1\t100\tA\tG\tMISSENSE\n")
    with pytest.raises(NotImplementedError, match="varimat"):
        detect.load_any(p)


def test_load_any_unknown_raises(tmp_path):
    p = tmp_path / "x.dat"
    p.write_text("nothing useful\n")
    with pytest.raises(ValueError, match="unrecognised"):
        detect.load_any(p)


# ================================================================== helpers ==

def assert_invariants(gm: GenotypeMatrix) -> None:
    assert gm.dosages.shape == (len(gm.variants), len(gm.sample_ids))
    assert gm.dosages.shape == (gm.n_variants, gm.n_samples)
    assert gm.dosages.dtype == np.int8
    assert set(np.unique(gm.dosages)) <= {MISSING, 0, 1, 2}
    assert gm.source_files and all(isinstance(f, str) for f in gm.source_files)


# ---------------------------------------------------- representative thinning --
def test_every_nth_samples_across_the_file_not_just_the_start(tmp_path):
    """max_variants alone takes the FIRST n records. On a chromosome-scale VCF
    that is the short arm and nothing else, so a thinned dataset has to sample
    the whole length to be representative."""
    from backend.app.research.formats import vcf as vcf_mod

    p = tmp_path / "long.vcf"
    header = ("##fileformat=VCFv4.2\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n")
    rows = "".join("1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n".format(1000 + i * 10)
                   for i in range(1000))
    p.write_text(header + rows)

    every = vcf_mod.read_vcf(p, every_nth=10)
    assert every.n_variants == 100
    positions = [v.pos for v in every.variants]
    # Spread across the file, and starting at the first record.
    assert positions[0] == 1000
    assert positions[-1] == 1000 + 990 * 10
    assert positions == sorted(positions)

    capped = vcf_mod.read_vcf(p, max_variants=100)
    assert [v.pos for v in capped.variants][-1] == 1000 + 99 * 10, (
        "capping should stop early — that is the behaviour every_nth exists to "
        "avoid")


def test_every_nth_and_max_variants_compose(tmp_path):
    from backend.app.research.formats import vcf as vcf_mod
    p = tmp_path / "c.vcf"
    p.write_text("##fileformat=VCFv4.2\n"
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                 + "".join("1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n".format(100 + i)
                           for i in range(500)))
    gm = vcf_mod.read_vcf(p, max_variants=20, every_nth=5)
    assert gm.n_variants == 20, "the cap applies to KEPT records, not records read"
    assert [v.pos for v in gm.variants] == [100 + i * 5 for i in range(20)]


def test_every_nth_of_one_is_the_unthinned_read(tmp_path):
    from backend.app.research.formats import vcf as vcf_mod
    p = tmp_path / "d.vcf"
    p.write_text("##fileformat=VCFv4.2\n"
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                 + "".join("1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\n".format(200 + i)
                           for i in range(50)))
    assert vcf_mod.read_vcf(p, every_nth=1).n_variants == 50
    assert vcf_mod.read_vcf(p).n_variants == 50


# ----------------------------------------------------------- density labels --
def test_a_whole_chromosome_is_not_called_a_single_gene():
    """The label a reader sees. Classifying on chromosome count alone put
    123,347 variants spanning chr22 — a thinned chromosome of 1000 Genomes —
    into "single_gene". A gene spans tens to hundreds of kilobases; no assay of
    one gene yields a hundred thousand sites."""
    from backend.app.research.profile import classify_density

    gene = classify_density(1_200, ["17"])
    assert gene["density_class"] == "single_gene"

    chrom = classify_density(123_347, ["22"])
    assert chrom["density_class"] == "single_chromosome", chrom
    assert "chromosome-scale" in chrom["density_rationale"]
    # And it must not claim to be genome-wide either, which is what gates GWAS.
    assert chrom["density_class"] not in ("genome_wide", "exome")


def test_density_classes_still_separate_panel_from_genome():
    from backend.app.research.profile import classify_density
    autosomes = [str(c) for c in range(1, 23)]
    assert classify_density(4_000, ["1", "2", "3", "7"])["density_class"] == "targeted"
    assert classify_density(500_000, autosomes)["density_class"] == "genome_wide"


def test_a_single_chromosome_does_not_unlock_genome_wide_analyses():
    """The point of the label: GWAS, burden and PRS need genome-wide data, and
    one chromosome is not that however many variants it carries."""
    from backend.app.research import capability
    from backend.app.research.profile import DataProfile

    p = DataProfile(n_samples=2504, n_variants=123_347,
                    density_class="single_chromosome", genome_build="GRCh37")
    p.phenotypes = {"status": {"kind": "binary", "cases": 1200, "controls": 1300}}
    p.n_binary_phenotypes = 1
    p.max_cases = 1200
    p.has_controls = True
    p.ancestry = {"n_pcs": 10}
    p.mean_call_rate = 0.99

    by_name = {c.analysis: c for c in capability.assess(p)}
    for locked in ("gwas", "prs"):
        assert not by_name[locked].available, (
            "{} unlocked on a single chromosome".format(locked))
