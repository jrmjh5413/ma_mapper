"""The VCF extractor, including the minus-strand defect fixed in 0.9.0."""

import numpy as np
import pytest

from ma_mapper import extract_vcf

CHROM = "chrT"


def test_plus_strand_places_values_at_their_own_base(position_vcf):
    window = extract_vcf.extract_vcf(position_vcf, CHROM, [100], [120], "+", "AF")
    expected = np.zeros(20)
    expected[0] = 0.5  # one-based 101 is the first base of the window
    np.testing.assert_allclose(window, expected)


def test_three_variants_in_one_minus_strand_window(position_vcf):
    """Each value lands at its own base, counting from the element's 5' end.

    copyB is chrT 200-220 on the minus strand, with SNPs at one-based 201, 203
    and 207, i.e. 0, 2 and 6 bases from the left end of the window. Reversed,
    they belong 19, 17 and 13 bases in.

    Before 0.9.0 the window was flipped once per variant rather than once per
    window, so with three variants the second one came back mirrored about the
    centre and the whole window was returned in the wrong orientation when the
    count was even.
    """
    window = extract_vcf.extract_vcf(position_vcf, CHROM, [200], [220], "-", "AF")
    expected = np.zeros(20)
    expected[19] = 0.1
    expected[17] = 0.2
    expected[13] = 0.25 + 0.15  # multiallelic: both ALT frequencies
    np.testing.assert_allclose(window, expected)


def test_even_variant_count_is_still_reversed(tmp_path):
    """Two variants in a minus-strand window: the old code returned it unflipped."""
    import pysam

    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={CHROM},length=400>",
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        f"{CHROM}\t101\t.\tA\tC\t.\tPASS\tAF=0.1",
        f"{CHROM}\t103\t.\tA\tG\t.\tPASS\tAF=0.2",
    ]
    plain = tmp_path / "two.vcf"
    plain.write_text("\n".join(lines) + "\n")
    path = str(tmp_path / "two.vcf.gz")
    pysam.tabix_compress(str(plain), path, force=True)
    pysam.tabix_index(path, preset="vcf", force=True)

    window = extract_vcf.extract_vcf(path, CHROM, [100], [110], "-", "AF")
    expected = np.zeros(10)
    expected[9] = 0.1
    expected[7] = 0.2
    np.testing.assert_allclose(window, expected)


def test_fragment_order_is_reversed_for_minus_strand(position_vcf):
    """Two fragments: the element's 5' end is the rightmost genomic base."""
    window = extract_vcf.extract_vcf(position_vcf, CHROM, [100, 300], [110, 310], "-", "AF")
    assert window.size == 20
    # chrT:305 (one-based) is 4 bases into the second fragment, which comes
    # first once the order is reversed: 10 - 1 - 4 = 5.
    assert window[5] == pytest.approx(0.3)
    # chrT:101 is the first base of the first fragment, now last overall.
    assert window[19] == pytest.approx(0.5)


def test_missing_info_key_counts_as_zero(position_vcf):
    window = extract_vcf.extract_vcf(position_vcf, CHROM, [100], [120], "+", "NOT_A_FIELD")
    np.testing.assert_allclose(window, np.zeros(20))


def test_io_rejects_a_missing_file(coordinate_table, tmp_path):
    with pytest.raises(FileNotFoundError):
        extract_vcf.vcf_io(
            coordinate_table, str(tmp_path / "nope.vcf.gz"), vcf_format="single"
        )


def test_io_returns_one_array_per_copy(coordinate_table, position_vcf):
    out = extract_vcf.vcf_io(coordinate_table, position_vcf, vcf_format="single")
    assert [len(a) for a in out] == [20, 20, 20, 15]


def test_io_saves_when_given_a_path(coordinate_table, position_vcf, tmp_path):
    """A string path used to be compared against True and silently ignored."""
    import compress_pickle

    destination = tmp_path / "saved.lzma"
    returned = extract_vcf.vcf_io(
        coordinate_table, position_vcf, vcf_format="single", save_to_file=str(destination)
    )
    assert destination.exists()
    assert returned == str(destination)
    assert len(compress_pickle.load(str(destination), compression='lzma')) == 4
