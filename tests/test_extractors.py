"""bigWig, BED and BAM extraction: orientation, fragment order, dtypes, errors."""

import numpy as np
import pytest

from ma_mapper import extract_bam, extract_bed, extract_bigwig

from .conftest import CHROM, expected_array


class TestBigwig:
    def test_plus_strand_fragment(self, position_bigwig):
        values = extract_bigwig.extract_bigwig(position_bigwig, CHROM, [100], [120], "+")
        np.testing.assert_allclose(values, expected_array("copyA"))

    def test_minus_strand_is_reversed(self, position_bigwig):
        values = extract_bigwig.extract_bigwig(position_bigwig, CHROM, [200], [220], "-")
        np.testing.assert_allclose(values, expected_array("copyB"))

    def test_fragments_are_concatenated_in_element_order(self, position_bigwig):
        values = extract_bigwig.extract_bigwig(
            position_bigwig, CHROM, [300, 320], [310, 330], "+"
        )
        np.testing.assert_allclose(values, expected_array("copyC"))

    def test_unmeasured_bases_are_nan_not_zero(self, position_bigwig):
        values = extract_bigwig.extract_bigwig(position_bigwig, CHROM, [248], [262], "+")
        assert np.isnan(values[2:12]).all()
        assert values[0] == pytest.approx(248)

    def test_unknown_chromosome_names_the_problem(self, position_bigwig):
        with pytest.raises(KeyError, match="chr-prefix"):
            extract_bigwig.extract_bigwig(position_bigwig, "1", [100], [120], "+")

    def test_out_of_range_fragment_is_rejected(self, position_bigwig):
        with pytest.raises(ValueError, match="outside the bigWig"):
            extract_bigwig.extract_bigwig(position_bigwig, CHROM, [390], [500], "+")

    def test_io_over_the_whole_family(self, coordinate_table, position_bigwig):
        out = extract_bigwig.bigwig_io(coordinate_table, position_bigwig)
        assert len(out) == 4
        for values, name in zip(out, ["copyA", "copyB", "copyC", "copyD"]):
            np.testing.assert_allclose(values, expected_array(name))

    def test_io_does_not_modify_the_caller_table(self, coordinate_table, position_bigwig):
        before = coordinate_table.copy()
        extract_bigwig.bigwig_io(coordinate_table, position_bigwig)
        assert list(coordinate_table.columns) == list(before.columns)

    def test_io_rejects_a_bad_type(self, coordinate_table):
        with pytest.raises(TypeError):
            extract_bigwig.extract_bigwig(42, CHROM, [100], [120], "+")


@pytest.mark.requires_bedtools
class TestBed:
    def test_integer_score_is_painted_over_the_overlap(self, coordinate_table, feature_bed):
        out = extract_bed.bed_io(coordinate_table, feature_bed)
        copya = out[0]
        # feat_int covers chrT 100-105, i.e. the first five bases of copyA
        np.testing.assert_allclose(copya[:5], [3, 3, 3, 3, 3])

    def test_float_score_is_not_truncated(self, coordinate_table, feature_bed):
        """int() used to turn a score of 0.75 into 0."""
        out = extract_bed.bed_io(coordinate_table, feature_bed)
        copya = out[0]
        np.testing.assert_allclose(copya[10:12], [0.75, 0.75])

    def test_dot_score_counts_as_one(self, coordinate_table, feature_bed):
        """int('.') used to raise ValueError."""
        out = extract_bed.bed_io(coordinate_table, feature_bed)
        copyb = out[1]
        # feat_dot covers chrT 205-210 on the minus strand; copyB is 200-220
        # minus, so those bases are 10..14 from the element's 5' end.
        np.testing.assert_allclose(copyb[10:15], [1, 1, 1, 1, 1])

    def test_uncovered_bases_are_zero(self, coordinate_table, feature_bed):
        out = extract_bed.bed_io(coordinate_table, feature_bed)
        assert out[3].sum() == 0  # copyD has no features

    def test_lengths_match_the_fragments(self, coordinate_table, feature_bed):
        out = extract_bed.bed_io(coordinate_table, feature_bed)
        assert [len(a) for a in out] == [20, 20, 20, 15]

    def test_dataframe_input_is_accepted(self, coordinate_table, feature_bed):
        """A DataFrame used to hit a leftover BedTool() call and raise TypeError."""
        import pandas as pd

        features = pd.read_csv(feature_bed, sep="\t", header=None)
        out = extract_bed.bed_io(coordinate_table, features)
        np.testing.assert_allclose(out[0][:5], [3, 3, 3, 3, 3])

    def test_rejects_a_bad_type(self, coordinate_table):
        with pytest.raises(TypeError, match="bed must be"):
            extract_bed.bed_io(coordinate_table, 42)


class TestBam:
    def test_read_counts_above_255_do_not_wrap(self, coordinate_table, read_bam):
        """The profile used to be uint8, which wraps to 44 at 300 reads."""
        out = extract_bam.bam_io(coordinate_table, read_bam, "read_forward", offset=0)
        copya = out[0]
        assert copya[5] == 300  # 300 reads start at chrT 105, 5 bases into copyA
        assert copya[10] == 2

    def test_default_offset_shifts_by_five(self, coordinate_table, read_bam):
        out = extract_bam.bam_io(coordinate_table, read_bam, "read_forward")
        assert out[0][10] == 300  # 105 + offset 5

    def test_forward_means_along_the_element(self, coordinate_table, read_bam):
        """On a minus-strand copy, a genome-forward read is a reverse read."""
        plus = extract_bam.bam_io(coordinate_table, read_bam, "read_forward", offset=0)
        assert plus[0].sum() == 302
        reverse = extract_bam.bam_io(coordinate_table, read_bam, "read_reverse", offset=0)
        assert reverse[0].sum() == 1  # the single genome-reverse read

    def test_read_sum_adds_both_profiles(self, coordinate_table, read_bam):
        out = extract_bam.bam_io(coordinate_table, read_bam, "read_sum", offset=0)
        assert out[0].sum() == 303

    def test_normal_profile_is_float32(self, coordinate_table, read_bam):
        out = extract_bam.bam_io(coordinate_table, read_bam, "normal_forward", offset=0)
        assert out[0].dtype == np.float32
        # float16 saturated at 65504; 300 summed unit Gaussians exceed that only
        # in aggregate, but the peak must at least be the read count.
        assert out[0].max() >= 300

    def test_unknown_format_is_rejected(self, coordinate_table, read_bam):
        with pytest.raises(ValueError, match="bam_format must be"):
            extract_bam.bam_io(coordinate_table, read_bam, "read_everything")

    def test_lengths_match_the_fragments(self, coordinate_table, read_bam):
        out = extract_bam.bam_io(coordinate_table, read_bam, "read_sum")
        assert [len(a) for a in out] == [20, 20, 20, 15]

    def test_missing_file_raises(self, coordinate_table, tmp_path):
        with pytest.raises(FileNotFoundError):
            extract_bam.bam_io(coordinate_table, str(tmp_path / "nope.bam"), "read_sum")
