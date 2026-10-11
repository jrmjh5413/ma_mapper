"""Parsing, filtering and placing values in alignment columns."""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ma_mapper import mapper

from .conftest import ALIGNMENT_ROWS

GAP, A, C, T, G, N = 0, 1, 2, 3, 4, 5


def _records(pairs):
    return [SeqRecord(Seq(seq), id=name, name=name) for name, seq in pairs]


class TestParseAlignment:
    def test_codes_bases_gaps_and_n(self):
        matrix = mapper.parse_alignment(_records([("r1", "ACTG-N")]))
        np.testing.assert_array_equal(matrix[0], [A, C, T, G, GAP, N])

    def test_dot_is_a_gap(self):
        """HMMER writes '.' in insert columns; it means gap, not an unknown base."""
        matrix = mapper.parse_alignment(_records([("r1", "AC..GT")]))
        np.testing.assert_array_equal(matrix[0], [A, C, GAP, GAP, G, T])

    def test_lowercase_is_accepted(self):
        matrix = mapper.parse_alignment(_records([("r1", "acgt")]))
        np.testing.assert_array_equal(matrix[0], [A, C, G, T])

    def test_unknown_characters_become_gaps_with_a_warning(self, caplog):
        with caplog.at_level("WARNING", logger="ma_mapper"):
            matrix = mapper.parse_alignment(_records([("r1", "ACXG")]))
        np.testing.assert_array_equal(matrix[0], [A, C, GAP, G])
        assert "treated as gaps" in caplog.text

    def test_is_deterministic_for_unknown_characters(self):
        """The matrix used to be allocated with np.ndarray, i.e. uninitialised."""
        for _ in range(5):
            matrix = mapper.parse_alignment(_records([("r1", "XXXX")]))
            np.testing.assert_array_equal(matrix[0], [GAP, GAP, GAP, GAP])

    def test_a_generator_is_not_consumed_before_it_is_read(self):
        """A one-shot iterator used to produce a silently all-gap matrix."""
        matrix = mapper.parse_alignment(r for r in _records([("r1", "ACGT"), ("r2", "ACGT")]))
        assert matrix.shape == (2, 4)
        assert matrix.sum() > 0

    def test_empty_input_raises(self):
        with pytest.raises(ValueError, match="no records"):
            mapper.parse_alignment([])

    def test_empty_file_raises(self, tmp_path):
        path = tmp_path / "empty.fasta"
        path.write_text("")
        with pytest.raises(ValueError, match="no records"):
            mapper.parse_alignment(str(path))

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            mapper.parse_alignment(str(tmp_path / "nope.fasta"))

    def test_unaligned_input_raises_naming_the_record(self):
        with pytest.raises(ValueError, match="not aligned.*r2"):
            mapper.parse_alignment(_records([("r1", "ACGT"), ("r2", "ACG")]))

    def test_h5_round_trip(self, tmp_path):
        """parse_alignment wrote one dataset name and the reader looked for another."""
        path = tmp_path / "parsed.h5"
        original = mapper.parse_alignment(_records(ALIGNMENT_ROWS))
        mapper.parse_alignment(_records(ALIGNMENT_ROWS), save_to_file=str(path))
        np.testing.assert_array_equal(mapper.import_parsed_alignment(str(path)), original)


class TestCreateFilter:
    def test_all_gap_row_does_not_crash(self):
        """np.nonzero(row)[0][0] used to raise IndexError on an empty row."""
        matrix = mapper.parse_alignment(_records([("r1", "ACGT"), ("r2", "----")]))
        row_filter, col_filter = mapper.create_filter(matrix, row_threshold=0.0)
        assert row_filter.shape == (2,)
        assert col_filter.shape == (4,)

    def test_second_pass_recounts_rather_than_accumulating(self):
        """The fixed and legacy second passes disagree once rows are dropped.

        Two long intact copies and two short gappy ones. With a row threshold
        that drops the gappy pair, the surviving columns should be judged on
        the two survivors alone.
        """
        rows = [
            ("intact1", "ACGTACGTACGT"),
            ("intact2", "ACGTACGTACGT"),
            ("gappy1", "A-G-A-G-A-G-"),
            ("gappy2", "-C-T-C-T-C-T"),
        ]
        matrix = mapper.parse_alignment(_records(rows))
        fixed_rows, fixed_cols = mapper.create_filter(
            matrix, col_threshold=0.9, col_content_threshold=0.4, row_threshold=0.9
        )
        legacy_rows, legacy_cols = mapper.create_filter(
            matrix,
            col_threshold=0.9,
            col_content_threshold=0.4,
            row_threshold=0.9,
            legacy_second_pass=True,
        )
        # The fixed pass keeps the two intact copies and, because they are
        # non-gap everywhere, every column.
        np.testing.assert_array_equal(fixed_rows, [True, True, False, False])
        assert fixed_cols.all()
        # The legacy arithmetic keeps a gappy copy and throws away nine of the
        # twelve columns, because the column counts are a mixture of the counts
        # from all four copies and the counts from the survivors.
        assert legacy_cols.sum() == 3
        assert legacy_rows.sum() == 3

    def test_last_covered_column_counts_towards_the_span(self):
        """A copy's final non-gap column was left out of the span denominator."""
        rows = [("r1", "ACGT"), ("r2", "AC--")]
        matrix = mapper.parse_alignment(_records(rows))
        _, col_filter = mapper.create_filter(
            matrix, col_threshold=0.99, col_content_threshold=0.0, row_threshold=0.0
        )
        # Column 1 is the last base of r2 and is non-gap in both copies, so its
        # occupancy is 2/2 and it cannot exceed a 0.99 threshold on a count it
        # was previously excluded from.
        assert col_filter[0] and col_filter[1]

    def test_returns_boolean_masks(self):
        matrix = mapper.parse_alignment(_records(ALIGNMENT_ROWS))
        row_filter, col_filter = mapper.create_filter(matrix)
        assert row_filter.dtype == bool and col_filter.dtype == bool

    def test_content_threshold_uses_the_surviving_copies(self):
        """After the row filter, "a fraction of all copies" means the survivors.

        Judging the second pass against the original copy count makes the
        threshold stricter in proportion to how many copies were dropped. On
        the real MamRep434 alignment, with a row threshold of 0.5 that keeps
        130 of 2,102 copies, that left no columns at all.

        Two intact copies, plus eighteen whose gaps alternate so that no single
        column is mostly gap and every gappy copy fails the row threshold.
        """
        width = 20
        rows = [("keep1", "A" * width), ("keep2", "A" * width)]
        for i in range(18):
            rows.append(
                (f"drop{i}", "".join("A" if (column + i) % 2 else "-" for column in range(width)))
            )
        matrix = mapper.parse_alignment(_records(rows))
        row_filter, col_filter = mapper.create_filter(
            matrix, col_threshold=0.5, col_content_threshold=0.5, row_threshold=0.9
        )
        assert row_filter.sum() == 2
        # Both survivors are non-gap everywhere, so every surviving column
        # passes a 50% content threshold measured on them.
        assert col_filter.sum() > 0
        legacy_rows, legacy_cols = mapper.create_filter(
            matrix,
            col_threshold=0.5,
            col_content_threshold=0.5,
            row_threshold=0.9,
            legacy_second_pass=True,
        )
        assert col_filter.sum() >= legacy_cols.sum()


class TestMapData:
    def test_values_land_in_the_columns_of_their_own_bases(self):
        matrix = mapper.parse_alignment(_records([("r1", "AC-GT")]))
        mapped = mapper.map_data([np.array([1.0, 2.0, 3.0, 4.0])], matrix, filter=False)
        np.testing.assert_allclose(mapped[0], [1.0, 2.0, 0.0, 3.0, 4.0])

    def test_length_mismatch_raises_a_useful_error(self):
        matrix = mapper.parse_alignment(_records([("r1", "AC-GT")]))
        with pytest.raises(ValueError, match="4 non-gap bases.*3 extracted"):
            mapper.map_data([np.array([1.0, 2.0, 3.0])], matrix, filter=False)

    def test_copy_count_mismatch_raises(self):
        matrix = mapper.parse_alignment(_records([("r1", "ACGT"), ("r2", "ACGT")]))
        with pytest.raises(ValueError, match="1 copies but the alignment has 2"):
            mapper.map_data([np.ones(4)], matrix, filter=False)

    def test_dict_payloads_go_on_an_object_canvas(self):
        matrix = mapper.parse_alignment(_records([("r1", "AC")]))
        mapped = mapper.map_data([[{"A": 1.0}, {"C": 2.0}]], matrix, filter=False)
        assert mapped.dtype == object
        assert mapped[0, 0] == {"A": 1.0}

    def test_numeric_payloads_go_on_a_float_canvas(self):
        matrix = mapper.parse_alignment(_records([("r1", "AC")]))
        mapped = mapper.map_data([np.array([1.5, 2.5])], matrix, filter=False)
        assert mapped.dtype == np.float32


class TestBaseCount:
    def test_counts_every_symbol(self):
        matrix = mapper.parse_alignment(_records([("r1", "AC"), ("r2", "A-")]))
        counts = mapper.base_count(matrix)
        assert counts.shape == (6, 2)
        assert counts[A, 0] == 2
        assert counts[C, 1] == 1
        assert counts[GAP, 1] == 1

    def test_n_in_a_later_column_does_not_break_the_shape(self):
        """np.bincount(minlength=5) sized the output on the first column alone."""
        matrix = mapper.parse_alignment(_records([("r1", "AN")]))
        counts = mapper.base_count(matrix)
        assert counts.shape == (6, 2)
        assert counts[N, 1] == 1


class TestNormalise:
    def test_average_divides_by_the_copies_with_a_base(self):
        alignment = np.array([[A, A], [A, GAP]])
        data = np.array([[2.0, 4.0], [4.0, 0.0]])
        np.testing.assert_allclose(
            mapper.normalise(alignment, data, method="average"), [3.0, 4.0]
        )

    def test_measured_denominator_ignores_unmeasured_bases(self):
        alignment = np.array([[A, A], [A, A]])
        data = np.array([[2.0, 4.0], [4.0, np.nan]])
        np.testing.assert_allclose(
            mapper.normalise(alignment, data, method="average", denominator="alignment"),
            [3.0, 2.0],
        )
        np.testing.assert_allclose(
            mapper.normalise(alignment, data, method="average", denominator="measured"),
            [3.0, 4.0],
        )

    def test_empty_column_is_nan_not_inf(self):
        alignment = np.array([[GAP], [GAP]])
        data = np.array([[1.0], [1.0]])
        assert np.isnan(mapper.normalise(alignment, data, method="average")[0])

    def test_unknown_method_raises(self):
        with pytest.raises(ValueError, match="method must be"):
            mapper.normalise(np.array([[A]]), np.array([[1.0]]), method="mean")


class TestParseAndFilter:
    def test_coordinates_come_back_with_the_alignment(self, alignment_fasta):
        matrix, coordinates = mapper.parse_and_filter(alignment_fasta, filter=None)
        assert matrix.shape == (4, 21)
        # copyC is two fragments, so it contributes two coordinate rows: the
        # table is fragment-level and the name column is what ties those rows
        # to one alignment row.
        assert coordinates.shape[0] == 5
        assert list(coordinates["name"].drop_duplicates()) == [
            "copyA", "copyB", "copyC", "copyD"
        ]

    def test_rejects_an_unknown_filter_name(self, alignment_fasta):
        with pytest.raises(ValueError, match="Invalid value for `filter`"):
            mapper.parse_and_filter(alignment_fasta, filter="aggressive")


class TestExtractCoordinateFromAlignment:
    def test_parses_multi_fragment_headers(self, alignment_fasta):
        table = mapper.extract_coordinate_from_alignment(alignment_fasta)
        copyc = table[table.name == "copyC"]
        assert list(copyc.start) == [300, 320]
        assert list(copyc.end) == [310, 330]
        assert set(copyc.strand) == {"+"}

    def test_rejects_a_header_without_coordinates(self):
        with pytest.raises(ValueError, match="Invalid format"):
            mapper.extract_coordinate_from_alignment(_records([("plain_name", "ACGT")]))
