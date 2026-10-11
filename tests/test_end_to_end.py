"""The whole chain, and the command line."""

import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from ma_mapper import mapper

from .conftest import CHROM, expected_array


def test_overlay_puts_every_value_in_its_own_column(alignment_fasta, position_bigwig):
    """Each copy's values land on its own bases, gaps stay empty."""
    matrix = mapper.map_and_overlay(
        alignment=alignment_fasta,
        data_file=position_bigwig,
        data_format="bigwig",
        filter=False,
    )
    assert matrix.shape == (4, 21)
    alignment_matrix = mapper.parse_alignment(alignment_fasta)
    for row, name in enumerate(["copyA", "copyB", "copyC", "copyD"]):
        non_gap = alignment_matrix[row] > 0
        np.testing.assert_allclose(matrix[row, non_gap], expected_array(name))
        assert (matrix[row, ~non_gap] == 0).all()


def test_overlay_returns_row_labels_on_request(alignment_fasta, position_bigwig):
    matrix, coordinates, alignment_matrix = mapper.map_and_overlay(
        alignment=alignment_fasta,
        data_file=position_bigwig,
        data_format="bigwig",
        filter=False,
        return_details=True,
    )
    assert list(coordinates["name"].drop_duplicates()) == ["copyA", "copyB", "copyC", "copyD"]
    assert alignment_matrix.shape == matrix.shape


def test_overlay_with_a_filter_keeps_the_tables_consistent(alignment_fasta, position_bigwig):
    """A fragmented copy used to break the coordinate filtering outright."""
    matrix, coordinates, alignment_matrix = mapper.map_and_overlay(
        alignment=alignment_fasta,
        data_file=position_bigwig,
        data_format="bigwig",
        filter=True,
        col_threshold=0.1,
        col_content_threshold=0.1,
        row_threshold=0.1,
        return_details=True,
    )
    assert matrix.shape[0] == coordinates["name"].nunique()
    assert alignment_matrix.shape == matrix.shape


def test_overlay_reorders_rows_to_the_coordinate_table(tmp_path, position_bigwig):
    """Row k of the matrix must be copy k of the coordinate table.

    Three copies with different numbers of bases, in a three-cycle against the
    alignment order. Before 0.9.0 the permutation was inverted, which paired
    each copy's alignment row with another copy's data.
    """
    rows = [
        ("cA::chrT:100-104(+)", "ACGT--"),
        ("cB::chrT:200-205(+)", "ACGTA-"),
        ("cC::chrT:300-306(+)", "ACGTAC"),
    ]
    alignment = tmp_path / "three.fasta"
    alignment.write_text("".join(f">{name}\n{seq}\n" for name, seq in rows))
    coordinates = pd.DataFrame(
        [
            [CHROM, 200, 205, "cB", 10, "+"],
            [CHROM, 300, 306, "cC", 10, "+"],
            [CHROM, 100, 104, "cA", 10, "+"],
        ]
    )
    coordinate_path = tmp_path / "coords.tsv"
    coordinates.to_csv(coordinate_path, sep="\t", header=False, index=False)

    matrix = mapper.map_and_overlay(
        alignment=str(alignment),
        data_file=position_bigwig,
        data_format="bigwig",
        filter=False,
        coordinate_file=str(coordinate_path),
    )
    assert matrix[0, 0] == pytest.approx(200)
    assert matrix[1, 0] == pytest.approx(300)
    assert matrix[2, 0] == pytest.approx(100)


def test_overlay_rejects_an_unknown_format(alignment_fasta, position_bigwig):
    with pytest.raises(ValueError, match="data_format must be"):
        mapper.map_and_overlay(
            alignment=alignment_fasta,
            data_file=position_bigwig,
            data_format="bam_max",  # the name the old type annotation advertised
            filter=False,
        )


def test_overlay_rejects_a_missing_coordinate_file(alignment_fasta, position_bigwig, tmp_path):
    with pytest.raises(FileNotFoundError):
        mapper.map_and_overlay(
            alignment=alignment_fasta,
            data_file=position_bigwig,
            data_format="bigwig",
            filter=False,
            coordinate_file=str(tmp_path / "nope.tsv"),
        )


def _cli(*args):
    return subprocess.run(
        [sys.executable, "-m", "ma_mapper.cli", *args],
        capture_output=True,
        text=True,
    )


class TestCommandLine:
    def test_version(self):
        result = _cli("--version")
        assert result.returncode == 0
        assert "ma_mapper" in result.stdout

    def test_help_lists_the_subcommands(self):
        result = _cli("--help")
        for command in ("sequence", "align", "overlay", "plot", "index-maf", "info"):
            assert command in result.stdout

    def test_info_reports_versions(self):
        result = _cli("info")
        assert result.returncode == 0
        assert "numpy" in result.stdout

    def test_overlay_writes_a_matrix_and_its_row_names(
        self, alignment_fasta, position_bigwig, tmp_path
    ):
        out = tmp_path / "matrix.npy"
        result = _cli(
            "overlay",
            "--alignment", alignment_fasta,
            "--data", position_bigwig,
            "--format", "bigwig",
            "--filter", "none",
            "--output", str(out),
        )
        assert result.returncode == 0, result.stderr
        matrix = np.load(out)
        assert matrix.shape == (4, 21)
        rows = pd.read_csv(tmp_path / "matrix.rows.tsv", sep="\t")
        assert set(rows["name"]) == {"copyA", "copyB", "copyC", "copyD"}

    def test_plot_writes_a_png(self, alignment_fasta, position_bigwig, tmp_path):
        matrix = tmp_path / "matrix.npy"
        _cli(
            "overlay", "--alignment", alignment_fasta, "--data", position_bigwig,
            "--format", "bigwig", "--filter", "none", "--output", str(matrix),
        )
        png = tmp_path / "profile.png"
        result = _cli(
            "plot", "--matrix", str(matrix), "--alignment", alignment_fasta,
            "--kind", "profile", "--output", str(png),
        )
        assert result.returncode == 0, result.stderr
        assert png.stat().st_size > 1000

    def test_a_missing_file_exits_cleanly(self, tmp_path):
        result = _cli(
            "overlay", "--alignment", str(tmp_path / "nope.fasta"),
            "--data", str(tmp_path / "nope.bw"), "--format", "bigwig",
            "--output", str(tmp_path / "out.npy"),
        )
        assert result.returncode == 2
        assert "not found" in result.stderr
        assert "Traceback" not in result.stderr

    def test_sequence_subcommand(self, coordinate_table, genome_fasta, tmp_path):
        coordinates = tmp_path / "coords.bed"
        coordinate_table.to_csv(coordinates, sep="\t", header=False, index=False)
        out = tmp_path / "seqs.fasta"
        result = _cli(
            "sequence", "--coordinates", str(coordinates),
            "--genome", genome_fasta, "--output", str(out),
        )
        assert result.returncode == 0, result.stderr
        assert out.read_text().count(">") == 4
