"""The command line, called in process.

`tests/test_end_to_end.py` runs the installed entry point as a subprocess,
which is the honest test of "does `ma-mapper` work". These call `main()`
directly, which is faster, puts the failure in the traceback where it belongs,
and lets coverage see the module.
"""

import json

import numpy as np
import pandas as pd
import pytest

from ma_mapper import cli


def _run(*args):
    return cli.main(list(args))


class TestInfoAndHelp:
    def test_info_lists_dependencies(self, capsys):
        assert _run("info") == 0
        out = capsys.readouterr().out
        assert "ma_mapper" in out and "numpy" in out
        # every dependency the package declares should resolve in a test env
        assert "MISSING" not in out.split("mafft")[0]

    def test_version_exits_zero(self, capsys):
        with pytest.raises(SystemExit) as exit_info:
            _run("--version")
        assert exit_info.value.code == 0

    def test_no_subcommand_is_an_error(self):
        with pytest.raises(SystemExit) as exit_info:
            _run()
        assert exit_info.value.code != 0


class TestSequence:
    def test_writes_one_record_per_copy(self, coordinate_table, genome_fasta, tmp_path):
        coordinates = tmp_path / "coords.bed"
        coordinate_table.to_csv(coordinates, sep="\t", header=False, index=False)
        out = tmp_path / "seqs.fasta"
        assert _run("sequence", "--coordinates", str(coordinates),
                    "--genome", genome_fasta, "--output", str(out)) == 0
        assert out.read_text().count(">") == 4

    def test_repeatmasker_without_subfamily_is_rejected(self, genome_fasta, tmp_path):
        with pytest.raises(SystemExit, match="needs --subfamily"):
            _run("sequence", "--repeatmasker", "x.tsv", "--genome", genome_fasta,
                 "--output", str(tmp_path / "o.fasta"))

    def test_neither_input_is_rejected(self, genome_fasta, tmp_path):
        with pytest.raises(SystemExit, match="either --coordinates or --repeatmasker"):
            _run("sequence", "--genome", genome_fasta, "--output", str(tmp_path / "o.fasta"))


class TestOverlay:
    def test_writes_four_files(self, alignment_fasta, position_bigwig, tmp_path, capsys):
        out = tmp_path / "m.npy"
        assert _run("overlay", "--alignment", alignment_fasta, "--data", position_bigwig,
                    "--format", "bigwig", "--filter", "none", "--output", str(out)) == 0
        printed = capsys.readouterr().out.split()
        assert len(printed) == 4
        assert np.load(out).shape == (4, 21)
        assert np.load(tmp_path / "m.alignment.npy").shape == (4, 21)
        assert len(pd.read_csv(tmp_path / "m.rows.tsv", sep="\t")) == 5  # copyC is two fragments
        record = json.loads((tmp_path / "m.run.json").read_text())
        assert record["parameters"]["filter"] == "none"

    def test_filter_thresholds_reach_the_filter(self, alignment_fasta, position_bigwig, tmp_path):
        out = tmp_path / "m.npy"
        assert _run("overlay", "--alignment", alignment_fasta, "--data", position_bigwig,
                    "--format", "bigwig", "--output", str(out),
                    "--col-threshold", "0.1", "--col-content-threshold", "0.1",
                    "--row-threshold", "0.1") == 0
        record = json.loads((tmp_path / "m.run.json").read_text())
        assert record["parameters"]["col_threshold"] == 0.1
        assert record["parameters"]["row_threshold"] == 0.1
        # The thresholds reached the filter, and at 0.1 every column of this
        # four-copy alignment is occupied enough to survive.
        matrix = np.load(out)
        assert matrix.shape == (4, 21)
        strict = tmp_path / "strict.npy"
        _run("overlay", "--alignment", str(record["inputs"]["alignment"]["path"]),
             "--data", str(record["inputs"]["data"]["path"]),
             "--format", "bigwig", "--output", str(strict),
             "--col-threshold", "0.9", "--col-content-threshold", "0.9",
             "--row-threshold", "0.1")
        assert np.load(strict).shape[1] < 21

    def test_iterative_filter_is_reachable(self, alignment_fasta, position_bigwig, tmp_path):
        out = tmp_path / "m.npy"
        assert _run("overlay", "--alignment", alignment_fasta, "--data", position_bigwig,
                    "--format", "bigwig", "--filter", "iterative",
                    "--col-threshold", "0.1", "--row-threshold", "0.1",
                    "--output", str(out)) == 0
        assert np.load(out).size > 0

    def test_a_missing_input_exits_two(self, tmp_path, capsys):
        with pytest.raises(SystemExit) as exit_info:
            _run("overlay", "--alignment", str(tmp_path / "nope.fasta"),
                 "--data", str(tmp_path / "nope.bw"), "--format", "bigwig",
                 "--output", str(tmp_path / "o.npy"))
        assert exit_info.value.code == 2

    def test_an_unknown_format_exits_two(self, alignment_fasta, position_bigwig, tmp_path):
        with pytest.raises(SystemExit) as exit_info:
            _run("overlay", "--alignment", alignment_fasta, "--data", position_bigwig,
                 "--format", "bigwug", "--output", str(tmp_path / "o.npy"))
        assert exit_info.value.code == 2


class TestPlot:
    def _matrix(self, alignment_fasta, position_bigwig, tmp_path):
        out = tmp_path / "m.npy"
        _run("overlay", "--alignment", alignment_fasta, "--data", position_bigwig,
             "--format", "bigwig", "--filter", "none", "--output", str(out))
        return out

    def test_profile_uses_the_sidecar_alignment(self, alignment_fasta, position_bigwig, tmp_path):
        matrix = self._matrix(alignment_fasta, position_bigwig, tmp_path)
        png = tmp_path / "profile.png"
        assert _run("plot", "--matrix", str(matrix), "--kind", "profile",
                    "--output", str(png)) == 0
        assert png.stat().st_size > 1000

    def test_heatmap(self, alignment_fasta, position_bigwig, tmp_path):
        matrix = self._matrix(alignment_fasta, position_bigwig, tmp_path)
        png = tmp_path / "heatmap.png"
        assert _run("plot", "--matrix", str(matrix), "--kind", "heatmap",
                    "--output", str(png)) == 0
        assert png.stat().st_size > 1000

    def test_profile_without_an_alignment_says_so(self, alignment_fasta, position_bigwig, tmp_path):
        matrix = self._matrix(alignment_fasta, position_bigwig, tmp_path)
        (tmp_path / "m.alignment.npy").unlink()
        with pytest.raises(SystemExit, match="needs the alignment"):
            _run("plot", "--matrix", str(matrix), "--kind", "profile",
                 "--output", str(tmp_path / "p.png"))

    def test_a_mismatched_alignment_is_refused(self, alignment_fasta, position_bigwig, tmp_path):
        matrix = self._matrix(alignment_fasta, position_bigwig, tmp_path)
        np.save(tmp_path / "m.alignment.npy", np.ones((2, 2), dtype=np.uint8))
        with pytest.raises(SystemExit, match="different filters"):
            _run("plot", "--matrix", str(matrix), "--kind", "profile",
                 "--output", str(tmp_path / "p.png"))


class TestIndexMaf:
    def test_builds_an_index(self, three_species_maf, tmp_path, capsys):
        assert _run("index-maf", "--maf", three_species_maf, "--chrom", "chrT") == 0
        assert capsys.readouterr().out.strip().endswith(".mafindex")


class TestLogging:
    def test_verbose_and_logfile(self, alignment_fasta, position_bigwig, tmp_path):
        logfile = tmp_path / "run.log"
        assert _run("--verbose", "--logfile", str(logfile), "overlay",
                    "--alignment", alignment_fasta, "--data", position_bigwig,
                    "--format", "bigwig", "--filter", "none",
                    "--output", str(tmp_path / "m.npy")) == 0
        assert logfile.exists() and logfile.stat().st_size > 0
