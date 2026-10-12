"""The run record written beside every matrix."""

import json

import numpy as np

from ma_mapper import __version__, provenance


def test_describes_a_small_file_with_a_checksum(tmp_path):
    path = tmp_path / "input.txt"
    path.write_text("hello")
    record = provenance.describe_file(str(path))
    assert record["bytes"] == 5
    assert len(record["sha256"]) == 64
    assert record["path"].startswith("/")


def test_large_files_are_identified_without_hashing(tmp_path):
    """Hashing a 182 GB gnomAD VCF would cost more than the analysis."""
    path = tmp_path / "big.bin"
    path.write_bytes(b"0" * 1024)
    record = provenance.describe_file(str(path), max_hash_bytes=512)
    assert "sha256" not in record
    assert "too large to hash" in record["identified_by"]
    assert record["bytes"] == 1024


def test_a_missing_file_is_recorded_not_raised(tmp_path):
    record = provenance.describe_file(str(tmp_path / "gone.txt"))
    assert "error" in record


def test_an_unused_input_role_stays_explicit():
    record = provenance.run_record({"coordinates": None}, {})
    assert record["inputs"]["coordinates"] is None


def test_record_round_trips(tmp_path):
    path = tmp_path / "run.json"
    provenance.write_run_record(
        str(path),
        inputs={"alignment": str(path)},
        parameters={"filter": "one_pass"},
        outputs={"matrix": {"shape": [2, 3]}},
        command="ma-mapper overlay ...",
    )
    record = provenance.read_run_record(str(path))
    assert record["ma_mapper_version"] == __version__
    assert record["parameters"]["filter"] == "one_pass"
    assert record["versions"]["numpy"]
    assert record["created_utc"].endswith("Z")


def test_overlay_writes_all_four_files(alignment_fasta, position_bigwig, tmp_path):
    from ma_mapper import cli

    out = tmp_path / "matrix.npy"
    assert cli.main([
        "overlay",
        "--alignment", alignment_fasta,
        "--data", position_bigwig,
        "--format", "bigwig",
        "--filter", "none",
        "--output", str(out),
    ]) == 0
    record = json.loads((tmp_path / "matrix.run.json").read_text())
    assert record["parameters"]["data_format"] == "bigwig"
    assert record["inputs"]["alignment"]["sha256"]
    assert record["outputs"]["matrix"]["shape"] == [4, 21]
    # the record must describe the matrix that is actually on disk
    assert list(np.load(out).shape) == record["outputs"]["matrix"]["shape"]
