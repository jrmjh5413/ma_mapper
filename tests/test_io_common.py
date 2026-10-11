"""The shared coordinate handling: grouping, ordering and saving."""

import os

import pandas as pd
import pytest

from ma_mapper import _io_common


def test_read_rejects_a_missing_path(tmp_path):
    with pytest.raises(FileNotFoundError, match="coordinate_table file not found"):
        _io_common.read_coordinate_table(str(tmp_path / "nope.bed"))


def test_read_rejects_a_bad_type():
    with pytest.raises(TypeError):
        _io_common.read_coordinate_table(42)


def test_groups_one_entry_per_copy(coordinate_table):
    meta_id, chrom, starts, ends, strands = _io_common.group_coordinates(coordinate_table)
    assert meta_id == ["copyA", "copyB", "copyC", "copyD"]
    assert starts[2] == [300, 320]
    assert strands == ["+", "-", "+", "-"]
    assert chrom == ["chrT"] * 4


def test_does_not_modify_the_caller_table(coordinate_table):
    before = coordinate_table.copy()
    _io_common.group_coordinates(coordinate_table)
    pd.testing.assert_frame_equal(coordinate_table, before)


def test_fragments_are_sorted_and_reported(caplog):
    table = pd.DataFrame(
        [
            ["chrT", 320, 330, "copyC", 10, "+"],
            ["chrT", 300, 310, "copyC", 10, "+"],
        ]
    )
    with caplog.at_level("WARNING", logger="ma_mapper"):
        _, _, starts, ends, _ = _io_common.group_coordinates(table)
    assert starts[0] == [300, 320]
    assert ends[0] == [310, 330]
    assert "ascending order" in caplog.text


def test_unsorted_fragments_can_be_left_alone():
    table = pd.DataFrame(
        [
            ["chrT", 320, 330, "copyC", 10, "+"],
            ["chrT", 300, 310, "copyC", 10, "+"],
        ]
    )
    _, _, starts, _, _ = _io_common.group_coordinates(table, sort_blocks=False)
    assert starts[0] == [320, 300]


def test_generate_new_id_makes_every_row_its_own_copy(coordinate_table):
    meta_id, _, starts, _, _ = _io_common.group_coordinates(
        coordinate_table, generate_new_id=True
    )
    assert len(meta_id) == 5
    assert all(len(s) == 1 for s in starts)


def test_too_few_columns_is_rejected():
    with pytest.raises(ValueError, match="at least 6 BED columns"):
        _io_common.group_coordinates(pd.DataFrame([["chrT", 1, 2]]))


def test_a_string_save_path_is_honoured(tmp_path):
    """`save_to_file == True` is False for a string, so a path was ignored."""
    destination = tmp_path / "out.p"
    result = _io_common.save_or_return([1, 2, 3], str(destination), None, "default.lzma", "thing")
    assert result == str(destination)
    assert destination.exists()


def test_true_saves_beside_the_coordinate_table(tmp_path):
    table = tmp_path / "coords.bed"
    table.write_text("chrT\t1\t2\tx\t10\t+\n")
    result = _io_common.save_or_return([1], True, str(table), "default.lzma", "thing")
    assert result == str(tmp_path / "default.lzma")
    assert os.path.isfile(result)


def test_falsy_returns_the_payload():
    payload = [1, 2]
    assert _io_common.save_or_return(payload, False, None, "x.lzma", "thing") is payload


def test_worker_count_reads_the_environment(monkeypatch):
    monkeypatch.setenv("MA_MAPPER_VCF_WORKERS", "3")
    assert _io_common.worker_count("MA_MAPPER_VCF_WORKERS") == 3
    monkeypatch.delenv("MA_MAPPER_VCF_WORKERS")
    monkeypatch.setenv("MA_MAPPER_WORKERS", "2")
    assert _io_common.worker_count("MA_MAPPER_VCF_WORKERS") == 2


def test_worker_count_defaults_to_the_allocation(monkeypatch):
    monkeypatch.delenv("MA_MAPPER_WORKERS", raising=False)
    monkeypatch.delenv("MA_MAPPER_VCF_WORKERS", raising=False)
    assert _io_common.worker_count() >= 1
    assert _io_common.worker_count(default=4) == 4
