from __future__ import annotations

import json
from pathlib import Path

import pytest

import pyslimmc
from pyslimmc.runs import SelectionError
from pyslimmc.tests.test_l2_5_moments_spectra import build


def _make_run(root: Path, run_id: str) -> None:
    root.mkdir(parents=True)
    build(root)
    metadata_path = root / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["run_id"] = run_id
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")


def _runs(tmp_path: Path):
    for run_id in ("f5a_DBI30_I3", "f5a_DBI50_I3", "f5a_DBI70_I3"):
        _make_run(tmp_path / run_id, run_id)
    return pyslimmc.scan(tmp_path)


def test_pack_defaults_to_full_run_id(tmp_path: Path):
    runs = _runs(tmp_path)
    packed = runs.pack()
    assert list(packed) == ["f5a_DBI30_I3", "f5a_DBI50_I3", "f5a_DBI70_I3"]
    assert packed["f5a_DBI30_I3"]["run"].run_id == "f5a_DBI30_I3"
    assert set(packed["f5a_DBI30_I3"]) == {"run"}


def test_pack_extracts_key_and_user_fields(tmp_path: Path):
    runs = _runs(tmp_path)
    packed = runs.pack(
        key="f5a_*_I3",
        label="f5a_*",
        color=["tab:blue", "tab:orange", "tab:green"],
        offset=(0.00, 0.05, 0.10),
        linewidth=2,
    )
    assert list(packed) == ["DBI30", "DBI50", "DBI70"]
    assert packed["DBI30"]["label"] == "DBI30_I3"
    assert packed["DBI50"]["color"] == "tab:orange"
    assert packed["DBI70"]["offset"] == 0.10
    assert all(item["linewidth"] == 2 for item in packed.values())


def test_pack_after_match_preserves_runs_order(tmp_path: Path):
    packed = _runs(tmp_path).match("f5a_DBI[35]*_I3").pack(key="f5a_*_I3")
    assert list(packed) == ["DBI30", "DBI50"]


def test_pack_rejects_invalid_key_patterns(tmp_path: Path):
    runs = _runs(tmp_path)
    with pytest.raises(ValueError, match="exactly one"):
        runs.pack(key="f5a_DBI_I3")
    with pytest.raises(ValueError, match="exactly one"):
        runs.pack(key="f5a_*_*_I3")
    with pytest.raises(SelectionError, match="does not match"):
        runs.pack(key="other_*_I3")


def test_pack_rejects_duplicate_extracted_keys(tmp_path: Path):
    _make_run(tmp_path / "a", "same_id")
    _make_run(tmp_path / "b", "same_id")
    runs = pyslimmc.scan(tmp_path)
    with pytest.raises(SelectionError, match="not unique"):
        runs.pack(key="*_id")


def test_pack_rejects_bad_field_values(tmp_path: Path):
    runs = _runs(tmp_path)
    with pytest.raises(ValueError, match="3 runs"):
        runs.pack(color=["blue", "orange"])
    with pytest.raises(ValueError, match="reserved"):
        runs.pack(run="anything")
    with pytest.raises(ValueError, match="exactly one"):
        runs.pack(label="f5a_*_DBI*_I3")
    with pytest.raises(SelectionError, match="does not match"):
        runs.pack(label="other_*")


def test_pack_empty_collection(tmp_path: Path):
    empty = pyslimmc.scan(tmp_path)
    assert empty.pack() == {}
    assert empty.pack(key="prefix_*", color=[]) == {}


def test_pack_has_pre_use_help():
    text = pyslimmc.Runs.pack.help()
    assert "ordered dictionary" in text
    assert 'key="prefix*suffix"' in text
    assert "reserved" in text
