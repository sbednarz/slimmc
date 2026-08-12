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


def test_direct_run_id_access_is_lazy_and_cached(tmp_path: Path):
    _make_run(tmp_path / "a", "case_A")
    _make_run(tmp_path / "b", "case_B")
    runs = pyslimmc.scan(tmp_path)

    assert runs._run_id_index_cache is None
    assert runs.case_A.run_id == "case_A"
    first_cache = runs._run_id_index_cache
    assert first_cache is not None
    assert runs.case_B.run_id == "case_B"
    assert runs._run_id_index_cache is first_cache


def test_dir_lists_only_unique_noncolliding_run_ids(tmp_path: Path):
    _make_run(tmp_path / "a", "case_A")
    _make_run(tmp_path / "dup1", "duplicate")
    _make_run(tmp_path / "dup2", "duplicate")
    _make_run(tmp_path / "api", "match")
    runs = pyslimmc.scan(tmp_path)

    names = dir(runs)
    assert "case_A" in names
    assert "duplicate" not in names
    assert "match" in names  # the real API method, never a dynamic run
    assert callable(runs.match)


def test_duplicate_direct_access_raises_selection_error(tmp_path: Path):
    _make_run(tmp_path / "x" / "one", "same_id")
    _make_run(tmp_path / "y" / "two", "same_id")
    runs = pyslimmc.scan(tmp_path)

    with pytest.raises(SelectionError, match="not unique"):
        _ = runs.same_id
    assert len(runs.run_id["same_id"]) == 2


def test_missing_direct_access_raises_attribute_error(tmp_path: Path):
    _make_run(tmp_path / "a", "case_A")
    runs = pyslimmc.scan(tmp_path)
    with pytest.raises(AttributeError, match="missing"):
        _ = runs.missing


def test_subcollections_have_independent_lazy_indexes(tmp_path: Path):
    _make_run(tmp_path / "a", "group_A")
    _make_run(tmp_path / "b", "group_B")
    runs = pyslimmc.scan(tmp_path)
    subset = runs.match("group_A")

    assert runs._run_id_index_cache is None
    assert subset._run_id_index_cache is None
    assert subset.group_A.run_id == "group_A"
    assert subset._run_id_index_cache is not None
    assert runs._run_id_index_cache is None


def test_thousand_run_dir_protocol_uses_metadata_only(tmp_path: Path):
    # Construct lightweight fake Runs directly: the protocol must only read run_id/path.
    class FakeRun:
        def __init__(self, index: int):
            self.run_id = f"run_{index:04d}"
            self.path = tmp_path / self.run_id

    fake = {str(tmp_path / f"run_{i:04d}"): FakeRun(i) for i in range(1000)}
    runs = pyslimmc.Runs(tmp_path, fake)  # type: ignore[arg-type]
    assert runs._run_id_index_cache is None
    names = dir(runs)
    assert "run_0000" in names
    assert "run_0999" in names
    assert len(runs._run_id_index_cache or {}) == 1000
