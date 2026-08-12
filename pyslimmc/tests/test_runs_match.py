from __future__ import annotations

import json
from pathlib import Path

import pytest

import pyslimmc
from pyslimmc.tests.test_l2_5_moments_spectra import build


def _make_run(root: Path, run_id: str) -> None:
    root.mkdir()
    build(root)
    metadata_path = root / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["run_id"] = run_id
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")


def _runs(tmp_path: Path):
    for run_id in (
        "f5a_DBI30_I3",
        "f5a_DBI50_I3",
        "f6b_DBI70_I2",
        "run_001",
        "run_012",
        "run_103",
    ):
        _make_run(tmp_path / run_id, run_id)
    return pyslimmc.scan(tmp_path)


def test_match_supports_full_glob_and_preserves_runs_type(tmp_path: Path):
    runs = _runs(tmp_path)
    selected = runs.match("f5a_DBI*_I3")
    assert isinstance(selected, pyslimmc.Runs)
    assert [run.run_id for run in selected] == ["f5a_DBI30_I3", "f5a_DBI50_I3"]
    assert len(runs) == 6


def test_match_supports_question_ranges_and_negated_ranges(tmp_path: Path):
    runs = _runs(tmp_path)
    assert [r.run_id for r in runs.match("run_00?")] == ["run_001"]
    assert [r.run_id for r in runs.match("run_0[0-2][0-9]")] == ["run_001", "run_012"]
    assert [r.run_id for r in runs.match("run_[!0]*")] == ["run_103"]


def test_match_is_case_sensitive_and_matches_complete_run_id(tmp_path: Path):
    runs = _runs(tmp_path)
    assert len(runs.match("F5A*")) == 0
    assert len(runs.match("DBI30")) == 0
    assert [r.run_id for r in runs.match("*DBI30*")] == ["f5a_DBI30_I3"]


def test_match_empty_result_and_chaining(tmp_path: Path):
    runs = _runs(tmp_path)
    assert len(runs.match("missing_*")) == 0
    selected = runs.match("f*_DBI*_I?").match("f5a_*")
    assert [r.run_id for r in selected] == ["f5a_DBI30_I3", "f5a_DBI50_I3"]


def test_match_rejects_non_string_pattern(tmp_path: Path):
    runs = _runs(tmp_path)
    with pytest.raises(TypeError, match="pattern must be a string"):
        runs.match(123)  # type: ignore[arg-type]


def test_match_has_pre_use_help():
    text = pyslimmc.Runs.match.help()
    assert "shell-style glob" in text
    assert "[!abc]" in text
