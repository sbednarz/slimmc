from __future__ import annotations

import json
from pathlib import Path

import pyslimmc
from pyslimmc.tests.test_l2_5_moments_spectra import build


def _make_run(root: Path, run_id: str) -> None:
    root.mkdir(parents=True)
    build(root)
    metadata_path = root / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["run_id"] = run_id
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")


def test_runs_help_covers_interactive_match_pack_sweep(tmp_path: Path):
    _make_run(tmp_path / "a", "f5a_DBI30_I3")
    runs = pyslimmc.scan(tmp_path)
    text = runs.help()
    assert "runs.<run_id>" in text
    assert "runs.match" in text
    assert "runs.pack" in text
    assert "runs.sweep" in text
    assert "runs.match.help()" in text
    assert "runs.pack.help()" in text
    assert "runs.sweep.help()" in text


def test_run_help_mentions_all_collection_paths_and_var(tmp_path: Path):
    _make_run(tmp_path / "a", "f5a_DBI30_I3")
    run = pyslimmc.scan(tmp_path).f5a_DBI30_I3
    text = run.help()
    assert "sl.open(path)" in text
    assert "runs.<run_id>" in text
    assert "runs.match" in text
    assert "runs.pack" in text
    assert 'run.var["name"].value' in text


def test_runs_info_reports_interactive_state_selection_and_collisions(tmp_path: Path):
    _make_run(tmp_path / "a", "case_A")
    _make_run(tmp_path / "b", "match")
    selected = pyslimmc.scan(tmp_path).match("*")
    text = selected.info_text(max_rows=0)
    assert "unique run_id: 2" in text
    assert "interactive attributes: 1" in text
    assert "selection: match('*')" in text
    assert "API-name collisions: match" in text
    assert "run-id index: lazy" in text


def test_info_does_not_eagerly_cache_interactive_index(tmp_path: Path):
    _make_run(tmp_path / "a", "case_A")
    runs = pyslimmc.scan(tmp_path)
    assert runs._run_id_index_cache is None
    runs.info_text(max_rows=0)
    assert runs._run_id_index_cache is None
