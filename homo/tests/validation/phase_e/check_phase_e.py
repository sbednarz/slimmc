from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent))
import pyslimmc as sl

MODELS = Path(__file__).resolve().parent / "models"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metadata(run_dir: Path) -> dict:
    return json.loads((run_dir / "run_metadata.json").read_text(encoding="utf-8"))


def run_model(engine: Path, model: Path, work: Path):
    work.mkdir(parents=True, exist_ok=True)
    local = work / model.name
    local.write_text(model.read_text(encoding="utf-8"), encoding="utf-8")
    cp = subprocess.run([str(engine), str(local)], capture_output=True, text=True)
    run_dir = work / "results" / model.stem
    return cp, run_dir


def assert_complete_storage(run_dir: Path) -> dict:
    md = metadata(run_dir)
    assert md["run_status"] == "completed", md
    assert md["exit_code"] == 0
    assert (run_dir / "RESULTS_COMPLETE").is_file()
    assert not (run_dir / ".work").exists()
    assert md["input_model_sha256"] == sha256(run_dir / "input.model")
    assert md["schema_sha256"] == sha256(run_dir / "schema.jsonl")
    finals = np.load(run_dir / "snapshots" / "is_final.npy", allow_pickle=False)
    assert np.count_nonzero(finals) == 1
    return md


def check_h20(engine: Path, work: Path) -> None:
    runs = {}
    for stem in ("H20_mass_repeat_units", "H20_mass_with_end_groups"):
        cp, run_dir = run_model(engine, MODELS / f"{stem}.model", work / stem)
        assert cp.returncode == 0, (cp.stdout, cp.stderr)
        assert_complete_storage(run_dir)
        run = sl.open(run_dir)
        assert run.is_complete
        chains = run.last.chains.all
        assert len(chains.dp) > 0
        runs[stem] = chains

    repeat = runs["H20_mass_repeat_units"]
    with_ends = runs["H20_mass_with_end_groups"]
    np.testing.assert_allclose(repeat.molar_mass, np.asarray(repeat.dp, float) * 100.0, rtol=0, atol=1e-10)
    np.testing.assert_allclose(with_ends.molar_mass, np.asarray(with_ends.dp, float) * 100.0 + 70.0, rtol=0, atol=1e-10)

    # Identical seed and chemistry: only the selected mass model may differ.
    np.testing.assert_array_equal(repeat.dp, with_ends.dp)
    np.testing.assert_array_equal(repeat.count, with_ends.count)
    np.testing.assert_allclose(with_ends.molar_mass - repeat.molar_mass, 70.0, rtol=0, atol=1e-10)


def check_h21(engine: Path, work: Path) -> None:
    stem = "H21_chain_aggregation"
    cp, run_dir = run_model(engine, MODELS / f"{stem}.model", work / stem)
    assert cp.returncode == 0, (cp.stdout, cp.stderr)
    assert_complete_storage(run_dir)
    run = sl.open(run_dir)
    chains = run.last.chains.all
    assert len(chains.dp) == 1, (chains.dp, chains.count)
    assert int(chains.dp[0]) == 1
    assert int(chains.count[0]) > 1
    assert str(chains.left_end[0]) == "R"
    assert str(chains.right_end[0]) == "ACTIVE"

    raw = run_dir / "chains"
    ids = np.load(raw / "chain_record_id.npy", allow_pickle=False)
    counts = np.load(raw / "count.npy", allow_pickle=False)
    assert np.array_equal(ids, np.arange(len(ids), dtype=np.uint64))
    assert np.sum(counts, dtype=np.uint64) == int(chains.count[0])

    # No duplicate aggregation key inside any snapshot.
    arrays = {name: np.load(raw / f"{name}.npy", allow_pickle=False) for name in (
        "snapshot_id", "population_id", "pool_id", "dp", "left_end_id", "right_end_id", "origin_id"
    )}
    for sid in np.unique(arrays["snapshot_id"]):
        rows = np.flatnonzero(arrays["snapshot_id"] == sid)
        keys = list(zip(*(arrays[name][rows].tolist() for name in (
            "population_id", "pool_id", "dp", "left_end_id", "right_end_id", "origin_id"
        ))))
        assert len(keys) == len(set(keys)), (sid, keys)


def check_h22(engine: Path, work: Path) -> None:
    cp, run_dir = run_model(engine, MODELS / "H22_mixed_snapshots.model", work / "mixed")
    assert cp.returncode == 0, (cp.stdout, cp.stderr)
    assert_complete_storage(run_dir)
    run = sl.open(run_dir)
    t = np.asarray(run.t, dtype=float)
    has_chains = np.asarray(run.snapshots.raw["has_chains"], dtype=bool)
    assert t[0] == 0.0, t
    assert np.all(np.diff(t) > 0.0), t
    assert np.count_nonzero(has_chains) == 2, has_chains
    chain_times = t[has_chains]
    np.testing.assert_allclose(chain_times, [0.30, 0.50], rtol=0, atol=1e-12)
    assert len(run.conc["M"]) == len(t)
    assert len(run.conv["M"]) == len(t)
    assert np.all(np.isfinite(np.asarray(run.conc["M"], float)))
    assert np.all(np.isfinite(np.asarray(run.conv["M"], float)))
    fcum = np.asarray(run.F.cum["M"], dtype=float)
    assert len(fcum) == len(t)
    assert np.all(np.isfinite(fcum)), fcum
    np.testing.assert_allclose(fcum, 1.0, rtol=0, atol=1e-12)

    cp, run_dir = run_model(engine, MODELS / "H22_same_time_snapshot.model", work / "same_time")
    assert cp.returncode == 0, (cp.stdout, cp.stderr)
    assert_complete_storage(run_dir)
    run = sl.open(run_dir)
    t = np.asarray(run.t, dtype=float)
    has_chains = np.asarray(run.snapshots.raw["has_chains"], dtype=bool)
    at_020 = np.flatnonzero(np.isclose(t, 0.20, rtol=0, atol=1e-12))
    assert at_020.size == 1, t
    assert has_chains[int(at_020[0])]
    assert len(np.unique(t)) == len(t), t


def check_h23_static(engine: Path, work: Path) -> None:
    for stem in ("H23_a0_zero", "H23_stop"):
        cp, run_dir = run_model(engine, MODELS / f"{stem}.model", work / stem)
        assert cp.returncode == 0, (stem, cp.stdout, cp.stderr)
        md = assert_complete_storage(run_dir)
        run = sl.open(run_dir)
        assert run.is_complete
        assert md["run_status"] == "completed"
        assert np.count_nonzero(np.asarray(run.snapshots.raw["is_final"], bool)) == 1
        assert np.all(np.diff(np.asarray(run.t, float)) >= 0.0)


def check_h23_failed(engine: Path, work: Path) -> None:
    base = (ROOT / "tests/regression/frp_all_channels_seeded/FRP_ALLCHAN01.model").read_text(encoding="utf-8")
    text = base.replace("rate kd const 1.0e-3", "rate kd const 1.0e308")
    model = work / "H23_failed.model"
    model.write_text(text, encoding="utf-8")
    cp = subprocess.run([str(engine), str(model)], capture_output=True, text=True)
    assert cp.returncode != 0, (cp.stdout, cp.stderr)
    run_dir = work / "results" / "H23_failed"
    md = metadata(run_dir)
    assert md["run_status"] == "failed" and md["exit_code"] == 1, md
    assert not (run_dir / "RESULTS_COMPLETE").exists()
    assert (run_dir / ".work").is_dir()
    finals = np.load(run_dir / "snapshots" / "is_final.npy", allow_pickle=False)
    assert np.count_nonzero(finals) == 0
    assert md["input_model_sha256"] == sha256(run_dir / "input.model")
    assert md["schema_sha256"] == sha256(run_dir / "schema.jsonl")


def check_h23_interrupted(engine: Path, work: Path) -> None:
    text = (ROOT / "tests/regression/frp_all_channels_seeded/FRP_ALLCHAN01.model").read_text(encoding="utf-8")
    text = text.replace("param t_end 4.0", "param t_end 1.0e12")
    text = text.replace("param max_steps 300000", "param max_steps 2000000000")
    text = text.replace("every 1.0 save", "every 1000000 save")
    text = text.replace("every 1.0 save_chains", "every 1000000 save_chains")
    model = work / "H23_interrupted.model"
    model.write_text(text, encoding="utf-8")
    proc = subprocess.Popen([str(engine), str(model)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    run_dir = work / "results" / "H23_interrupted"
    deadline = time.time() + 8.0
    while time.time() < deadline and not (run_dir / "run_metadata.json").exists() and proc.poll() is None:
        time.sleep(0.01)
    assert proc.poll() is None, "run ended before SIGINT could be tested"
    proc.send_signal(signal.SIGINT)
    out, err = proc.communicate(timeout=20)
    assert proc.returncode == 0, (proc.returncode, out, err)
    md = metadata(run_dir)
    assert md["run_status"] == "interrupted" and md["exit_code"] == 130, md
    assert not (run_dir / "RESULTS_COMPLETE").exists()
    assert (run_dir / ".work").is_dir()
    finals = np.load(run_dir / "snapshots" / "is_final.npy", allow_pickle=False)
    assert np.count_nonzero(finals) == 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", type=Path, required=True)
    args = parser.parse_args()
    engine = args.engine.resolve()
    assert engine.is_file(), engine

    with tempfile.TemporaryDirectory(prefix="slimmc_phase_e_") as raw:
        work = Path(raw)
        for name in ("H20", "H21", "H22", "H23"):
            (work / name).mkdir(parents=True)
        check_h20(engine, work / "H20")
        check_h21(engine, work / "H21")
        check_h22(engine, work / "H22")
        check_h23_static(engine, work / "H23")
        check_h23_failed(engine, work / "H23")
        check_h23_interrupted(engine, work / "H23")

    print("Homo phase E: mass/storage/snapshot/finalization validation: PASS")


if __name__ == "__main__":
    main()
