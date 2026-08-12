from __future__ import annotations

import argparse
import hashlib
import json
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]  # copo/
REPO = ROOT.parent
sys.path.insert(0, str(REPO))
import pyslimmc as sl

MODELS = Path(__file__).resolve().parent / "models"
REGRESSION_MODEL = REPO / "tests/validation/copo/engine/binary_terminal_seeded/COP_REG_BINARY01.model"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metadata(run_dir: Path) -> dict:
    return json.loads((run_dir / "run_metadata.json").read_text(encoding="utf-8"))


def run_model(engine: Path, model: Path, work: Path):
    work.mkdir(parents=True, exist_ok=True)
    local = work / model.name
    local.write_text(model.read_text(encoding="utf-8"), encoding="utf-8")
    cp = subprocess.run([str(engine), str(local)], capture_output=True, text=True)
    return cp, work / "results" / model.stem


def assert_complete_storage(run_dir: Path) -> dict:
    md = metadata(run_dir)
    assert md["run_status"] == "completed", md
    assert md["exit_code"] == 0
    assert (run_dir / "RESULTS_COMPLETE").is_file()
    assert not (run_dir / ".work").exists()
    assert md["input_model_sha256"] == sha256(run_dir / "input.model")
    assert md["schema_sha256"] == sha256(run_dir / "schema.jsonl")
    finals = np.load(run_dir / "snapshots/is_final.npy", allow_pickle=False)
    assert np.count_nonzero(finals) == 1
    return md


def check_c20(engine: Path, work: Path) -> None:
    runs = {}
    for stem in ("C20_mass_repeat_units", "C20_mass_with_end_groups"):
        cp, run_dir = run_model(engine, MODELS / f"{stem}.model", work / stem)
        assert cp.returncode == 0, (cp.stdout, cp.stderr)
        assert_complete_storage(run_dir)
        run = sl.open(run_dir)
        chains = run.last.chains.all
        assert len(chains.dp) > 0
        runs[stem] = chains

    repeat = runs["C20_mass_repeat_units"]
    with_ends = runs["C20_mass_with_end_groups"]
    expected_repeat = np.asarray(repeat.counts["A"], float) * 100.0 + np.asarray(repeat.counts["B"], float) * 128.0
    np.testing.assert_allclose(repeat.molar_mass, expected_repeat, rtol=0, atol=1e-10)
    np.testing.assert_allclose(with_ends.molar_mass, expected_repeat + 70.0, rtol=0, atol=1e-10)
    np.testing.assert_array_equal(repeat.dp, with_ends.dp)
    np.testing.assert_array_equal(repeat.count, with_ends.count)
    np.testing.assert_array_equal(repeat.counts["A"], with_ends.counts["A"])
    np.testing.assert_array_equal(repeat.counts["B"], with_ends.counts["B"])
    np.testing.assert_allclose(with_ends.molar_mass - repeat.molar_mass, 70.0, rtol=0, atol=1e-10)


def check_c21(engine: Path, work: Path) -> None:
    cp, run_dir = run_model(engine, MODELS / "C21_chain_aggregation.model", work / "aggregation")
    assert cp.returncode == 0, (cp.stdout, cp.stderr)
    assert_complete_storage(run_dir)
    run = sl.open(run_dir)
    chains = run.last.chains.all
    # Copo currently writes one physical row per chain, unlike homo's structural
    # aggregation. Validate that this explicit-row representation is lossless.
    assert len(chains.dp) > 1
    np.testing.assert_array_equal(chains.dp, np.ones(len(chains.dp), dtype=chains.dp.dtype))
    np.testing.assert_array_equal(chains.count, np.ones(len(chains.count), dtype=chains.count.dtype))
    np.testing.assert_array_equal(chains.counts["A"], np.ones(len(chains.dp), dtype=chains.counts["A"].dtype))
    np.testing.assert_array_equal(chains.counts["B"], np.zeros(len(chains.dp), dtype=chains.counts["B"].dtype))
    assert set(map(str, chains.left_end)) == {"R"}
    assert set(map(str, chains.right_end)) == {"ACTIVE"}

    raw = run_dir / "chains"
    ids = np.load(raw / "chain_record_id.npy", allow_pickle=False)
    counts = np.load(raw / "count.npy", allow_pickle=False)
    assert np.array_equal(ids, np.arange(len(ids), dtype=np.uint64))
    assert np.all(counts == 1)
    assert int(np.sum(counts, dtype=np.uint64)) == int(chains.total_chains)

    # The structural key is deliberately repeated in current copo storage;
    # chain_record_id still remains unique and composition rows remain dense.
    comp = run_dir / "chain_composition"
    comp_ids = np.load(comp / "chain_record_id.npy", allow_pickle=False)
    mids = np.load(comp / "monomer_id.npy", allow_pickle=False)
    units = np.load(comp / "unit_count.npy", allow_pickle=False)
    assert len(comp_ids) == 2 * len(ids)
    for rid in ids:
        rows = np.flatnonzero(comp_ids == rid)
        assert rows.size == 2
        assert set(mids[rows].tolist()) == {0, 1}
        by_mid = {int(m): int(u) for m, u in zip(mids[rows], units[rows])}
        assert by_mid == {0: 1, 1: 0}


def check_c22(engine: Path, work: Path) -> None:
    cp, run_dir = run_model(engine, MODELS / "C22_mixed_snapshots.model", work / "mixed")
    assert cp.returncode == 0, (cp.stdout, cp.stderr)
    assert_complete_storage(run_dir)
    run = sl.open(run_dir)
    t = np.asarray(run.t, dtype=float)
    has_chains = np.asarray(run.snapshots.raw["has_chains"], dtype=bool)
    assert t[0] == 0.0, t
    assert np.all(np.diff(t) > 0.0), t
    assert np.count_nonzero(has_chains) == 2, has_chains
    np.testing.assert_allclose(t[has_chains], [0.30, 0.50], rtol=0, atol=1e-12)
    for monomer in ("A", "B"):
        conc = np.asarray(run.conc[monomer], float)
        conv = np.asarray(run.conv[monomer], float)
        fcum = np.asarray(run.F.cum[monomer], float)
        assert len(conc) == len(conv) == len(fcum) == len(t)
        assert np.all(np.isfinite(conc)), conc
        assert np.all(np.isfinite(conv)), conv
        assert np.all(np.isfinite(fcum)), fcum
    np.testing.assert_allclose(np.asarray(run.F.cum["A"], float) + np.asarray(run.F.cum["B"], float), 1.0, rtol=0, atol=1e-12)

    cp, run_dir = run_model(engine, MODELS / "C22_same_time_snapshot.model", work / "same_time")
    assert cp.returncode == 0, (cp.stdout, cp.stderr)
    assert_complete_storage(run_dir)
    run = sl.open(run_dir)
    t = np.asarray(run.t, dtype=float)
    has_chains = np.asarray(run.snapshots.raw["has_chains"], dtype=bool)
    at_020 = np.flatnonzero(np.isclose(t, 0.20, rtol=0, atol=1e-12))
    assert at_020.size == 1, t
    assert has_chains[int(at_020[0])]
    assert len(np.unique(t)) == len(t), t


def check_c23_static(engine: Path, work: Path) -> None:
    for stem in ("C23_a0_zero", "C23_stop"):
        cp, run_dir = run_model(engine, MODELS / f"{stem}.model", work / stem)
        assert cp.returncode == 0, (stem, cp.stdout, cp.stderr)
        md = assert_complete_storage(run_dir)
        run = sl.open(run_dir)
        assert run.is_complete and md["run_status"] == "completed"
        assert np.count_nonzero(np.asarray(run.snapshots.raw["is_final"], bool)) == 1
        assert np.all(np.diff(np.asarray(run.t, float)) >= 0.0)


def check_c23_failed(engine: Path, work: Path) -> None:
    # Fail immediately. A scheduled failure at t=0.10 is not reliable for a
    # deliberately step-limited regression model, which may stop earlier.
    text = REGRESSION_MODEL.read_text(encoding="utf-8").replace(
        "rate ki_a 2.0e6", "rate ki_a 1.0e308"
    )
    model = work / "C23_failed.model"
    model.write_text(text, encoding="utf-8")
    cp = subprocess.run([str(engine), str(model)], capture_output=True, text=True)
    assert cp.returncode != 0, (cp.stdout, cp.stderr)
    run_dir = work / "results/C23_failed"
    md = metadata(run_dir)
    assert md["run_status"] == "failed" and md["exit_code"] == 1, md
    assert not (run_dir / "RESULTS_COMPLETE").exists()
    assert (run_dir / ".work").is_dir()
    final_path = run_dir / "snapshots/is_final.npy"
    if final_path.exists():
        finals = np.load(final_path, allow_pickle=False)
        assert np.count_nonzero(finals) == 0
    assert md["input_model_sha256"] == sha256(run_dir / "input.model")
    assert md["schema_sha256"] == sha256(run_dir / "schema.jsonl")


def check_c23_interrupted(engine: Path, work: Path) -> None:
    text = REGRESSION_MODEL.read_text(encoding="utf-8")
    text = text.replace("param t_end 0.5", "param t_end 1.0e12")
    text = text.replace("param max_steps 20000", "param max_steps 2000000000")
    text = text.replace("every 0.1 save", "every 1000000 save")
    text = text.replace("every 0.25 save_chains", "every 1000000 save_chains")
    model = work / "C23_interrupted.model"
    model.write_text(text, encoding="utf-8")
    proc = subprocess.Popen([str(engine), str(model)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    run_dir = work / "results/C23_interrupted"
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
    finals = np.load(run_dir / "snapshots/is_final.npy", allow_pickle=False)
    assert np.count_nonzero(finals) == 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", type=Path, required=True)
    ns = ap.parse_args()
    engine = ns.engine.resolve()
    assert engine.is_file(), engine
    with tempfile.TemporaryDirectory(prefix="slimmc_copo_phase_e_") as raw:
        work = Path(raw)
        check_c20(engine, work / "C20")
        check_c21(engine, work / "C21")
        check_c22(engine, work / "C22")
        check_c23_static(engine, work / "C23")
        check_c23_failed(engine, work / "C23")
        check_c23_interrupted(engine, work / "C23")
    print("Copo phase E: mass/storage/snapshot/finalization validation: PASS")


if __name__ == "__main__":
    main()
