from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent))
import pyslimmc as sl

HERE = Path(__file__).resolve().parent


def output_for(model: Path) -> Path:
    return model.parent / "results" / model.stem


def action_types(run):
    return [a.type for a in run.actions]


def actions_of(run, type_name: str):
    return [a for a in run.actions if a.type == type_name]


def close(a, b, atol=5e-12):
    return abs(float(a) - float(b)) <= atol * max(1.0, abs(float(b)))


def run_success(model: Path, engine: Path):
    out = output_for(model)
    if out.exists():
        shutil.rmtree(out)
    subprocess.run([str(engine), str(model)], check=True)
    run = sl.open(out)
    assert run.is_complete, model.stem
    return run


def run_failure(model: Path, engine: Path, message_fragment: str):
    out = output_for(model)
    if out.exists():
        shutil.rmtree(out)
    proc = subprocess.run([str(engine), str(model)], text=True, capture_output=True)
    assert proc.returncode != 0, model.stem
    combined = proc.stdout + "\n" + proc.stderr
    assert message_fragment in combined, (model.stem, combined)
    metadata = json.loads((out / "run_metadata.json").read_text())
    assert metadata["run_status"] == "failed", metadata


def check_c11_no_t0(run):
    expected = np.array([0.0, 0.10, 0.20, 0.30, 0.35])
    assert np.allclose(run.t, expected, atol=1e-12), run.t
    saves = actions_of(run, "save")
    assert len(saves) == 4
    assert all(a.trigger == "every" for a in saves)
    assert np.allclose([a.t for a in saves], [0.0, 0.1, 0.2, 0.3])


def check_c11_coalesce(run):
    at_idx = np.flatnonzero(np.isclose(np.asarray(run.t), 0.10))
    assert at_idx.size == 1, run.t
    i = int(at_idx[0])
    assert bool(run.snapshots.raw["has_chains"][i])
    actions = list(run.actions)
    assert [a.type for a in actions] == ["save", "save_chains"]
    assert all(close(a.t, 0.10) for a in actions)
    assert actions[0].snapshot is not None
    assert actions[1].snapshot is not None
    assert actions[0].snapshot.index == actions[1].snapshot.index == i


def check_c11_tend(run):
    acts = actions_of(run, "set_c")
    assert len(acts) == 1
    a = acts[0]
    assert close(a.t, 0.30)
    assert close(a.before_value, 0.10, 2e-5)
    assert close(a.after_value, 0.25, 2e-5)
    assert close(run.conc["Q"][-1], 0.25, 2e-5)
    assert close(run.t[-1], 0.30)


def check_c12_rate_actions(run):
    acts = list(run.actions)
    assert [a.type for a in acts] == ["set_k", "add_k", "set_temp", "add_temp"]
    assert np.allclose([a.t for a in acts], [0.1, 0.2, 0.3, 0.4])
    assert close(acts[0].after_value, 2.5)
    assert close(acts[1].before_value, 2.5)
    assert close(acts[1].after_value, 4.0)
    assert close(acts[2].before_value, 300.0)
    assert close(acts[2].after_value, 320.0)
    assert close(acts[3].before_value, 320.0)
    assert close(acts[3].after_value, 325.0)
    assert [a.kinetic_parameter_set_id for a in acts] == [1, 2, 3, 4]
    assert all(a.snapshot is not None for a in acts)
    assert np.all(np.diff(np.asarray(run.snapshots.raw["kinetic_parameter_set_id"], dtype=int)) >= 0)
    # set_k converts the Arrhenius rate to fixed; later temperature changes do not alter it.
    assert close(run.last.k["k"], 4.0)
    assert close(run.last.temp, 325.0)


def check_c12_switch_off(run):
    a = actions_of(run, "set_k")[0]
    assert close(a.after_value, 0.0)
    name = run.channels.event_count.names[0]
    counts = np.asarray(run.channels.event_count[name], dtype=np.int64)
    before = int(np.flatnonzero(np.asarray(run.t) <= 0.20 + 1e-12)[-1])
    assert counts[before] > 0
    assert np.all(counts[before:] == counts[before]), (run.t, counts)


def check_c13(run):
    conc_actions = [a for a in run.actions if a.type in {"set_c", "add_c"}]
    assert len(conc_actions) == 6
    # Declaration order at equal times is preserved: Q change, save, M change.
    assert np.allclose([a.t for a in conc_actions], [0.1, 0.1, 0.2, 0.2, 0.3, 0.3])
    q = conc_actions[0::2]
    m = conc_actions[1::2]
    for seq in (q, m):
        assert close(seq[0].before_value, 0.10, 2e-5)
        assert close(seq[0].after_value, 0.20, 2e-5)
        assert close(seq[1].before_value, 0.20, 2e-5)
        assert close(seq[1].after_value, 0.30, 2e-5)
        assert close(seq[2].before_value, 0.30, 2e-5)
        assert close(seq[2].after_value, 0.25, 2e-5)
        assert all(a.kinetic_parameter_set_id is None for a in seq)
        assert all(a.snapshot is None for a in seq)
    assert close(run.conc["Q"][-1], 0.25, 2e-5)
    assert close(run.conc["M"][-1], 0.25, 2e-5)
    assert np.all(np.asarray(run.count["M"], dtype=np.int64) >= 0)
    assert np.all(np.asarray(run.snapshots.raw["kinetic_parameter_set_id"], dtype=int) == 0)


def check_c14_initial(run):
    acts = list(run.actions)
    assert [a.type for a in acts] == ["set_c", "set_c", "print"]
    assert all(a.trigger == "when" for a in acts)
    assert all(close(a.t, 0.0) for a in acts)
    assert len(acts[0].conditions) == 1
    assert len(acts[1].conditions) == 2
    assert acts[2].message == "cascade complete"
    assert close(run.conc["R"][-1], 0.10, 2e-5)
    assert close(run.conc["S"][-1], 0.20, 2e-5)


def check_c14_scheduled(run):
    acts = list(run.actions)
    assert [a.type for a in acts] == ["set_c", "set_c", "set_c", "print"]
    assert acts[0].trigger == "at"
    assert all(a.trigger == "when" for a in acts[1:])
    assert all(close(a.t, 0.10) for a in acts)
    assert acts[-1].message == "scheduled cascade complete"
    # Each line fires once even though the cascade scanner repeats.
    assert len(acts) == 4


def check_c14_stop(run):
    stop = actions_of(run, "stop")
    assert len(stop) == 1
    a = stop[0]
    assert a.trigger == "when"
    assert a.event % 10 == 0 and a.event > 0
    assert len(a.conditions) == 1
    c = a.conditions[0]
    assert c.met and c.observed_value < c.threshold
    assert run.t[-1] < 10.0
    assert run.status == "completed"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", type=Path, default=ROOT.parent / "bin" / "slimmc")
    args = parser.parse_args()
    engine = args.engine.resolve()
    models = HERE / "models"

    checks = {
        "C11_no_t0_save": check_c11_no_t0,
        "C11_same_time_coalesce": check_c11_coalesce,
        "C11_tend_action": check_c11_tend,
        "C12_rate_actions": check_c12_rate_actions,
        "C12_rate_switch_off": check_c12_switch_off,
        "C13_concentration_actions": check_c13,
        "C14_initial_and_independent": check_c14_initial,
        "C14_scheduled_cascade": check_c14_scheduled,
        "C14_stop_cadence": check_c14_stop,
    }

    for name, check in checks.items():
        run = run_success(models / f"{name}.model", engine)
        check(run)
        print(f"[PASS] {name}")

    run_failure(models / "C12_negative_rate_fail.model", engine, "add_k would make negative rate constant")
    print("[PASS] C12_negative_rate_fail")
    run_failure(models / "C13_negative_count_fail.model", engine, "add_c would make negative molecule count")
    print("[PASS] C13_negative_count_fail")
    print(f"[PASS] copo phase C integration ({len(checks) + 2} models)")


if __name__ == "__main__":
    main()
