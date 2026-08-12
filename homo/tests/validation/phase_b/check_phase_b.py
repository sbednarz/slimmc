from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT.parent))
import pyslimmc as sl


def output_for(model: Path) -> Path:
    return model.parent / "results" / model.stem


def run_model(engine: Path, model: Path, trace_limit: int) -> sl.StorageRun:
    out = output_for(model)
    if out.exists():
        shutil.rmtree(out)
    subprocess.run(
        [str(engine), "--trace-channels", str(trace_limit), str(model)],
        check=True,
    )
    run = sl.open(out)
    assert run.is_complete or run.stop_reason == "max_steps"
    trace = run.diagnostics.channel_trace
    assert trace.enabled
    assert trace.complete
    assert not trace.truncated
    assert len(trace) > 1000
    return run


def check_h08_h09(engine: Path, model: Path) -> None:
    run = run_model(engine, model, 9000)
    trace = run.diagnostics.channel_trace
    names = np.asarray(trace.channel, dtype=object)

    n1 = int(np.count_nonzero(names == "rxn_A_B_k1"))
    n2 = int(np.count_nonzero(names == "rxn_A_C_k2"))
    total = n1 + n2
    assert total == len(trace), (n1, n2, len(trace))

    observed = n1 / total
    expected = 100.0 / (100.0 + 300.0)
    # Binomial 6-sigma envelope plus a small numerical guard.
    sigma = np.sqrt(expected * (1.0 - expected) / total)
    assert abs(observed - expected) <= 6.0 * sigma + 0.002, (observed, expected, sigma)

    dt = np.asarray(trace.dt, dtype=float)
    a0 = np.asarray(trace.total_propensity, dtype=float)
    z = dt * a0
    assert np.all(np.isfinite(z)) and np.all(z > 0.0)

    # For exact SSA, z=a0*dt follows Exp(1), even while a0 changes by event.
    mean = float(np.mean(z))
    variance = float(np.var(z, ddof=1))
    second_moment = float(np.mean(z * z))
    assert abs(mean - 1.0) < 0.045, mean
    assert abs(variance - 1.0) < 0.10, variance
    assert abs(second_moment - 2.0) < 0.14, second_moment

    # Coarse empirical CDF checks avoid depending on scipy.
    for x, expected_cdf in [(0.25, 1 - np.exp(-0.25)), (1.0, 1 - np.exp(-1.0)), (2.0, 1 - np.exp(-2.0))]:
        observed_cdf = float(np.mean(z <= x))
        se = np.sqrt(expected_cdf * (1.0 - expected_cdf) / len(z))
        assert abs(observed_cdf - expected_cdf) <= 6.0 * se + 0.006, (x, observed_cdf, expected_cdf)

    print(f"[PASS] H08 channel selection: n1={n1}, n2={n2}, share={observed:.5f}")
    print(f"[PASS] H09 SSA time law: mean={mean:.5f}, var={variance:.5f}, E[z^2]={second_moment:.5f}")


def snapshot(run: sl.StorageRun) -> dict[str, np.ndarray]:
    trace = run.diagnostics.channel_trace
    return {
        "event": np.asarray(trace.kmc_event).copy(),
        "time": np.asarray(trace.t).copy(),
        "dt": np.asarray(trace.dt).copy(),
        "channel_id": np.asarray(trace.channel_id).copy(),
        "propensity": np.asarray(trace.propensity).copy(),
        "a0": np.asarray(trace.total_propensity).copy(),
        "A": np.asarray(run.count["A"]).copy(),
        "B": np.asarray(run.count["B"]).copy(),
        "C": np.asarray(run.count["C"]).copy(),
        "snapshot_t": np.asarray(run.t).copy(),
    }


def check_h10(engine: Path, model: Path) -> None:
    first = snapshot(run_model(engine, model, 4000))
    second = snapshot(run_model(engine, model, 4000))
    for key in first:
        np.testing.assert_array_equal(first[key], second[key], err_msg=f"fixed-seed mismatch in {key}")
    print(f"[PASS] H10 fixed seed reproducibility ({len(first['event'])} traced events)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", type=Path, default=ROOT.parent / "bin" / "slimmc")
    args = parser.parse_args()
    engine = args.engine.resolve()
    models = Path(__file__).parent / "models"
    check_h08_h09(engine, models / "H08_H09_competing_channels.model")
    check_h10(engine, models / "H10_seed_reproducibility.model")
    print("[PASS] phase B integration")


if __name__ == "__main__":
    main()
