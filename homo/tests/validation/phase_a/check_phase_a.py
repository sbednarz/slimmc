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


def channel_name(run, needle: str) -> str:
    names = [n for n in run.channels.event_count.names if needle in n]
    assert len(names) == 1, (needle, names)
    return names[0]


def events(run, needle: str, idx: int = -1) -> int:
    return int(run.channels.event_count[channel_name(run, needle)][idx])


def weighted_sum(values, counts) -> int:
    return int(np.sum(np.asarray(values, dtype=np.int64) * np.asarray(counts, dtype=np.int64)))


def total_units(chains) -> int:
    return weighted_sum(chains.dp, chains.count)


def total_chains(chains) -> int:
    return int(np.sum(np.asarray(chains.count, dtype=np.int64)))


def mean_dp(chains) -> float:
    n = total_chains(chains)
    return float(total_units(chains) / n) if n else float("nan")


def count_end(chains, end: str) -> int:
    return int(sum(int(n) for e, n in zip(chains.right_end, chains.count) if str(e) == end))


def phase_start_index(run, switch_time: float) -> int:
    idx = np.flatnonzero(np.asarray(run.t) <= switch_time + 1.0e-12)
    assert idx.size
    return int(idx[-1])


def assert_common(run, name: str) -> None:
    assert run.is_complete
    assert len(run.t) > 1
    assert np.all(np.diff(np.asarray(run.t)) >= 0.0)
    assert np.all(np.asarray(run.count["M"], dtype=np.int64) >= 0)

    totals = []
    has_chains = np.asarray(run.snapshots.raw["has_chains"], dtype=bool)
    for idx, present in enumerate(has_chains):
        if not present:
            continue
        snap = run.at_snapshot(idx)
        all_chains = snap.chains.all
        assert np.all(np.asarray(all_chains.dp, dtype=np.int64) >= 1)
        assert np.all(np.asarray(all_chains.count, dtype=np.int64) >= 1)
        totals.append(int(run.count["M"][idx]) + total_units(all_chains))
    assert totals and len(set(totals)) == 1, (name, totals)


def run_model(model: Path, engine: Path):
    out = model.parent / "results" / model.stem
    if out.exists():
        shutil.rmtree(out)
    subprocess.run([str(engine), str(model)], check=True)
    run = sl.open(out)
    assert_common(run, model.stem)
    return run


def check_single(name: str, run) -> None:
    if name == "H01_prop":
        n_init = events(run, "macro_init")
        n_prop = events(run, "macro_prop")
        assert n_init > 0 and n_prop > 0
        assert total_chains(run.last.chains.live) == n_init
        assert total_chains(run.last.chains.dead) == 0
        assert total_units(run.last.chains.live) == n_init + n_prop

    elif name == "H01_prop_zero":
        assert events(run, "macro_prop") == 0
        assert set(np.asarray(run.last.chains.live.dp, dtype=int)) == {1}
        assert total_units(run.last.chains.live) == events(run, "macro_init")

    elif name == "H01_prop_dpmax":
        assert events(run, "macro_prop") > 0
        assert np.max(np.asarray(run.last.chains.live.dp, dtype=int)) == 4
        assert np.min(np.asarray(run.last.chains.live.dp, dtype=int)) == 4
        assert total_units(run.last.chains.live) == events(run, "macro_init") + events(run, "macro_prop")

    elif name == "H02_deprop":
        i0 = phase_start_index(run, 0.20)
        delta_e = events(run, "macro_deprop") - events(run, "macro_deprop", i0)
        delta_m = int(run.count["M"][-1]) - int(run.count["M"][i0])
        unit_loss = total_units(run.at_snapshot(i0).chains.live) - total_units(run.last.chains.live)
        assert delta_e > 0
        assert delta_e == delta_m == unit_loss
        assert total_chains(run.last.chains.live) == total_chains(run.at_snapshot(i0).chains.live)

    elif name == "H02_D02_control":
        i0 = phase_start_index(run, 0.20)
        assert events(run, "macro_deprop") == 0
        assert int(run.count["M"][-1]) == int(run.count["M"][i0])
        assert total_units(run.last.chains.live) == total_units(run.at_snapshot(i0).chains.live)

    elif name.startswith("H02_D03_kdp_"):
        i0 = phase_start_index(run, 0.20)
        delta_e = events(run, "macro_deprop") - events(run, "macro_deprop", i0)
        delta_m = int(run.count["M"][-1]) - int(run.count["M"][i0])
        assert delta_e == delta_m

    elif name.startswith("H02_D04_"):
        i0 = phase_start_index(run, 0.20)
        p = events(run, "macro_prop") - events(run, "macro_prop", i0)
        d = events(run, "macro_deprop") - events(run, "macro_deprop", i0)
        dp0 = mean_dp(run.at_snapshot(i0).chains.live)
        dpf = mean_dp(run.last.chains.live)
        if "deprop_dominant" in name:
            assert d > p and dpf < dp0
        elif "prop_dominant" in name:
            assert p > d and dpf > dp0
        else:
            assert p + d > 0
            assert abs(p - d) / (p + d) < 0.12, (p, d)
            assert abs(dpf - dp0) < 5.0, (dp0, dpf)

    elif name == "H03_term_c":
        i0 = phase_start_index(run, 0.15)
        e = events(run, "macro_term_c") - events(run, "macro_term_c", i0)
        live0 = total_chains(run.at_snapshot(i0).chains.live)
        dead0 = total_chains(run.at_snapshot(i0).chains.dead)
        assert e > 0
        assert live0 - total_chains(run.last.chains.live) == 2 * e
        assert total_chains(run.last.chains.dead) - dead0 == e
        assert int(run.count["M"][-1]) == int(run.count["M"][i0])

    elif name == "H03_term_c_control":
        i0 = phase_start_index(run, 0.15)
        assert events(run, "macro_term_c") == 0
        assert total_chains(run.last.chains.live) == total_chains(run.at_snapshot(i0).chains.live)
        assert total_chains(run.last.chains.dead) == 0

    elif name == "H04_term_d":
        i0 = phase_start_index(run, 0.15)
        e = events(run, "macro_term_d") - events(run, "macro_term_d", i0)
        live0 = total_chains(run.at_snapshot(i0).chains.live)
        assert e > 0
        assert live0 - total_chains(run.last.chains.live) == 2 * e
        assert total_chains(run.last.chains.dead) == 2 * e
        assert count_end(run.last.chains.dead, "H") == e
        assert count_end(run.last.chains.dead, "U") == e
        assert int(run.count["M"][-1]) == int(run.count["M"][i0])

    elif name == "H04_term_d_control":
        i0 = phase_start_index(run, 0.15)
        assert events(run, "macro_term_d") == 0
        assert total_chains(run.last.chains.live) == total_chains(run.at_snapshot(i0).chains.live)
        assert total_chains(run.last.chains.dead) == 0

    elif name == "H05_transfer_m":
        i0 = phase_start_index(run, 0.15)
        e = events(run, "macro_transfer_m") - events(run, "macro_transfer_m", i0)
        live0 = total_chains(run.at_snapshot(i0).chains.live)
        dead0 = total_chains(run.at_snapshot(i0).chains.dead)
        assert e > 0
        assert total_chains(run.last.chains.live) == live0
        assert total_chains(run.last.chains.dead) - dead0 == e
        assert int(run.count["M"][i0]) - int(run.count["M"][-1]) == e
        assert total_chains(run.last.chains.all) - total_chains(run.at_snapshot(i0).chains.all) == e
        assert set(np.asarray(run.last.chains.live.dp, dtype=int)) == {1}
        assert count_end(run.last.chains.dead, "H") == e

    elif name == "H05_transfer_m_control":
        i0 = phase_start_index(run, 0.15)
        assert events(run, "macro_transfer_m") == 0
        assert total_chains(run.last.chains.live) == total_chains(run.at_snapshot(i0).chains.live)
        assert total_chains(run.last.chains.dead) == 0
        assert int(run.count["M"][-1]) == int(run.count["M"][i0])


def summary_for_cross_model(name: str, run) -> dict[str, float | int]:
    if name.startswith("H02_D03_kdp_"):
        return {"deprop_events": events(run, "macro_deprop")}
    pairs = {
        "H02_deprop": 0.20,
        "H02_D02_control": 0.20,
        "H03_term_c": 0.15,
        "H03_term_c_control": 0.15,
        "H04_term_d": 0.15,
        "H04_term_d_control": 0.15,
        "H05_transfer_m": 0.15,
        "H05_transfer_m_control": 0.15,
    }
    if name in pairs:
        i = phase_start_index(run, pairs[name])
        return {
            "M": int(run.count["M"][i]),
            "units": total_units(run.at_snapshot(i).chains.all),
            "live": total_chains(run.at_snapshot(i).chains.live),
        }
    return {}


def check_cross_model(summaries: dict[str, dict[str, float | int]]) -> None:
    counts = [int(summaries[f"H02_D03_kdp_{k}"]["deprop_events"]) for k in (10, 25, 50)]
    assert counts[0] < counts[1] < counts[2], counts
    assert counts[1] / counts[0] > 1.3
    assert counts[2] / counts[1] > 1.2

    for active, control in [
        ("H02_deprop", "H02_D02_control"),
        ("H03_term_c", "H03_term_c_control"),
        ("H04_term_d", "H04_term_d_control"),
        ("H05_transfer_m", "H05_transfer_m_control"),
    ]:
        assert summaries[active] == summaries[control], (active, control, summaries[active], summaries[control])

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", type=Path, default=ROOT.parent / "bin" / "slimmc")
    args = parser.parse_args()

    models = sorted((Path(__file__).parent / "models").glob("H*.model"))
    import gc
    summaries = {}
    for model in models:
        run = run_model(model, args.engine.resolve())
        check_single(model.stem, run)
        summaries[model.stem] = summary_for_cross_model(model.stem, run)
        del run
        gc.collect()
        print(f"[PASS] {model.stem}")

    check_cross_model(summaries)
    print(f"[PASS] detailed phase A integration ({len(models)} models)")


if __name__ == "__main__":
    main()
