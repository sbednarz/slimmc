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


def channel(run, needle: str | None = None) -> str:
    names = list(run.channels.event_count.names)
    if needle is None:
        assert len(names) == 1, names
        return names[0]
    hits = [n for n in names if needle in n]
    assert len(hits) == 1, (needle, names)
    return hits[0]


def ev(run, needle: str | None = None) -> int:
    return int(run.channels.event_count[channel(run, needle)][-1])


def productive(run, needle: str | None = None) -> int:
    return int(run.channels.productive_event_count[channel(run, needle)][-1])


def nonproductive(run, needle: str | None = None) -> int:
    return int(run.channels.nonproductive_event_count[channel(run, needle)][-1])


def total_chains(chains) -> int:
    return int(np.sum(np.asarray(chains.count, dtype=np.int64)))


def count_end(chains, end: str) -> int:
    return int(sum(int(n) for e, n in zip(chains.right_end, chains.count) if str(e) == end))


def phase_index(run, t: float) -> int:
    idx = np.flatnonzero(np.asarray(run.t) <= t + 1e-12)
    assert idx.size
    return int(idx[-1])


def run_model(model: Path, engine: Path):
    out = model.parent / "results" / model.stem
    if out.exists():
        shutil.rmtree(out)

    cp = subprocess.run(
        [str(engine), str(model)],
        capture_output=True,
        text=True,
    )
    if cp.returncode != 0:
        raise RuntimeError(
            f"engine failed for {model.name} with exit code {cp.returncode}\n"
            f"stdout:\n{cp.stdout}\nstderr:\n{cp.stderr}"
        )

    # H19 intentionally ends at max_steps after exactly 500 efficiency trials.
    # Suppress only that expected engine warning; preserve every other warning.
    expected_max_steps = model.stem.startswith("H19_eff_")
    stderr_lines = []
    for line in cp.stderr.splitlines():
        if expected_max_steps and line.strip() == "WARNING: max_steps reached before t_end":
            continue
        stderr_lines.append(line)
    if stderr_lines:
        print("\n".join(stderr_lines), file=sys.stderr)

    run = sl.open(out)
    assert run.is_complete
    assert np.all(np.diff(np.asarray(run.t)) >= 0)
    return run


def check(name: str, run) -> None:
    if name == "H15_init":
        e = ev(run, "macro_init")
        assert e > 0
        ce = np.asarray(run.channels.event_count[channel(run, "macro_init")], dtype=np.int64)
        assert np.all(np.asarray(run.count["R"], dtype=np.int64) + ce == int(run.count["R"][0]) + int(ce[0]))
        assert np.all(np.asarray(run.count["M"], dtype=np.int64) + ce == int(run.count["M"][0]) + int(ce[0]))
        assert total_chains(run.last.chains.live) == e
        assert set(np.asarray(run.last.chains.live.dp, dtype=int)) == {1}
    elif name == "H15_init_control":
        assert ev(run, "macro_init") == 0
        assert total_chains(run.last.chains.live) == 0

    elif name == "H16_transfer_h":
        i = phase_index(run, 0.15)
        e = ev(run, "macro_transfer") - int(run.channels.event_count[channel(run,"macro_transfer")][i])
        assert e > 0
        assert int(run.count["CTA"][i]) - int(run.count["CTA"][-1]) == e
        assert int(run.count["Rcta"][-1]) - int(run.count["Rcta"][i]) == e
        assert total_chains(run.at_snapshot(i).chains.live) - total_chains(run.last.chains.live) == e
        assert total_chains(run.last.chains.dead) - total_chains(run.at_snapshot(i).chains.dead) == e
        assert count_end(run.last.chains.dead, "H") == e
    elif name == "H16_transfer_h_control":
        i = phase_index(run, 0.15)
        assert ev(run, "macro_transfer") == 0
        assert total_chains(run.last.chains.live) == total_chains(run.at_snapshot(i).chains.live)

    elif name == "H17_term_x":
        i = phase_index(run, 0.15)
        e = ev(run, "macro_term_x") - int(run.channels.event_count[channel(run,"macro_term_x")][i])
        assert e > 0
        assert int(run.count["CAP"][i]) - int(run.count["CAP"][-1]) == e
        assert total_chains(run.at_snapshot(i).chains.live) - total_chains(run.last.chains.live) == e
        assert total_chains(run.last.chains.dead) - total_chains(run.at_snapshot(i).chains.dead) == e
        assert count_end(run.last.chains.dead, "CAP") == e
    elif name == "H17_term_x_control":
        i = phase_index(run, 0.15)
        assert ev(run, "macro_term_x") == 0
        assert total_chains(run.last.chains.live) == total_chains(run.at_snapshot(i).chains.live)

    elif name == "H18_rxn_uni":
        e = ev(run)
        assert e > 0
        ce = np.asarray(run.channels.event_count[channel(run)], dtype=np.int64)
        assert np.all(np.asarray(run.count["A"], dtype=np.int64) + ce == int(run.count["A"][0]) + int(ce[0]))
        assert np.all(np.asarray(run.count["B"], dtype=np.int64) - ce == int(run.count["B"][0]) - int(ce[0]))
    elif name == "H18_rxn_bidiff":
        e = ev(run)
        assert e > 0
        ce = np.asarray(run.channels.event_count[channel(run)], dtype=np.int64)
        assert np.all(np.asarray(run.count["A"], dtype=np.int64) + ce == int(run.count["A"][0]) + int(ce[0]))
        assert np.all(np.asarray(run.count["B"], dtype=np.int64) + ce == int(run.count["B"][0]) + int(ce[0]))
        assert np.all(np.asarray(run.count["C"], dtype=np.int64) - ce == int(run.count["C"][0]) - int(ce[0]))
    elif name == "H18_rxn_bisame":
        e = ev(run)
        assert e > 0
        ce = np.asarray(run.channels.event_count[channel(run)], dtype=np.int64)
        assert np.all(np.asarray(run.count["A"], dtype=np.int64) + 2 * ce == int(run.count["A"][0]) + 2 * int(ce[0]))
        assert np.all(np.asarray(run.count["C"], dtype=np.int64) - ce == int(run.count["C"][0]) - int(ce[0]))

    elif name.startswith("H19_eff_"):
        e = ev(run)
        p = productive(run)
        q = nonproductive(run)
        assert e == p + q == 500
        ce = np.asarray(run.channels.event_count[channel(run)], dtype=np.int64)
        cp = np.asarray(run.channels.productive_event_count[channel(run)], dtype=np.int64)
        assert np.all(np.asarray(run.count["A"], dtype=np.int64) + ce == int(run.count["A"][0]) + int(ce[0]))
        assert np.all(np.asarray(run.count["B"], dtype=np.int64) - cp == int(run.count["B"][0]) - int(cp[0]))
        if name.endswith("_0"):
            assert p == 0 and q == e
        elif name.endswith("_1"):
            assert p == e and q == 0
        else:
            frac = p / e
            assert abs(frac - 0.25) < 0.06, (p, e, frac)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", type=Path, default=ROOT.parent / "bin" / "slimmc")
    args = ap.parse_args()
    models = sorted((Path(__file__).parent / "models").glob("H*.model"))
    for model in models:
        run = run_model(model, args.engine.resolve())
        check(model.stem, run)
        print(f"[PASS] {model.stem}")
    print(f"[PASS] phase D integration ({len(models)} models)")

if __name__ == "__main__":
    main()
