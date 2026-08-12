from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import pyslimmc as sl

HERE = Path(__file__).resolve().parent
MODELS = HERE / "models"
CASES = ("EQ01_init", "EQ02_prop", "EQ03_deprop", "EQ04_term_c", "EQ05_term_d")


def replace_seed(text: str, seed: int) -> str:
    out, n = re.subn(r"(?m)^param seed\s+\d+\s*$", f"param seed {seed}", text)
    assert n == 1
    return out


def run_model(engine: Path, model: Path, seed: int, work: Path):
    work.mkdir(parents=True, exist_ok=True)
    local = work / model.name
    local.write_text(replace_seed(model.read_text(encoding="utf-8"), seed), encoding="utf-8")
    cp = subprocess.run([str(engine), str(local)], cwd=work, capture_output=True, text=True)
    assert cp.returncode == 0, (model, cp.stdout, cp.stderr)
    return sl.open(work / "results" / model.stem)


def weighted_total(chains, field: str = "dp") -> int:
    counts = np.asarray(chains.count, dtype=np.int64)
    values = np.asarray(getattr(chains, field), dtype=np.int64)
    return int(np.sum(counts * values))


def chain_count(chains) -> int:
    return int(np.sum(np.asarray(chains.count, dtype=np.int64)))


def events(run, needle: str) -> int:
    names = [name for name in run.channels.event_count.names if needle in name]
    assert names, (needle, run.channels.event_count.names)
    return int(sum(int(run.channels.event_count[name][-1]) for name in names))


def snapshot_at(run, t: float):
    idx = np.flatnonzero(np.isclose(np.asarray(run.t, dtype=float), t, atol=1e-12, rtol=0.0))
    assert idx.size, (t, run.t)
    return run.at_snapshot(int(idx[-1]))


def monomer_count(run, engine_kind: str, idx: int = -1) -> int:
    key = "M" if engine_kind == "homo" else "monomer_A"
    return int(run.count[key][idx])


def assert_one_monomer_copo(run) -> None:
    assert np.all(np.asarray(run.count["monomer_B"], dtype=np.int64) == 0)
    chains = run.last.chains.all
    if len(chains.dp):
        assert np.all(np.asarray(chains.counts["B"], dtype=np.int64) == 0)
        assert np.array_equal(
            np.asarray(chains.dp, dtype=np.int64),
            np.asarray(chains.counts["A"], dtype=np.int64),
        )


def metrics(case: str, run, engine_kind: str) -> dict[str, float]:
    if engine_kind == "copo":
        assert_one_monomer_copo(run)
    assert run.is_complete
    all_chains = run.last.chains.all
    live = run.last.chains.live
    dead = run.last.chains.dead
    out = {
        "chains": float(chain_count(all_chains)),
        "live": float(chain_count(live)),
        "dead": float(chain_count(dead)),
        "units": float(weighted_total(all_chains)),
        "monomer": float(monomer_count(run, engine_kind)),
    }

    if case == "EQ01_init":
        n_init = events(run, "init")
        assert chain_count(all_chains) == n_init
        assert weighted_total(all_chains) == n_init
        out["init"] = float(n_init)

    elif case == "EQ02_prop":
        n_init = events(run, "init")
        n_prop = events(run, "prop")
        assert chain_count(live) == n_init
        assert chain_count(dead) == 0
        assert weighted_total(live) == n_init + n_prop
        out.update(init=float(n_init), prop=float(n_prop))

    elif case == "EQ03_deprop":
        before = snapshot_at(run, 0.10)
        before_units = weighted_total(before.chains.live)
        before_m = monomer_count(run, engine_kind, int(np.flatnonzero(np.isclose(run.t, 0.10))[-1]))
        n_deprop = events(run, "deprop")
        assert before_units - weighted_total(live) == n_deprop
        assert monomer_count(run, engine_kind) - before_m == n_deprop
        assert chain_count(before.chains.live) == chain_count(live)
        out.update(deprop=float(n_deprop), units_before=float(before_units))

    elif case == "EQ04_term_c":
        n_init = events(run, "init")
        n_term = events(run, "term_c")
        assert chain_count(live) in (0, 1)
        assert chain_count(dead) == n_term
        assert chain_count(live) + 2 * n_term == n_init
        assert n_term > 0
        out["term"] = float(n_term)

    elif case == "EQ05_term_d":
        n_init = events(run, "init")
        n_term = events(run, "term_d")
        assert chain_count(live) in (0, 1)
        assert chain_count(dead) == 2 * n_term
        assert chain_count(live) + 2 * n_term == n_init
        assert n_term > 0
        out["term"] = float(n_term)

    return out


def rel_diff(a: float, b: float) -> float:
    scale = max(abs(a), abs(b), 1.0)
    return abs(a - b) / scale


def compare_ensemble(case: str, homo: list[dict[str, float]], copo: list[dict[str, float]]) -> None:
    hmean = {k: float(np.mean([row[k] for row in homo])) for k in homo[0]}
    cmean = {k: float(np.mean([row[k] for row in copo])) for k in copo[0]}

    # Exact structural quantities should agree for these deliberately saturated cases.
    for key in ("chains", "live", "dead"):
        assert hmean[key] == cmean[key], (case, key, hmean[key], cmean[key])

    # Stochastic trajectories use separate engine implementations/RNG call streams;
    # compare ensemble observables, not event-by-event identity.
    assert rel_diff(hmean["units"], cmean["units"]) < 0.06, (case, "units", hmean, cmean)

    if case == "EQ01_init":
        assert hmean["init"] == cmean["init"]
    elif case == "EQ02_prop":
        assert hmean["init"] == cmean["init"]
        assert rel_diff(hmean["prop"], cmean["prop"]) < 0.06
    elif case == "EQ03_deprop":
        assert rel_diff(hmean["units_before"], cmean["units_before"]) < 0.06
        assert rel_diff(hmean["deprop"], cmean["deprop"]) < 0.30
    elif case in ("EQ04_term_c", "EQ05_term_d"):
        assert hmean["term"] == cmean["term"]

    print(f"{case}: PASS")
    print("  homo means:", hmean)
    print("  copo means:", cmean)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--homo-engine", type=Path, required=True)
    ap.add_argument("--copo-engine", type=Path, required=True)
    ap.add_argument("--seeds", type=int, default=8)
    ns = ap.parse_args()
    assert ns.seeds >= 4

    with tempfile.TemporaryDirectory(prefix="slimmc_homo_copo_equiv_") as td:
        work = Path(td)
        for ci, case in enumerate(CASES):
            hm, cm = [], []
            for j in range(ns.seeds):
                seed = 9400 + 100 * ci + j
                hr = run_model(ns.homo_engine.resolve(), MODELS / "homo" / f"{case}.model", seed, work / case / f"h{j}")
                cr = run_model(ns.copo_engine.resolve(), MODELS / "copo" / f"{case}.model", seed, work / case / f"c{j}")
                hm.append(metrics(case, hr, "homo"))
                cm.append(metrics(case, cr, "copo"))
            compare_ensemble(case, hm, cm)

    print("Homo vs effective one-monomer copo equivalence: PASS")


if __name__ == "__main__":
    main()
