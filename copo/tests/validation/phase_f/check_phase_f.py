"""Copo phase F -- Mayo-Lewis and Alfrey-Goldfinger analytic validation.

The existing copo suites verify that every propagation channel fires, that
pool metadata is consistent, and that engine dyad counters agree with the
counters reconstructed from literal sequences. All of those are internal
consistency checks: transposing the kp matrix would leave every one of them
passing.

This group compares the simulated composition and dyad statistics with the
closed-form terminal-model results. `pyslimmc.copolymerization.mayo_lewis()`
is exercised here as a check on the engine rather than as a user-facing
analysis, which is the one use it currently does not have.

All models run to low conversion, where composition drift is negligible and
the instantaneous forms apply directly.
"""

from __future__ import annotations
import argparse, subprocess, sys, tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
REPO = ROOT.parent
sys.path.insert(0, str(REPO))
import pyslimmc as sl

MODELS = Path(__file__).resolve().parent / 'models'


def run_model(engine: Path, model: Path, work: Path):
    work.mkdir(parents=True, exist_ok=True)
    local = work / model.name
    local.write_text(model.read_text())
    cp = subprocess.run([str(engine), str(local)], capture_output=True, text=True)
    assert cp.returncode == 0, (model.name, cp.stdout, cp.stderr)
    return sl.open(work / 'results' / model.stem)


def band(name, value, expected, half_width):
    ok = abs(value - expected) <= half_width
    print(f'  {name:<32} {value:8.4f}   expected {expected:.4f} '
          f'+/- {half_width:.4f}   {"ok" if ok else "FAIL"}')
    assert ok, (name, value, expected, half_width)


def free_fractions(run, index=-1):
    a = float(np.asarray(run.conc['A'], dtype=float)[index])
    b = float(np.asarray(run.conc['B'], dtype=float)[index])
    return a / (a + b), a, b


def polymer_fraction_A(run):
    c = run.final.chains.all
    n = np.asarray(c.count, dtype=np.int64)
    ca = int((n * np.asarray(c.counts['A'], dtype=np.int64)).sum())
    cb = int((n * np.asarray(c.counts['B'], dtype=np.int64)).sum())
    return ca / (ca + cb), ca + cb


def conversion(run):
    a = np.asarray(run.conc['A'], dtype=float)
    b = np.asarray(run.conc['B'], dtype=float)
    return 1.0 - (a[-1] + b[-1]) / (a[0] + b[0])


def mayo_lewis_F_A(f_a, r1, r2):
    """Instantaneous mole fraction of A in the copolymer, terminal model."""
    f_b = 1.0 - f_a
    num = r1 * f_a ** 2 + f_a * f_b
    den = r1 * f_a ** 2 + 2.0 * f_a * f_b + r2 * f_b ** 2
    return num / den


def dyad_fractions(run):
    """Normalized AA, AB+BA, BB dyad fractions from the engine counters."""
    tab = run.microstructure.dyads()
    d = {str(k).replace('|', ''): float(v) for k, v in zip(tab['motif'], tab['count'])}
    tot = sum(d.values())
    return {k: v / tot for k, v in d.items()}, tot


# --------------------------------------------------------------------------

def check_CML01(engine, work):
    """Terminal model with r1 = 0.5, r2 = 2.0 against Mayo-Lewis.

    A transposed kp matrix swaps r1 and r2 and moves F_A to the other side of
    the azeotrope; the asymmetric ratios here make that unmissable.
    """
    r1, r2 = 0.5, 2.0
    run = run_model(engine, MODELS / 'CML01_terminal_r05_r20.model', work)
    x = conversion(run)
    f_a, _, _ = free_fractions(run, 0)
    F_obs, units = polymer_fraction_A(run)
    F_pred = mayo_lewis_F_A(f_a, r1, r2)
    print(f'CML01 terminal r1={r1} r2={r2}: conversion {x*100:.2f} %, '
          f'{units} repeat units, f_A = {f_a:.4f}')
    assert x < 0.05, x                       # instantaneous form must apply
    assert units >= 200000, units
    band('F_A (Mayo-Lewis)', F_obs, F_pred, 0.02)

    # dyads: Alfrey-Goldfinger conditional probabilities
    dy, tot = dyad_fractions(run)
    p11 = r1 * f_a / (r1 * f_a + (1.0 - f_a))
    p22 = r2 * (1.0 - f_a) / (r2 * (1.0 - f_a) + f_a)
    # stationary chain statistics: P(AA) = F_A * p11 etc.
    band('AA dyad fraction', dy.get('AA', 0.0), F_pred * p11, 0.03)
    band('BB dyad fraction', dy.get('BB', 0.0), (1 - F_pred) * p22, 0.03)

    # cross-check pyslimmc's own Mayo-Lewis series against the engine
    ml = run.copolymerization.compare_mayo_lewis()
    diff = np.abs(np.asarray(ml.composition_difference[ml.monomers[0]], dtype=float))
    diff = diff[np.asarray(ml.is_defined, dtype=bool)]
    print(f'  pyslimmc compare_mayo_lewis max |F_obs - F_pred| = {diff.max():.4f}')
    assert diff.max() < 0.03, diff.max()


def check_CML02(engine, work):
    """r1 = r2 = 1: the terminal model degenerates to Bernoulli statistics.

    Composition must equal the feed and dyads must factorize, independently
    of the pool machinery. Any leakage between terminal pools shows up here
    as a dyad imbalance even though the composition stays correct.
    """
    run = run_model(engine, MODELS / 'CML02_ideal_r1_r1.model', work)
    x = conversion(run)
    f_a, _, _ = free_fractions(run, 0)
    F_obs, units = polymer_fraction_A(run)
    print(f'CML02 ideal r1=r2=1: conversion {x*100:.2f} %, f_A = {f_a:.4f}')
    assert x < 0.05, x
    band('F_A (= f_A)', F_obs, f_a, 0.015)
    dy, _ = dyad_fractions(run)
    band('AA dyad (= f_A^2)', dy.get('AA', 0.0), f_a ** 2, 0.02)
    band('BB dyad (= f_B^2)', dy.get('BB', 0.0), (1 - f_a) ** 2, 0.02)


def check_CML03(engine, work):
    """r1 = r2 = 0: strictly alternating.

    F_A = 0.5 whatever the feed, and homo-dyads must be exactly absent. This
    is an exact integer assertion, not a statistical one.
    """
    run = run_model(engine, MODELS / 'CML03_alternating_r0_r0.model', work)
    f_a, _, _ = free_fractions(run, 0)
    F_obs, units = polymer_fraction_A(run)
    print(f'CML03 alternating: f_A = {f_a:.4f}, {units} repeat units')
    band('F_A (= 0.5 regardless of feed)', F_obs, 0.5, 0.01)
    tab = run.microstructure.dyads()
    d = {str(k).replace('|', ''): float(v) for k, v in zip(tab['motif'], tab['count'])}
    print(f'  homo-dyad counts: AA = {d.get("AA", 0):.0f}, BB = {d.get("BB", 0):.0f}')
    assert d.get('AA', 0.0) == 0.0, d
    assert d.get('BB', 0.0) == 0.0, d


CHECKS = {'CML01': check_CML01, 'CML02': check_CML02, 'CML03': check_CML03}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--engine', required=True, type=Path)
    ap.add_argument('--only', nargs='*', choices=sorted(CHECKS))
    args = ap.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for key in (args.only or sorted(CHECKS)):
            CHECKS[key](args.engine, work / key)
            print()
    print('copo phase F: all analytic checks passed')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
