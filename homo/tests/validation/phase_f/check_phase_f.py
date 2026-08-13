"""Homo phase F -- comparison against closed-form kinetic results.

Phases A-E verify that the engine does what the model declares: exact event
counts, conservation identities, controls, reproducibility. They pass under any
self-consistent kinetics, including wrong ones. This phase asks the separate
question of whether the declared kinetics reproduce laws derived outside
slimmc.

The models here deliberately do NOT use the build-then-switch pattern of
phase A. Most-probable and Poisson distributions are steady-state results;
switching termination on after propagation stops destroys the very statistics
being tested.

Tolerances are stated as absolute bands chosen to be several standard errors
wide at the configured chain counts, following the phase B convention: these
are deterministic seeded regression tests, not hypothesis tests.
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


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def dead_population(run):
    """(dp, multiplicity) of dead chains at the final chain snapshot."""
    c = run.final.chains.dead
    dp = np.asarray(c.dp, dtype=np.int64)
    n = np.asarray(c.count, dtype=np.int64)
    assert dp.size, 'no dead chains stored'
    return dp, n


def moments(dp, n):
    N = int(n.sum())
    s1 = float((n * dp).sum())
    s2 = float((n * dp.astype(float) ** 2).sum())
    dpn = s1 / N
    dpw = s2 / s1
    return N, dpn, dpw, dpw / dpn


def band(name, value, expected, half_width, extra=''):
    lo, hi = expected - half_width, expected + half_width
    ok = lo <= value <= hi
    print(f'  {name:<34} {value:10.4f}   expected {expected:.4f} '
          f'+/- {half_width:.4f}   {"ok" if ok else "FAIL"}  {extra}')
    assert ok, (name, value, expected, half_width)


# --------------------------------------------------------------------------
# F01 / F02 -- dispersity limits
# --------------------------------------------------------------------------

def check_F01(engine, work):
    """Termination exclusively by disproportionation, no transfer.

    The instantaneous CLD is the most-probable (geometric) distribution.
    For a geometric DP distribution on n >= 1, D = 2 - 1/DPn; the familiar
    D -> 2 result is its long-chain limit.
    """
    run = run_model(engine, MODELS / 'F01_disproportionation_mwd.model', work)
    dp, n = dead_population(run)
    N, dpn, dpw, D = moments(dp, n)
    print(f'F01 disproportionation: {N} dead chains, DPn = {dpn:.1f}')
    assert N >= 5000, N
    assert dpn >= 50, dpn                       # long-chain limit must hold
    expected = 2.0 - 1.0 / dpn
    band('D (geometric)', D, expected, 0.10)
    check_geometric_shape(dp, n, dpn)


def check_geometric_shape(dp, n, dpn):
    """Binned comparison against P(n) = (1-p) p^(n-1), p = 1 - 1/DPn.

    Bins are made wide enough that Poisson counting noise stays below ~5 %,
    otherwise the test measures noise rather than shape.
    """
    p = 1.0 - 1.0 / dpn
    N = int(n.sum())
    width = max(1, int(round(dpn / 4)))
    edges = np.arange(1, int(dp.max()) + width + 1, width)
    if edges[-1] <= int(dp.max()):
        edges = np.append(edges, int(dp.max()) + 1)
    obs, _ = np.histogram(dp, bins=edges, weights=n)
    # Exact mass of P(n)=(1-p)p^(n-1) over each integer histogram bin
    # [lo, hi): sum_{n=lo}^{hi-1} P(n) = p^(lo-1)-p^(hi-1).
    pred = N * (p ** (edges[:-1] - 1) - p ** (edges[1:] - 1))
    keep = obs >= 400                            # ~5 % counting noise
    dev = np.abs(obs[keep] - pred[keep]) / pred[keep]
    print(f'  geometric shape: max deviation {dev.max()*100:.1f} % '
          f'over {keep.sum()} bins')
    assert dev.max() < 0.15, dev.max()


def check_F02(engine, work):
    """Termination exclusively by combination.

    Each dead chain is the sum of two independent geometric chains. For the
    finite-DP distribution used here, D = 1.5 - 1/DPn; the familiar D -> 1.5
    result is the long-chain limit. This distinguishes a correct combination
    event from one that merely retires two radicals.
    """
    run = run_model(engine, MODELS / 'F02_combination_mwd.model', work)
    dp, n = dead_population(run)
    N, dpn, dpw, D = moments(dp, n)
    print(f'F02 combination: {N} dead chains, DPn = {dpn:.1f}')
    assert N >= 2500, N
    expected = 1.5 - 1.0 / dpn
    band('D (combination)', D, expected, 0.10)


# --------------------------------------------------------------------------
# F03 -- Mayo equation
# --------------------------------------------------------------------------

def check_F03(engine, work):
    """1/DPn is linear in C_M with unit slope.

    Absolute DPn depends on the termination convention; the slope does not.
    This test is therefore insensitive to the kt factor-of-two ambiguity and
    isolates the transfer channel.
    """
    cms = [0.000, 0.010, 0.020]
    names = ['F03_mayo_CM0.model', 'F03_mayo_CM10.model', 'F03_mayo_CM20.model']
    inv = []
    for cm, name in zip(cms, names):
        run = run_model(engine, MODELS / name, work)
        dp, n = dead_population(run)
        N, dpn, _, _ = moments(dp, n)
        inv.append(1.0 / dpn)
        print(f'F03 C_M = {cm:.3f}: {N} chains, DPn = {dpn:.1f}, 1/DPn = {1/dpn:.5f}')
    slope, intercept = np.polyfit(cms, inv, 1)
    resid = np.asarray(inv) - (slope * np.asarray(cms) + intercept)
    band('Mayo slope', slope, 1.0, 0.15)
    print(f'  max residual {np.abs(resid).max():.6f}')
    assert np.abs(resid).max() < 0.1 * max(inv)


# --------------------------------------------------------------------------
# F04 -- depropagation equilibrium
# --------------------------------------------------------------------------

def check_F04(engine, work):
    """Free monomer settles at [M]eq = kdp/kp regardless of starting side.

    Two runs reach the plateau from opposite directions. Agreement of both
    with kdp/kp is a single-number test of the depropagation rate law that
    the existing event-count identities cannot provide.
    """
    kp, kdp = 1.0e3, 100.0
    m_eq = kdp / kp
    finals = []
    for name in ('F04_deprop_equilibrium_high.model',
                 'F04_deprop_equilibrium_low.model'):
        run = run_model(engine, MODELS / name, work)
        m = np.asarray(run.conc['M'], dtype=float)
        tail = m[-3:]                            # last three saved points
        drift = abs(tail[-1] - tail[0]) / tail[-1]
        print(f'F04 {name.split("_")[-1][:-6]:>5}: [M]final = {m[-1]:.4f} '
              f'mol/L, tail drift {drift*100:.2f} %')
        assert drift < 0.05, (name, drift)       # plateau actually reached
        finals.append(float(m[-1]))
    for v in finals:
        band('[M]eq', v, m_eq, 0.10 * m_eq)
    rel = abs(finals[0] - finals[1]) / m_eq
    print(f'  two-sided agreement: {rel*100:.2f} %')
    assert rel < 0.10, rel


# --------------------------------------------------------------------------
# F05 -- Poisson limit
# --------------------------------------------------------------------------

def check_F05(engine, work):
    """Instantaneous initiation, no termination, no transfer.

    All chains are initiated at DP=1 and then accumulate propagation events
    with Poisson statistics. Thus DP = 1 + Poisson(lambda), so
    Var(DP) = DPn - 1 and D = 1 + (DPn - 1)/DPn^2. This exercises the engine
    in a regime with no dead-chain bookkeeping at all.
    """
    run = run_model(engine, MODELS / 'F05_poisson_living.model', work)
    c = run.final.chains.live
    dp = np.asarray(c.dp, dtype=np.int64)
    n = np.asarray(c.count, dtype=np.int64)
    N, dpn, dpw, D = moments(dp, n)
    expected = 1.0 + (dpn - 1.0) / (dpn ** 2)
    print(f'F05 Poisson: {N} live chains, DPn = {dpn:.2f}')
    assert N >= 5000, N
    band('D (Poisson)', D, expected, 0.02)
    var = float((n * (dp - dpn) ** 2).sum()) / N
    band('Var (shifted Poisson)', var, dpn - 1.0, 0.10 * max(1.0, dpn - 1.0))


# --------------------------------------------------------------------------

CHECKS = {'F01': check_F01, 'F02': check_F02, 'F03': check_F03,
          'F04': check_F04, 'F05': check_F05}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--engine', required=True, type=Path)
    ap.add_argument('--only', nargs='*', choices=sorted(CHECKS))
    args = ap.parse_args()
    selected = args.only or sorted(CHECKS)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for key in selected:
            CHECKS[key](args.engine, work / key)
            print()
    print('phase F: all analytic checks passed')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
