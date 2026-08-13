# Copo phase F — Mayo–Lewis and Alfrey–Goldfinger

The existing copo suites verify that every propagation channel fires, that
pool metadata is consistent, and that engine dyad counters agree with counters
reconstructed from literal sequences. All are internal consistency checks:
**transposing the `kp` matrix leaves every one of them passing.** This group
closes that gap.

## Coverage

- **CML01** terminal model, `r1 = 0.5`, `r2 = 2.0`: instantaneous `F_A`
  against Mayo–Lewis; `AA` and `BB` dyad fractions against Alfrey–Goldfinger
  conditional probabilities; and `pyslimmc.compare_mayo_lewis()` used as a
  check on the engine rather than as user-facing analysis.
- **CML02** `r1 = r2 = 1`: Bernoulli limit. `F_A = f_A` and dyads factorize.
  Terminal-pool leakage shows up as a dyad imbalance even when composition
  stays correct.
- **CML03** `r1 = r2 = 0`: strictly alternating. `F_A = 0.5` regardless of
  feed, and homo-dyad counts are exactly zero — an integer assertion, not a
  statistical one.

The asymmetric ratios in CML01 are deliberate: a transposed matrix swaps `r1`
and `r2` and moves `F_A` to the other side of the azeotrope.

All models stop below 5 % conversion, where the instantaneous forms apply and
composition drift is negligible; the checker asserts this rather than assuming
it.

## Run

    python check_phase_f.py --engine ../../../../bin/slimmc
