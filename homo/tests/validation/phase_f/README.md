# Homo phase F — analytic laws

Phases A–E verify that the engine does what the model declares: exact event
counts, conservation identities, seeded controls, reproducibility. Those
checks pass under any self-consistent kinetics, including incorrect ones.
Phase F asks the separate question of whether the declared kinetics reproduce
results derived outside slimmc.

## Coverage

- **F01** disproportionation-only steady state → most-probable geometric CLD,
  `D = 2 - 1/DPn`, plus a binned comparison of `P(n)` with `(1-p)p^(n-1)`.
- **F02** combination-only steady state → convolution of two geometric chains,
  `D = 1.5 - 1/DPn` (approaching 1.5 at long chain length).
  Distinguishes a correct combination event from one that merely retires two
  radicals.
- **F03** transfer to monomer → Mayo equation: `1/DPn` linear in `C_M` with
  unit slope. Slope-based, therefore insensitive to the `kt` convention.
- **F04** depropagation equilibrium: free monomer settles at `[M]eq = kdp/kp`,
  approached from both sides in two runs.
- **F05** instantaneous initiation without termination → shifted-Poisson CLD
  (`DP = 1 + Poisson(lambda)`), with `Var(DP) = DPn - 1` and the corresponding
  exact dispersity relation.

## Why these models do not use the phase A pattern

Phase A builds chains, then switches a channel on. That pattern is required
for exact event-count identities and makes them unambiguous. It is *wrong*
here: most-probable and Poisson distributions are steady-state results, and
switching termination on after propagation has stopped destroys the statistic
under test. F01–F03 and F05 therefore run all channels concurrently from
`t = 0`. F04 keeps the switch, because an equilibrium plateau is a
steady-state result of the post-switch chemistry alone.

## Conventions

Tolerances are absolute bands, several standard errors wide at the configured
chain counts, in the phase B style: deterministic seeded regression tests, not
flaky hypothesis tests. Each check asserts a minimum chain count first, so a
model that silently produces too little statistics fails loudly instead of
passing on noise.

## Run

    python check_phase_f.py --engine ../../../release/slimmc-homo
    python check_phase_f.py --engine ... --only F01 F04
