# Homo validation phase B: SSA reliability

This phase complements chemistry phase A with tests of the stochastic engine.

- **H07** (Nim unit tests): exact propensity scaling with rate constants,
  populations, volume, `dp_max`, and eligibility for depropagation.
- **H08** (full model -> storage -> pyslimmc): competing channels are selected
  in the theoretical propensity ratio.
- **H09** (full model -> channel trace -> pyslimmc): the normalized waiting
  time `z = a0 * dt` follows the unit exponential law.
- **H10** (full model -> storage -> pyslimmc): two runs with the same model and
  seed have byte-identical numerical event traces and state series.

Run from `homo/`:

```bash
make test-phase-b
```

The statistical checks use fixed seeds, thousands of events, and conservative
six-standard-error tolerances. They are deterministic regression tests, not
flaky hypothesis tests.
