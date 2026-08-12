# Homo chemistry phase A — detailed validation

This suite complements the fast Nim unit tests in `tests/chemistry/test_phase_a.nim`.
It runs complete `.model -> slimmc-storage -> pyslimmc` workflows without adding
or changing model-language syntax.

## Coverage

- H01 propagation: positive run, `kp=0` control, exact event/unit relation, and
  `dp_max` saturation.
- H02 depropagation:
  - D01 isolated depropagation after a chain-building phase,
  - D02 same-seed inactive control,
  - D03 ordered rate scaling for `kdp = 10, 25, 50 s^-1`,
  - D04 propagation-dominant, near-balanced, and depropagation-dominant regimes.
  Exact checks use `deprop events = released monomers = lost repeat units`.
- H03 combination: active model plus same-seed inactive control; verifies
  `2 live -> 1 dead` per event and exact repeat-unit conservation.
- H04 disproportionation: active model plus control; verifies
  `2 live -> 2 dead`, exactly one H and one U product per event.
- H05 transfer to monomer: active model plus control; verifies constant radical
  count, one new dead chain and one consumed monomer per event, and growth of
  total chain count by one per event.
- H06 common invariants at every chain snapshot: nonnegative counts, `DP >= 1`,
  and constant `free M + repeat units in chains` for all closed models.

Controls use the same seed and identical pre-switch chemistry. The checker also
verifies that paired active/control models have identical state at the switch.

## Run

```bash
cd homo
make test-phase-a
```

The command runs parser tests, all Nim phase-A/phase-B unit tests, builds the
family CLI, runs 17 validation models, and reads every result through the public
`pyslimmc` API.
