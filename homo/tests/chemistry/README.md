# Homo chemistry unit tests

`test_phase_a.nim` covers the first critical chemistry phase:

- H01 propagation, including `dp_max` eligibility;
- H02 depropagation and the DP=1 guard;
- H03 termination by combination;
- H04 termination by disproportionation;
- H05 transfer to monomer;
- H06 shared monomer-balance and pool invariants.

Run from `homo/`:

```bash
make test
```

These are fast deterministic unit tests of `computePropensities` and
`applyChannel`. Full model/storage/pyslimmc validation remains a separate
integration layer.
