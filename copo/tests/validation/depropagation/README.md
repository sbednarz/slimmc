# Detailed copolymer depropagation validation

This suite complements `tests/chemistry/test_depropagation*.nim` with black-box
models executed through the family `slimmc` dispatcher and read only through the
public `pyslimmc` API.

It checks:

- all terminal transitions AA, BA, AB and BB;
- independent A/B repeat-unit balance;
- exact equality between depropagation event counts and lost polymer units;
- a matched `kdeprop=0` control;
- rate scaling for `kd = 10, 20, 40 s-1`;
- clear propagation-dominant and depropagation-dominant regimes;
- consistency of DP, composition and terminal pools after depropagation.

Run with:

```sh
make -C copo test-depropagation
```
