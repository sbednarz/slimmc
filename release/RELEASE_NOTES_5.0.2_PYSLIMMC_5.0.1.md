# Slimmc 5.0.2 / pyslimmc 5.0.1

This corrective pyslimmc release repairs the CLD/MWD mathematical contract introduced
in pyslimmc 5.0.0. The Slimmc 5.0.2 simulator, model language, Storage format, and kMC
chemistry semantics are unchanged.

## Important

pyslimmc 5.0.0 should not be used for logarithmic CLD/MWD analysis. Its `form="log"`
representation changed the abscissa to `log10(M)` while retaining discrete mass fractions
on the ordinate, which is not `dW/dlog10(M)`.

## pyslimmc 5.0.1

- restores a strict separation between exact discrete source distributions and derived
  density representations;
- keeps `dp_counts()` and `mass_counts()` as the exact raw source projections;
- defines `cld(weighting="number"|"mass"|"z")` as a normalized exact discrete
  chain-length distribution on integer DP support;
- adds `mass_distribution(weighting="number"|"mass"|"z")` as the normalized exact
  discrete distribution on exact molar-mass support;
- redefines `mwd()` as the normalized logarithmic mass-weighted density
  `dW/dlog10(M)`;
- reconstructs `mwd()` using an mcPolymer-style `N M^2` transform followed by
  piecewise-linear interpolation in `log10(M)` and area normalization;
- zero-fills missing integer-DP states for homopolymer MWD reconstruction when the
  declared mass model gives a unique mass for each DP;
- uses occupied exact-mass support without zero-filling for general/copolymer MWD
  reconstruction;
- keeps `sec()` as a separate instrument-response transform applied directly to the
  exact mass measure, independent of `mwd()` reconstruction;
- removes the `form="number"|"mass"|"z"|"log"` distribution contract;
- updates plotting, API references, quick-start material, cookbook examples, glossary,
  concepts, and generated signatures to the corrected semantics;
- replaces regression tests that encoded `log.y == mass_fraction.y` with normalization,
  exact-distribution, reconstruction, and SEC-independence tests.

## Public API

```python
run.dp_counts()
run.cld(weighting="number")

run.mass_counts()
run.mass_distribution(weighting="mass")

run.mwd()
run.sec(sigma_log10M=0.05)
```

`mass_distribution()` is the exact discrete representation and is the preferred output
when discrete oligomer species matter. `mwd()` is a derived density representation for
polymer-distribution analysis. `sec()` remains a separate apparent/instrument-response
representation.

## Compatibility

- Slimmc engine: unchanged at 5.0.2.
- Storage: unchanged.
- Model language: unchanged.
- pyslimmc 5.0.0 distribution calls using `form=` require migration.
- pyslimmc-opt: unchanged.

## Validation

The corrected pyslimmc unit/API suite passes 127 tests. The MWD contract is additionally
checked against exact discrete distributions, sparse/heavy-tail edge cases, homopolymer
zero-filling behavior, copolymer exact-mass support, area normalization, and independent
SEC transformation.
