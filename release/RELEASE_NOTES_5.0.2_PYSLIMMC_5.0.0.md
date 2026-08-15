# Slimmc 5.0.2 / pyslimmc 5.0.0

This family update keeps the Slimmc 5.0.x simulator and Storage contracts stable
while releasing a breaking redesign of the pyslimmc distribution-analysis API.

## Slimmc 5.0.2

- no change to the model-language or kMC chemistry semantics;
- integration tests were updated to the new pyslimmc distribution contract;
- build/version metadata now reports pyslimmc 5.0.0.

## pyslimmc 5.0.0

- exact `dp_counts()` and `mass_counts()` replace the old chain-count / mass-spectrum split;
- CLD/MWD use `form="number"|"mass"|"z"|"log"`;
- histogram, KDE and generic Gaussian distribution methods are removed;
- logarithmic CLD/MWD are exact discrete distributions, not axis-scaling aliases;
- `sec()` provides explicit Gaussian SEC broadening in `log10(M)` using required `sigma_log10M`;
- MWD is derived from actual chain masses, independently of CLD;
- canonical mass-model resolution is shared by counts, MWD, moments and SEC;
- `PopulationMoments` distinguishes DP and molar-mass dispersity;
- multi-series `per_series` and `combined` normalization is supported with overlap validation.

### Migration highlights

```text
chain_counts()              -> dp_counts()
chain_mass_spectrum()       -> mass_counts()
basis/method/coordinate/output -> form
histogram/KDE/Gaussian MWD  -> removed
instrumental broadening     -> sec(sigma_log10M=...)
```

Exact chain amounts belong to `dp_counts()` / `mass_counts()`; CLD/MWD are normalized
mathematical distributions. SEC is a separate apparent instrumental response.
