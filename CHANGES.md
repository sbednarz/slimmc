# Changelog

This file records public, user-visible changes. Detailed implementation history
is available from the Git history.

## 5.0.0 — final unified release

- unified the homo and copo engines behind the public `slimmc` dispatcher;
- released pyslimmc 4.0.0 and pyslimmc-opt 1.0.0 alongside the engine family;
- finalized the common model-language contracts, canonical `with_end_groups`
  spelling, and explicit homo conversion-condition syntax;
- added native release artifacts for Linux glibc 2.28+, static Linux musl, and
  Windows x86-64;
- finalized build provenance and reproducibility verification, including Git and
  binary/storage hashes when available;
- completed the independent integration suite with 92 end-to-end tests;
- moved worked examples and literature studies out of the core repository so
  release validation depends only on core-owned fixtures;
- reorganized the documentation into task-oriented, exact-reference, and
  developer/AI layers, including quick start, concepts, limitations, Storage
  specification, testing contracts, and generated API signatures.

## 5.0.0-rc.3

- simplified component Makefiles and centralized integration/release validation
  in the root Makefile;
- added release-configuration and literal-path checks;
- completed independent integration coverage for CLI, both engines, Storage,
  pyslimmc, and pyslimmc-opt.

## 5.0.0-rc.2

- made `param mass_model with_end_groups` the sole end-group-aware spelling;
- made homo conversion conditions use the explicit
  `when X MONOMER >|< VALUE ...` form.

## 5.0.0-rc.1

- established the unified family CLI contract `slimmc [options] model.model`;
- introduced one public engine-family version with independently versioned
  pyslimmc and pyslimmc-opt packages;
- added concise family help/version output and build metadata reporting.

## Selected pre-5.0 history

### slimmc 3.3.1

- synchronized component/package versions and current reference documentation;
- made omitted copolymer propagation transitions explicit zero-rate warnings;
- improved chain-snapshot, chain-count, feed-event, and CLI conveniences.

### slimmc 3.3.0 — semibatch

- added portion-wise semibatch feeds to both engines;
- added dynamic physical/KMC volumes, feed events, and physical species
  balances;
- added read-only semibatch analysis in pyslimmc;
- renamed the stochastic model volume to `kmc_volume`;
- advanced Slimmc Storage to 1.2.0.

### slimmc 3.1.0

- added canonical `run_id`, multidimensional `var` declarations, `Runs.match()`,
  lazy run access, `Runs.pack()`, and multidimensional `Runs.sweep()`.

### pyslimmc 3.7.1

- stabilized cached full-sequence statistics and filtered views;
- kept cache publication atomic and clarified multi-character n-gram labels.

### pyslimmc 3.7.0

- added full-sequence motif/ngram, relative-position, and DP × microstructure
  analyses with matching plotting shortcuts.

### pyslimmc 3.6.0

- unified common plotting conveniences under `run.plot` while preserving the
  data-first analysis objects and plain Matplotlib workflow.

### pyslimmc 3.5.0

- added full-sequence transition/block statistics, filters, maps, shared cache,
  and optional progress reporting.

### pyslimmc 3.4.0

- added chain-composition filters, composition-by-DP statistics, 2D
  composition maps, and homo/copolymer/terpolymer component classes.

### pyslimmc 3.0.0

- cleaned up snapshot/distribution APIs, made run-level analyses default to the
  final snapshot, renamed `chain_spectrum()` to `chain_mass_spectrum()`, and
  made `.x` physical with `.log10_x` as the explicit transformed coordinate.

### pyslimmc 2.1.1

- standardized `run.desc`;
- removed artificial zero-conversion insertion when no `t=0` snapshot exists;
- allowed cumulative/interval copolymer composition reconstruction from monomer
  depletion when chain snapshots are sparse.
