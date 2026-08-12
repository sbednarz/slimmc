# slimmc architecture

The public executable is one in-process dispatcher. It counts non-comment
`monomer` declarations and routes models with 0–1 monomers to homo and models
with 2–3 monomers to copo. Each engine parses and validates the complete model.

```text
cli/                         unified dispatcher and summary binary
common/                      shared model, Storage, hashing, and numeric code
homo/src/                    homopolymer parser, SSA engine, IO, and writer
copo/src/                    copolymer parser, SSA engine, IO, and writer
pyslimmc/                    read-only Python analysis package and tests
pyslimmc_opt/                optimization layer and tests
tests/integration/           independent CLI/engine/API/optimization fixtures
tests/validation/            shared technical black-box cases
homo/tests/                  homo unit, chemistry, regression, validation
copo/tests/                  copo unit, chemistry, validation
validation/regressions/      seeded cross-component regression workflows
docs/                        canonical user and developer documentation
scripts/                     repository, release, docs, and CI checks
.github/workflows/           continuous integration and release workflows
```

## Data flow

![Slimmc architecture](../assets/slimmc-arch.svg)

```text
.model -> dispatcher -> engine/parser -> SSA state -> Storage writer
       -> slimmc Storage (.npy + JSON/JSONL) -> pyslimmc -> user exports
```

The Storage directory is the only results backend. TSV, JSON summaries,
figures, and PDF reports produced by pyslimmc are derived user exports, not
inputs to the reader and not part of a run.

## Ownership boundaries

- engines own model parsing, reaction execution, state and chain bookkeeping,
  snapshot creation, and Storage writing;
- `common/` owns contracts shared by both engines;
- pyslimmc is read-only and never repairs or mutates a Storage run;
- pyslimmc-opt builds on public model/run interfaces rather than owning an
  alternative result format;
- validation diagnostics do not replace scientific source tables;
- user examples live outside the core repository; technical validation cases
  in the core tree pin invariants.

## Eligibility counters in the SSA hot path

The engines maintain incremental counters so common chain-selection operations
remain O(1) when an entire pool is eligible. The counters are an optimization,
not a second source of chemical truth: exact pool contents remain canonical.

### homo

`State` keeps `poolPropagatableCounts`, `poolDepropableCounts`, and
`poolEligibilityTrackedLen`. Engine-owned mutations use `trackAdd`,
`trackRemove`, and `trackDpChange`, which update the counters together with the
pool. `poolEligibilityTrackedLen[poolId] == pools[poolId].len` means the cached
counts correspond to the current pool. If they are not current, propensity and
selection helpers rebuild or scan exactly rather than trusting stale values.

In debug mode `debugCheckState` recomputes propagatable and depropagatable
counts by a full pool scan and reports a `counter mismatch` if the stored
counter disagrees. A change that mutates a chain or pool without going through
an engine-owned tracking path is therefore a correctness risk even when normal
release output still looks plausible.

### copo

Copo additionally tracks pool invariant compatibility and terminal-dependent
depropagation eligibility through `poolEligibleCounts`,
`poolPropagatableCounts`, `poolEligibilityTrackedLen`, and
`channelDepropEligibleCounts`. `trackAdd` and `trackRemove` update all related
counters. A negative tracked length is a validity sentinel for a pool whose
state was manually assembled or mutated: that pool stays scan-backed instead
of allowing a partial mutation to make stale counters appear current.

Fast selection is used only when the counters prove that the whole pool is
eligible for the requested operation. Otherwise the engine falls back to the
older exact eligible-index scan. This fallback is intentional hardening and
must not be removed merely to improve benchmark numbers.

### Mutation rule

When adding a new engine operation that changes a live pool or chain state:

1. use or extend the engine-owned tracking helpers rather than editing the pool
   directly;
2. update every eligibility class affected by the mutation;
3. preserve the scan fallback for untracked/manually mutated states;
4. add a focused test that compares the optimized path with exact eligibility;
5. run the affected engine validation suites and `make test`.

## Why these boundaries exist

A few architectural choices are deliberate and should not drift silently:

- **One dispatcher, two engines.** Homo and copo share the public CLI and
  Storage contract, but retain separate parsers/state machines where their
  kinetic representations differ materially.
- **Storage is the durable boundary.** Engines write it; pyslimmc reads it.
  Analysis code never silently repairs a run, so a corrupted or incomplete run
  remains observable as such.
- **No pickle in canonical results.** The format uses typed NumPy columns and
  JSON/JSONL metadata so it can be inspected independently of pyslimmc and does
  not depend on Python object serialization.
- **Scientific validation is separate from examples.** User-facing examples
  are allowed to evolve pedagogically; core validation fixtures pin numerical,
  chemical, Storage, and API invariants.
- **Optimization counters are disposable caches.** Exact state wins whenever a
  cache cannot prove that it is current.

The exact Storage schema is specified in
[`../reference/STORAGE.md`](../reference/STORAGE.md); test obligations are in
[`TESTING.md`](TESTING.md).

## See also

- [`../reference/STORAGE.md`](../reference/STORAGE.md) — Storage format
- [`TESTING.md`](TESTING.md) — Testing and validation
- [`DEVELOPMENT.md`](DEVELOPMENT.md) — Development workflow
