# slimmc testing and validation

Stochastic output can look plausible while violating chemical, scheduling,
storage, or API invariants. The suite therefore combines direct Nim tests,
complete `.model` runs, Storage checks, and public pyslimmc assertions.

## Canonical commands

```bash
make test-fast                     # pyslimmc only
make test                          # normal local suite
make test-phase-a
make test-phase-b
make test-phase-c
make test-phase-d
make test-phase-e
make test-validation
make test-depropagation
make test-terminal-microstructure
make test-homo-copo-equivalence
make test-devel                    # phases and targeted technical suites
make test-integration              # 92 real end-to-end API tests
make test-full                     # standard plus development regression
make test-release                  # release metadata and packaging gates
```

The scenario-to-contract map is maintained in
[`INTEGRATION_COVERAGE.md`](INTEGRATION_COVERAGE.md).

Component Makefiles are intentionally limited to standalone engine work:

```bash
make -C homo build
make -C homo test
make -C copo build
make -C copo test
```

All black-box, phase, cross-engine, integration, examples-documentation, and
release checks are owned by the root Makefile. `make check-makefiles` rejects
literal script, Nim source, or model paths that do not exist in the distributed
tree.

## Test layers

1. **Nim unit tests** pin parser behavior, propensities, reaction application,
   chain bookkeeping, action order, safe arithmetic, and deterministic helpers.
2. **Black-box models** execute small `.model` files through the normal family
   dispatcher and assert observable invariants.
3. **Storage validation** checks schema, dtypes, row alignment, identity,
   hashes, finalization, balances, and chain/composition consistency.
4. **Public pyslimmc tests** read results only through supported APIs unless
   the Storage contract itself is under test.
5. **Statistical tests** check SSA selection and waiting-time laws with explicit
   tolerances; exact replay is required only for a fixed seed and engine build.
6. **Independent integration tests** run compact models from
   `tests/integration/models/` through the public CLI, both engines, real
   Storage, pyslimmc, and pyslimmc-opt. Session fixtures share ten compact
   engine scenarios, real completed/failed/interrupted outputs, a three-trial
   optimization, and focused optimization-surface runs across 92 tests.

## Chemistry and behavior phases

- **Phase A:** initiation, propagation, depropagation, combination,
  disproportionation, transfer-to-monomer, and material/radical balances.
- **Phase B:** exact propensity scaling, channel probability, waiting-time law,
  and fixed-seed replay.
- **Phase C:** `at`, `every`, `from`, condition cadence/order, time barriers,
  output actions, parameter/concentration actions, feeds, and `stop`.
- **Phase D:** regulator/H-donor transfer, reinitiation, capping, elementary
  reactions, and productive/nonproductive accounting.
- **Phase E:** both mass bases, chain identity, mixed `save`/`save_chains`,
  dense series, semibatch volumes/balances, finalization, failure, and SIGINT.

## Targeted suites

The homo/copo equivalence suite compares initiation, propagation,
depropagation, combination, and disproportionation using reactive monomer A and
zero-concentration spectator B in copo. It compares invariants and ensemble
summaries rather than forcing identical RNG trajectories.

The copo depropagation suite covers AA, AB, BA, and BB terminal cases, returned
monomer, pool movement, DP=1 protection, event balance, and rate scaling. The
terminal/penultimate/microstructure suite covers four terminal and eight
penultimate propagation transitions, dyads, triads, blocks, combination
boundaries, and parity between sequence modes where information overlaps.

## Required integration invariants

Where applicable, an integration test verifies:

- monotonic continuous snapshot IDs and aligned `(snapshot_id, event, time)`;
- no artificial initial snapshot or initial zero in a series;
- a true final snapshot only after successful finalization;
- chain DP equals summed monomer composition;
- live chains occupy compatible active pools;
- sequences, composition, terminal fields, masses, and moments agree;
- cumulative firing counters are monotonic and sum to the snapshot event;
- physical volume, KMC volume, feed amounts, and balances agree;
- advertised optional outputs exist and unavailable outputs raise explicitly;
- `RESULTS_COMPLETE`, `.work`, status, exit code, and hashes are consistent.


## Gate hierarchy and change policy

Not every failing check means the same thing. A developer or coding agent must
first identify the contract protected by the failure rather than editing the
expected value until the suite turns green.

| Gate / suite | What it protects | When expectations may change |
|---|---|---|
| `check-versions` | family/component version consistency | only with an intentional version change |
| `check-release-config` / provenance checks | reproducible release configuration, clean Git provenance, target platform contract | only with an intentional release-policy change |
| documentation contract | live links, parser-token coverage, public API/signature coverage, checked model blocks | together with the corresponding public contract |
| Nim unit tests | local parser, propensity, mutation, arithmetic, and bookkeeping invariants | only when the underlying contract deliberately changes |
| phases A–E / validation suites | scientific, SSA, scheduling, Storage, balance, and finalization invariants | only after a reviewed scientific/format contract change |
| homo–copo equivalence | shared behavior represented by both engines | never by changing only one side's expected result without explaining the intended divergence |
| pyslimmc / pyslimmc-opt tests | public Python behavior | together with the public API and documentation |
| integration tests | end-to-end CLI → engine → Storage → Python contracts | only with an intentional public contract change |

Release/provenance checks must not be weakened to make a dirty or untraceable
build pass. Scientific validation fixtures must not be rewritten solely because
an implementation refactor changes their result. If a deliberate scientific
change really requires a new expected value, the commit must state the changed
contract and update its theory/reference documentation at the same time.

### Before changing an expected value

Ask, in order:

1. Is the implementation wrong while the documented contract is unchanged?
   Fix the implementation.
2. Is the test asserting an obsolete public contract?
   Update the implementation, test, and canonical documentation together.
3. Is the change numerical only?
   Demonstrate why the old tolerance/fixture was invalid rather than widening
   it reflexively.
4. Does the failure involve a seeded validation case, balance, propensity,
   Storage invariant, or homo/copo parity?
   Treat it as a scientific/format contract change and review it explicitly.

## Adding or changing functionality

1. Add or update a focused unit test.
2. Add a black-box model if the invariant crosses module boundaries.
3. Check it through public pyslimmc if users can observe the result.
4. Update a Storage fixture only for a deliberate format/contract change.
5. Run the affected phase, documentation check, and `make test-full`.
6. Run `make clean` before packaging.

## See also

- [`INTEGRATION_COVERAGE.md`](INTEGRATION_COVERAGE.md) — Integration coverage
- [`ARCHITECTURE.md`](ARCHITECTURE.md) — Architecture
- [`DEVELOPMENT.md`](DEVELOPMENT.md) — Development workflow
