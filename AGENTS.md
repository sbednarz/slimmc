# AGENTS.md

This file is a compact change-safety guide for coding agents and contributors.
It does not replace the canonical documentation linked below.

## Start here

Before changing code, identify which contract the change touches:

- model syntax: `docs/MODEL_SYNTAX.md` and `docs/reference/HOMO.md` / `COPO.md`
- scientific equations and conventions: `docs/THEORY.md`
- supported scope: `docs/LIMITATIONS.md`
- Storage: `docs/reference/STORAGE.md`
- public Python API: `docs/reference/PYSLIMMC_API.md` and signature inventories
- architecture/ownership: `docs/development/ARCHITECTURE.md`
- change workflow: `docs/development/DEVELOPMENT.md`
- tests and protected invariants: `docs/development/TESTING.md`

## Repository ownership

- `cli/`: unified public dispatcher and summary executable
- `common/`: contracts shared by engines
- `homo/`: homo parser, state, SSA, and writer
- `copo/`: copo parser, pools, SSA, microstructure, and writer
- `pyslimmc/`: read-only analysis of Storage runs
- `pyslimmc_opt/`: optimization built on public model/run interfaces
- `tests/` plus engine test trees: regression, validation, and integration
- `docs/`: canonical documentation

Do not move ownership across these boundaries casually. In particular,
pyslimmc must not silently repair or mutate canonical Storage results.

## Scientific contracts

Do not change a scientific expected value merely to make a test pass. Treat the
following as protected contracts unless the task explicitly changes them:

- SSA propensity definitions and waiting-time/channel-selection laws;
- factor-of-two conventions for bimolecular termination;
- material, radical, chain, and composition balances;
- homo/copo shared-mechanism equivalence where applicable;
- terminal/penultimate pool compatibility and sequence invariants;
- fixed-seed reproducibility within the documented build contract;
- Storage row alignment, identity, finalization, hashes, and status semantics.

A deliberate change to one of these requires a focused test and an update to
the relevant theory/reference documentation in the same change.

## Eligibility-counter safety

The homo and copo engines maintain incremental eligibility counters for the SSA
hot path. They are caches over canonical chain pools. Live-pool mutations must
go through the engine tracking helpers described in
`docs/development/ARCHITECTURE.md`; preserve the exact-scan fallback for stale
or manually mutated states. Never optimize by assuming a counter is current
unconditionally.

## Public API and format rules

- Public model syntax must be present in the exact engine reference and parser
  coverage check.
- Public Python exports must appear in the API reference and generated
  signatures inventory.
- Do not hand-edit generated signature inventories; run
  `python scripts/update_api_signatures.py`.
- A Storage schema change requires compatible writer/reader/validator/test and
  documentation changes; do not reuse the existing format contract for an
  incompatible layout.
- Canonical Storage uses typed `.npy` plus JSON/JSONL; do not introduce pickle
  into the durable results format.

## Tests

For ordinary changes run at least:

```bash
make test
```

Use the targeted phases/suites listed in `docs/development/TESTING.md` while
developing. Engine/scientific/Storage changes normally require `make test-full`
before release. Run `make clean` before packaging.

If a validation or provenance gate fails, diagnose the protected contract
before modifying the gate or fixture. Release provenance checks must not be
weakened to accept dirty or untraceable builds.

Concrete red flags:

- a refactor only passes after changing an expected value in `phase_*`;
- a numerical regression only passes after widening a tolerance without a
  quantified reason;
- homo/copo equivalence is restored by changing only one side's fixture;
- an optimization removes the exact-scan fallback or bypasses a chain-pool
  tracking helper;
- a Storage/API change makes an old reader silently reinterpret existing data.

Any of these should trigger contract-level review before the change is merged.

## Documentation

Public documentation is English. Keep one canonical home for each concept:
user guides explain tasks, reference files specify exact syntax/API, THEORY
defines equations, and development documents describe implementation/testing.
Examples and literature studies are maintained outside the core repository and
must not become core regression dependencies.
