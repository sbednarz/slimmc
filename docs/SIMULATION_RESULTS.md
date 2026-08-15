# Simulation results

A slimmc result is more than a final conversion value. A run is a structured,
auditable record containing scientific state, selected chain populations,
process history, model/build metadata, validation information and checksums.
`pyslimmc` provides the user-facing read-only view; the exact on-disk contract
is [`reference/STORAGE.md`](reference/STORAGE.md).

## Three layers of a run

### 1. Scientific result

Depending on what was saved, a run can contain:

- simulation time and KMC event count;
- concentrations and discrete counts of small species and monomers;
- monomer conversion;
- temperature and kinetic parameter states;
- physical reactor volume, moles and feed balances when `init_volume` is used;
- live, dead and total chain counts;
- chain populations at snapshots with `save_chains`;
- DP and molar-mass moments;
- copolymer composition and sequence-derived microstructure where applicable;
- channel propensities and realized firings where enabled by the model/output.

The final scientific value is not always synonymous with the last stored row.
Use the snapshot semantics below.

### 2. Run state and history

A run records its lifecycle and chronology:

- `run_id` and path identity;
- `completed`, `failed`, or `interrupted` status;
- snapshots and their reasons;
- scheduled and conditional actions;
- feed history and cumulative balances;
- the last stored snapshot;
- a true final snapshot when the engine completed/finalized one.

This lets an interrupted run remain inspectable without pretending that it is a
completed result.

### 3. Reproducibility and provenance

The Storage metadata can identify the model/run and the software that produced
it. Through `pyslimmc`, reproducibility information includes input/model,
binary and storage hashes where available, Git commit/dirty metadata embedded
by official builds, engine/CLI versions and the RNG seed.

```python
run.reproducibility.input_hash
run.reproducibility.model_hash
run.reproducibility.binary_hash
run.reproducibility.storage_hash
run.reproducibility.git_commit
run.reproducibility.git_dirty
run.reproducibility.verify()
```

Opening a run validates the structure needed to read it; full checksum and
provenance verification is deliberately opt-in through
`run.reproducibility.verify()`. Runtime analysis does not infer provenance from
the user's current Git working tree.

This provenance model is a deliberate strength of slimmc: a simulation result
can be tied to a discrete model, seed, software version/build and stored data,
rather than being only an anonymous table exported after a calculation.

## Runs and status

`pyslimmc.open(path)` opens one Storage run. `pyslimmc.scan(path)` finds run
directories recursively and can return completed, failed and interrupted runs.

Useful identity/lifecycle fields include:

```python
run.run_id
run.desc
run.status
run.termination_reason
run.engine
run.engine_version
run.cli_version
run.storage_format_version
```

Direct opening of failed/interrupted data is strict unless
`allow_incomplete=True` is explicit.

## Snapshots: first, last and final

A snapshot is a stored state at a particular simulated time and KMC event.
Common navigation is:

```python
run.first
run.last
run.final
run.at_time(...)
run.at_event(...)
run.at_conversion(...)
```

`run.last` means the last stored snapshot. `run.final` means a snapshot marked
as truly final by the engine and is unavailable when the lifecycle does not
provide one. Keeping these concepts separate prevents an interrupted run from
being silently treated as a completed calculation.

Each snapshot exposes identity/flags such as time, KMC event, reason,
`is_final`, `has_chains` and `has_sequences`.

## State, concentration and conversion

For a run-wide time series:

```python
run.t
run.conc["Sty"]
run.conv["Sty"]
run.conv.total
run.temp
```

Storage keeps discrete KMC counts as the scientific source state. Concentration
is related to count through the stochastic simulation volume. Initial declared
concentrations can require integer rounding; accepted runs report the realized
discrete state.

## KMC volume and physical reactor volume

Two volumes may appear and must not be confused:

- `run.kmc_volume` is the stochastic simulation volume used to convert between
  molecular counts and concentration;
- `run.volume` is the physical reactor volume and requires `init_volume`.

A batch model can use only `kmc_volume`. Semibatch/feed models use
`init_volume` to track physical moles and cumulative feed volume while the
stochastic representation remains controlled by `kmc_volume`.

See [`THEORY.md`](THEORY.md) and [`LIMITATIONS.md`](LIMITATIONS.md) for interpretation and convergence issues.

## Chain counts and saved chains

Every ordinary `save` stores chain-count summary columns, including live, dead
and total chain counts. It does **not** necessarily store the molecular chain
population.

`save_chains` marks a snapshot with chain data. Chain-resolved operations such
as MWD, CLD, exact chain-mass counts, SEC broadening, and many microstructure analyses require an
appropriate snapshot with chains (and, for some analyses, sequence data).

Useful navigation includes:

```python
run.snapshots_with_chains
run.first_with_chains
run.last_with_chains
run.final.chains
```

Do not interpret a missing chain table as zero chains; it can simply mean that
chains were not saved at that snapshot.

## Moments and distributions

Molar-mass and DP moments are weighted directly from the discrete stored chain
population. Typical exact quantities are:

- `Mn`, `Mw`, `Mz`, dispersity;
- `DPn`, `DPw`, `DPz`.

MWD and CLD objects provide plot/distribution representations. Histogram,
Gaussian or KDE smoothing does not redefine the exact moments. This is an
important distinction when comparing figures with reported numeric values.

For a final MWD:

```python
mwd = run.final.mwd()
print(mwd.mn, mwd.mw, mwd.mz, mwd.dispersity)
```

See [`PYSLIMMC.md`](PYSLIMMC.md) for analysis workflows and
[`THEORY.md`](THEORY.md) for definitions.

## Process history

Scheduled/conditional actions and semibatch feeds are stored as part of the
run. This permits reconstruction of important process changes such as
temperature changes, concentration/rate changes and doses. With physical
volume enabled, cumulative component amounts and feed volume support mass and
material-balance analyses.

## Validation, audit and integrity

Before relying on a result, distinguish several checks:

- Storage/readability validation: can the run be read consistently?
- scientific/mass audits: do stored populations and declared mass information
  satisfy the implemented balances/checks?
- reproducibility verification: do hashes and available provenance match?
- statistical convergence: is the stochastic population/ensemble sufficient
  for the quantity being reported?

The first three have explicit software support; the fourth is a scientific
modeling responsibility and should be assessed by sensitivity to
`kmc_volume`, seed/replicate behavior, and the observable of interest rather
than by a universal hard threshold.

## Exact storage format

Users normally work through `pyslimmc`. Developers and independent readers can
use [`reference/STORAGE.md`](reference/STORAGE.md), which is the canonical
Storage schema and source-of-truth specification.

## See also

- [`PYSLIMMC.md`](PYSLIMMC.md) — Analyse results with pyslimmc
- [`reference/STORAGE.md`](reference/STORAGE.md) — Storage format
- [`CONCEPTS.md`](CONCEPTS.md) — Core concepts
