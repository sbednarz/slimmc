# pyslimmc API tree

This page is a user-oriented map of the public pyslimmc API. It is organized
by **the object the user has in hand and what can be accessed from it**, rather
than by Python source modules or implementation classes.

The tables use four labels:

- **function** — call it with `(...)`;
- **method** — call it on an object with `(...)`;
- **attribute** — read it without `(...)`;
- **namespace/object** — an attribute that opens another part of the API.

`Returns` gives the practical result type. Public user roles such as `Run`,
`Snapshot`, `Chains`, and `ChainRecord` are used in the tree; concrete
Storage-backed implementation class names are mentioned only where useful.
Exact parameters, defaults, and Python annotations are in
[`PYSLIMMC_SIGNATURES.md`](PYSLIMMC_SIGNATURES.md). Definitions and semantics
are in [`PYSLIMMC_API.md`](PYSLIMMC_API.md).

The central object path is:

```text
sl.scan(...) -> Runs -> Run -> Snapshots -> Snapshot -> Chains -> ChainRecord
sl.open(...) ----------^
```

Not every analysis requires descending to `ChainRecord`; most run-wide work is
done directly from `Run`, and chain-resolved analyses usually start from a
chain-bearing `Snapshot` or `Chains`.

## 1. Start here: `pyslimmc`

The normal entry points are `open()` for one run and `scan()` for a collection.

```text
pyslimmc
├── open(...)               -> Run
├── scan(...)               -> Runs
├── help()                  -> str
├── report(...)             -> Report
├── available_styles()      -> tuple[str, ...]
├── get_style(...)          -> PlotStyle
├── figure_size(...)        -> tuple[float, float]
└── options                 -> global options object
```

| Name | Kind | Returns | Purpose |
|---|---|---|---|
| `sl.open(path, allow_incomplete=False)` | function | `Run` (`StorageRun`) | Open one Slimmc Storage run. |
| `sl.scan(path=".", recursive=True, skip_bad=False)` | function | `Runs` | Find runs below a directory. |
| `sl.help()` | function | `str` | Print and return a compact workflow reminder. |
| `sl.report(...)` | function | `Report` | Create an optional PDF report builder. |
| `sl.available_styles()` | function | `tuple[str, ...]` | List plotting styles. |
| `sl.get_style(name="screen")` | function | `PlotStyle` | Return one plotting-style definition. |
| `sl.figure_size(...)` | function | `tuple[float, float]` | Return a figure width and height. |
| `sl.options` | object | options object | Global analysis options, including progress display. |
| `sl.__version__` | attribute | `str` | Installed pyslimmc version. |

Most workflows therefore branch immediately into either:

```python
run = sl.open("results/run_000001")
```

or:

```python
runs = sl.scan("results")
```

## 2. Collections: `Runs`

`Runs` is the object returned by `scan()` and by selectors that intentionally
preserve a collection.

```text
Runs
├── selection
│   ├── run_id[...]          -> Run
│   ├── by_path(...)         -> Run
│   ├── one(...)             -> Run
│   ├── first()              -> Run
│   ├── match(...)           -> Runs
│   ├── filter(...)          -> Runs
│   ├── completed            -> Runs
│   ├── failed               -> Runs
│   └── interrupted          -> Runs
├── study helpers
│   ├── sweep(...)           -> Runs
│   ├── sweep_variables      -> tuple[str, ...]
│   ├── pack(...)            -> dict
│   ├── as_table()           -> Table
│   └── model_diff(...)      -> Table
└── indexes
    ├── run_id               -> RunIdIndex
    ├── prefix               -> PrefixIndex
    └── var                  -> VarIndex
```

| Name | Kind | Returns | Purpose |
|---|---|---|---|
| `runs.run_id["..."]` | namespace/index | `Run` | Exact run-ID lookup. |
| `runs.<run_id>` | attribute lookup | `Run` | Convenient run-ID lookup when the ID is a valid attribute name. |
| `runs.by_path(path)` | method | `Run` | Explicit path lookup. |
| `runs.one(**filters)` | method | `Run` | Require exactly one matching run. |
| `runs.first()` | method | `Run` | First run in deterministic collection order. |
| `runs.match(pattern)` | method | `Runs` | Shell-style glob selection on the complete `run_id`. |
| `runs.filter(...)` | method | `Runs` | Filter by engine, status, variables, path, output, and related metadata. |
| `runs.completed` | namespace/object | `Runs` | Completed runs only. |
| `runs.failed` | namespace/object | `Runs` | Failed runs only. |
| `runs.interrupted` | namespace/object | `Runs` | Interrupted runs only. |
| `runs.sweep(*variables)` | method | `Runs` | Declare variables used to describe a study/sweep. |
| `runs.sweep_variables` | attribute | `tuple[str, ...]` | Declared sweep-variable names. |
| `runs.pack(...)` | method | `dict[str, dict]` | Build a lightweight user metadata mapping. |
| `runs.as_table()` | method | `Table` | Tabular summary of the collection. |
| `runs.model_diff(...)` | method | `Table` | Compare model text across runs. |
| `runs.paths` | attribute | collection | Run paths. |
| `runs.engines` | attribute | collection | Engines represented in the collection. |
| `runs.schemas` | attribute | collection | Storage schemas represented in the collection. |

String indexing is run-ID lookup; use `by_path()` when path lookup is intended.
Selectors such as `match()`, `filter()`, and lifecycle namespaces return
`Runs`, so further selection remains composable.

### Shell-style glob selection

`runs.match(pattern)` applies a **case-sensitive shell-style glob to the full
`run_id`** and returns another `Runs` collection. It supports `*`, `?`, character
sets/ranges such as `[0-3]`, and negated sets such as `[!0]`. For example:

```python
runs.match("feed_*_T3")
runs.match("run_00?")
runs.match("run_[0-3]*")
```

This is run-ID selection after discovery; `sl.scan()` itself scans the directory
tree and does not take a glob pattern.

## 3. One simulation: `Run`

A `Run` is the central pyslimmc object. `sl.open()` returns the Storage-backed
implementation `StorageRun`, but users normally work with it simply as the run.

```text
Run
├── identity and metadata
├── chronology and snapshots
├── simulation state and conversion
├── polymer quantities
├── feeds and balances
├── chains
├── distributions
├── kinetics and event diagnostics
├── copolymerization and microstructure
├── plotting
└── validation, raw data, and summary
```

### 3.1 Identity, model, and variables

| Name | Kind | Returns | Purpose |
|---|---|---|---|
| `run.run_id` | attribute | `str` | Run identifier. |
| `run.path` | attribute | `Path` | Run directory. |
| `run.relative_dir` | attribute | `str` | Path relative to scan root when available. |
| `run.prefix` | attribute | `str` | Run-ID prefix. |
| `run.desc` | attribute | `str | None` | Model/run description. |
| `run.status` | attribute | `str` | Lifecycle status. |
| `run.termination_reason` | attribute | `str | None` | Why execution ended. |
| `run.engine` | attribute | `str` | Engine identifier. |
| `run.engine_family` | attribute | `str` | Engine family. |
| `run.engine_version` | attribute | `str` | Engine version. |
| `run.cli_version` | attribute | `str` | CLI version. |
| `run.storage_format` | attribute | `str` | Storage format name. |
| `run.storage_format_version` | attribute | `str` | Storage format version. |
| `run.kinetic_model` | attribute | model metadata | Kinetic-model metadata. |
| `run.metadata` | namespace/object | mapping | Run metadata. |
| `run.model` | namespace/object | model view | Stored model information. |
| `run.input` | namespace/object | input view | Input/provenance information. |
| `run.execution` | namespace/object | execution view | Execution metadata. |
| `run.storage` | namespace/object | storage view | Storage metadata. |
| `run.var` | namespace/object | `Variables` | Model variables. |
| `run.monomers` | namespace/object | mapping/view | Monomer definitions. |
| `run.endgroups` | namespace/object | mapping/view | End-group definitions. |
| `run.monomer_mw(name)` | method | `float` | Monomer molar mass. |
| `run.endgroup_mw(name)` | method | `float` | End-group mass contribution. |

`run.var` exposes declared model variables. Each variable has at least its
name, value, unit, and kind; for example:

```python
run.var.temperature.value
run.var["temperature"].unit
```

### 3.2 Reproducibility

```text
run.reproducibility
├── input_hash              -> str | None
├── model_hash              -> str | None
├── binary_hash             -> str | None
├── storage_hash            -> str | None
├── git_commit              -> str | None
├── git_dirty               -> bool | None
├── verify(...)             -> ReproducibilityReport
└── compare(other)          -> dict[str, str]
```

Full checksum verification is opt-in through `verify()`.

### 3.3 Snapshots and snapshot selection

`run.snapshots` is the complete `Snapshots` collection for the run. Individual
`Snapshot` objects are selected from that collection or through the convenience
selectors on `Run`.

```text
Run
├── snapshots               -> Snapshots
│   ├── ids                 -> ndarray
│   ├── time                -> ndarray
│   ├── kmc_event           -> ndarray
│   ├── first               -> Snapshot
│   ├── last                -> Snapshot
│   └── final               -> Snapshot
├── first                   -> Snapshot
├── last                    -> Snapshot
├── final                   -> Snapshot
├── first_with_chains       -> Snapshot
├── last_with_chains        -> Snapshot
├── snapshots_with_chains   -> tuple[Snapshot, ...]
├── at_snapshot(id)         -> Snapshot
├── at_time(...)            -> Snapshot
├── at_event(...)           -> Snapshot
├── at_conversion(...)      -> Snapshot
└── at_temperature(...)     -> tuple[Snapshot, ...]
```

| Name | Kind | Returns | Purpose |
|---|---|---|---|
| `run.first` | attribute | `Snapshot` | First stored snapshot. |
| `run.last` | attribute | `Snapshot` | Last stored row, regardless of final status. |
| `run.final` | attribute | `Snapshot` | True final snapshot; unavailable if none was written. |
| `run.first_with_chains` | attribute | `Snapshot` | First chain-bearing snapshot. |
| `run.last_with_chains` | attribute | `Snapshot` | Last chain-bearing snapshot. |
| `run.snapshots` | namespace/object | `Snapshots` (`StorageSnapshots`) | Snapshot collection and chronology. |
| `run.at_snapshot(id)` | method | `Snapshot` | Exact snapshot-ID lookup. |
| `run.at_time(time, method=...)` | method | `Snapshot` | Select by simulation time. |
| `run.at_event(event, method=...)` | method | `Snapshot` | Select by KMC event number. |
| `run.at_conversion(x, monomer=..., method=...)` | method | `Snapshot` | Select by conversion. |
| `run.at_temperature(T)` | method | `tuple[Snapshot, ...]` | Return all equally nearest temperature snapshots. |

Run-wide chronology arrays are available directly as `run.sid`, `run.t`, and
`run.event`; `run.sid` contains the snapshot IDs (`np.ndarray`).

### 3.4 Simulation state and conversion

```text
Run
├── state                   -> StateSeries
│   ├── names               -> tuple/list[str]
│   ├── counts[...]         -> ndarray
│   ├── moles[...]          -> ndarray
│   └── concentrations[...] -> ndarray
├── count[...]              -> ndarray
├── conc[...]               -> ndarray
├── moles[...]              -> ndarray
├── c0[...]                 -> scalar mapping/view
├── count0[...]             -> scalar mapping/view
├── moles0[...]             -> scalar mapping/view
├── kmc_volume              -> float
├── volume                  -> ndarray
├── temp                    -> ndarray
├── conv[...]               -> ndarray
├── conv.total              -> ndarray
├── f[...]                  -> ndarray
├── f0[...]                 -> scalar/array view
└── F
    ├── ins[...]            -> ndarray
    ├── int[...]            -> ndarray
    └── cum[...]            -> ndarray
```

`run.state.moles` means amounts represented in the KMC volume. Top-level
`run.moles` means physical reactor amounts and therefore requires physical
volume information. The same distinction applies at snapshot level.

### 3.5 Polymer moments and chain counts

| Name | Kind | Returns | Purpose |
|---|---|---|---|
| `run.chain_count` | namespace/object | chain-count series | Living, dead, and total chain counts at ordinary snapshots. |
| `run.chain_count.live` | attribute | `ndarray` | Living-chain count series. |
| `run.chain_count.dead` | attribute | `ndarray` | Dead-chain count series. |
| `run.chain_count.total` | attribute | `ndarray` | Total-chain count series. |
| `run.moments` | namespace/object | moments series | Stored polymer moments by population and mass model. |
| `run.dpn` | attribute | `ndarray` | Number-average DP series. |
| `run.dpw` | attribute | `ndarray` | Weight-average DP series. |
| `run.mn` | attribute | `ndarray` | Number-average molar mass series. |
| `run.mw` | attribute | `ndarray` | Weight-average molar mass series. |
| `run.mz` | attribute | `ndarray` | z-average molar mass series. |
| `run.dispersity` | attribute | `ndarray` | Molar-mass dispersity series. |

### 3.6 Feeds and physical balances

```text
Run
├── feeds                   -> Feeds
│   └── feeds[name]         -> Feed
│       ├── concentration
│       ├── fraction
│       ├── events
│       ├── volume_cum
│       └── moles_cum
├── feed_events             -> FeedEvents
└── balance                 -> Balance
    ├── initial
    ├── dosed
    ├── total
    ├── free
    ├── consumed
    └── incorporated
```

These physical amount/balance views require the corresponding physical-volume
data. A technical `set_c` can make a physical balance inapplicable rather than
turning it into a false zero.

### 3.7 Chain-bearing data

There are two useful routes:

```python
chains = run.last_with_chains.chains
# or, when appropriate:
chains = run.chains
```

Both lead to the `Chains` API described in section 5.

### 3.8 Exact counts, distributions, and SEC

| Name | Kind | Returns | Purpose |
|---|---|---|---|
| `run.mwd(...)` | method | `MolarMassDistribution` | Molecular-weight distribution from a selected snapshot. |
| `run.cld(...)` | method | `ChainLengthDistribution` | Chain-length distribution. |
| `run.dp_counts(...)` | method | `DPCounts` | Exact chain counts grouped by DP. |
| `run.mass_counts(...)` | method | `MassCounts` | Exact chain counts grouped by actual neutral molar mass. |
| `run.sec(...)` | method | `SECDistribution` | Continuous apparent SEC response in `log10(M)`. |
| `run.mwd_series(...)` / `run.cld_series(...)` | method | `DistributionGroup` | Multi-series composition with `per_series` or `combined` normalization. |

These methods normally default to `snapshot="final"`; a snapshot ID,
`"last"`, or a `Snapshot` can be selected when supported.

### 3.9 Reaction channels, firings, kinetics, and actions

```text
Run
├── channels                -> ChannelsSeries
├── firings                 -> Firings
├── kinetics                -> KineticsSeries
├── k[...]                  -> ndarray
├── actions                 -> Actions
├── event_counts            -> event-count view
└── channel_events          -> channel-event view
```

Typical access patterns include:

```python
run.channels.event_count["prop_AA"]
run.channels.interval_event_counts()
run.firings.total_fires()
run.firings.fire_shares()
run.firings.rate_shares()
run.kinetics.names
run.kinetics.rate_constants
run.k["kp_aa"]
```

`actions` is iterable; each action record exposes its trigger, conditions,
before/after values, source line, and links to generated output/snapshots.

### 3.10 Copolymerization and microstructure

```text
Run
├── copolymerization        -> Copolymerization
│   ├── capabilities
│   ├── monomer_composition()
│   ├── incremental_composition()
│   ├── cumulative_composition()
│   ├── polymer_composition()
│   ├── reactivity_ratios()
│   ├── mayo_lewis()
│   ├── compare_mayo_lewis()
│   ├── composition_drift()
│   ├── terminal_diagnostics()
│   ├── penultimate_parameters()
│   ├── penultimate_composition()
│   ├── compare_penultimate()
│   └── penultimate_diagnostics()
└── microstructure          -> Microstructure
    ├── dyads()
    ├── triads()
    ├── run_lengths(...)
    ├── transition_fraction()
    ├── homodyad_fraction()
    ├── blockiness()
    └── check_sequence_consistency()
```

These namespaces are capability-dependent. Chemically inapplicable analyses
raise a specific exception rather than returning invented data.

### 3.11 Plotting shortcuts

`run.plot` is a discoverable plotting namespace:

```text
run.plot
├── conversion()
├── concentrations()
├── counts()
├── moles()
├── temperature()
├── volume()
├── mwd()
├── cld()
├── dp_counts()
├── mass_counts()
├── monomer_composition()
├── incremental_composition()
├── cumulative_composition()
├── composition_drift()
├── mayo_lewis()
├── compare_mayo_lewis()
├── composition_by_dp()
├── composition_dp_map(...)
├── composition_mass_map(...)
├── composition_map(...)
├── component_classes()
└── sequence/microstructure plotting helpers
```

Plot methods return Matplotlib-compatible plotting results according to the
individual method and accept the common style/axes conventions documented in
[`PYSLIMMC_API.md`](PYSLIMMC_API.md).

### 3.12 Validation, diagnostics, raw Storage, and summaries

```text
Run
├── validation              -> Validation
├── validate(...)           -> ValidationReport
├── mass_audit(...)         -> MassAuditResult
├── diagnostics             -> Diagnostics
│   ├── validation
│   ├── memory
│   ├── run_log
│   ├── debug_log
│   └── channel_trace
├── raw                     -> Raw
│   ├── metadata
│   ├── schema
│   ├── tables
│   ├── dictionaries
│   ├── table(...)          -> StorageTable
│   └── dictionary(...)     -> mapping
├── table(...)              -> StorageTable
├── dictionary(...)         -> dict
├── summary()               -> RunSummary
├── refresh()               -> Run
└── info()/help()           -> str
```

`run.raw` is intentionally low-level. Prefer the structured namespaces above
unless direct Storage access is needed.

## 4. One saved state: `Snapshot`

A snapshot is the scalar/single-state counterpart of many run-wide series.
Concrete objects are `StorageSnapshot` instances.

```text
Snapshot
├── identity
│   ├── id                  -> int
│   ├── t, time             -> float
│   ├── event, kmc_event    -> int
│   ├── reason_id           -> int
│   ├── reason              -> str
│   ├── is_final            -> bool
│   ├── has_chains          -> bool
│   └── has_sequences       -> bool
├── state and composition
│   ├── state               -> StateSnapshot
│   ├── count, conc, and moles
│   ├── conv
│   ├── f, f0, and F
│   ├── kmc_volume
│   ├── volume
│   └── temp
├── polymer
│   ├── moments
│   ├── dpn, dpw
│   ├── mn, mw, mz
│   └── dispersity
├── chains                  -> Chains
├── channels                -> ChannelsSnapshot
├── kinetics                -> KineticsSnapshot
├── mwd()                   -> distribution
├── cld()                   -> distribution
├── dp_counts()             -> DPCounts
├── mass_counts()           -> MassCounts
├── sec(...)                -> SECDistribution
└── validate(...)           -> ValidationReport
```

A snapshot that was not saved with chain records can still expose ordinary
state and moment information; chain-level analyses require a chain-bearing
snapshot.

## 5. Chain population: `Chains`

`Chains` is the compressed chain-record population associated with saved chain
data. The concrete implementation is `StorageChains`. A single `ChainRecord`
(`StorageChainRecord`) can represent many identical physical chains;
`record.count` is that multiplicity. There is no separate public `Chain` class.

```text
Chains
├── population selectors
│   ├── all                 -> Chains
│   ├── live                -> Chains
│   ├── dead                -> Chains
│   ├── population_activity(...) -> Chains
│   ├── pool(...)           -> Chains
│   ├── origin(...)         -> Chains
│   └── select(...)         -> Chains
├── record access
│   ├── [index]             -> ChainRecord
│   ├── iteration           -> ChainRecord
│   ├── record(id)          -> ChainRecord
│   ├── at_snapshot(id)     -> Chains
│   └── last                -> Chains
├── arrays and metadata
├── filters
├── distributions
├── composition analyses
└── full-sequence analyses
```

### 5.1 Core selectors and arrays

| Name | Kind | Returns | Purpose |
|---|---|---|---|
| `chains.all` | attribute | `Chains` | Full available population. |
| `chains.live` | attribute | `Chains` | Living chains. |
| `chains.dead` | attribute | `Chains` | Dead chains. |
| `chains.population_activity(name)` | method | `Chains` | Select by activity classification. |
| `chains.pool(name)` | method | `Chains` | Select by kinetic pool. |
| `chains.origin(name)` | namespace/object selector | `Chains` | Select by chain-record origin. |
| `chains[index]` | indexing | `ChainRecord` (`StorageChainRecord`) | Positional record lookup. |
| iteration over `chains` | iteration | `ChainRecord` (`StorageChainRecord`) | Iterate over compressed chain records. |
| `chains.record(id)` | method | `ChainRecord` (`StorageChainRecord`) | Record-ID lookup. |
| `chains.n_records` | attribute | `int` | Number of compressed records. |
| `chains.n_chains` | attribute | numeric | Multiplicity-weighted number of chains. |
| `chains.dp` | attribute | `ndarray` | Degree of polymerization per record. |
| `chains.molar_mass` | attribute | `ndarray` | Molar mass per record where available. |
| `chains.count` | attribute | `ndarray` | Record multiplicity. |
| `chains.moles` | attribute | `ndarray` | Record amount. |
| `chains.conc` | attribute | `ndarray` | Record concentration. |
| `chains.composition` | namespace/object | chain-composition view | Counts/fractions by monomer. |

`pool` and `origin` answer different questions: pool describes kinetic/state
membership, while origin describes how the stored chain record arose.

### 5.2 One compressed record: `ChainRecord`

Indexing or iterating over `Chains`, or calling `chains.record(id)`, returns a
`ChainRecord` (`StorageChainRecord`):

```text
ChainRecord
├── chain_record_id         -> int
├── snapshot_id             -> int
├── dp                      -> int
├── molar_mass              -> float
├── count                   -> int
├── moles                   -> float
├── conc                    -> float
├── population              -> str
├── pool                    -> str
├── origin                  -> str
├── left_end                -> str
├── right_end               -> str
├── composition             -> counts and fractions
└── sequence                -> tuple[str, ...]
```

`count` is essential: a `ChainRecord` is a compressed record and may stand for
more than one identical physical polymer chain.

### 5.3 Filters

All basic filters return another `Chains` object, so they can be composed:

```python
chains.where(dp_min=10, dp_max=100)
chains.where_count("A", min=10, max=30)
chains.where_fraction("A", min=0.4, max=0.6)
chains.where_component_count(min=2, max=2)
chains.where_components(("A", "B"), exact=True)
```

### 5.4 Distributions and composition analyses

```text
Chains
├── dp_counts()                 -> DPCounts
├── mass_counts(...)            -> MassCounts
├── mwd(...)                    -> MolarMassDistribution
├── cld(...)                    -> ChainLengthDistribution
├── sec(...)                    -> SECDistribution
├── dp_counts()                 -> DPCounts
├── mass_counts(...)            -> MassCounts
├── sec(...)                    -> SECDistribution
├── mwd_series(...)             -> DistributionGroup
├── cld_series(...)             -> DistributionGroup
├── composition_by_dp(...)      -> CompositionByDP
├── composition_dp_map(...)     -> CompositionMap
├── composition_mass_map(...)   -> CompositionMap
├── composition_map(...)        -> CompositionMap
└── component_classes()         -> ComponentClasses
```

These result objects normally expose their numerical arrays plus export and/or
plot methods. See their entries in the public API reference for the supported
operations.

### 5.5 Full-sequence analyses

When complete sequences were stored, `Chains` additionally exposes methods for
motifs, blocks, transitions, n-grams, positional profiles, dyads/triads, and
microstructure maps, including:

```python
chains.sequence_stats()
chains.contains_motif("A|B")
chains.block_lengths("A")
chains.block_count("A")
chains.junction_positions("A", "B")
chains.transition_matrix()
chains.ngrams(n=4)
chains.position_profile()
chains.microstructure_map("blockiness", monomer="A")
```

These methods are data-dependent and may require `sequence_mode full`.

## 6. Returned analysis objects

Many methods return another public result object instead of a bare NumPy array.
For example:

```python
result = run.mwd()
```

Common result families are:

| Result family | Typical source | Typical user-facing contents |
|---|---|---|
| `ChainLengthDistribution` | `run.mwd()`, `run.cld()`, chain equivalents | `x`, `y`, exact moments, export, plot. |
| `DPCounts` | `run.dp_counts()` | exact DP/count projection, totals, export/plot. |
| `MassCounts` | `run.mass_counts()` | exact mass/count projection, mass model, export/plot. |
| `SECDistribution` | `run.sec()` | continuous apparent `log10(M)` SEC density plus exact source moments. |
| `DistributionGroup` | `run.mwd_series()` / `run.cld_series()` | named exact distributions with explicit normalization. |
| `CompositionByDP` | `chains.composition_by_dp()` | composition statistics versus DP, plot. |
| `CompositionMap` | composition-map methods | 2-D map data, plot. |
| `ComponentClasses` | `chains.component_classes()` | chain component-class summary. |
| composition-series objects | `run.copolymerization.*()` | arrays plus snapshot/index selection, export/plot where applicable. |
| diagnostics objects | copolymerization diagnostics | structured diagnostic values plus `info()`. |
| `MassAuditResult` | `run.mass_audit()` | pass/fail status, entries, warnings, table, failure check. |
| `ValidationReport` | `run.validate()` | structured validation result. |
| `RunSummary` | `run.summary()` | text/dict/JSON serialization and file output. |
| `Table` / `StorageTable` | collection/raw/table routes | columns, shape, rows, NumPy conversion. |

The exact result class can depend on whether one or several series are
requested. The callable-signature inventory is the authoritative list of
concrete result classes.

## 7. Reports and plotting styles

These are useful, but they are not the core data-access path.

### `Report`

```text
Report
├── text(...)               -> Report
├── text_raw(...)           -> Report
├── math(...)               -> Report
├── vspace(...)             -> Report
├── plot(...)               -> Report
├── page_break()            -> Report
└── save(path)              -> Path
```

`Report` is an optional PDF report builder. It is not a GUI and is not required
for analysis.

### `PlotStyle`

Create/read styles through the package functions:

```python
sl.available_styles()
style = sl.get_style("screen")
size = sl.figure_size("screen", span="double")
```

## 8. Uniform object inspection

A large part of the public API follows a common inspection contract:

```text
object.info_text()   -> str    # where implemented
object.info()        -> str    # prints and returns information
object.help()        -> str    # supported operations / usage help
```

This is especially useful when an object is data-dependent and the available
analysis is not obvious from its name alone.

## 9. Recommended navigation path

For a new user, the shortest mental model is:

```text
sl.scan(...) -> Runs
                 │
                 ├── filter/match/status -> Runs
                 │
                 └── first/one/run_id    -> Run
                                            │
                    ┌───────────────────────┼────────────────────────┐
                    │                       │                        │
                 Snapshot                Chains                 analysis
              run.final / at_*      snap.chains / run.chains   mwd/cld/...
                    │                       │                        │
                 scalar state          filters / sequence       result object
                                                                    │
                                                               plot/export
```

Or, for one known run:

```text
sl.open(path) -> Run -> Snapshots -> Snapshot -> Chains -> ChainRecord
                  └──────────────────────────────-> analysis results
                  └──────────────────────────────-> run-wide series and diagnostics
```

This tree is intentionally organized by **what the user has in hand**, not by
Python source module or implementation class.

## See also

- [pyslimmc guide](../PYSLIMMC.md) — practical analysis workflow and examples.
- [pyslimmc public API](PYSLIMMC_API.md) — semantics of the supported API.
- [pyslimmc callable signatures](PYSLIMMC_SIGNATURES.md) — exact parameters, defaults, and return annotations.
- [Slimmc Storage](STORAGE.md) — files and schema underlying the read-only API.
