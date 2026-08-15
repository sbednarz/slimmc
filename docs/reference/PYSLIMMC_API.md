# pyslimmc 5.0.1 public API

This is the canonical public reference for the read-only Python interface to
Slimmc Storage 1.2.0. It documents supported names rather than implementation
classes. Numerical leaves are read-only NumPy arrays; methods that export TSV,
JSON, text, plots, or PDF reports write user products outside the Storage run.
Exact signatures for every callable returned to users are generated in
[`PYSLIMMC_SIGNATURES.md`](PYSLIMMC_SIGNATURES.md). This file defines meaning;
the generated inventory defines names, parameters, defaults, and return
annotations.

## Package entry points

```python
import pyslimmc as sl

sl.__version__
run = sl.open("results/run_000001", allow_incomplete=False)
runs = sl.scan("results", recursive=True, skip_bad=False)
sl.help()
```

`allow_incomplete=True` explicitly permits a run without
`RESULTS_COMPLETE`. `skip_bad=True` skips invalid run directories during a
scan; it does not weaken validation of runs that are returned. `scan()`
intentionally includes valid `completed`, `failed`, and `interrupted` Storage
runs so its lifecycle namespaces are useful. Direct `open()` remains strict
for a failed/interrupted path unless `allow_incomplete=True` is explicit.

Plotting and reports:

```python
sl.available_styles()
sl.get_style("screen")
sl.figure_size("screen", span="double")
report = sl.report("Title")
```

## Run identity, metadata, and reproducibility

```python
run.run_id
run.path
run.relative_dir
run.prefix
run.desc
run.status
run.termination_reason

run.engine
run.engine_family
run.engine_version
run.cli_version
run.storage_format
run.storage_format_version
run.schema
run.kinetic_model

run.metadata
run.model
run.input
run.execution
run.storage
run.var
run.monomers
run.endgroups
run.monomer_mw("A")
run.endgroup_mw("R")
```

`run.var` is a `Variables` collection. It supports iteration, string and
attribute lookup, for example `run.var.temperature.value` and
`run.var["temperature"].unit`.

Reproducibility:

```python
run.reproducibility.input_hash
run.reproducibility.model_hash
run.reproducibility.binary_hash
run.reproducibility.storage_hash
run.reproducibility.git_commit
run.reproducibility.git_dirty
run.reproducibility.verify(binary=None)
run.reproducibility.compare(other_run)
```

`pyslimmc.open()` validates the Storage structure needed to read the run, but
it does **not** hash every result file automatically. Full checksum/provenance
verification is opt-in through `run.reproducibility.verify()`. Official Slimmc
release builds embed their Git commit/tag/dirty state at compile time; runtime
analysis never infers the Slimmc source commit from the user's current working
Git repository.

## Snapshot navigation

```python
run.first
run.last
run.final

run.snapshots_with_chains
run.first_with_chains
run.last_with_chains

run.at_snapshot(snapshot_id)
run.at_time(time, method="before")
run.at_event(event, method="before")
run.at_conversion(0.5, monomer=None, method="before")
run.at_temperature(333.15)
```

`method` for time, event, and conversion selection is `before`, `after`, or
`nearest`. `run.last` is the last stored row; `run.final` exists only when the
engine marked a true final snapshot. `at_temperature()` deliberately returns a
tuple because the same or equally near temperature can occur more than once.
Its only selection policy is nearest temperature; all equally nearest
snapshots are returned in chronological order.

Snapshot identity and flags:

```python
snap.id
snap.t                 # alias: snap.time
snap.event             # alias: snap.kmc_event
snap.reason_id
snap.reason
snap.is_final
snap.has_chains
snap.has_sequences
snap.kinetic_parameter_set_id
```

Run-wide chronology:

```python
run.sid
run.t
run.event
run.snapshots.ids
run.snapshots.time
run.snapshots.kmc_event
```

## State, volume, amounts, and conversion

```python
run.state.names
run.state.counts["A"]       # alias: run.count["A"]
run.state.moles["A"]        # KMC-representation mol
run.state.concentrations["A"]
run.conc["A"]               # mol/L

run.kmc_volume               # L
run.volume                   # physical reactor L; requires init_volume
run.moles["A"]              # physical reactor mol; requires init_volume

run.c0["A"]
run.count0["A"]
run.moles0["A"]             # requires init_volume
```

Structured state views additionally expose `names`, `entity_id`, `matrix`, and
`raw`. For a snapshot, `row_values()` returns the matching dense state row.
`run.monomer_names`, `run.free_monomer_composition`, and
`run.initial_monomer_composition` are the long names behind `f` and `f0`.
`run.is_complete` and `run.is_ok` are convenience lifecycle/validation
booleans; neither substitutes for inspecting `run.status` and
`run.output_status` when diagnosing a run.

The corresponding `snap.*` accessors return scalars or mappings for one
snapshot. `state.moles` always means `count / N_A` in the KMC volume, whereas
top-level `run.moles` and `snap.moles` mean physical reactor amounts. This
distinction matters in semibatch runs.

Chain population counts are available at every ordinary `save`:

```python
run.chain_count.live
run.chain_count.dead
run.chain_count.total
```

Conversion and composition:

```python
run.conv["A"]
run.conv.total
run.f["A"]
run.f0["A"]
run.F.ins["A"]
run.F.int["A"]
run.F.cum["A"]
```

`f` is current free-monomer composition; `f0` is initial free-monomer
composition; `F.ins`, `F.int`, and `F.cum` are theoretical instantaneous,
actual interval, and cumulative incorporated-polymer compositions.

## Semibatch feeds and balances

```python
run.feeds.names
run.feeds["solution"].concentration
run.feeds["solution"].fraction
run.feeds["solution"].events
run.feeds["solution"].volume_cum
run.feeds["solution"].moles_cum

run.feed_events.time
run.feed_events.dose_mL
run.feed_events.cumulative_volume
run.feed_events.cumulative_amount
run.feed_events.volume_before
run.feed_events.volume_after
run.feed_events.kmc_volume_before
run.feed_events.kmc_volume_after
```

Physical balances, in mol:

```python
run.balance.names
run.balance.initial["A"]
run.balance.dosed["A"]
run.balance.total["A"]
run.balance.free["A"]
run.balance.consumed["A"]
run.balance.incorporated["A"]  # monomers only
```

Physical amounts and balances require `param init_volume`. A `set_c` action
makes the affected physical balance inapplicable and raises
`AnalysisNotApplicableError`; missing data are not represented by a false zero.

## Chains

Chain records exist only at `save_chains` snapshots:

```python
chains = run.last_with_chains.chains
chains.all
chains.live
chains.dead
chains.population_activity("dead")
chains.pool("PA")
chains.origin("init")
chains.record(chain_record_id)
chains.row(index)
chains.rows()
chains.at_snapshot(snapshot_id)
chains.last
```

`population_scope`, `population_activity`, kinetic `pool`, and `origin` are
distinct selectors. A **pool** is current kinetic classification/eligibility
(e.g. a terminal live pool or a dead pool). **Origin** is provenance: the
mechanism that formed the stored chain record. A chain can therefore be in one
pool while carrying an origin such as `init`, `transfer_m`, or `term_c`. Origin
does not mean “last reaction”. One row is a compressed structural record;
`count` is its multiplicity.

```python
chains.n_records
chains.n_chains
chains.total_chains
chains.total_repeat_units
chains.chain_record_id
chains.dp
chains.molar_mass
chains.count
chains.moles
chains.conc
chains.population_activity_names
chains.pool_names
chains.origin_names
chains.left_end
chains.right_end
chains.first_monomer
chains.penultimate_monomer
chains.last_monomer
chains.composition.counts["A"]
chains.composition.fractions["A"]
chains.composition.matrix
chains.component_count
chains.has_sequences
chains.sequences
chains.masses(mass_model="with_end_groups")
```

Basic filters:

```python
chains.where(dp_min=10, dp_max=100)
chains.where_count("A", min=10, max=30)
chains.where_fraction("A", min=0.4, max=0.6)
chains.where_component_count(min=2, max=2)
chains.where_components(("A", "B"), exact=True)
chains.dead
chains.pool("P_A")
```

Population activity and kinetic pool are intentionally separate concepts:
use `chains.live` / `chains.dead` for activity and `chains.pool(name)` for a
kinetic pool. `at_snapshot()` selects the chain block
for a chain-bearing snapshot; `last` selects the last such block.

Composition analyses:

```python
chains.composition_by_dp(bins=None)
chains.composition_dp_map("A", dp_bins=None, fraction_bins=None)
chains.composition_mass_map("A", mass_model="with_end_groups")
chains.composition_map("A", "B", bins=None)
chains.component_classes()
```

Full-sequence methods require `sequence_mode full` and complete sequences for
the selected records:

```python
chains.sequence_stats(progress=None)
chains.contains_motif("A|B", min_occurrences=1)
chains.starts_with("A|B")
chains.ends_with("B|A")
chains.where_transition_count(min=1, max=None)
chains.where_transition_fraction(min=0.1, max=None)
chains.where_block_count("A", min=2)
chains.where_max_block("A", min=5)
chains.block_lengths("A")
chains.block_count("A")
chains.junction_positions("A", "B")
chains.junction_position("A", "B")
chains.transition_matrix(normalize=None)
chains.dyads_by_dp(bins=16)
chains.triads_by_composition("A", bins=12)
chains.microstructure_by_dp("blockiness", monomer="A")
chains.motif_counts("A|B")
chains.ngrams(n=4, min_count=1)
chains.position_profile(bins=20)
chains.microstructure_map("blockiness", monomer="A")
```

All aggregations use compressed-row multiplicity. Expensive full-sequence
methods accept `progress=None|True|False`; global defaults are controlled by
`pyslimmc.options`.

## Moments, exact counts, distributions, and SEC

Simple moment series remain directly available:

```python
run.dpn
run.dpw
run.dp_dispersity
run.mn
run.mw
run.mz
run.dispersity
```

Explicit population moments use one callable:

```python
m = run.moments(snapshot="final", population="dead", mass_model="repeat_units")
```

`PopulationMoments` exposes `dpn`, `dpw`, `dpz`, `mn`, `mw`, `mz`,
`dp_dispersity`, `mass_dispersity`, `total_chains`, `mass_model`, `source`, and
`has_dpz`. Exact moments come from source chains or stored aggregate moments,
never from displayed curves.

Exact unnormalized source projections are:

```python
dp = run.dp_counts(snapshot="final", pool="dead")
mass = run.mass_counts(snapshot="final", pool="dead", mass_model=None)
```

`DPCounts` exposes `dp`, `count`, `total_chains`, `total_repeat_units`,
`min_dp`, `max_dp`, `to_tsv()`, and `plot()`. `MassCounts` exposes `mass`,
`count`, `total_chains`, `min_mass`, `max_mass`, `mass_model`, `to_tsv()`, and
`plot()`. These objects intentionally do not provide generic `x`/`y` aliases.

Single-population exact discrete distributions are:

```python
cld = run.cld(snapshot="final", pool="dead", weighting="number", mass_model=None)
mass_dist = run.mass_distribution(
    snapshot="final", pool="dead", weighting="mass", mass_model=None
)
```

`cld()` accepts `weighting="number"`, `"mass"`, or `"z"` and remains discrete
on integer DP support. `mass_distribution()` accepts the same weightings and
remains discrete on exact neutral-mass support. Both normalize by summation.

The reconstructed polymer MWD is:

```python
mwd = run.mwd(snapshot="final", pool="dead", mass_model=None)
```

`MolarMassDistribution.x` is `log10(M)` and `y` is the normalized density
`dW/dlog10(M)`. The reconstruction uses the documented mcPolymer-style linear
interpolation. Exact source `mn`, `mw`, `mz`, and dispersity remain attached to
the result and are not estimated from the displayed curve.

`MassDistribution` exposes `x`, `y`, `mass`, `weighting`, `mass_model`, exact
source moments, export, plot, `info()`, and `help()`.
`ChainLengthDistribution` analogously exposes `x`, `y`, `dp`, `weighting`,
`dpn`, `dpw`, `dpz`, and DP dispersity.

Multi-series composition is explicit:

```python
g = run.mwd_series(
    snapshot="final",
    series=("live", "dead"),
    normalization="per_series",
)
```

`mwd_series()` contains reconstructed logarithmic MWD densities. `cld_series()`
contains exact discrete CLDs and additionally accepts `weighting=`. Only
`per_series` and `combined` normalizations are accepted. `combined` requires
pairwise-disjoint populations.

SEC instrumental broadening is separate and acts directly on the exact mass
measure:

```python
sec = run.sec(
    snapshot="final",
    pool="dead",
    sigma_log10M=0.05,
    mass_model=None,
    step_log10M=None,
)
```

`SECDistribution.y` is the continuous apparent density
`dW_app/dlog10(M)`. `sigma_log10M` is required. SEC does not depend on the MWD
reconstruction grid. Exact source `mn`, `mw`, `mz`, and dispersity remain
attached to the result.

## Channels, firings, kinetics, and actions

```python
run.channels.event_count["prop_AA"]
run.channels.productive_event_count["init_A"]
run.channels.nonproductive_event_count["init_A"]
run.channels.interval_event_counts()

run.firings.total_fires()
run.firings.channel_fires("prop_AA")
run.firings.final_row()
run.firings.final_fires()
run.firings.delta_fires("prop_AA")
run.firings.delta_fires_series("prop_AA")
run.firings.fire_shares()
run.firings.fire_shares_series()
run.firings.rate_shares()
run.firings.rate_shares_series()
run.firings.validate()

run.temp
run.k["kp_aa"]
run.kinetics.names
run.kinetics.rate_constants
run.kinetics.arrhenius_A
run.kinetics.arrhenius_Ea
run.kinetics.efficiency
run.kinetics.by_kind("rate_constant")

run.actions
```

Iterating `run.actions` yields action records with `trigger`, `conditions`,
before/after values, output/snapshot links, and the source line. Condition
records expose `observable`, `operator`, `threshold`, `observed_value`, and
`met`. `run.event_counts` and `run.channel_events` are convenience routes to
cumulative event information.

Fire shares describe realized intervals. Rate/propensity shares describe
instantaneous competition at stored states.

## Copolymerization and microstructure

```python
cp = run.copolymerization
cp.capabilities
cp.monomer_composition()
cp.incremental_composition()
cp.cumulative_composition()
cp.polymer_composition()
cp.reactivity_ratios()
cp.mayo_lewis()
cp.compare_mayo_lewis()
cp.composition_drift()
cp.terminal_diagnostics()
cp.penultimate_parameters()
cp.penultimate_composition()
cp.compare_penultimate()
cp.penultimate_diagnostics()

run.microstructure.dyads()
run.microstructure.triads()
run.microstructure.run_lengths("A")
run.microstructure.transition_fraction()
run.microstructure.homodyad_fraction()
run.microstructure.blockiness()
run.microstructure.check_sequence_consistency()
```

`composition_drift(monomer_reference="start")` accepts `start`, `end`, or
`midpoint`. `compare_mayo_lewis()`, `terminal_diagnostics()`,
`compare_penultimate()`, and `penultimate_diagnostics()` accept both
`monomer_reference` and `parameter_reference`, each with the same three
values. The reference selects where free-monomer composition or kinetic
parameters are sampled for a saved interval.

Chemically inapplicable calls raise `ChemicalAnalysisNotApplicableError`.
Composition-series result objects support `at_index()`, `at_snapshot()`, and,
for interval-defined quantities, `ending_at_snapshot()`. Their explicit arrays
include `fraction_array`, `mole_fractions`, or `repeat_unit_fractions` as
applicable.

## Plot namespace

Every plotting shortcut accepts `style`, optional Matplotlib `ax`, `span`, and
normal plotting keywords where applicable:

```python
run.plot.conversion()
run.plot.concentrations()
run.plot.counts()
run.plot.moles()
run.plot.temperature()
run.plot.volume()
run.plot.mwd()
run.plot.cld()
run.plot.dp_counts()
run.plot.mass_counts()
run.plot.monomer_composition()
run.plot.incremental_composition()
run.plot.cumulative_composition()
run.plot.composition_drift()
run.plot.mayo_lewis()
run.plot.compare_mayo_lewis()
run.plot.composition_by_dp()
run.plot.composition_dp_map("A")
run.plot.composition_mass_map("A")
run.plot.composition_map("A", "B")
run.plot.component_classes()
run.plot.block_lengths("A")
run.plot.transition_matrix()
run.plot.microstructure_by_dp("blockiness", monomer="A")
run.plot.ngrams()
run.plot.position_profile()
run.plot.microstructure_map("blockiness", monomer="A")
```

## Collections of runs

```python
runs = sl.scan("results")
runs.run_id["run_000001"]
runs.run_000001
runs.by_path("path/to/run")
runs.one(run_id="run_000001")
runs.first()
```

Selection and status namespaces:

```python
runs.match("feed_*_T3")
runs.filter(engine="copo", status="completed", var_name="T", var_value=333.15)
runs.completed
runs.failed
runs.interrupted
```

Study helpers:

```python
runs.sweep("temperature", "feed_fraction")
runs.sweep_variables
runs.pack(key="feed_*_T3", color=[...], label=[...])
runs.as_table()
runs.model_diff(include_same=False)
```

String indexing is run-ID lookup, while `by_path()` is explicit path lookup.
Selection preserves `Runs`, declared sweep variables, and deterministic order.
`one()` requires exactly one match and raises `SelectionError` otherwise.

The complete filter signature is:

```python
runs.filter(
    engine=None,
    schema=None,
    version=None,
    model_class=None,
    has_output=None,
    path=None,
    prefix=None,
    var_name=None,
    var_value=None,
    run_id=None,
    status=None,
)
```

String filters are exact except `path` and `prefix`, which are
case-insensitive substring matches. `status` accepts one string or a
collection. `var_value` requires `var_name` and performs the library's numeric
value match. `has_output` is one name from `run.available_outputs()`.

## Validation, raw data, summaries, and reports

```python
run.output_status
run.validation
run.validate(strict=False)
run.mass_audit(snapshot="final")
run.mass_audit(snapshot="final", mass_model="with_end_groups")
run.refresh()

run.raw.metadata
run.raw.schema
run.raw.tables
run.raw.dictionaries
run.raw.table("chains")
run.raw.dictionary("monomers")
run.table("chains")
run.dictionary("monomers")
run.column_unit("state", "concentration")

summary = run.summary()
summary.to_dict()
summary.to_json()
summary.to_text()
summary.write("summary.json")
```

Diagnostics:

```python
run.diagnostics.validation
run.diagnostics.memory
run.diagnostics.run_log
run.diagnostics.debug_log
run.diagnostics.channel_trace
```

The validation view exposes `status`, `warning_count`, `error_count`, `passed`,
`failed`, and individual records. Channel trace exposes `complete`,
`truncated`, `limit`, `by_channel()`, and `channel_counts()`. Feed events also
expose `n_events`.

Tables provide `columns`, `shape`, `row()`, `rows()`, `head()`, `tail()`,
`equals()`, `to_numpy()`, `info_text()`, `info()`, and `help()`. Columns provide
`iloc`, `tolist()`, and `to_numpy()`. Storage tables also expose `n_rows` and a
read-only `filtered()` view. Storage tables are read-only; analysis tables may
be exported by the owning analysis object.

All major run, snapshot, collection, analysis, diagnostics, and table objects
follow the same inspection convention: `info_text()` returns text without
printing, `info()` prints and returns it, and `help()` describes supported
operations before data-dependent execution.

Optional PDF reports:

`Report` is a presentation helper layered on top of the analysis API. It is not a GUI and is not required for normal pyslimmc workflows.


```python
report = sl.report("Run report")
report.text(run.summary().to_text())
report.text_raw("literal text")
report.math(r"M_n")
report.vspace(6)
report.plot(run.final.mwd())
report.page_break()
report.save("report.pdf")
```

Constructor and formatting parameters:

```python
report = sl.Report(title=None, orientation="portrait",
                   font="DejaVu Sans", font_size=10.0,
                   math_font="dejavusans")
report.text(value, size=None, font=None, align="left", weight="normal")
report.text_raw(value, size=None, font="DejaVu Sans Mono",
                align="left", weight="normal")
report.math(expression, size=None, align="center", font=None)
report.vspace(lines=1.0)
report.plot(value, height=None, span=None, align="center", **plot_kwargs)
```

`orientation` is `portrait` or `landscape`; alignment is `left`, `center`, or
`right`. `plot()` accepts a callable taking axes, an object exposing `plot()`,
or a Matplotlib figure. `height` and `span` are mutually exclusive.

Global progress control is ordinary mutable configuration:

```python
import pyslimmc

pyslimmc.options.progress = True   # force progress
pyslimmc.options.progress = False  # suppress progress
```

## Top-level exports and exceptions

The public root exports are `open`, `scan`, `help`, `Run`, `Variable`,
`Variables`, `Runs`, `SelectionError`, `Report`, `report`, `PlotStyle`,
`available_styles`, `get_style`, `figure_size`, and `MassAuditResult`, plus the
documented exception hierarchy:

```text
PyslimmcError
FeatureUnavailableError
ChemicalAnalysisNotApplicableError
AnalysisNotApplicableError
ChemicalModelIncompatibleError
DataUnavailableError
IncompleteSequenceDataError
InvalidOutputError
ValidationFailedError
NumericalAnalysisError
DataConsistencyError
UnknownColumnError
UnknownMonomerError
UnsupportedChainSchema
SnapshotUnavailableError
FinalSnapshotUnavailableError
MassModelUnavailableError
InvalidDistributionConfigurationError
```

Unavailable or undefined analyses raise a specific exception; public methods do
not fabricate zeros or silently substitute `run.last` for a missing
`run.final`.

## See also

- [`../PYSLIMMC.md`](../PYSLIMMC.md) — pyslimmc guide
- [`PYSLIMMC_API_TREE.md`](PYSLIMMC_API_TREE.md) — User-oriented API tree
- [`PYSLIMMC_SIGNATURES.md`](PYSLIMMC_SIGNATURES.md) — Exact callable signatures
- [`STORAGE.md`](STORAGE.md) — Storage format
