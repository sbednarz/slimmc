# pyslimmc callable signatures

This is the exhaustive callable inventory for objects intentionally returned
to users by `pyslimmc`. It is generated from the installed source by
`scripts/update_api_signatures.py`; CI rejects a stale inventory. The tutorial
and semantic reference remain in [`../PYSLIMMC.md`](../PYSLIMMC.md) and
[`PYSLIMMC_API.md`](PYSLIMMC_API.md).

Classes exported at package root are also typing/inspection names. Construct
runs and collections with `open()` and `scan()`; do not instantiate their
Storage implementation classes directly.

## Shared parameter contract

| Parameter | Accepted values and meaning |
|---|---|
| `snapshot` | `"final"` (default), `"last"`, an integer snapshot ID, or a `StorageSnapshot`. |
| `pool` | `"all"`, `"live"`, `"dead"`, or a kinetic pool name; a sequence requests a grouped result where supported. |
| `series` | Mapping/name-to-population selectors for several distributions; mutually exclusive with a non-`all` `pool`. |
| `mass_model` | `"repeat_units"` or `"with_end_groups"`; `None` uses stored/default mass semantics. |
| `progress` | `None` uses `pyslimmc.options.progress`; `True` forces and `False` suppresses progress. |
| `method` | Distribution representation: `"sticks"`, `"hist"`, `"gaussian"`, or `"kde"`. |
| `basis` | `"number"` or `"mass"`. |
| `coordinate` | `"linear"` or `"log10"`. |
| `output` | `"amount"`, `"fraction"`, or `"density"`. |
| `normalization` | `"absolute"`, `"per_series"`, `"combined"`, or `"reference"`. |
| `bins`, `bin_width` | Alternative grid controls; they are mutually exclusive and must be positive. |
| `sigma` | Positive smoothing width for Gaussian/KDE methods; units follow `coordinate`. |
| `grid_step` | Positive output-grid step; units follow `coordinate`. |
| `reference` | Reference series name when `normalization="reference"`. |
| `style` | Plot style from `available_styles()`; default `"screen"`. |
| `span` | `None`, `"column"`, or `"double"`; controls owned figure geometry. |
| `ax` | Existing Matplotlib axes. If supplied, `span` must be omitted. |
| `path` | Optional output path. Methods return the written `Path` where applicable. |
| `dpi` | Positive raster resolution, normally `300`. |
| `metadata` | Export metadata mode; `"comments"` writes a commented header. |
| `layout` | Export layout `"wide"` or `"long"` where supported. |

For MWD the defaults are `basis="mass"`, `method="gaussian"`,
`coordinate="log10"`, `output="density"`, and
`normalization="per_series"`. For CLD they are `basis="number"`,
`method="sticks"`, `coordinate="linear"`, `output="fraction"`, and
`normalization="per_series"`. A neutral chain-mass spectrum uses exact sticks;
`normalize` is `"count"`, `"fraction"`, or `"base_peak"`.

## Top-level callables
- `pyslimmc.open(path: 'str | Path', *, allow_incomplete: 'bool' = False) -> 'Run'`
- `pyslimmc.scan(path: 'str | Path' = '.', *, recursive: 'bool' = True, skip_bad: 'bool' = False) -> 'Runs'`
- `pyslimmc.help() -> 'str'`
- `pyslimmc.Run(path: 'Path', _metadata: 'dict[str, Any]' = <factory>, run_id: 'str' = '', relative_dir: 'str' = '.', _prefix: 'str' = '') -> None`
- `pyslimmc.Variable(kind: 'str', name: 'str', value: 'float', unit: 'str') -> None`
- `pyslimmc.Variables(records: 'list[dict[str, Any]] | tuple[dict[str, Any], ...] | None' = None)`
- `pyslimmc.Runs(root: 'Path', _runs: 'dict[str, Run]', _sweep_variables: 'tuple[str, ...]' = <factory>, _order_paths: 'tuple[str, ...] | None' = None, _selection_note: 'str | None' = None) -> None`
- `pyslimmc.SelectionError(signature unavailable)`
- `pyslimmc.Report(title: 'str | None' = None, *, orientation: 'str' = 'portrait', font: 'str' = 'DejaVu Sans', font_size: 'float' = 10.0, math_font: 'str' = 'dejavusans')`
- `pyslimmc.report(*args: 'Any', **kwargs: 'Any') -> 'Report'`
- `pyslimmc.PlotStyle(name: 'str', palette: 'tuple[str, ...]', foreground: 'str', background: 'str', axes_background: 'str', grid_color: 'str', font_family: 'str', font_size: 'float', title_size: 'float', line_width: 'float', grid: 'bool', column_size: 'tuple[float, float]', double_size: 'tuple[float, float]', default_span: 'str') -> None`
- `pyslimmc.available_styles() -> 'tuple[str, ...]'`
- `pyslimmc.get_style(name: 'str' = 'screen') -> 'PlotStyle'`
- `pyslimmc.figure_size(name: 'str' = 'screen', *, span: 'str | None' = None) -> 'tuple[float, float]'`
- `pyslimmc.MassAuditResult(ok: 'bool', mass_model: 'str', checked_records: 'int' = 0, checked_chains: 'int' = 0, missing_monomers: 'tuple[str, ...]' = (), missing_endgroups: 'tuple[str, ...]' = (), implicit_zero_monomers: 'tuple[str, ...]' = (), implicit_zero_endgroups: 'tuple[str, ...]' = (), warnings: 'tuple[str, ...]' = (), details: 'Any' = None, entries: 'tuple[MassEntry, ...]' = ()) -> None`
- `pyslimmc.PyslimmcError(signature unavailable)`
- `pyslimmc.FeatureUnavailableError(signature unavailable)`
- `pyslimmc.ChemicalAnalysisNotApplicableError(signature unavailable)`
- `pyslimmc.AnalysisNotApplicableError(signature unavailable)`
- `pyslimmc.ChemicalModelIncompatibleError(signature unavailable)`
- `pyslimmc.DataUnavailableError(signature unavailable)`
- `pyslimmc.IncompleteSequenceDataError(signature unavailable)`
- `pyslimmc.InvalidOutputError(signature unavailable)`
- `pyslimmc.ValidationFailedError(signature unavailable)`
- `pyslimmc.NumericalAnalysisError(signature unavailable)`
- `pyslimmc.DataConsistencyError(signature unavailable)`
- `pyslimmc.UnknownColumnError(name: 'str', context: 'str', available: 'list[str]')`
- `pyslimmc.UnknownMonomerError(name: 'str', available: 'list[str]')`
- `pyslimmc.UnsupportedChainSchema(signature unavailable)`
- `pyslimmc.SnapshotUnavailableError(signature unavailable)`
- `pyslimmc.FinalSnapshotUnavailableError(signature unavailable)`
- `pyslimmc.MassModelUnavailableError(signature unavailable)`
- `pyslimmc.InvalidDistributionConfigurationError(signature unavailable)`

## Object callables

### `MassAuditResult`

- `MassAuditResult.as_table()`
- `MassAuditResult.entry(name: 'str') -> 'MassEntry | None'`
- `MassAuditResult.info() -> 'str'`
- `MassAuditResult.raise_if_failed() -> 'None'`

### `Reproducibility`

Properties: `binary_hash`, `git_commit`, `git_dirty`, `input_hash`, `model_hash`, `storage_hash`.

- `Reproducibility.compare(other) -> 'dict[str, str]'`
- `Reproducibility.info() -> 'str'`
- `Reproducibility.verify(binary: 'str | Path | None' = None) -> 'ReproducibilityReport'`

### `ReproducibilityReport`

Properties: `ok`.

- `ReproducibilityReport.info() -> 'str'`
- `ReproducibilityReport.info_text() -> 'str'`

### `Run`

Properties: `endgroups`, `engine`, `engine_family`, `execution`, `input`, `metadata`, `model`, `monomers`, `prefix`, `reproducibility`, `schema`, `status`, `storage`, `termination_reason`, `var`, `version`.

- `Run.endgroup_mw(name: 'str') -> 'float'`
- `Run.help() -> 'str'`
- `Run.monomer_mw(name: 'str') -> 'float'`

### `Variables`

- `Variables.help() -> 'str'`
- `Variables.info() -> 'str'`
- `Variables.info_text() -> 'str'`

### `PrefixIndex`

- `PrefixIndex.keys() -> 'tuple[str, ...]'`

### `RunIdIndex`

- `RunIdIndex.keys() -> 'tuple[str, ...]'`

### `Runs`

Properties: `completed`, `engines`, `failed`, `interrupted`, `paths`, `prefix`, `run_id`, `schemas`, `sweep_variables`, `var`.

- `Runs.as_table() -> 'Table'`
- `Runs.by_path(path: 'str | Path') -> 'Run'`
- `Runs.filter(*, engine: 'str | None' = None, schema: 'str | None' = None, version: 'str | None' = None, model_class: 'str | None' = None, has_output: 'str | None' = None, path: 'str | None' = None, prefix: 'str | None' = None, var_name: 'str | None' = None, var_value: 'Any' = None, run_id: 'str | None' = None, status: 'str | set[str] | tuple[str, ...] | list[str] | None' = None) -> "'Runs'"`
- `Runs.first() -> 'Run'`
- `Runs.help() -> 'str'`
- `Runs.info(max_rows: 'int' = 10) -> 'str'`
- `Runs.info_text(max_rows: 'int' = 10) -> 'str'`
- `Runs.match(pattern: 'str') -> "'Runs'"`
- `Runs.model_diff(*, include_same: 'bool' = False) -> 'Table'`
- `Runs.one(**filters: 'Any') -> 'Run'`
- `Runs.pack(*, key: 'str | None' = None, **fields: 'Any') -> 'dict[str, dict[str, Any]]'`
- `Runs.sweep(*variables: 'str') -> "'Runs'"`

### `VarIndex`

- `VarIndex.keys() -> 'tuple[str, ...]'`

### `ConversionSeries`

- `ConversionSeries.help() -> 'str'`
- `ConversionSeries.info() -> 'str'`

### `PolymerCompositionSeries`

Properties: `cum`, `cumulative`, `ins`, `instantaneous`, `int`, `interval`.

- `PolymerCompositionSeries.help() -> 'str'`
- `PolymerCompositionSeries.info() -> 'str'`

### `SeriesView`

Properties: `dtype`, `shape`.

- `SeriesView.help() -> 'str'`
- `SeriesView.info() -> 'str'`
- `SeriesView.info_text() -> 'str'`
- `SeriesView.to_numpy(*, copy: 'bool' = False) -> 'np.ndarray'`

### `StorageAction`

Properties: `after_value`, `before_value`, `conditions`, `event`, `id`, `kinetic_parameter_set_id`, `message`, `output_written`, `snapshot`, `source_line`, `state_changed`, `t`, `trigger`, `type`.

- `StorageAction.help() -> 'str'`
- `StorageAction.info() -> 'str'`

### `StorageActions`

Properties: `raw`.

- `StorageActions.help() -> 'str'`
- `StorageActions.info() -> 'str'`
- `StorageActions.info_text() -> 'str'`

### `StorageBalance`

Properties: `consumed`, `dosed`, `free`, `incorporated`, `initial`, `names`, `total`.


### `StorageChainComposition`

Properties: `counts`, `fractions`, `matrix`, `names`.

- `StorageChainComposition.monomer_id(name: 'str') -> 'int'`

### `StorageChainCountSeries`

Properties: `dead`, `live`, `total`.


### `StorageChainOrigin`

- `StorageChainOrigin.summary()`

### `StorageChainRecord`

Properties: `chain_record_id`, `composition`, `conc`, `count`, `dp`, `first_monomer`, `last_monomer`, `left_end`, `molar_mass`, `moles`, `origin`, `penultimate_monomer`, `pool`, `population`, `right_end`, `sequence`, `snapshot_id`.


### `StorageChains`

Properties: `all`, `chain_record_id`, `component_count`, `composition`, `conc`, `count`, `dead`, `dp`, `first_monomer`, `has_sequences`, `kmc_event`, `last`, `last_monomer`, `left_end`, `live`, `molar_mass`, `moles`, `n_chains`, `n_records`, `origin`, `origin_names`, `penultimate_monomer`, `pool_names`, `population_activity_names`, `right_end`, `sequences`, `snapshot_id`, `t`.

- `StorageChains.at_snapshot(snapshot_id: 'int') -> "'StorageChains'"`
- `StorageChains.block_count(monomer=None, *, progress=None)`
- `StorageChains.block_lengths(monomer=None, *, progress=None)`
- `StorageChains.chain_mass_spectrum(*, mass_model: 'str | None' = None, series=None, normalize: 'str' = 'count', **kwargs)`
- `StorageChains.cld(*, mass_model: 'str | None' = None, series=None, **kwargs)`
- `StorageChains.component_classes()`
- `StorageChains.composition_by_dp(*, bins=None)`
- `StorageChains.composition_dp_map(monomer: 'str', *, dp_bins=None, fraction_bins=None)`
- `StorageChains.composition_map(x: 'str', y: 'str', *, bins=None)`
- `StorageChains.composition_mass_map(monomer: 'str', *, mass_model='with_end_groups', mass_bins=None, fraction_bins=None)`
- `StorageChains.contains_motif(motif, *, min_occurrences=1)`
- `StorageChains.dyads_by_dp(*, bins=16)`
- `StorageChains.ends_with(motif)`
- `StorageChains.help() -> 'str'`
- `StorageChains.info() -> 'str'`
- `StorageChains.info_text() -> 'str'`
- `StorageChains.junction_position(left: 'str', right: 'str') -> 'np.ndarray'`
- `StorageChains.junction_positions(left: 'str', right: 'str') -> 'tuple[tuple[int, ...], ...]'`
- `StorageChains.masses(*, mass_model: 'str' = 'with_end_groups') -> 'np.ndarray'`
- `StorageChains.microstructure_by_dp(statistic, *, monomer=None, bins=None, progress=None)`
- `StorageChains.microstructure_map(statistic, *, monomer=None, dp_bins=None, value_bins=None, progress=None)`
- `StorageChains.motif_counts(motif, *, progress=None)`
- `StorageChains.mwd(*, mass_model: 'str | None' = None, series=None, **kwargs)`
- `StorageChains.ngrams(n=4, *, min_count=1, progress=None)`
- `StorageChains.pool(name: 'str') -> "'StorageChains'"`
- `StorageChains.population_activity(name: 'str') -> "'StorageChains'"`
- `StorageChains.position_profile(*, bins=20, progress=None)`
- `StorageChains.record(chain_record_id: 'int') -> 'StorageChainRecord'`
- `StorageChains.select(*, pool: 'str') -> "'StorageChains'"`
- `StorageChains.sequence_stats(*, progress=None)`
- `StorageChains.starts_with(motif)`
- `StorageChains.transition_matrix(*, normalize=None, progress=None)`
- `StorageChains.triads_by_composition(monomer: 'str', *, bins=12)`
- `StorageChains.where(*, dp_min: 'int | None' = None, dp_max: 'int | None' = None) -> "'StorageChains'"`
- `StorageChains.where_block_count(monomer, *, min=None, max=None, progress=None)`
- `StorageChains.where_component_count(*, min: 'int | None' = None, max: 'int | None' = None) -> "'StorageChains'"`
- `StorageChains.where_components(components, *, exact: 'bool' = True) -> "'StorageChains'"`
- `StorageChains.where_count(monomer: 'str', *, min: 'int | None' = None, max: 'int | None' = None) -> "'StorageChains'"`
- `StorageChains.where_fraction(monomer: 'str', *, min: 'float | None' = None, max: 'float | None' = None) -> "'StorageChains'"`
- `StorageChains.where_max_block(monomer, *, min=None, max=None, progress=None)`
- `StorageChains.where_transition_count(*, min=None, max=None, progress=None)`
- `StorageChains.where_transition_fraction(*, min=None, max=None, progress=None)`

### `StorageChannelTrace`

Properties: `channel`, `channel_id`, `complete`, `dt`, `enabled`, `kmc_event`, `limit`, `propensity`, `rate`, `raw`, `t`, `total_propensity`, `truncated`.

- `StorageChannelTrace.by_channel(name: 'str') -> 'StorageTable'`
- `StorageChannelTrace.channel_counts() -> 'dict[str, int]'`
- `StorageChannelTrace.info() -> 'str'`
- `StorageChannelTrace.info_text() -> 'str'`

### `StorageChannelsSeries`

Properties: `event_count`, `nonproductive`, `nonproductive_event_count`, `productive`, `productive_event_count`, `raw`.

- `StorageChannelsSeries.fire_shares()`
- `StorageChannelsSeries.help() -> 'str'`
- `StorageChannelsSeries.info() -> 'str'`
- `StorageChannelsSeries.info_text() -> 'str'`
- `StorageChannelsSeries.interval_event_counts()`

### `StorageChannelsSnapshot`

Properties: `event_count`, `nonproductive`, `productive`.

- `StorageChannelsSnapshot.help() -> 'str'`
- `StorageChannelsSnapshot.info() -> 'str'`

### `StorageCondition`

Properties: `met`, `observable`, `observed_value`, `operator`, `threshold`.


### `StorageDiagnostics`

Properties: `channel_trace`, `debug_log`, `memory`, `run_log`, `validation`.

- `StorageDiagnostics.help() -> 'str'`
- `StorageDiagnostics.info() -> 'str'`
- `StorageDiagnostics.info_text() -> 'str'`

### `StorageFeed`

Properties: `concentration`, `events`, `fraction`, `moles_cum`, `volume_cum`.


### `StorageFeedEvents`

Properties: `cumulative_amount`, `cumulative_volume`, `dose`, `dose_mL`, `kmc_volume_after`, `kmc_volume_before`, `n_events`, `raw`, `time`, `volume_after`, `volume_before`.


### `StorageFeeds`

Properties: `names`.


### `StorageKineticsSeries`

Properties: `arrhenius_A`, `arrhenius_Ea`, `definitions`, `efficiency`, `names`, `rate_constants`, `temperature`.

- `StorageKineticsSeries.by_kind(kind)`
- `StorageKineticsSeries.help() -> 'str'`
- `StorageKineticsSeries.info() -> 'str'`
- `StorageKineticsSeries.info_text() -> 'str'`
- `StorageKineticsSeries.keys()`

### `StorageKineticsSnapshot`

Properties: `k`, `raw`, `temperature`.

- `StorageKineticsSnapshot.help() -> 'str'`
- `StorageKineticsSnapshot.info() -> 'str'`

### `StorageMemory`

Properties: `columns`.


### `StorageMomentsSeries`

Properties: `all`, `dead`, `default`, `live`, `raw`.

- `StorageMomentsSeries.help() -> 'str'`
- `StorageMomentsSeries.info() -> 'str'`
- `StorageMomentsSeries.info_text() -> 'str'`
- `StorageMomentsSeries.select(*, population_scope='all', mass_basis='with_end_groups')`

### `StorageMomentsSnapshot`

Properties: `all`, `dead`, `default`, `live`, `raw`.

- `StorageMomentsSnapshot.help() -> 'str'`
- `StorageMomentsSnapshot.info() -> 'str'`
- `StorageMomentsSnapshot.select(*, population_scope='all', mass_basis='with_end_groups')`

### `StorageRaw`

Properties: `dictionaries`, `metadata`, `schema`, `tables`.

- `StorageRaw.dictionary(name: 'str') -> 'Mapping[int, Mapping[str, Any]]'`
- `StorageRaw.help() -> 'str'`
- `StorageRaw.info() -> 'str'`
- `StorageRaw.info_text() -> 'str'`
- `StorageRaw.table(name: 'str') -> 'StorageTable'`

### `StorageRun`

Properties: `F`, `actions`, `balance`, `c0`, `chain_count`, `chains`, `channel_events`, `channels`, `cli_version`, `conc`, `conv`, `copolymerization`, `count`, `count0`, `desc`, `diagnostics`, `dispersity`, `dpn`, `dpw`, `endgroups`, `engine`, `engine_version`, `event`, `event_counts`, `f`, `f0`, `feed_events`, `feeds`, `final`, `firings`, `first`, `first_with_chains`, `free_monomer_composition`, `initial_monomer_composition`, `is_complete`, `is_ok`, `k`, `kinetic_model`, `kinetics`, `kmc_volume`, `last`, `last_with_chains`, `microstructure`, `mn`, `moles`, `moles0`, `moments`, `monomer_names`, `monomers`, `mw`, `mz`, `output_status`, `plot`, `polymer_composition`, `raw`, `schema`, `sid`, `snapshots`, `snapshots_with_chains`, `state`, `status`, `storage_format`, `storage_format_version`, `t`, `tables`, `temp`, `validation`, `var`, `version`, `volume`.

- `StorageRun.at_conversion(conversion: 'float', *, monomer: 'str | None' = None, method: 'str' = 'before') -> 'StorageSnapshot'`
- `StorageRun.at_event(event: 'int', *, method: 'str' = 'before') -> 'StorageSnapshot'`
- `StorageRun.at_snapshot(snapshot_id: 'int') -> 'StorageSnapshot'`
- `StorageRun.at_temperature(temperature: 'float') -> 'tuple[StorageSnapshot, ...]'`
- `StorageRun.at_time(time: 'float', *, method: 'str' = 'before') -> 'StorageSnapshot'`
- `StorageRun.chain_counts(*, snapshot='final', pool='all', grouping: 'str' = 'dp')`
- `StorageRun.chain_mass_spectrum(*, snapshot='final', **kwargs)`
- `StorageRun.cld(*, snapshot='final', **kwargs)`
- `StorageRun.column_unit(table: 'str', column: 'str') -> 'str | None'`
- `StorageRun.dictionary(name: 'str') -> 'dict[int, dict]'`
- `StorageRun.help() -> 'str'`
- `StorageRun.info() -> 'str'`
- `StorageRun.info_text() -> 'str'`
- `StorageRun.mass_audit(*, tolerance: 'float' = 1e-09, snapshot='final', mass_model: 'str | None' = None) -> 'MassAuditResult'`
- `StorageRun.mwd(*, snapshot='final', **kwargs)`
- `StorageRun.refresh()`
- `StorageRun.summary(path: 'str | Path | None' = None)`
- `StorageRun.table(name: 'str') -> 'StorageTable'`
- `StorageRun.validate(*, strict: 'bool' = False) -> 'ValidationReport'`

### `StorageSnapshot`

Properties: `F`, `chains`, `channels`, `conc`, `conv`, `count`, `dispersity`, `dpn`, `dpw`, `endgroups`, `event`, `f`, `f0`, `free_monomer_composition`, `has_chains`, `has_sequences`, `id`, `initial_monomer_composition`, `is_final`, `k`, `kinetic_parameter_set_id`, `kinetics`, `kmc_event`, `kmc_volume`, `mn`, `moles`, `moments`, `monomer_names`, `monomers`, `mw`, `mz`, `output_status`, `polymer_composition`, `reason`, `reason_id`, `state`, `t`, `temp`, `time`, `validation`, `volume`.

- `StorageSnapshot.chain_mass_spectrum(**kwargs)`
- `StorageSnapshot.cld(**kwargs)`
- `StorageSnapshot.help() -> 'str'`
- `StorageSnapshot.info() -> 'str'`
- `StorageSnapshot.info_text() -> 'str'`
- `StorageSnapshot.mwd(**kwargs)`
- `StorageSnapshot.refresh()`
- `StorageSnapshot.validate(*, strict: 'bool' = False) -> 'ValidationReport'`

### `StorageSnapshots`

Properties: `final`, `first`, `ids`, `kmc_event`, `last`, `raw`, `time`.

- `StorageSnapshots.at_event(event: 'int', *, method: 'str' = 'before') -> 'StorageSnapshot'`
- `StorageSnapshots.at_time(time: 'float', *, method: 'str' = 'before') -> 'StorageSnapshot'`
- `StorageSnapshots.help() -> 'str'`
- `StorageSnapshots.info() -> 'str'`

### `StorageStateSeries`

Properties: `conc`, `concentrations`, `count`, `counts`, `moles`, `names`, `raw`.

- `StorageStateSeries.entity_id(name: 'str') -> 'int'`
- `StorageStateSeries.help() -> 'str'`
- `StorageStateSeries.info() -> 'str'`
- `StorageStateSeries.matrix(column: 'str') -> 'np.ndarray'`

### `StorageStateSnapshot`

Properties: `conc`, `concentrations`, `count`, `counts`, `moles`, `names`, `raw`, `run`.

- `StorageStateSnapshot.entity_id(name: 'str') -> 'int'`
- `StorageStateSnapshot.row_values(column: 'str') -> 'np.ndarray'`

### `StorageTable`

Properties: `columns`, `n_rows`.

- `StorageTable.filtered(mask: 'np.ndarray') -> "'StorageTable'"`

### `StorageValidation`

Properties: `error_count`, `failed`, `passed`, `records`, `status`, `warning_count`.

- `StorageValidation.help() -> 'str'`
- `StorageValidation.info() -> 'str'`

### `TextLog`

Properties: `exists`, `lines`, `text`.


### `ValidationCheck`

Properties: `passed`.


### `ChainCountsView`

Properties: `array`.


### `ChainEndgroups`

- `ChainEndgroups.summary()`

### `ChainPopulation`

Properties: `compressed_rows`, `count`, `counts`, `endgroups`, `total_chains`, `total_repeat_units`.

- `ChainPopulation.all() -> "'ChainPopulation'"`
- `ChainPopulation.chain_mass_spectrum(*, mass_model: 'str | None' = None, series=None, normalize: 'str' = 'count', **kwargs)`
- `ChainPopulation.cld(*, mass_model: 'str | None' = None, series=None, **kwargs)`
- `ChainPopulation.counts_total() -> 'np.ndarray'`
- `ChainPopulation.mwd(*, mass_model: 'str | None' = None, series=None, **kwargs)`
- `ChainPopulation.row(index: 'int') -> 'ChainRow'`
- `ChainPopulation.rows() -> 'tuple[ChainRow, ...]'`

### `ChainRow`

Properties: `count`, `counts`, `dp`, `endgroups`.

- `ChainRow.as_dict() -> 'dict[str, Any]'`
- `ChainRow.counts_total() -> 'int'`

### `ChainCounts`

Properties: `max_dp`, `min_dp`, `total_chains`, `total_repeat_units`, `x`, `y`.

- `ChainCounts.as_table() -> 'Table'`
- `ChainCounts.plot(*, ax=None, path: 'str | Path | None' = None, dpi: 'int' = 300, style: 'str' = 'screen', span: 'str | None' = None, **plot_kwargs)`
- `ChainCounts.to_tsv(path: 'str | Path') -> 'Path'`

### `ChainCountsGroup`

- `ChainCountsGroup.plot(*, ax=None, path: 'str | Path | None' = None, dpi: 'int' = 300, style: 'str' = 'screen', span: 'str | None' = None, **plot_kwargs)`

### `ChainLengthDistribution`

Properties: `dp_n`, `dp_w`, `dp_z`.


### `ChainMassSpectrum`

Properties: `base_peak_intensity`, `base_peak_mass`, `intensity`, `mass`.

- `ChainMassSpectrum.info() -> 'str'`

### `MultiChainMassSpectrum`

- `MultiChainMassSpectrum.plot(*args, mode: 'str' = 'overlay', **kwargs)`

### `MultiDistribution`

Properties: `basis`, `coordinate`, `data`, `is_empty`, `kmc_event`, `log10_x`, `meta`, `method`, `normalization`, `series_names`, `snapshot_id`, `time`, `x`, `y`.

- `MultiDistribution.as_table() -> 'Table'`
- `MultiDistribution.info() -> 'str'`
- `MultiDistribution.info_text() -> 'str'`
- `MultiDistribution.plot(ax=None, *, path: 'str | Path | None' = None, dpi: 'int' = 300, mode: 'str' = 'overlay', xscale: 'str | None' = None, yscale: 'str' = 'linear', display_normalization: 'str | None' = None, title: 'str | None' = None, styles: 'Mapping[str, Mapping[str, Any]] | None' = None, style: 'str' = 'screen', span: 'str | None' = None)`
- `MultiDistribution.to_tsv(path: 'str | Path', *, metadata: 'str' = 'comments', layout: 'str' = 'wide') -> 'Path'`

### `ComponentClasses`

- `ComponentClasses.plot(path=None, *, value='number_fraction', style='screen', ax=None, span=None, dpi=300, title=None)`

### `CompositionByDP`

- `CompositionByDP.plot(path=None, *, statistic='mean', interval=None, style='screen', ax=None, span=None, dpi=300, title=None)`

### `CompositionMap`

- `CompositionMap.plot(path=None, *, log=False, style='screen', ax=None, span=None, dpi=300, title=None)`

### `RunPlotNamespace`

- `RunPlotNamespace.block_lengths(monomer=None, *, snapshot='final', progress=None, **plot_kwargs)`
- `RunPlotNamespace.chain_counts(*, snapshot='final', pool='all', grouping='dp', **plot_kwargs)`
- `RunPlotNamespace.chain_mass_spectrum(*, snapshot='final', **kwargs)`
- `RunPlotNamespace.cld(*, snapshot='final', **kwargs)`
- `RunPlotNamespace.compare_mayo_lewis(*, monomer_reference='start', parameter_reference='start', **plot_kwargs)`
- `RunPlotNamespace.component_classes(*, snapshot='final', **plot_kwargs)`
- `RunPlotNamespace.composition_by_dp(*, snapshot='final', bins=None, **plot_kwargs)`
- `RunPlotNamespace.composition_dp_map(monomer: 'str', *, snapshot='final', dp_bins=None, fraction_bins=None, **plot_kwargs)`
- `RunPlotNamespace.composition_drift(*, monomer_reference='start', **plot_kwargs)`
- `RunPlotNamespace.composition_map(x: 'str', y: 'str', *, snapshot='final', bins=None, **plot_kwargs)`
- `RunPlotNamespace.composition_mass_map(monomer: 'str', *, snapshot='final', mass_model='with_end_groups', mass_bins=None, fraction_bins=None, **plot_kwargs)`
- `RunPlotNamespace.concentrations(*, x='time', entities=None, **plot_kwargs)`
- `RunPlotNamespace.conversion(*, x='time', monomers=None, total=True, **plot_kwargs)`
- `RunPlotNamespace.counts(*, x='time', entities=None, **plot_kwargs)`
- `RunPlotNamespace.cumulative_composition(*, x='conversion', **plot_kwargs)`
- `RunPlotNamespace.incremental_composition(*, x='conversion', **plot_kwargs)`
- `RunPlotNamespace.mayo_lewis(**plot_kwargs)`
- `RunPlotNamespace.microstructure_by_dp(statistic, *, snapshot='final', monomer=None, bins=None, progress=None, **plot_kwargs)`
- `RunPlotNamespace.microstructure_map(statistic, *, snapshot='final', monomer=None, dp_bins=None, value_bins=None, progress=None, **plot_kwargs)`
- `RunPlotNamespace.moles(*, x='time', entities=None, **plot_kwargs)`
- `RunPlotNamespace.monomer_composition(*, x='conversion', **plot_kwargs)`
- `RunPlotNamespace.mwd(*, snapshot='final', **kwargs)`
- `RunPlotNamespace.ngrams(n=4, *, snapshot='final', min_count=1, progress=None, **plot_kwargs)`
- `RunPlotNamespace.position_profile(*, snapshot='final', bins=20, progress=None, **plot_kwargs)`
- `RunPlotNamespace.temperature(*, x='time', **plot_kwargs)`
- `RunPlotNamespace.transition_matrix(*, snapshot='final', normalize=None, progress=None, **plot_kwargs)`
- `RunPlotNamespace.volume(*, x='time', **plot_kwargs)`

### `Capabilities`

- `Capabilities.info() -> 'str'`
- `Capabilities.items()`
- `Capabilities.keys()`

### `Capability`

Properties: `available`.

- `Capability.info() -> 'str'`

### `CompositionDrift`

- `CompositionDrift.at_index(index: 'int') -> '_Row'`
- `CompositionDrift.ending_at_snapshot(snapshot_id: 'int') -> '_Row'`
- `CompositionDrift.info() -> 'str'`
- `CompositionDrift.plot(path=None, *, style: 'str' = 'screen', ax=None, span: 'str | None' = None)`
- `CompositionDrift.to_tsv(path) -> 'Path'`

### `CompositionResult`

- `CompositionResult.info() -> 'str'`

### `IntervalCompositionSeries`

Properties: `fraction_array`, `repeat_unit_fractions`.

- `IntervalCompositionSeries.at_index(index: 'int') -> '_Row'`
- `IntervalCompositionSeries.ending_at_snapshot(snapshot_id: 'int') -> '_Row'`
- `IntervalCompositionSeries.info() -> 'str'`
- `IntervalCompositionSeries.plot(path=None, *, x: 'str' = 'conversion', style: 'str' = 'screen', ax=None, span: 'str | None' = None)`
- `IntervalCompositionSeries.to_tsv(path) -> 'Path'`

### `MayoLewisComparison`

- `MayoLewisComparison.at_index(index)`
- `MayoLewisComparison.ending_at_snapshot(snapshot_id)`
- `MayoLewisComparison.info() -> 'str'`
- `MayoLewisComparison.plot(path=None, *, style='screen', ax=None, span=None)`
- `MayoLewisComparison.to_tsv(path)`

### `MayoLewisSeries`

Properties: `final`, `fraction_array`.

- `MayoLewisSeries.at_index(index: 'int')`
- `MayoLewisSeries.at_snapshot(snapshot_id: 'int')`
- `MayoLewisSeries.info() -> 'str'`
- `MayoLewisSeries.plot(path=None, *, x='conversion', style='screen', ax=None, span=None)`
- `MayoLewisSeries.to_tsv(path)`

### `PairValues`

Properties: `array`.

- `PairValues.info() -> 'str'`
- `PairValues.items()`
- `PairValues.keys()`

### `PenultimateComparison`

- `PenultimateComparison.info() -> 'str'`

### `PenultimateCompositionSeries`

Properties: `final`, `fraction_array`.

- `PenultimateCompositionSeries.at_index(index)`
- `PenultimateCompositionSeries.info() -> 'str'`

### `PenultimateDiagnostics`

- `PenultimateDiagnostics.info() -> 'str'`

### `PenultimateParameterSeries`

Properties: `final`.

- `PenultimateParameterSeries.at_index(index)`
- `PenultimateParameterSeries.at_snapshot(snapshot_id)`
- `PenultimateParameterSeries.info() -> 'str'`

### `ReactivityRatioSeries`

Properties: `array`, `final`.

- `ReactivityRatioSeries.at_index(index: 'int')`
- `ReactivityRatioSeries.at_snapshot(snapshot_id: 'int')`
- `ReactivityRatioSeries.info() -> 'str'`

### `SnapshotCompositionSeries`

Properties: `fraction_array`, `mole_fractions`, `repeat_unit_fractions`.

- `SnapshotCompositionSeries.at_index(index: 'int') -> '_Row'`
- `SnapshotCompositionSeries.at_snapshot(snapshot_id: 'int') -> '_Row'`
- `SnapshotCompositionSeries.info() -> 'str'`
- `SnapshotCompositionSeries.plot(path=None, *, x: 'str' = 'conversion', style: 'str' = 'screen', ax=None, span: 'str | None' = None)`
- `SnapshotCompositionSeries.to_tsv(path) -> 'Path'`

### `TerminalBlockDiagnostics`

- `TerminalBlockDiagnostics.info() -> 'str'`

### `TerminalDiagnostics`

- `TerminalDiagnostics.info() -> 'str'`

### `TerminalTransitionDiagnostics`

- `TerminalTransitionDiagnostics.info() -> 'str'`

### `TripleValues`

- `TripleValues.info() -> 'str'`
- `TripleValues.items()`
- `TripleValues.keys()`

### `StorageCopolymerization`

Properties: `capabilities`.

- `StorageCopolymerization.compare_mayo_lewis(monomer_reference='start', parameter_reference='start')`
- `StorageCopolymerization.compare_penultimate(monomer_reference='start', parameter_reference='start')`
- `StorageCopolymerization.composition()`
- `StorageCopolymerization.composition_drift(monomer_reference='start')`
- `StorageCopolymerization.cumulative_composition()`
- `StorageCopolymerization.help() -> 'str'`
- `StorageCopolymerization.incremental_composition()`
- `StorageCopolymerization.info() -> 'str'`
- `StorageCopolymerization.mayo_lewis()`
- `StorageCopolymerization.monomer_composition()`
- `StorageCopolymerization.penultimate_composition()`
- `StorageCopolymerization.penultimate_diagnostics(monomer_reference='start', parameter_reference='start')`
- `StorageCopolymerization.penultimate_parameters()`
- `StorageCopolymerization.polymer_composition()`
- `StorageCopolymerization.reactivity_ratios()`
- `StorageCopolymerization.terminal_diagnostics(monomer_reference='start', parameter_reference='start')`

### `StorageFirings`

- `StorageFirings.channel_fires(channel: 'str | None' = None)`
- `StorageFirings.channels() -> 'list[str]'`
- `StorageFirings.delta_fires(channel: 'str | None' = None)`
- `StorageFirings.delta_fires_series(channel: 'str') -> 'np.ndarray'`
- `StorageFirings.final_fires() -> 'dict[str, int]'`
- `StorageFirings.final_row() -> 'dict[str, Any]'`
- `StorageFirings.fire_shares() -> 'dict[str, float]'`
- `StorageFirings.fire_shares_series() -> 'dict[str, np.ndarray]'`
- `StorageFirings.help() -> 'str'`
- `StorageFirings.info()`
- `StorageFirings.info_text()`
- `StorageFirings.propensity_shares() -> 'dict[str, float]'`
- `StorageFirings.propensity_shares_series() -> 'dict[str, np.ndarray]'`
- `StorageFirings.rate_shares() -> 'dict[str, float]'`
- `StorageFirings.rate_shares_series() -> 'dict[str, np.ndarray]'`
- `StorageFirings.rows() -> 'list[dict[str, Any]]'`
- `StorageFirings.total_fires() -> 'int'`
- `StorageFirings.validate()`

### `StorageMicrostructure`

- `StorageMicrostructure.blockiness(*, source='engine') -> 'dict[str, float]'`
- `StorageMicrostructure.check_sequence_consistency(*, snapshot=None) -> 'dict[str, bool]'`
- `StorageMicrostructure.dyads(*, source='engine', snapshot=None)`
- `StorageMicrostructure.help() -> 'str'`
- `StorageMicrostructure.homodyad_fraction(*, source='engine') -> 'float'`
- `StorageMicrostructure.info()`
- `StorageMicrostructure.info_text()`
- `StorageMicrostructure.run_lengths(monomer: 'str | None' = None, *, snapshot=None)`
- `StorageMicrostructure.transition_fraction(*, source='engine') -> 'float'`
- `StorageMicrostructure.triads(*, source='engine', snapshot=None)`

### `Report`

- `Report.help() -> 'str'`
- `Report.info() -> 'str'`
- `Report.math(expression: 'str', *, size: 'float | None' = None, align: 'str' = 'center', font: 'str | None' = None) -> "'Report'"`
- `Report.page_break() -> "'Report'"`
- `Report.plot(value: 'Any', *, height: 'float | None' = None, span: 'str | None' = None, align: 'str' = 'center', **kwargs: 'Any') -> "'Report'"`
- `Report.save(path: 'str | Path') -> 'Path'`
- `Report.text(value: 'Any', *, size: 'float | None' = None, font: 'str | None' = None, align: 'str' = 'left', weight: 'str' = 'normal') -> "'Report'"`
- `Report.text_raw(value: 'Any', *, size: 'float | None' = None, font: 'str' = 'DejaVu Sans Mono', align: 'str' = 'left', weight: 'str' = 'normal') -> "'Report'"`
- `Report.vspace(lines: 'float' = 1.0) -> "'Report'"`

### `RunSummary`

- `RunSummary.help() -> 'str'`
- `RunSummary.info() -> 'str'`
- `RunSummary.to_dict() -> 'dict[str, Any]'`
- `RunSummary.to_json(*, indent: 'int' = 2) -> 'str'`
- `RunSummary.to_text() -> 'str'`
- `RunSummary.write(path: 'str | Path') -> 'Path'`

### `Column`

Properties: `iloc`.

- `Column.to_numpy()`
- `Column.tolist() -> 'list[Any]'`

### `ColumnNames`

- `ColumnNames.tolist() -> 'list[str]'`

### `Table`

Properties: `shape`.

- `Table.equals(other: 'object') -> 'bool'`
- `Table.head(n: 'int' = 5) -> "'Table'"`
- `Table.help() -> 'str'`
- `Table.info() -> 'str'`
- `Table.info_text() -> 'str'`
- `Table.row(idx: 'int') -> 'dict[str, Any]'`
- `Table.rows() -> 'list[dict[str, Any]]'`
- `Table.tail(n: 'int' = 5) -> "'Table'"`
- `Table.to_numpy()`


## See also

- [`../PYSLIMMC.md`](../PYSLIMMC.md) — task-oriented analysis guide
- [`PYSLIMMC_API.md`](PYSLIMMC_API.md) — semantic API reference
