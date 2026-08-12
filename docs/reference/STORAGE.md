# slimmc Storage 1.2.0

Current on-disk storage format version: `1.2.0`.

Status: implemented format and developer specification.

## 1. Scope

Slimmc Storage defines the canonical on-disk results format for Slimmc homo, copo, and pyslimmc.

Architecture:

```text
Engine -> SnapshotBuilder -> ResultWriter (.npy + JSON/JSONL + text logs) -> pyslimmc
```

Core rules:

- all declared scientific state must remain unambiguously recoverable from the normalized Storage tables;
- one `.model` file describes exactly one run, including seed and actions;
- one run is stored in one results directory;
- directory = table;
- `.npy` file = column;
- shared row index across columns of a table;
- `run_metadata.json` is a human-readable, multi-line JSON document; record streams use JSONL;
- pickle is forbidden;
- column table is not the primary results format;
- little-endian data is canonical;
- UTF-8 without BOM is canonical for text, JSON, and JSONL.

## 2. Run directory

Example:

```text
run_000001.model -> results/run_000001/
```

The union of current homo and copo tables is shown below. Tables marked
engine-specific are declared only by the engine that writes them; every table
declared in that run's `schema.jsonl` is required even when it has zero rows.

```text
results/
└── run_000001/
    ├── input.model
    ├── run_metadata.json
    ├── schema.jsonl
    ├── snapshots/
    ├── state/
    ├── chains/
    ├── chain_composition/          # copo
    ├── sequences/
    ├── memory/
    ├── moments/
    ├── channel_events/
    ├── channel_propensities/       # copo
    ├── microstructure_motifs/      # copo
    ├── block_statistics/           # copo
    ├── actions/
    ├── action_conditions/
    ├── feed_events/
    ├── monomer_balance/
    ├── species_balance/
    ├── kinetic_parameters/
    │   ├── sets/
    │   └── values/
    └── diagnostics/
        ├── validation.jsonl
        ├── channel_trace/
        ├── run.log
        └── debug.log          # optional
```

Declared table directories exist even when empty and contain valid empty
`.npy` columns of the declared dtype. The schema, not this union listing, is
authoritative for one run.

## 3. Scalar and identifier conventions

- record identifiers, large counters, offsets and lengths: `uint64`;
- dictionary identifiers: `uint32`;
- physical continuous values: `float64`;
- flags: NumPy `bool`;
- required dictionary IDs are normally continuous from zero;
- optional dictionary fields may use `0 = not_applicable`, with real values starting at `1`;
- optional `uint64` foreign keys use a separate `has_*` boolean and store `0` when absent;
- undefined floating-point values may use `NaN` only where explicitly permitted;
- `Inf` is invalid everywhere.

## 4. JSON and JSONL rules

`run_metadata.json` is the only multi-line JSON document in the core results format. It contains one object and is written with two-space indentation for convenient human inspection. Record streams use JSONL, including `schema.jsonl`, `actions/messages.jsonl`, and `diagnostics/validation.jsonl`.

Common rules:

- UTF-8 without BOM;
- comments are forbidden;
- key order is semantically irrelevant, but writers use stable ordering;
- `uint64` values that may exceed `2^53-1` are written as decimal strings;
- JSON and JSONL files are written atomically through temporary file + rename.

Additional JSONL rules:

- each non-empty line is one complete JSON object;
- blank lines are forbidden.

## 5. `run_metadata.json`

`run_metadata.json` contains exactly one multi-line JSON object describing the concrete run. It is formatted with two-space indentation and rewritten atomically as the run progresses and finalizes.

Required or conditionally required fields include:

```text
run_id
input_model_file
source_model_name
input_model_sha256
results_format
results_format_version
schema_sha256
seed
engine
cli_version
engine_version
git_commit? 
git_dirty?
started_at_utc
finished_at_utc?
wall_time_s?
run_status
exit_code
kmc_volume_L
initial_kmc_volume_L
initial_volume_mL?
volume_mode
avogadro_constant_mol_inv
initial_temperature_K
kinetic_model
n_monomers
platform
threads
compiler
build_mode
validation_status
validation_warning_count
validation_error_count
channel_trace_enabled
channel_trace_complete
model
variables
feeds
hashes
```

Rules:

- `run_id` must match the run directory name;
- the copied model is always `input.model`;
- full original filesystem paths are not stored;
- `seed` is logically `uint64` and is written as a decimal string if needed;
- `run_status`: `running`, `completed`, `failed`, or `interrupted`;
- `volume_mode` is `constant` or `variable`;
- `initial_kmc_volume_L` is the starting KMC representation volume;
- `initial_volume_mL` is present/positive when `param init_volume` was defined;
- `kmc_volume_L` metadata is the initial KMC volume; the snapshot column is
  authoritative after feeds;
- KMC volumes are positive and finite;
- `avogadro_constant_mol_inv = 6.02214076e23` unless explicitly versioned otherwise;
- resolved model declarations, variables, and feed definitions are stored in
  structured metadata so pyslimmc need not infer them from comments or paths;
- table contents, not cached snapshot counts in metadata, are authoritative;
- validation summary may be repeated here, while detailed checks live in `diagnostics/validation.jsonl`.

## 6. `schema.jsonl`

`schema.jsonl` describes tables, columns, dtypes, units, dictionaries, foreign keys, and named validation rules. It is frozen before the first `.npy` data write.

Recommended record order:

```text
schema_header
table
column
dictionary
dictionary_entry
rule
```

Example record types:

```json
{"record_type":"schema_header","schema_name":"slimmc-storage","schema_version":"1.2.0","byte_order":"little","column_storage":"npy","npy_format_version":"1.0"}
{"record_type":"table","name":"snapshots","required":true}
{"record_type":"column","table":"snapshots","name":"time","file":"time.npy","dtype":"float64","unit":"s","required":true}
{"record_type":"dictionary_entry","dictionary":"population_scope","id":0,"name":"all"}
{"record_type":"rule","name":"state_dense","description":"Each snapshot contains all state entities in ascending entity_id order."}
```

Unknown JSONL fields are ignored by readers for forward compatibility, optionally with a warning. Missing required fields are fatal.

## 7. `snapshots/`

Columns:

```text
snapshot_id.npy                 uint64
time.npy                        float64, s
kmc_event.npy                   uint64
snapshot_reason_id.npy          uint32
is_final.npy                    bool
has_chains.npy                  bool
has_sequences.npy               bool
kinetic_parameter_set_id.npy    uint64
volume_mL.npy                   float64, mL
kmc_volume_L.npy                float64, L
chain_count_live.npy            uint64
chain_count_dead.npy            uint64
chain_count_total.npy           uint64
```

Rules:

- `snapshot_id = row_index`;
- IDs are continuous from zero;
- there is no implicit snapshot at `t=0`; the first stored snapshot is the
  first requested output or the final state;
- a snapshot at `t=0` exists only when an action requests it or the run itself
  finishes there;
- a completed run contains exactly one `is_final=true` snapshot;
- no duplicate final snapshot is written if a scheduled snapshot already occurs at the end;
- an initial snapshot may also be final for a zero-event run;
- failed or interrupted runs may contain no final snapshot;
- `time` and `kmc_event` are nondecreasing;
- `(time, kmc_event)` need not be unique;
- every snapshot always has state, so `has_state` is not stored;
- `has_sequences=true` implies `has_chains=true`;
- `snapshot_reason_id` records the actual cause, while `is_final` independently marks finality.

Current `snapshot_reasons` names are:

```text
0 initial
1 scheduled
2 action
3 manual
4 final
```

The exact active dictionary is declared in `schema.jsonl`.

## 8. `state/`

Columns:

```text
snapshot_id.npy       uint64
entity_id.npy         uint32
count.npy             uint64
moles.npy             float64
concentration.npy     float64, mol/L
```

`state/` is dense over `snapshot_id × entity_id`.

Rules:

- each snapshot contains every state entity, including zero-valued entities;
- row order: `snapshot_id` ascending, then `entity_id` ascending;
- within each snapshot, entity IDs are exactly `0..n_entities-1`;
- total rows = `n_snapshots × n_state_entities`;
- `count` is the actual number of physical objects represented in `V_KMC`, not a scaled pseudo-particle count;
- `moles = count / N_A`;
- `concentration = moles / snapshots.kmc_volume_L` for the matching snapshot;
- all three values are stored;
- `count=0` requires exact `moles=0.0` and `concentration=0.0`;
- all values are nonnegative and finite;
- validation uses `rtol=1e-12`, `atol=0` for the two conversion identities.

`state_entities` covers monomers, ordinary species, and useful chain aggregates
such as `live_chains` and `dead_chains`. Homo additionally exposes declared
kinetic pools as state entities.

Each dictionary entry has at least:

```text
id
name
kind
```

and may have:

```text
label
description
```

Machine names use `snake_case`; chemical symbolic capitals are allowed, e.g. `monomer_A`, `terminal_B`.

`state/` does not contain derived analysis such as conversion, rates, fractions, Mn, Mw, or Mz.

## 9. `chains/`

Definition: `chains` is a table of compressed chain records.

One row represents `count` identical chains.

The shared columns are:

```text
chain_record_id.npy    uint64
snapshot_id.npy        uint64
population_id.npy      uint32
pool_id.npy            uint32
origin_id.npy          uint32
dp.npy                 uint64
molar_mass.npy         float64, g/mol
count.npy              uint64
moles.npy              float64
concentration.npy      float64, mol/L
left_end_id.npy        uint32
right_end_id.npy             uint32
```

Copo additionally writes terminal summaries and explicit sequence presence:

```text
has_first_monomer.npy        bool
first_monomer_id.npy         uint32
has_penultimate_monomer.npy  bool
penultimate_monomer_id.npy   uint32
has_last_monomer.npy         bool
last_monomer_id.npy          uint32
has_sequence.npy             bool
sequence_offset.npy          uint64
sequence_length.npy          uint64
```

Rules:

- `chain_record_id = row_index`, continuous from zero, unique within the run;
- it identifies a stored row, not a persistent chain tracked across snapshots;
- the same logical chain type in two snapshots receives two different IDs;
- `dp >= 1`;
- `count >= 1`;
- `molar_mass > 0`, finite;
- `moles = count / N_A`;
- `concentration = moles / snapshots.kmc_volume_L` for the matching snapshot;
- homo merges identical records within a snapshot and represents multiplicity
  with `count`;
- copo currently writes one physical chain per row (`count=1`), preserving a
  lossless but not compressed representation;
- records of one snapshot form one contiguous block;
- primary order: `snapshot_id`, `population_id`, `pool_id`, `dp`;
- further order is stable and deterministic, but lexicographic sequence sorting is not required;
- `population_id` is required; initial dictionary: `0=live`, `1=dead`;
- `pool_id` is optional; `0=not_applicable`, real pools start at `1`;
- `origin_id` is required; `0=unknown`, real origins start at `1`;
- `origin_id` records the provenance mechanism that created the current chain record, not the last reaction that modified it;
- `left` and `right` are representational conventions, not spatial orientation;
- `left_end_id` and `right_end_id` use the same `chain_end_types` dictionary;
- first, penultimate, and last monomer summaries remain available even when literal sequence storage is absent;
- optional monomer summaries use explicit `has_*` flags because monomer ID zero is valid;
- chain end codes begin with `0=not_applicable`, `1=unknown`;
- both ends may have the same type.

Suggested origin names include `init`, `transfer_m`, `term_c`, `term_d_H`, `term_d_U`, `term_x`, and `transfer_h`, as supported by the engine. Here `transfer_h` is a result/storage origin label; the only valid model keyword for the corresponding reaction is `macro transfer`.

## 10. `chain_composition/`

`chain_composition/` is a copo table that preserves complete per-monomer
composition without model-specific `n_A`, `n_B`, ... columns. Homo composition
is completely determined by `dp` and its single monomer and therefore has no
separate table.

Columns:

```text
chain_record_id.npy    uint64
monomer_id.npy         uint32
unit_count.npy         uint64
```

Rules:

- rows are dense by `chain_record_id × monomer_id`; zero counts are stored;
- row order is chain record ID, then monomer ID;
- for every chain record, `sum(unit_count) = dp`;
- composition participates in compressed-record identity;
- the table provides per-monomer unit counts and repeat-unit molar mass;

## 11. `sequences/`

Column:

```text
symbols.npy    uint32
```

All stored linear sequences share one symbol array. A chain record references:

```text
symbols[sequence_offset : sequence_offset + sequence_length]
```

Rules:

- symbol IDs refer to `sequence_symbols` in `schema.jsonl`;
- symbol IDs are continuous from zero and `0` may be a real symbol;
- no symbol value encodes absence;
- absent sequence: `has_sequence=false`, `sequence_offset=0`, `sequence_length=0`;
- a present sequence has `has_sequence=true` and is complete: `sequence_length=dp`;
- partial sequence (`0 < sequence_length < dp`) is forbidden in v1;
- `snapshots.has_sequences=true` means at least one chain record in that snapshot has a complete stored sequence;
- mixed sequence coverage is permitted because copo may retain all live sequences but only short dead-chain sequences; per-record `has_sequence` is authoritative;
- `snapshots.has_sequences=false` means every corresponding record has `has_sequence=false` and zero offset/length;
- symbols are ordered from `left_end_id` to `right_end_id`;
- chain-end symbols are not included in `symbols.npy`;
- identical full ranges may be shared by multiple records, including across snapshots;
- ranges may be identical or disjoint, but partial overlap is forbidden;
- deduplication is allowed but not required;
- an empty valid `symbols.npy` exists when no sequences are stored;
- v1 supports linear sequences only; general branching, cycles, and networks require a future graph format.

Each `sequence_symbols` entry has:

```text
id
name
kind
label?
description?
monomer_id?
molar_mass_increment?
```

`molar_mass_increment` uses g/mol and may be omitted when unknown. `chain_end_types` may include optional `molar_mass_contribution` and `has_known_molar_mass_contribution`. Chain molar-mass reconstruction is validated only when the full sequence is stored and every symbol and both ends have known mass contributions. Missing end-group masses are permitted and cause that validation check to be skipped, not failed.

## 12. `moments/`

One row represents:

```text
snapshot_id × population_scope_id × mass_basis_id
```

Columns:

```text
snapshot_id.npy                 uint64
population_scope_id.npy         uint32
mass_basis_id.npy               uint32
chain_count.npy                 uint64
sum_dp.npy                      float64
sum_dp2.npy                     float64
dp_n.npy                        float64
dp_w.npy                        float64
sum_molar_mass.npy              float64
sum_molar_mass2.npy             float64
sum_molar_mass3.npy             float64
mn.npy                          float64
mw.npy                          float64
mz.npy                          float64
dispersity.npy                  float64
```

Initial `population_scope` dictionary:

```text
0 all
1 live
2 dead
```

The dictionary selects the chain population included in an aggregate. Future activity states such as active or dormant should be represented by a separate dimension rather than overloaded into `population_scope`.

Initial `mass_bases` dictionary:

```text
0 repeat_units
1 with_end_groups
```

Rules:

- row order: snapshot, population scope, mass basis;
- `chain_count = Σ count_i`, not number of compressed records;
- `sum_dp = Σ(count_i × dp_i)`;
- `sum_dp2 = Σ(count_i × dp_i²)`;
- `sum_molar_mass = Σ(count_i × molar_mass_i)`;
- second and third mass sums use powers 2 and 3;
- `moments/` exists only for snapshots with `has_chains=true`;
- for an empty selected population, raw sums are exactly zero and derived averages (`dp_n`, `dp_w`, `mn`, `mw`, `mz`, `dispersity`) are `NaN`;
- `NaN` is allowed only for mathematically undefined derived values;
- `Inf` is invalid;
- `moments/` is derived; `chains/` remains the source of truth;
- only canonical mass bases are stored.

MWD means molecular-weight distribution. CLD means chain-length distribution.

## 13. `channel_events/`

One row represents one channel at one snapshot.

Columns:

```text
snapshot_id.npy                    uint64
channel_id.npy                     uint32
event_count.npy                    uint64
productive_event_count.npy         uint64
nonproductive_event_count.npy      uint64
```

Semantics:

- all three counts are cumulative from the beginning of the run;
- for every snapshot, `snapshots.kmc_event` equals the sum of `event_count` over all channels at that snapshot;
- each selected KMC/SSA channel increments `kmc_event` by one, including nonproductive events; actions do not increment `kmc_event`;
- `event_count = productive_event_count + nonproductive_event_count`;
- for channels without an efficiency/nonproductive mechanism, `productive_event_count = event_count` and `nonproductive_event_count = 0`;
- a productive event consumes substrates and creates the defined products;
- a nonproductive event consumes substrates but does not create the products;
- counts are nondecreasing across snapshots;
- `total_fires` and `share` are derived and therefore not stored.

## 13a. `channel_propensities/`

In the current format this is a copo table. One row represents one SSA channel
at one snapshot. It preserves the
instantaneous event competition needed for propensity/rate shares; it is not a
copy of realized firing counts.

Columns:

```text
snapshot_id.npy        uint64
channel_id.npy         uint32
propensity.npy         float64, s^-1
total_propensity.npy   float64, s^-1
```

Rules:

- rows are dense by `snapshot_id × channel_id`;
- row order is snapshot, then channel ID;
- `propensity` is the exact channel propensity evaluated for the saved state;
- `total_propensity` equals the sum over all channels at that snapshot and is
  repeated on each row for direct validation;
- values are finite and nonnegative;
- when total propensity is zero, propensity shares are mathematically undefined
  and readers return `NaN`, not zero;
- firing shares remain derived from `channel_events/` and must not be presented
  as propensity shares.

## 13b. `microstructure_motifs/`

This table stores engine-counted sequence motifs so dyads and triads remain
available in `sequence_mode=composition`, where literal chain order is not
written. One row represents one motif at one snapshot.

Columns:

```text
snapshot_id.npy   uint64
motif_order.npy   uint32       # 2 = dyad, 3 = triad
motif_id.npy      uint32
count.npy         uint64
```

Motif IDs refer to the order-specific dictionaries
`microstructure_dyads` and `microstructure_triads`. Counts are abundance-weighted
whole-population counts produced by the engine and include live and dead chains.
Zero counts are stored, making the table dense over declared motifs.

## 13c. `block_statistics/`

This table preserves the exact block-length histogram accumulated by the copo
engine, independently of literal sequence storage.

Columns:

```text
snapshot_id.npy   uint64
monomer_id.npy    uint32
block_length.npy  uint64
block_count.npy   uint64
```

Rules:

- rows are ordered by snapshot, monomer, then block length;
- `block_count` is abundance weighted;
- only positive block lengths are valid;
- the table is the source for `run_lengths()` and block-derived diagnostics;
- in `sequence_mode=full`, readers may independently recompute motifs and blocks
  from sequences and compare both sources.

## 14. `actions/`

`actions/` stores executed-action history. Action definitions remain authoritative in `input.model`.

Columns:

```text
action_id.npy                       uint64
kmc_event.npy                       uint64
time.npy                            float64
source_line.npy                     uint32
action_type_id.npy                  uint32
trigger_type_id.npy                 uint32
scheduled_time.npy                  float64
target_id.npy                       uint32
requested_value.npy                 float64
before_value.npy                    float64
after_value.npy                     float64
state_changed.npy                   bool
output_written.npy                  bool
has_snapshot.npy                    bool
snapshot_id.npy                     uint64
has_kinetic_parameter_set.npy       bool
kinetic_parameter_set_id.npy        uint64
```

## 14a. `action_conditions/`

One row represents one atomic condition belonging to an executed `when` action.

Columns:

```text
condition_record_id.npy              uint64
action_id.npy                        uint64
condition_index.npy                  uint32
observable_id.npy                    uint32
operator_id.npy                      uint32
threshold.npy                        float64
observed_value.npy                   float64
condition_met.npy                    bool
```

Rules:

- `condition_record_id = row_index`, continuous from zero;
- rows are ordered by `action_id`, then `condition_index`;
- `condition_index` starts at zero for each action;
- all rows with one `action_id` form an `AND` conjunction from one physical `when` line;
- separate `when` lines are separate actions and are evaluated independently;
- `at` and `every` actions have no rows in this table;
- because `actions/` stores only executed actions, every stored condition has `condition_met=true`;
- `observed_value` is the value evaluated on the common state that triggered the action;
- concrete species concentrations have distinct observable dictionary entries;
- general `OR`, parentheses, and expression trees are outside v1.

Text messages are stored in `actions/messages.jsonl`, with one JSONL record per
action that has a message, for example:

```json
{"action_id":"7","message":"Temperature changed."}
```

For `actions/` and its messages:

- `action_id = row_index`, continuous from zero;
- `source_line` is the physical 1-based line in `input.model`;
- triggers are `at`, `every`, and `when`;
- `scheduled_time` is `NaN` for `when` actions;
- absent numeric values use `NaN` and absent foreign keys use `has_*=false`
  with ID zero;
- actions without messages have no JSONL record.

## 14b. `feed_events/`

One row represents one executed feed action:

```text
action_id.npy                    uint64
feed_id.npy                      uint32
dose_mL.npy                      float64, mL
volume_before_mL.npy             float64, mL
volume_after_mL.npy              float64, mL
kmc_volume_before_L.npy          float64, L
kmc_volume_after_L.npy           float64, L
```

`feed_id` refers to the `feeds` dictionary. The row is linked to the executed
action, and the physical and KMC volume changes must have the same scale
factor. Cumulative feed volume and component amounts are derived by pyslimmc
from this table and the structured feed definition in metadata.

## 14c. `monomer_balance/` and `species_balance/`

`monomer_balance/` is dense over snapshot × monomer:

```text
snapshot_id.npy          uint64
monomer_id.npy           uint32
initial_moles.npy        float64, mol
introduced_moles.npy     float64, mol
free_moles.npy           float64, mol
incorporated_moles.npy   float64, mol
conversion.npy           float64
```

`species_balance/` is dense over snapshot × balance entity:

```text
snapshot_id.npy       uint64
entity_id.npy         uint32
initial_moles.npy     float64, mol
dosed_moles.npy       float64, mol
total_moles.npy       float64, mol
free_moles.npy        float64, mol
consumed_moles.npy    float64, mol
```

The balance tables describe the physical reactor and require an initial
physical volume. A `set_c` action is a technical override rather than a
physical addition/removal; readers mark the affected balance inapplicable from
that action onward rather than pretending it closes.

## 15. `kinetic_parameters/`

The directory stores complete kinetic parameter sets active during the run.

Structure:

```text
kinetic_parameters/
├── sets/
│   ├── kinetic_parameter_set_id.npy
│   ├── start_kmc_event.npy
│   ├── start_time.npy
│   ├── has_source_action.npy
│   └── source_action_id.npy
└── values/
    ├── kinetic_parameter_set_id.npy
    ├── kinetic_parameter_id.npy
    └── value.npy
```

Dtypes:

```text
set identifiers and event/action IDs  uint64
kinetic_parameter_id                  uint32
start_time and value                  float64
has_source_action                     bool
```

Rules:

- set IDs are continuous from zero;
- set zero is the initial set active from event zero and time zero;
- each snapshot points to one complete set;
- each set contains exactly one value for every model-specific entry in `kinetic_parameter_definitions`;
- values are dense by `set_id × parameter_id`;
- model-specific definitions include only parameters that exist for that model: a constant rate contributes `k`, while an Arrhenius rate contributes current `k(T)`, `A`, and `Ea`;
- parameter-set values are finite and do not use `NaN`;
- row order: set ID, then parameter ID;
- after an action changes one value, a full new set is stored;
- temperature is a kinetic parameter with `kind=temperature`;
- current `k(T)`, Arrhenius `A`, and `Ea` are ordinary parameter records with kinds `rate_constant`, `arrhenius_A`, and `arrhenius_Ea`;
- no `effective_value` column and no `k_*_effective` naming convention are used.

Definitions are stored in the `kinetic_parameter_definitions` dictionary with:

```text
id
name
kind
unit
label?
description?
```

Initial kinds include `temperature`, `rate_constant`, `efficiency`, `arrhenius_A`, and `arrhenius_Ea`.

## 16. `diagnostics/`

`diagnostics/` is a technical container, not a source-of-truth results table.

```text
diagnostics/
├── validation.jsonl
├── channel_trace/
├── run.log
└── debug.log
```

Rules:

- `run.log` is required;
- `debug.log` is optional;
- empty `validation.jsonl` means validation was not run;
- each validation line is one check object;
- validation summary is duplicated in `run_metadata.json`;
- absence or incompleteness of debug trace does not invalidate scientific results.

## 17. `diagnostics/channel_trace/`

Columns:

```text
kmc_event.npy          uint64
time.npy               float64, s
dt.npy                 float64, s
channel_id.npy         uint32
rate.npy               float64, engine-native
propensity.npy         float64, s^-1
total_propensity.npy   float64, s^-1
```

Rules:

- detailed event trace is optional and diagnostic only;
- directory and declared empty columns exist even when trace is disabled;
- first recorded event normally has `kmc_event=1`;
- `time` is the time after the event;
- `dt` is the sampled waiting time;
- `rate`, selected-channel `propensity`, and `total_propensity` record the
  values used for that event;
- full state is not duplicated in the trace;
- RNG draws are not stored in the basic trace;
- partial trace by event range or channel is allowed;
- metadata records whether trace is enabled and complete.

## 18. Chain DP, sequence mode, and stored memory metadata

The concrete run configuration is preserved in `run_metadata.json`. At minimum metadata records the following fields. The `memory_limit_*` names below are **Storage metadata keys**, not model-language commands; the user-facing model statement is [`at_memory`](../MODEL_SYNTAX.md#memory-policy).

```text
t_end_s
max_steps
when_check_events
dp_max
sequence_mode
memory_limit_enabled
memory_limit_bytes?
memory_snapshot_on_limit
memory_stop_on_limit
engine_chain_dp_dtype
engine_chain_dp_max
max_dp_observed
termination_reason?
```

Rules:

- `dp_max` is the single user-facing maximum chain DP for both homo and copo;
- `dp_max > 0`; propagation that would create `dp > dp_max` is not eligible;
- reaching `dp_max` does not by itself kill or remove a chain: other applicable channels may still act on it;
- `sequence_mode` is exactly `composition` or `full`;
- `composition` stores DP, complete monomer composition, first/penultimate/last monomer, ends, mass, pool, population, and origin, but not the full linear order;
- `full` stores the same information plus the complete linear sequence for every stored chain record;
- `dead_sequence_max_dp`, `sequence_storage_max_dp`, and `oligomer_max_dp` are not part of Slimmc Storage or the new canonical model language;
- oligomer classification is reader-side analysis with an explicit query cutoff, not an engine output mode;
- `engine_chain_dp_max` is the technical bound of the in-memory DP representation and is distinct from `dp_max`;
- `max_dp_observed` is an observed result, not configuration;
- memory-policy configuration is preserved independently of whether the threshold was reached;
- readers should treat `memory_limit_*` as historical/internal Storage field names, not as syntax to put in a `.model` file;
- `termination_reason` distinguishes `t_end`, `max_steps`, `memory_limit`, user stop, interruption, validation failure, and runtime failure where applicable.

The optional `memory/` table preserves implementation-dependent memory estimates; it is diagnostic rather than scientific state.

## 19. Source-of-truth hierarchy

- `input.model`: source of model definitions and action definitions;
- `snapshots/`: source of chronology, current physical/KMC volume, chain
  population counts, and available content;
- `state/`: source of engine counters and their direct unit conversions;
- `chains/`: source of chain counts, DP, molecular weights, ends, terminal summaries, populations, pools, and origins;
- `chain_composition/`: source of complete per-chain monomer composition;
- `sequences/`: source of stored linear sequence data;
- `moments/`, CLD, and MWD: derived aggregates;
- `channel_events/`: source of cumulative per-channel event statistics;
- `channel_propensities/`: source of instantaneous per-channel propensities at snapshots;
- `microstructure_motifs/`: source of engine-counted dyads and triads;
- `block_statistics/`: source of engine-counted block-length histograms;
- `actions/` and `action_conditions/`: source of executed-action history and evaluated `when` atoms;
- `feed_events/`: source of executed feed doses and volume changes;
- `monomer_balance/` and `species_balance/`: source of physical amount and
  conversion ledgers;
- `kinetic_parameters/`: source of parameter sets active during the run;
- `memory/`: source of engine-specific memory estimates, but not a scientific state table;
- `diagnostics/`: non-authoritative validation and debug material.

## 20. Atomicity, interruption, and finalization

Writers must prevent a partially finalized run from appearing complete.

Run status rules:

- `running`: writer has initialized the run and no final snapshot is required;
- `completed`: all required canonical files are finalized and exactly one snapshot has `is_final=true`;
- `failed`: execution ended because of an engine or writer error; no final snapshot is required;
- `interrupted`: execution was stopped externally or intentionally before normal completion; no final snapshot is required.

Working-data model:

```text
results/run_000001/.work/
```

- while the simulation runs, append-only working columns are stored under `.work/`;
- canonical `.npy` columns are produced from complete working records during finalization;
- each canonical column is written to a temporary path, flushed, closed, and atomically renamed;
- finalization may be attempted after failure or interruption to publish all complete snapshots already present in `.work/`;
- the last available snapshot of a failed or interrupted run remains an ordinary snapshot and is never relabeled as final.

Completion marker:

```text
RESULTS_COMPLETE
```

- the marker is created atomically only after every required canonical file, hash, and completed-run metadata field has been finalized;
- `RESULTS_COMPLETE` may exist only when `run_status=completed`;
- absence of the marker means the run is running, failed, interrupted, or not fully finalized;
- readers may open such results only in an explicit partial/recovery mode;
- `.work/` is implementation data and is not part of the scientific source-of-truth hierarchy.

Additional requirements:

- `run_metadata.json` is created at initialization with `run_status=running` and rewritten atomically on status changes;
- JSON and JSONL files are written atomically;
- `run_status` is not set to `completed` before all required columns and hashes are finalized;
- schema is frozen before the first working-data append;
- all columns of one table must have equal row count;
- a failed or interrupted run may retain valid canonical `.npy` files for complete snapshots, but never receives `RESULTS_COMPLETE`.

## 21. Validation summary

A conforming validator checks at least:

- JSONL syntax and required schema records;
- required directory and column presence;
- exact dtype and little-endian profile;
- equal column lengths per table;
- continuous record IDs where required;
- foreign keys and dictionary bounds;
- dense and ordered `state/`;
- `moles` and concentration conversions;
- snapshot physical/KMC volumes, feed-event scaling, and balance identities;
- finite/nonnegative constraints;
- snapshot ordering and finality;
- chain merging uniqueness;
- full-or-absent sequence rule;
- nonoverlapping sequence storage ranges except identical sharing;
- moment reconstruction from chains;
- cumulative channel-event identities, including `kmc_event = Σ event_count`;
- complete, dense, finite kinetic parameter sets;
- action optional-key masks;
- optional chain-mass reconstruction only when all mass contributions are known;
- consistency of `run_status`, `is_final`, and `RESULTS_COMPLETE`;
- hashes of `input.model` and `schema.jsonl`.


## 22. Minimal independent reader

The format is intentionally readable without pyslimmc. For example, the
following code reads run metadata and the canonical snapshot axis directly:

```python
from pathlib import Path
import json
import numpy as np

run = Path("results/main")
metadata = json.loads((run / "run_metadata.json").read_text())
snapshot_id = np.load(run / "snapshots" / "snapshot_id.npy", allow_pickle=False)
time_s = np.load(run / "snapshots" / "time.npy", allow_pickle=False)
kmc_event = np.load(run / "snapshots" / "kmc_event.npy", allow_pickle=False)

assert len(snapshot_id) == len(time_s) == len(kmc_event)
print(metadata["run_status"], time_s[-1], kmc_event[-1])
```

An independent reader should use `schema.jsonl` rather than assuming optional
tables exist, keep `allow_pickle=False`, preserve row alignment within each
table, and honor the completion/status rules in section 20. pyslimmc remains
the supported high-level analysis API; this example demonstrates that the
Storage specification is independently consumable rather than Python-object
serialization.

## 23. Deferred beyond v1

- general graph representation for branching, cycles, and networks;
- partial stored sequences;
- a global schema shared by all models;
- mandatory sequence deduplication;
- persisted CLD/MWD as authoritative data;
- restart/checkpoint format.

### Implementation note: homo validator

The homo writer runs named validation checks after writing canonical columns and before publishing `RESULTS_COMPLETE`. Detailed results are written to `diagnostics/validation.jsonl`; summary status and counts are copied to `run_metadata.json`. A requested completed run with any validation error is finalized as failed and receives no completion marker.


## 24. Canonical output

Slimmc writes only Slimmc Storage v1. `pyslimmc` reads only this format. User-requested analysis exports such as `to_tsv()` are external products and are not part of the run backend.

## See also

- [`../SIMULATION_RESULTS.md`](../SIMULATION_RESULTS.md) — Meaning of simulation results
- [`../PYSLIMMC.md`](../PYSLIMMC.md) — pyslimmc guide
- [`../development/ARCHITECTURE.md`](../development/ARCHITECTURE.md) — Architecture
