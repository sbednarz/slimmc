# Model syntax

A slimmc `.model` file is a text description of a stochastic kinetic model and
its process/storage actions. This guide organizes the syntax by what a chemist
usually wants to declare. The exact engine contracts remain canonical in
[`reference/HOMO.md`](reference/HOMO.md) and
[`reference/COPO.md`](reference/COPO.md).

## Command groups

| Group | Main statements | Purpose |
|---|---|---|
| Metadata | `desc`, `var` | description and sweep metadata |
| Global parameters | `param` | stochastic scale, reactor/process settings, output |
| Components | `monomer`, `species`, `endgroup`, `feed` | chemical components and masses |
| Chain pools | `polymer` | active/dead polymer populations |
| Kinetics | `rate`, `rxn`, `macro ...` | rate constants and reaction channels |
| Scheduled actions | `at`, `every`, `from ... repeat ... every ...` | time-based process/output changes |
| Conditional actions | `when ... [and ...] ...` | one-shot state/conversion conditions |
| Memory policy | `at_memory` | save and/or stop when estimated memory reaches a threshold |

Keywords are case-sensitive. Names use `[A-Za-z_][A-Za-z0-9_]*`; names that
differ only by ASCII letter case are rejected. Comments begin with `#`.

## Global parameters and defaults

The common parameter contract is:

| Parameter | Required | Default | Meaning / unit |
|---|---:|---:|---|
| `kmc_volume` | yes | none | stochastic simulation volume, L |
| `init_volume` | with `feed` | absent | initial physical reactor volume; L by default, `mL/ml/ML` accepted |
| `temperature` | no | `298.15` | K |
| `t_end` | yes | none | simulation end time, s |
| `max_steps` | no | `10000000000` | maximum SSA events |
| `when_check_events` | no | `1` | event cadence for conditional checks |
| `seed` | no | `12345` | signed 64-bit RNG seed |
| `output_dir` | no | `results/<model_stem>/` | Storage output path |
| `dp_max` | no | `2147483647` | maximum chain DP |
| `sequence_mode` | no | `composition` | `composition` or `full` |
| `mass_model` | no | `repeat_units` | `repeat_units` or `with_end_groups` |

Example:

```text
param kmc_volume 1.0e-15
param init_volume 100 mL
param temperature 353.15
param t_end 7200
param seed 12345
param output_dir "results/run_01"
```

`kmc_volume` controls the discrete KMC population. `init_volume`, when present,
represents the physical reactor volume and enables physical moles, feed-volume
tracking and reactor-volume analyses. They are different quantities.

Relative `output_dir` paths are resolved from the model directory. Parent
segments (`..`), empty segments, extensions, spaces and glob characters are
rejected by the model-path contract.

## Chemical declarations

### Monomers

```text
monomer Sty 1.0 104.15
monomer MMA 0.50 100.12
```

The fields are `name`, initial concentration in mol/L, and repeat-unit molar
mass in g/mol. Homo accepts zero or one monomer; copo accepts two or three.

### Small species

```text
species AIBN 0.010
species R 0.0
species CTA 0.020
```

Initial values are mol/L.

### End groups and chain mass

```text
endgroup R 68.0
endgroup H 1.008
param mass_model with_end_groups
```

End-group masses are in g/mol and contribute only with
`mass_model with_end_groups`. With the default `repeat_units`, chain molar mass
uses repeat-unit masses only.

### Polymer pools

```text
polymer P active
polymer D dead
```

Homo normally uses active and dead chain pools. Copo uses terminal or explicit
penultimate active pools plus exactly one dead pool; see the exact copo
reference for topology rules.

### Feed mixtures

```text
param init_volume 50 mL
feed feed1 Sty 2.0
feed feed1 CTA 0.02
```

Repeating a feed name builds one constant-composition mixture. Feed component
concentrations are mol/L. A feed declaration requires `init_volume`. The feed
mixture is immutable; define another feed name for another composition.

A feed dose is an explicit volume. Without a suffix it is in L; `mL`, `ml`,
and `ML` are accepted millilitre forms:

```text
at 400 feed feed1 0.001
at 400 feed feed1 1 mL
```

These two actions are equivalent. A dose increases physical reactor volume and
scales the represented KMC volume consistently; existing polymer chains are not
removed. Dilution follows from the increased volume.

## Rate constants

Three forms are accepted:

```text
rate kp 500.0
rate kt const 1.0e7
rate kd arr 2.89e15 130230
```

`rate NAME VALUE` and `rate NAME const VALUE` are fixed rates. Arrhenius rates
use pre-exponential factor `A` and activation energy `Ea` in J/mol:

$$
k(T)=A\exp\left(-\frac{E_a}{RT}\right).
$$

The dimensional unit of a rate constant follows its reaction molecularity.
The exact propensity convention, including the factor-of-two convention for
same-population bimolecular channels, is defined in [`THEORY.md`](THEORY.md).

### Where kinetic parameters come from

slimmc does not provide a hidden kinetic-parameter database. Use
conditions-appropriate primary literature or critically evaluated datasets and
check the defining rate equation, units, solvent/composition range, pressure,
and temperature before transferring a number into a model. Useful source types
include:

- IUPAC/PLP critically evaluated propagation-rate datasets where available;
- manufacturer kinetic/half-life data for commercial radical initiators,
  checked against primary literature when the application requires it;
- primary kinetic studies for termination, transfer, depropagation, and other
  mechanism-specific channels;
- measured reactivity ratios or elementary rate coefficients for
  copolymerization.

A handbook or database can be a useful index, but the defining equation and
experimental conditions remain part of the parameter. The quick-start model
shows explicit provenance for its styrene `kp` and AIBN `kd` values.

## Reaction channels

### Elementary reactions

`rxn` represents small-species channels. Examples:

```text
rxn I -> R + R kd
rxn A + B -> C k
rxn A + A -> B k
```

An optional final efficiency factor is supported where documented, for example
initiator decomposition:

```text
rxn AIBN -> R + R kd 0.6
```

### Polymer reactions

Shared macro channel families include:

| Statement | Chemical role |
|---|---|
| `macro init ...` | radical addition / chain initiation |
| `macro prop ...` | propagation |
| `macro deprop ...` | depropagation where supported |
| `macro term_c ...` | termination by combination |
| `macro term_d ...` | termination by disproportionation |
| `macro transfer ...` | transfer to a small species and radical product |
| `macro transfer_m ...` | transfer to monomer |
| `macro term_x ...` | capping / chain-ending reaction with a small species |

The exact left/right pools and products are engine-specific. Copy exact forms
from [`reference/HOMO.md`](reference/HOMO.md) or
[`reference/COPO.md`](reference/COPO.md), especially for terminal and
penultimate copolymer propagation.

## Process and storage actions

### Scheduled actions

```text
at 0 save
at 3600 set_temp 373.15
every 300 save
from 0 repeat 12 every 300 feed feed1 0.50 mL
```

`at` fires at one simulated time. `every` repeats on a time cadence. The
bounded `from ... repeat ... every ...` form is useful for semibatch dosing and
other finite schedules. In

```text
from START repeat COUNT every PERIOD ACTION
```

the action fires at `START + n*PERIOD` for `n = 0..COUNT-1`; occurrences after
`t_end` are not executed. `repeat 1` is valid.

For a simple constant feed with no finite count, the engine-specific references
also document the supported repeating forms. Complex linear/nonlinear dosing
profiles are represented by explicit dose blocks rather than a hidden profile
expression evaluator.

Available action kinds are documented exhaustively in the engine references;
common actions include:

```text
save
save_chains
stop
print
print_info
set_c
add_c
set_rate
add_rate
set_temp
feed
```

`save` writes a state snapshot. `save_chains` retains the chain population
needed for chain-resolved distributions and microstructure at that snapshot.
See [`SIMULATION_RESULTS.md`](SIMULATION_RESULTS.md) for the consequence of this
distinction.

`add_c` is a physical inventory adjustment at constant volume and may be
positive or negative as long as the resulting population remains nonnegative.
`set_c` is a technical concentration override, not a balance-preserving physical
operation; the engine warns and the physical balance for that species becomes
inapplicable.

### Conditional actions

Conditions are one-shot and use `<` or `>`:

```text
when X Sty > 0.80 save
when X Sty > 0.80 save_chains
when X Sty > 0.80 stop
```

Multiple conditions can be combined with `and`:

```text
when X A > 0.90 and X B > 0.90 stop
```

Homo and copo use the same explicit monomer form `X MONOMER`.

## From laboratory description to model statement

| Laboratory/modeling idea | Typical slimmc declaration |
|---|---|
| styrene at 1 mol/L | `monomer Sty 1.0 104.15` |
| AIBN at 10 mmol/L | `species AIBN 0.010` |
| reactor at 80 °C | `param temperature 353.15` |
| simulate for 2 h | `param t_end 7200` |
| Arrhenius initiator decay | `rate kd arr A Ea` + `rxn ...` |
| propagation | `macro prop ... kp` |
| combination | `macro term_c ... ktc` |
| disproportionation | `macro term_d ... ktd` |
| chain-transfer agent | `macro transfer ... ktr` |
| dose a feed periodically | `feed ...` + `from ... repeat ... every ... feed ...` |
| save concentrations/summary state | `save` |
| retain chains for MWD/CLD | `save_chains` |
| stop after conversion threshold | `when X MONOMER > ... stop` |

This table maps concepts to syntax only; it does not provide kinetic constants.
Use kinetic data from appropriate primary literature or validated data sources
and check the defining rate equation before transferring a numeric constant.

## Sweep metadata

`var` identifies declarations varied across a model series for later scanning:

```text
var rate kp L_mol_s
var param temperature K
var species CTA mol_L
var monomer Sty mol_L
```

It is metadata: the actual value remains in the normal `rate`, `param`,
`species`, `monomer`, or `endgroup` declaration.

## Memory policy

`at_memory SIZE save|stop` defines what Slimmc should do when its estimated
memory use reaches `SIZE`. Byte units `B`, `KB`, `MB`, and `GB` are accepted.
The actions may be combined:

```text
at_memory 3GB save stop
```

`save` writes one full logical snapshot at the threshold, including the chain
population (equivalent in data content to a chain-bearing snapshot rather than
a normal state-only `save`). `stop` then terminates the run cleanly with a
memory-policy termination reason. With `save` alone, the threshold-triggered
snapshot is written only once and the simulation continues.

## What is intentionally not repeated here

This page is the shared user-facing syntax map. The formal references contain
all parser tokens, cardinality rules, engine-specific channel topology,
validation constraints, engine-specific actions and omitted transitions:

- [`reference/HOMO.md`](reference/HOMO.md)
- [`reference/COPO.md`](reference/COPO.md)

## See also

- [`reference/HOMO.md`](reference/HOMO.md) — Homo engine reference
- [`reference/COPO.md`](reference/COPO.md) — Copo engine reference
- [`QUICKSTART.md`](QUICKSTART.md) — Quick start
