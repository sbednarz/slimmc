# slimmc homo engine



This is the user and language reference for the single-monomer `slimmc`
engine. It describes the current implementation, not older model dialects.

## Purpose and scope

The homo engine performs stochastic simulation of homogeneous free-radical
polymerization using a compact chain representation. A chain stores its degree
of polymerization, kinetic pool, end groups and formation channel, but not a
per-mer sequence. This makes it the preferred engine for large
single-monomer simulations.

Implemented mechanisms include elementary reactions, initiation,
propagation, combination, disproportionation, capping, transfer,
transfer-to-monomer, terminal depropagation, scheduled/conditional actions,
fixed and Arrhenius rates, and repeat-unit or end-group-aware masses.

Batch and semibatch operation are supported. In semibatch mode the engine
updates both the physical reactor volume and the KMC representation volume
after every feed action. Not represented: copolymer sequences,
branching/crosslink topology, spatial effects and diffusion.

## Running the engine

From the family root, the normal interface is the dispatcher:

```bash
make build
bin/slimmc path/to/model.model
```

The component binary can also be built and used directly:

```bash
make -C homo build
homo/slimmc path/to/example.model
slimmc --check path/to/example.model
```

`--check` parses and validates a model without simulating it. See
`homo/slimmc --help` for diagnostic and trace options supported by the built
binary.

## Model file

Comments begin with `#`. Keywords are case-sensitive. Values use ordinary
decimal or scientific notation. Paths or text containing whitespace or `#`
must be quoted.

Declaration names follow `[A-Za-z_][A-Za-z0-9_]*`. Keywords are reserved, and
two names that differ only by ASCII letter case are rejected to keep stored
dictionaries unambiguous. An `output_dir` may be relative or absolute, but
each non-root path segment follows the same identifier grammar; `.`, `..`,
spaces, extensions, and glob characters are not accepted in model paths.

The dispatcher also routes zero-monomer pure-kinetics models to homo. Such
models may use elementary `rxn` channels but cannot use polymer mechanisms.

### Metadata and parameters

```text
desc "Free-radical polymerization"

param kmc_volume 1.0e-15
param init_volume 0.100 L
param temperature 333.15
param t_end 10.0
param max_steps 1000000
param when_check_events 1000
param seed 12345
param output_dir "results/FRP01"
param dp_max 1000000
param sequence_mode composition
param mass_model repeat_units
```

The parameter contract is:

| Parameter | Required | Default | Accepted value and unit | Eligible for `var param` |
|---|---:|---:|---|---:|
| `kmc_volume` | yes | none | finite float `>0`, L | yes |
| `init_volume` | with `feed` | absent | finite float `>0`; no suffix or `L/l` means L, `mL/ml/ML` means mL | yes |
| `temperature` | no | `298.15` | finite float `>0`, K | yes |
| `t_end` | yes | none | finite float `>=0`, s | yes |
| `max_steps` | no | `10000000000` | integer `>0` | yes |
| `when_check_events` | no | `1` | integer `>0` | yes |
| `seed` | no | `12345` | signed 64-bit integer | yes |
| `output_dir` | no | `results/<model_stem>/` | quoted validated path | no |
| `dp_max` | no | `2147483647` | integer `>0` | yes |
| `sequence_mode` | no | `composition` | `composition` or `full` | no |
| `mass_model` | no | `repeat_units` | `repeat_units` or `with_end_groups` | no |

Relative `output_dir` paths are resolved from the directory containing the
`.model` file. Parent-directory segments (`..`), empty segments, extensions,
spaces and glob characters are rejected. `init_volume` enables physical moles,
feed balances and reactor-volume analyses. The homo chain state has no
monomer-order sequence, so `sequence_mode` mainly preserves the shared model
contract.

### Complete declaration grammar

| Declaration | Syntax | Cardinality and constraints |
|---|---|---|
| Description | `desc "TEXT"` | optional, at most once |
| Sweep value | `var KIND NAME UNIT` | zero or more; `KIND` is `rate`, `param`, `species`, `monomer`, or `endgroup`; target must exist and be unique |
| Parameter | `param NAME VALUE [UNIT]` | only `init_volume` accepts a unit token |
| Species | `species NAME C0` | `C0>=0` mol/L; zero or more |
| Monomer | `monomer NAME C0 MW` | `C0>=0`, `MW>0`; zero or one for this engine |
| Feed component | `feed FEED NAME C` | repeat `FEED` for a mixture; each component once; `C>=0` mol/L |
| End group | `endgroup NAME MASS` | zero or more; contribution in g/mol |
| Polymer pool | `polymer NAME active\|dead` | required by any referenced macro channel |
| Rate | `rate NAME VALUE`, `rate NAME const VALUE`, or `rate NAME arr A EA` | fixed value and `A` are `>=0`; `EA` is J/mol |
| Elementary reaction | `rxn ...` | forms listed below |
| Polymer reaction | `macro KIND ...` | forms listed below |
| Scheduled action | `at`, `every`, or `from` | forms listed under Actions |
| Conditional action | `when CONDITION [and CONDITION ...] ACTION` | one-shot; `>` and `<` only |
| Memory policy | `at_memory SIZE save [stop]` | threshold with `B`, `KB`, `MB`, or `GB`; `save`, `stop`, or both |

### Species, monomer, end groups and pools

```text
species I 0.01
species R 0.0
species CTA_H 0.02
species CTA_ 0.0

monomer M 1.0 100.12

endgroup R 68.0
endgroup CTA 121.2
endgroup H 1.008
endgroup U -1.008
endgroup ACTIVE 0.0

polymer P active
polymer D dead
```

Initial small-species and monomer values are concentrations. `monomer` also
declares repeat-unit molar mass in g/mol. `endgroup` masses contribute only
under `with_end_groups`. `H` and `U` are the conventional right-end labels for
the two products of disproportionation: hydrogen-terminated and unsaturated,
respectively. The homo engine already defines `H=+1.008`, `U=-1.008`, and
`ACTIVE=0.0` g/mol relative to the active-end baseline; explicit `endgroup`
lines override those built-in contributions. Polymer pools are kinetic
populations used to determine reaction eligibility; they are distinct from a
stored chain `origin`, which records the mechanism that formed a chain record.

Feed mixtures are declared one component per line. Repeating the feed name
builds one constant-composition mixture:

```text
feed solution M 2.0
feed solution CTA_H 0.02
```

The concentration is in mol/L. A declared feed requires `param init_volume`.

### Rate constants

```text
rate kp 1.0e3
rate kt const 1.0e7
rate kd arr 1.0e15 125000.0
```

The first two forms are fixed rates. For an Arrhenius rate,

$$
k(T)=A\exp\left(-\frac{E_a}{RT}\right),
\qquad R=8.31446261815324\;\mathrm{J\,mol^{-1}\,K^{-1}}.
$$

`Ea` is in J/mol and `T` in K.

### Sweep metadata

One declaration may identify the value varied across a model series:

```text
var rate kp L_mol_s
var param temperature K
var species CTA_H mol_L
var monomer M mol_L
var endgroup CTA g_mol
```

Zero or more `var` lines are allowed. Each target name may appear only once. It describes a real declaration elsewhere in
the model; it does not assign or substitute the value.

For `var param`, only numeric built-ins are valid: `kmc_volume`,
`init_volume`, `temperature`, `t_end`, `max_steps`, `when_check_events`,
`seed`, and `dp_max`. `output_dir`, `sequence_mode`, and `mass_model` cannot be
optimization variables.

## Reaction syntax

### Elementary reactions

```text
rxn A -> B k
rxn A -> B + C k
rxn A -> 2B k
rxn A + B -> C k
rxn 2A -> C k
rxn A -> B k f
```

Elementary channels operate on small species and monomer reservoirs, not
polymer pools. Optional `f` is in `[0,1]`. It is the probability that a firing
forms products; reactants are consumed on every firing and successful plus
failed counts equal the total firing count.

### Polymer mechanisms

```text
macro init R + M -> P ki
macro prop P + M -> P kp
macro term_c P + P -> D ktc
macro term_d P + P -> D + D ktd
macro term_x P + X -> D kx
macro transfer P + CTA_H -> D + CTA_ ktr
macro transfer_m P + M -> D + P ktrm
macro deprop P -> P + M kdeprop
```

- `init` consumes a radical and monomer and creates a live DP=1 chain.
- `prop` consumes one monomer and increases DP by one.
- `term_c` combines two live chains into one dead chain.
- `term_d` produces two dead chains, conventionally with `H` and `U` ends.
- `term_x` caps a live chain with a small species.
- `transfer` terminates a chain and produces a new radical reservoir.
  `transfer_h` is not a valid model keyword; it is retained only as an
  internal/result origin label for chains formed by this mechanism.
- `transfer_m` terminates the old chain and creates a new monomer-born live
  chain.
- `deprop` returns one terminal monomer and never shortens below DP=1.

Reinitiation is ordinary `macro init` from a radical created by transfer.

## Actions

All schedule forms are:

```text
at TIME ACTION [ARGS...]
every PERIOD ACTION [ARGS...]
from START step PERIOD ACTION [ARGS...]
from START repeat COUNT every PERIOD ACTION [ARGS...]
```

`every PERIOD` starts at `t=0`. `from ... step ...` repeats without a count;
the bounded `repeat` form executes exactly `COUNT` times, except occurrences
after `t_end` are ignored.

Examples:

```text
every 0.1 print_info
every 0.5 save
every 1.0 save_chains
at 5.0 set_k kp 800.0
at 6.0 add_k kp 100.0
at 7.0 set_temp 353.15
at 8.0 add_temp 5.0
at 2.0 set_c CTA_H 0.01
at 3.0 add_c M 0.10
at 4.0 feed solution 1.0 mL
from 0.0 repeat 10 every 0.5 feed solution 0.10 mL
at 1.0 print "one second"
at 1.0 print_memory
```

The complete action vocabulary is:

| Action | Arguments | Effect |
|---|---|---|
| `print` | one quoted message | print and store a message |
| `print_info` | none | print progress information |
| `save` | none | save state and summary data |
| `save_chains` | none | save state plus chain-level data |
| `print_memory` | none | print memory diagnostics |
| `set_k`, `add_k` | `RATE VALUE` | replace or increment a rate |
| `set_temp`, `add_temp` | `VALUE` | replace or increment temperature |
| `set_c`, `add_c` | `NAME VALUE` | replace or increment concentration |
| `feed` | `FEED VOLUME [UNIT]` | add a declared mixture and volume |
| `stop` | none | stop cleanly; valid only after `when` |

Feed volume defaults to litres and accepts the same `L`/`mL` spellings as
`init_volume`. `add_c` is a physical amount change at the current volume;
`set_c` is a technical override and invalidates the physical balance of the
target entity.

Conditional actions fire once. Multiple predicates on one line are joined
with `and` and are evaluated on the same state:

```text
when X M > 0.50 print "half conversion"
when X M > 0.80 and c CTA_H < 1.0e-4 save_chains
```

The `stop` action is valid only after `when`. It requests clean termination
after the complete current conditional-action scan. It does not write a state
or chain snapshot implicitly; request those outputs explicitly:

```text
when X M > 0.80 save
when X M > 0.80 save_chains
when X M > 0.80 stop
```

The lines above are independent one-shot actions stored in source order. Since
`save`, `save_chains` and `stop` do not change chemical state, they observe the
same conversion, event and simulated time. A successful stop gives
`status=completed` and `termination_reason=stop_condition`.

Parameter-changing actions are time barriers: the SSA event loop does not
step past a scheduled boundary before executing the action. They are recorded
in `actions/` and `kinetic_parameter_sets/` and create a reconstructible state
snapshot.

A memory policy uses byte units `B`, `KB`, `MB` or `GB`:

```text
at_memory 3GB save stop
```

When the estimated memory use first reaches the threshold, `save` writes one
full logical snapshot **including chains**. This is intentionally stronger than
a normal state-only `save`. `stop` terminates the run after the threshold is
reached. With `save` but without `stop`, the threshold snapshot is written only
once and the simulation continues.

## SSA/kMC theory

For state `x`, every channel `j` has propensity `a_j(x)`. With
`a_0=Σ_j a_j`, two independent uniform random numbers select waiting time and
channel:

$$
\tau=-\frac{\ln u_1}{a_0},
\qquad
\sum_{j<\mu}a_j < u_2a_0 \leq \sum_{j\leq\mu}a_j.
$$

For molecule counts `n` in volume `V`, unimolecular propensity is
proportional to `k n`. A bimolecular cross reaction is proportional to
`k n_A n_B/(N_A V)`, while a same-species channel uses
`k n(n-1)/(N_A V)`. In the engine's rate convention the same-species symmetry
factor is carried by the declared rate constant rather than an additional
`1/2` in the propensity. The engine converts declared concentrations to
integer molecule populations and reports discretization diagnostics before a
run.

Polymer channels use the eligible live-chain population and relevant small
species/monomer population in the same SSA selection. Pool and DP eligibility
are checked before firing.

## Chain and molar-mass theory

For compressed rows `i` with multiplicity `c_i`, degree `d_i` and molar mass
`M_i`:

$$
DP_n=\frac{\sum_i c_i d_i}{\sum_i c_i},\qquad
DP_w=\frac{\sum_i c_i d_i^2}{\sum_i c_i d_i},
$$

$$
M_n=\frac{\sum_i c_i M_i}{\sum_i c_i},\qquad
M_w=\frac{\sum_i c_i M_i^2}{\sum_i c_i M_i},
$$

$$
M_z=\frac{\sum_i c_i M_i^3}{\sum_i c_i M_i^2},\qquad
Đ=\frac{M_w}{M_n}.
$$

Mass models are

$$
M_{\mathrm{repeat}}=DP\,M_0,
$$

and

$$
M_{\mathrm{ends}}=DP\,M_0+M_{\mathrm{left}}+M_{\mathrm{right}}.
$$

Undeclared end-group masses make a structural mass audit incomplete; declare
every observed end label when `with_end_groups` is used in the model.

## Output and analysis

The default output directory is `results/<model_stem>/`. The common schemas,
snapshot identity and tables are documented in [`../development/ARCHITECTURE.md`](../development/ARCHITECTURE.md).

```python
import pyslimmc as sl

run = sl.open("results/FRP01")
print(run.info())
print(run.t)
snap = run.final
print(snap.chains)
mwd = snap.mwd()
print(mwd.mn, mwd.mw, mwd.dispersity)
```

See [`../THEORY.md`](../THEORY.md) for derivations, [`../PYSLIMMC.md`](../PYSLIMMC.md) for the analysis tutorial and
[`PYSLIMMC_API.md`](PYSLIMMC_API.md) for the complete public reference.

## See also

- [`../THEORY.md`](../THEORY.md) — Theory and derivations
- [`../PYSLIMMC.md`](../PYSLIMMC.md) — pyslimmc guide
- [`PYSLIMMC_API.md`](PYSLIMMC_API.md) — pyslimmc API
