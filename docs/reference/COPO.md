# slimmc copolymer engine



This is the user and language reference for `slimmc-copo`. It describes the
current implementation for linear binary copolymerization and terminal
terpolymerization.

## Purpose and scope

The copolymer engine stores the exact sequence of every live chain and keeps
composition, terminal/penultimate units, dyads, triads and block statistics.
Dead chains are compact summaries; their full sequence may be retained below
a configurable DP limit.

Implemented mechanisms:

- two or three monomers;
- terminal and explicit penultimate propagation pools;
- elementary reactions, initiation and propagation;
- combination, disproportionation and small-species capping;
- transfer, transfer-to-monomer and reinitiation;
- terminal depropagation;
- fixed and Arrhenius rates;
- scheduled and conditional actions;
- repeat-unit and end-group-aware masses;
- optional oligomer tables and detailed diagnostics.

Batch and constant-composition, portion-wise semibatch feeds are supported,
including changing reactor and KMC volumes. Not represented:
branching/crosslinks, gel topology, spatial/diffusion effects, penultimate
depropagation and general N-monomer kinetics above three monomers.

## Running the engine

The family dispatcher chooses copo for a model with two or three `monomer`
declarations:

```bash
make build
bin/slimmc path/to/model.model
```

Direct component use is intended mainly for engine development:

```bash
make -C copo build
copo/slimmc-copo path/to/model.model
```

For normal use and diagnostics, prefer the family dispatcher:

```bash
slimmc --check path/to/model.model
slimmc --debug path/to/model.model
slimmc --trace-channels 100000 path/to/model.model
```

`--check` validates without simulating. Completed-run metadata and advertised
diagnostic files are recorded in `run_metadata.json`.

## Model file

A minimal complete binary terminal model is shown below. It is deliberately
small and uses abstract monomers so the syntax remains the focus:

```model
desc "Minimal binary terminal copolymerization"

param output_dir "results/COP_MIN"
param kmc_volume 1.0e-18
param temperature 343.15
param t_end 100
param seed 101
param sequence_mode composition

monomer A 2.0 100.12
monomer B 2.0 128.17
species R 1.0e-4

polymer PA active
polymer PB active
polymer D dead

rate ki_a 1.0e4
rate ki_b 1.0e4
rate kp_aa 500.0
rate kp_ab 450.0
rate kp_ba 450.0
rate kp_bb 500.0
rate ktc 1.0e7

macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
macro term_c PA + PA -> D ktc
macro term_c PA + PB -> D ktc
macro term_c PB + PB -> D ktc

at 0 save
at 100 save_chains
```

Comments begin with `#`; keywords are case-sensitive. Quote free text and
paths containing whitespace or `#`.

Declaration names follow `[A-Za-z_][A-Za-z0-9_]*`. Keywords are reserved, and
names differing only by ASCII letter case are rejected. An `output_dir` may be
relative or absolute, but each non-root segment follows the same identifier
grammar; `.`, `..`, spaces, extensions, and glob characters are rejected in
model paths.

### Metadata and global parameters

```text
desc "terminal binary copolymerization"

param kmc_volume 1.0e-19
param init_volume 100 mL
param temperature 333.15
param t_end 2.0
param max_steps 200000
param seed 101
param output_dir "results/COP01"
param sequence_mode full
param dp_max 1000000
param mass_model with_end_groups
param when_check_events 1000
```

The parameter contract is identical in both engines except where noted:

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
| `dp_max` | no | `2147483647` | integer in `1..2147483647` | yes |
| `sequence_mode` | no | `composition` | `composition` or `full` | no |
| `mass_model` | no | `repeat_units` | `repeat_units` or `with_end_groups` | no |

`sequence_mode full` retains complete order. `composition` retains composition,
terminal summaries, dyads, triads and block statistics without literal
sequences. Copolymer `deprop` requires `full`. Oligomer subsets are selected in
`pyslimmc` with explicit DP filters rather than by a model parameter.

Relative `output_dir` paths are resolved from the model directory. The same
path validation and volume-unit rules as homo apply.

### Complete declaration grammar

| Declaration | Syntax | Cardinality and constraints |
|---|---|---|
| Description | `desc "TEXT"` | optional, at most once |
| Sweep value | `var KIND NAME UNIT` | zero or more; `KIND` is `rate`, `param`, `species`, `monomer`, or `endgroup`; target must exist and be unique |
| Parameter | `param NAME VALUE [UNIT]` | only `init_volume` accepts a unit token |
| Monomer | `monomer NAME C0 MW` | exactly two or three; `C0>=0`, `MW>0` |
| Species | `species NAME C0` | `C0>=0` mol/L; zero or more |
| Feed component | `feed FEED NAME C` | repeat `FEED` for a mixture; component once; `C>=0` mol/L |
| End group | `endgroup NAME MASS` | zero or more; contribution in g/mol |
| Polymer pool | `polymer NAME active\|dead` | exactly one dead pool is required; live pools follow the kinetic topology |
| Rate | `rate NAME VALUE`, `rate NAME const VALUE`, or `rate NAME arr A EA` | fixed value and `A` are `>=0`; `EA` is J/mol |
| Elementary reaction | `rxn ...` | forms listed below |
| Polymer reaction | `macro KIND ...` | forms listed below |
| Scheduled action | `at`, `every`, or `from` | forms listed under Actions |
| Conditional action | `when CONDITION [and CONDITION ...] ACTION` | one-shot; `>` and `<` only |
| Memory policy | `at_memory SIZE save [stop]` | threshold with `B`, `KB`, `MB`, or `GB`; `save`, `stop`, or both |

### Monomers and species

```text
monomer A 0.20 100.12
monomer B 0.20 128.17
monomer C 0.10 142.20

species R 1.0e-4
species CTA 1.0e-3
species Rcta 0.0
species SOLV 1.0
```

`monomer NAME c0 MW` declares initial concentration and repeat-unit molar mass.
The engine accepts two or three monomers. Identifiers may have multiple
characters; canonical stored sequences use `|` separators, for example
`Sty|MMA|Sty`, so tokenization is unambiguous.

Feed mixtures use one line per component. Repeating a feed name extends that
mixture:

```text
feed solution A 1.5
feed solution B 0.5
feed solution CTA 0.01
```

Concentrations are in mol/L. The same feed can later be dosed by scheduled or
conditional `feed` actions.

### End groups and mass model

```text
endgroup R 68.0
endgroup CTA 121.2
endgroup H 1.008
endgroup U -1.008
endgroup ACTIVE 0.0
```

`H` and `U` are the conventional right-end labels assigned to the two
disproportionation products: hydrogen-terminated and unsaturated, respectively.
For `mass_model with_end_groups`, declare their molar-mass contributions
explicitly; `+1.008` and `-1.008` g/mol preserve the zero net end-group mass
change of a disproportionation pair relative to `ACTIVE=0.0`. Polymer pools
represent current kinetic state/eligibility. They are independent of stored
chain `origin`, which records the mechanism that formed the chain record.

`repeat_units` uses only monomer masses. Model value `with_end_groups` also adds declared
left and right end-group masses. Missing declarations are visible in the mass
audit and should be fixed before quantitative MWD work.

### Polymer pools

```text
polymer PA active
polymer PB active
polymer D dead
```

Terminal pools accept chains with the matching last unit. Explicit
penultimate pools can be declared:

```text
polymer PAA active
polymer PAB active
polymer PBA active
polymer PBB active
```

Their terminal and penultimate requirements are inferred from macro channels.
Pool compatibility is validated before propensity calculation and firing.

### Rates

```text
rate kp_ab 120.0
rate kt const 1.0e7
rate kd arr 1.0e15 125000.0
```

Arrhenius rates use

$$
k(T)=A\exp\left(-\frac{E_a}{RT}\right),
\qquad R=8.31446261815324\;\mathrm{J\,mol^{-1}\,K^{-1}}.
$$

### Sweep metadata

Zero or more varied declarations may be identified for `pyslimmc.scan()`:

```text
var rate kp_ab L_mol_s
var param temperature K
var species CTA mol_L
var monomer A mol_L
var endgroup CTA g_mol
```

The `var` line is metadata; the actual value remains in its normal
declaration.

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

These operate on small species and monomers, never polymer pools. Optional
`f∈[0,1]` is a success probability: reactants are consumed per firing and
products form on successful firings.

### Initiation

```text
macro init R + A -> PA ki_a
macro init R + B -> PB ki_b
```

Each firing creates a live DP=1 chain whose sequence and composition contain
the consumed monomer.

### Terminal propagation

```text
macro prop PA + A -> PA kp_aa
macro prop PA + B -> PB kp_ab
macro prop PB + A -> PA kp_ba
macro prop PB + B -> PB kp_bb
```

The input pool selects the old terminal state; the output pool represents the
new terminal state after adding the incoming monomer.

### Explicit penultimate propagation

For binary monomers, first-step terminal channels feed penultimate pools:

```text
macro prop PA + A -> PAA kp_a_a
macro prop PA + B -> PAB kp_a_b
macro prop PB + A -> PBA kp_b_a
macro prop PB + B -> PBB kp_b_b
```

Then eight channels describe `(penultimate,terminal,incoming)` combinations:

```text
macro prop PAA + A -> PAA kp_aa_a
macro prop PAA + B -> PAB kp_aa_b
macro prop PAB + A -> PBA kp_ab_a
macro prop PAB + B -> PBB kp_ab_b
macro prop PBA + A -> PAA kp_ba_a
macro prop PBA + B -> PAB kp_ba_b
macro prop PBB + A -> PBA kp_bb_a
macro prop PBB + B -> PBB kp_bb_b
```

### Termination

Combination:

```text
macro term_c PA + PB -> D ktc
```

stores `sequence1 + reverse(sequence2)` as one dead chain and accounts for the
new block boundary.

Disproportionation:

```text
macro term_d PA + PB -> D ktd
macro term_d PA + PB -> D + D ktd
```

creates two dead summaries, conventionally with `H` and `U` right ends. The
two-output spelling requires the same dead pool on both sides.

Small-species capping:

```text
macro term_x PA + Cap -> D kcap
```

consumes one live chain and one small species; the species label becomes the
right end.

### Transfer and reinitiation

```text
macro transfer PA + CTA -> D + Rcta ktr_a
macro transfer PB + CTA -> D + Rcta ktr_b
macro init Rcta + A -> PA ki_cta_a
macro init Rcta + B -> PB ki_cta_b
```

Transfer preserves the old chain's DP, mass, composition and microstructure,
consumes the transfer species and produces a radical reservoir. Reinitiation
is an ordinary `macro init`; omit those channels to model a non-reinitiating
radical.

Solvent transfer uses the same pattern with a solvent species.

### Transfer to monomer

```text
macro transfer_m PA + A -> D + PA ktrm_a
```

The old chain becomes dead with right end `H`; one monomer is consumed and a
new DP=1 chain is created with a monomer-derived left end such as `A_tr`.

### Terminal depropagation

```text
macro deprop PA -> PB + A kdeprop_ba
```

This represents

$$
\ldots-B-A^\bullet \longrightarrow \ldots-B^\bullet + A.
$$

Eligibility requires `DP≥2`, terminal `A`, and penultimate `B`. The operation
returns one monomer, shortens the sequence, updates composition, mass, dyads,
triads, blocks and terminal metadata, and moves the chain to `PB`.

## Actions and boundaries

All schedule forms are:

```text
at TIME ACTION [ARGS...]
every PERIOD ACTION [ARGS...]
from START step PERIOD ACTION [ARGS...]
from START repeat COUNT every PERIOD ACTION [ARGS...]
```

`every` begins at `t=0`; the bounded form executes exactly `COUNT` times unless
later occurrences lie beyond `t_end`.

```text
every 0.2 print_info
every 0.5 save
every 1.0 save_chains
at 1.0 set_k kp_ab 200.0
at 1.2 add_k kp_ab 20.0
at 1.4 set_temp 353.15
at 1.6 add_c CTA 1.0e-4
at 1.8 feed solution 0.50 mL
from 0.0 repeat 10 every 0.2 feed solution 0.10 mL
at 0.0 print "start"
at 0.0 print_memory
```

Available actions are `print`, `print_info`, `save`, `save_chains`,
`print_memory`, `set_k`, `add_k`, `set_temp`, `add_temp`, `set_c`, `add_c`,
`feed`, and conditional-only `stop`. `print` requires one quoted message;
`print_info`, `save`, `save_chains`, `print_memory`, and `stop` take no
arguments. Feed volume defaults to litres and accepts `L`/`mL` case variants.
`set_c` is a technical override and invalidates the physical balance of its
target; `add_c` and `feed` preserve explicit amount accounting.

Conditional one-shot actions:

```text
when X A > 0.20 save_chains
when X A > 0.50 print_info
when c CTA < 1.0e-5 save
```

Multiple predicates on one line are joined with `and`. There is no two-term
limit, so ternary and higher declared-monomer models can use one joint
condition:

```text
when X A > 0.80 and X B > 0.90 save
when X A > 0.80 and X B > 0.90 save_chains
when X A > 0.80 and X B > 0.90 stop

when X A > 0.80 and X B > 0.90 and X C > 0.70 stop
```

All atoms on one line are evaluated on the same state. Separate `when` lines
remain independent and therefore provide OR semantics. No parentheses or
inline `or` are part of this contract.

The `stop` action is valid only after `when`. It requests clean termination
after the complete current conditional-action scan and never performs `save`
or `save_chains` implicitly. Users must place explicit output actions before
(or elsewhere in) the matching set of one-shot lines. A successful stop gives
`status=completed` and `termination_reason=stop_condition`.

`X` is global conversion; `X A` is monomer-specific conversion. Scheduled
actions are time barriers. If the next SSA event would exceed `t_end`, it is
not fired; the final state is written exactly at `t_end`.

Memory policy:

```text
at_memory 6GB save stop
```

When the estimated memory use first reaches the threshold, `save` writes one
full logical snapshot **including chains**. A normal scheduled/conditional
`save` is state-only; use `save_chains` when you explicitly want chains at such
checkpoints. `stop` terminates after the threshold is reached. With `save` but
without `stop`, the threshold snapshot is written only once and the simulation
continues.

## SSA and copolymerization theory

The engine uses the direct SSA. For total propensity `a0=Σa_j`, waiting time
and selected channel are

$$
\tau=-\frac{\ln u_1}{a_0},
\qquad
\sum_{j<\mu}a_j < u_2a_0 \leq \sum_{j\leq\mu}a_j.
$$

For molecule counts in volume `V`, a cross-species bimolecular channel uses
`k n_A n_B/(N_A V)`. A same-species channel uses
`k n(n-1)/(N_A V)`; the declared rate convention carries the symmetry factor.

For a terminal binary model, a chain ending in `i` adds monomer `j` with a
propensity proportional to `k_{ij}[P_i^\bullet][M_j]`. The instantaneous
probability of adding A rather than B to terminal `i` is therefore

$$
p(A\mid i)=\frac{k_{iA}[A]}{k_{iA}[A]+k_{iB}[B]}.
$$

Reactivity ratios are

$$
r_A=\frac{k_{AA}}{k_{AB}},\qquad
r_B=\frac{k_{BB}}{k_{BA}}.
$$

With monomer mole fractions `f_A` and `f_B`, the terminal Mayo–Lewis
instantaneous copolymer composition is

$$
F_A=\frac{r_Af_A^2+f_Af_B}
{r_Af_A^2+2f_Af_B+r_Bf_B^2},
\qquad F_B=1-F_A.
$$

The explicit penultimate model replaces `k_{ij}` by `k_{hij}`, where `h` is
the penultimate unit. pyslimmc can compare simulated transition/composition
statistics with terminal and penultimate predictions.

## Sequence and memory model

Live chains retain exact token sequences. Dead chains retain composition,
mass, terminal information, dyads, triads and block statistics regardless of
whether full text is dropped. With

```text
param sequence_mode full
```

complete sequence order is retained for saved live and dead chains. In
`composition` mode, complete composition, terminal summaries, dyads, triads,
and block statistics remain available without literal sequence order. The
optional `memory/` table reports live and dead storage estimates.

For sequence `s_1…s_N`, dyad count `N_{ij}` counts adjacent pairs
`(s_k,s_{k+1})=(i,j)`; triads analogously use triples. A block is a maximal
contiguous run of one monomer. These counts survive compact dead summaries.

## Omitted propagation transitions

A missing `macro prop` declaration is a legal zero-rate transition. The engine
does not infer a channel. `--check` reports every omitted active-pool/monomer
pair as a warning so the omission is explicit.

## Instantaneous composition in pyslimmc

For binary terminal models, `run.F.ins` identifies the four propagation pairs
from the declared `macro prop` topology. Kinetic parameter names are arbitrary;
`kp_aa`, `kp_ab`, `kp_ba`, and `kp_bb` are conventions, not requirements.
`F.ins` is the terminal Mayo--Lewis prediction, while `F.int` and `F.cum` are
derived from the simulated trajectory.

## Output and analysis

See [`../development/ARCHITECTURE.md`](../development/ARCHITECTURE.md) for common tables and copolymer-specific Storage outputs.

```python
import pyslimmc as sl

run = sl.open("results/COP01")
print(run.copolymerization.cumulative_composition())
print(run.microstructure.dyads())
print(run.final.chains.where(dp_max=10))
print(run.final.mwd().dispersity)
```

See [`../THEORY.md`](../THEORY.md) for derivations, [`../PYSLIMMC.md`](../PYSLIMMC.md) for the analysis tutorial and
[`PYSLIMMC_API.md`](PYSLIMMC_API.md) for the complete public reference.


## See also

- [`../THEORY.md`](../THEORY.md) — Theory and derivations
- [`../PYSLIMMC.md`](../PYSLIMMC.md) — pyslimmc guide
- [`PYSLIMMC_API.md`](PYSLIMMC_API.md) — pyslimmc API
