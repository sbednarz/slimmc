# Core concepts

This page explains the slimmc ideas that most strongly affect how a model is
built and how its results should be interpreted. It is intentionally less
formal than [`THEORY.md`](THEORY.md) and less exhaustive than the engine
references.

For each topic, the emphasis is: **what it is, why it matters, how slimmc
represents it, and what can go wrong when it is misunderstood.**

## Kinetic Monte Carlo and the SSA

### Kinetic Monte Carlo / direct SSA

**What it is.** slimmc follows a discrete stochastic reaction trajectory. At
each step, the direct Gillespie stochastic simulation algorithm (SSA) chooses
both the next reaction channel and the waiting time to that event.

**Why it matters.** A run is one stochastic realization, not a deterministic
ODE curve. Two runs with different seeds can differ even when every chemical
parameter is identical.

**In slimmc.** The engines maintain the current state, calculate channel
propensities, draw the next event, mutate the state, and advance simulation
time. The exact equations and molecularity conventions are in
[`THEORY.md`](THEORY.md).

**Caveat.** A reproducible seeded trajectory is not the same thing as an
ensemble average. Use replicate runs when the scientific question depends on
stochastic variability.

### Event, channel, and propensity

A **channel** is one possible reaction pathway. An **event** or **firing** is
one realized occurrence of that channel. A **propensity** is the instantaneous
SSA rate controlling how likely that channel is to fire next.

This distinction is important in analysis:

- firing shares describe what actually happened over a trajectory;
- rate/propensity shares describe competition between channels at saved states.

Do not interpret a firing share as an instantaneous kinetic rate.

## Simulation volume `V_kMC`

### What it is

`param kmc_volume` is the stochastic volume represented explicitly by the KMC
state. For a species at concentration $c_i$, the represented molecule count is
approximately

$$
n_i = c_i V_{kMC} N_A.
$$

The parser converts continuous concentrations into integer molecular
populations because the SSA operates on discrete objects.

### Why it matters

`V_kMC` controls two things at once:

1. **population size and statistical resolution** — a larger represented
   volume contains more molecules and usually more polymer chains;
2. **computational cost** — more represented molecules generally require more
   reaction events to reach the same conversion or physical time.

A very small `V_kMC` can therefore make rare species disappear through
rounding, leave only a few polymer chains in an important population, or make
distribution observables noisy. A larger value improves discrete resolution at
the cost of more work.

### How to choose it

There is no universal "correct" value. Choose `V_kMC` by **convergence of the
observable you care about**:

1. start with a practical value;
2. repeat the same model with a larger `V_kMC`;
3. compare conversion, moments, composition, MWD/CLD, rare populations, or the
   specific quantity used in your study;
4. continue until the scientific conclusion is insensitive to further volume
   increases within the precision you require.

A numerical scale example helps interpret the parameter. In the styrene
quick-start model, `[Sty] = 1 mol/L` and `V_kMC = 1e-15 L`, so the initial
represented styrene population is approximately

$$
1\times10^{-15} N_A \approx 6.0\times10^8
$$

molecules. This is an example, **not a recommended universal minimum**. A
different observable, chemistry, or rare population may require a substantially
different volume. `run.chain_count.live`, `.dead`, and `.total` are useful
diagnostics when judging whether the chain populations relevant to your
observable are well represented, but slimmc deliberately does not turn them
into universal pass/fail thresholds.

`slimmc --check` performs a discretization preflight and can warn about poor
initial integer representation, but passing the preflight is not a proof that a
chain-distribution statistic has converged.

### `kmc_volume` is not reactor volume

`kmc_volume` describes the stochastic representation. `init_volume` describes
the initial **physical reactor volume** used for physical amounts and semibatch
volume accounting. In feed runs both evolve consistently, but they remain
different concepts.

See also: [`MODEL_SYNTAX.md`](MODEL_SYNTAX.md),
[`SIMULATION_RESULTS.md`](SIMULATION_RESULTS.md), and
[`LIMITATIONS.md`](LIMITATIONS.md).

## Discretization

Concentrations are continuous quantities; SSA populations are integers.
Discretization is the mapping between them at `V_kMC`.

It matters most for small concentrations, small simulation volumes, rare
radicals, and small chain populations. A concentration that is chemically
nonzero can map to zero represented molecules if the chosen KMC population is
too small.

Discretization error and stochastic sampling noise are related but not
identical. Increasing `V_kMC` usually improves both, but a scientific convergence
check should still target the observable of interest.

## Trajectory, seed, and ensemble

A **trajectory** is one complete stochastic run. `param seed` fixes the random
number sequence and makes the same build/model/seed reproducible when the
reproducibility contract is satisfied.

An **ensemble** is a collection of trajectories, normally using different
seeds while keeping the scientific model fixed. Ensembles are useful for:

- estimating stochastic spread;
- confidence intervals or uncertainty summaries;
- checking whether a conclusion is dominated by one realization;
- rare-event studies.

Do not vary both kinetic parameters and seed and then call the resulting spread
"stochastic uncertainty" unless those sources of variation are separated in
the analysis.

## Chain populations and compressed rows

### Chain population

A chain population is the set of polymer chains stored at a `save_chains`
snapshot. slimmc can represent many chemically identical chains by one
**compressed row** plus an integer multiplicity (`abundance` / chain `count`).

This is an implementation-efficient representation of an exact discrete
population; it is not a histogram bin.

Consequences for analysis:

- the number of stored rows is not the number of physical chains;
- all moments and population aggregations must weight rows by multiplicity;
- pyslimmc does this automatically.

### Live and dead chains

A **live** chain belongs to an active polymer population and can participate in
channels allowed for that population. A **dead** chain no longer propagates in
the current engine model.

"Dead" is a model state, not a universal statement about every possible
post-polymerization chemistry. Reactions of dead chains that are outside the
engine scope remain outside the model; see [`LIMITATIONS.md`](LIMITATIONS.md).

## Pools

A **pool** is a chain population with a defined kinetic identity. In homo this
is often simply active versus dead. In copo, active pools encode terminal or
penultimate information required by the selected kinetic model.

Pool identity determines which propagation, depropagation, transfer, and
termination channels are eligible. Do not treat a pool name as a cosmetic
label.

## Terminal and penultimate copolymer models

### Terminal model

In a terminal model, the propagation rate depends on the identity of the last
incorporated monomer. For a binary system this naturally gives terminal-A and
terminal-B active populations.

### Penultimate model

In a penultimate model, propagation depends on the last **two** monomer
identities. slimmc therefore uses explicit penultimate active pools so that the
necessary kinetic memory is part of the state.

These models are kinetic assumptions, not merely alternative analysis modes.
The exact supported transitions are specified in
[`reference/COPO.md`](reference/COPO.md), and the equations are in
[`THEORY.md`](THEORY.md).

## Snapshot, `save`, and `save_chains`

A **snapshot** is the state recorded at a particular simulation time/event.

`save` records run-level state such as concentrations, conversion, temperature,
volume, channel information, and live/dead/total chain counts.

`save_chains` additionally stores the chain-resolved population required for
MWD, CLD, exact mass counts, SEC broadening, composition-resolved chain analysis, and other
chain-level calculations.

Therefore:

- a run may have many ordinary snapshots but only a few chain snapshots;
- chain distributions are only available where chain data were stored;
- `run.final` and `run.last_with_chains` can refer to different times.

See [`SIMULATION_RESULTS.md`](SIMULATION_RESULTS.md) for navigation rules.

## Final versus last

`final` means a successfully finalized terminal snapshot according to the run
contract. `last` means the last snapshot physically present in Storage.

For a completed run they may coincide. For failed or interrupted output they
need not. pyslimmc keeps the distinction explicit so that incomplete output is
not silently presented as a successful final result.

## Physical reactor volume and semibatch feed

A `feed` is a named, constant-composition mixture. A feed action adds an
explicit dose and increases physical reactor volume. The represented KMC volume
is scaled consistently so concentrations continue to map to discrete counts.

`init_volume` is required for feed models because physical dose accounting
needs a physical starting volume.

A scheduled semibatch profile is deliberately expressed as explicit portions:

```text
from 0 repeat 60 every 10 feed F 0.10 mL
```

The engine does not hide an arbitrary ramp/expression evaluator behind this
statement. Complex profiles can be written as multiple explicit portions.

## `add_c` versus `set_c`

`add_c` is a balance-preserving instantaneous addition/removal at constant
volume. It represents a physical inventory change and cannot drive the species
population below zero.

`set_c` forcibly replaces a concentration. It is useful for technical model
experiments but is **not** interpreted as a physical material operation. The
engine warns, and the physical balance for that species becomes inapplicable;
pyslimmc reports this explicitly rather than fabricating a balance.

## Mass model

A chain contains repeat units and may also carry stored end groups. The model
parameter `mass_model` determines the canonical chain-mass interpretation:

- `repeat_units` — sum repeat-unit masses only;
- `with_end_groups` — include stored end-group contributions.

Always state which interpretation is used when comparing calculated masses to
an experiment or external model. `mass_audit()` checks whether the requested
interpretation is chemically supported by the stored metadata.

## Exact source populations, counts, and distributions

A stored chain population is the source of truth. pyslimmc separates four
operations that should not be conflated:

```text
selected ChainPopulation
        |
        +-- dp_counts()   -> exact DP -> N projection
        +-- mass_counts() -> exact M  -> N projection
        +-- cld()/mwd()   -> normalized discrete mathematical forms
        +-- sec()         -> continuous apparent instrumental response
```

Compressed rows are always weighted by their stored multiplicity.

## Exact moments versus represented curves

`DPn`, `DPw`, `DPz`, `Mn`, `Mw`, and `Mz` are source-population statistics.
They are independent of CLD/MWD form, plotting, SEC broadening, and numerical
sampling. If a required source aggregate was not stored and chains are not
available, pyslimmc reports the quantity as unavailable rather than estimating
it from a displayed curve.

## MWD and CLD forms

Both analyses accept:

```text
form = number | mass | z | log
```

`number`, `mass`, and `z` specify the discrete weighting. `log` is the canonical
mass-weighted logarithmic form: the exact mass fractions are placed on
`log10(DP)` or `log10(M)` support.

For MWD, actual neutral chain masses are used. General MWD is therefore an
independent projection of the selected chain population and is not obtained by
multiplying or transforming an aggregated CLD.

For mass-weighted CLD, pyslimmc sums the actual mass carried by all chains in a
DP class. Only in a simple homopolymer with a linear DP-to-mass relation does
this reduce to the familiar `DP * count` weighting.

## Discrete logarithmic support and the Jacobian

For a continuous mass density `w(M) = dW/dM`, changing variable to
`u = log10(M)` gives

```text
dW/dlog10(M) = M ln(10) w(M).
```

That Jacobian is required for continuous densities. It is **not** applied to
an exact discrete atom: a point carrying mass fraction `w_i` still carries
`w_i` after its location is re-expressed as `log10(M_i)`.

Therefore `mwd(form="log")` is an exact discrete measure on logarithmic
support, not a finite continuous density. `plt.xscale("log")` only changes
visual presentation.

## Multi-series normalization

Multi-series analyses are composed from already valid single-population
results. `per_series` independently normalizes each series and permits overlap.
`combined` uses one exact common source denominator and therefore requires
pairwise-disjoint populations. No common numerical grid is introduced.

## SEC broadening

SEC is a separate forward representation of instrumental broadening. With
`u = log10(M)`, pyslimmc uses

```text
S(u) = sum_i w_i G_sigma(u - u_i),
```

where `w_i` are exact polymer-mass fractions and `G_sigma` is a normalized
Gaussian response. Consequently `S(u)` is a continuous apparent
`dW_app/dlog10(M)` curve. The public width parameter is `sigma_log10M`; under a
linear SEC calibration `log10(M) = a - b v`, it corresponds to `b * sigma_v`.

Histogramming, KDE, and generic Gaussian smoothing are intentionally not part
of the CLD/MWD API.

## Sequence modes and microstructure

`sequence_mode composition` stores compact composition/microstructure counters
without retaining every complete monomer sequence. It supports aggregated
composition, dyad/triad and block statistics that the engine explicitly stores.

`sequence_mode full` additionally retains complete sequences and enables
position-resolved, motif, n-gram and full-sequence consistency analyses.

Full sequences cost more memory. Choose them because the scientific question
requires sequence-resolved information, not by default.

## Validation, audit, and convergence

These are different checks:

- `run.validate(strict=True)` checks Storage and run invariants;
- `run.mass_audit()` checks chain-mass accounting under a selected mass model;
- `run.reproducibility.verify()` checks stored reproducibility/provenance data;
- varying `V_kMC` checks numerical/statistical convergence of an observable;
- replicate seeds assess trajectory-to-trajectory stochastic variability.

A run can pass structural validation and still be scientifically under-resolved
for a particular distribution statistic. Validation is necessary; convergence
and uncertainty analysis answer different questions.

## See also

- equations and rate conventions: [`THEORY.md`](THEORY.md)
- model statements: [`MODEL_SYNTAX.md`](MODEL_SYNTAX.md)
- meaning of stored output: [`SIMULATION_RESULTS.md`](SIMULATION_RESULTS.md)
- practical analysis: [`PYSLIMMC.md`](PYSLIMMC.md) and
  [`COOKBOOK.md`](COOKBOOK.md)
- scope boundaries: [`LIMITATIONS.md`](LIMITATIONS.md)
- exact engine/API contracts: [`reference/`](reference/)
