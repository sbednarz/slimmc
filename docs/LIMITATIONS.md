# Scope and limitations

slimmc intentionally models a defined subset of polymer reaction kinetics. This
page separates deliberate scope boundaries from numerical/statistical
limitations so that unsupported chemistry is not mistaken for a bug or an
unverified approximation.

For equations and implemented mechanisms see [`THEORY.md`](THEORY.md); for the
formal per-engine contract see [`reference/HOMO.md`](reference/HOMO.md) and
[`reference/COPO.md`](reference/COPO.md).

## Deliberate engine scope

### Linear-chain state representation

The current engines represent linear chains. They do not provide a general
polymer-topology graph engine. Consequently, the current core scope excludes:

- branching as a general topology operation;
- crosslinking and network formation;
- gelation/topology-percolation calculations;
- general intramolecular cyclization.

These features require a different chain/topology representation rather than a
small extension of the current linear-chain model.

### Controlled radical polymerization mechanisms

slimmc can represent generic elementary and macro kinetic channels that fit the
implemented state model, but it is **not** a dedicated mechanistic ATRP, RAFT,
or other reversible-deactivation radical-polymerization engine. A model should
not be described as mechanistic ATRP/RAFT merely because some effective
activation/deactivation or transfer-like reactions can be approximated with
available channels.

### Chain-length-dependent kinetics

Rate constants are not general functions of chain DP. Diffusion-controlled or
chain-length-dependent propagation/termination laws are therefore outside the
current kinetic contract unless reduced to supported time/state changes by the
user.

### Energy balance

Temperature may be set initially and changed by scheduled/conditional actions.
slimmc does not solve a non-isothermal reactor energy balance that couples heat
release, heat transfer and temperature dynamically.

### Post-polymerization chemistry of dead chains

Dead chains are product populations rather than a general reactive topology.
Arbitrary subsequent chemistry of dead polymer molecules is outside the current
engine scope.

### Copolymer depropagation and memory

The copo engine supports terminal depropagation under its documented sequence
requirements. General penultimate depropagation is not part of the current
contract. Sequence-dependent calculations are limited by the selected
`sequence_mode` and the information retained in Storage.

## Numerical and stochastic limitations

### Finite KMC population

A stochastic run represents a finite discrete population controlled by
`kmc_volume`. If the represented population is too small, rare species may be
poorly sampled and distributional quantities may be noisy. Increasing
`kmc_volume` increases the number of represented molecules and usually the
computational cost.

There is no universal chain-count threshold that guarantees accuracy for every
observable. Convergence should be demonstrated for the quantity of interest.

### Initial discretization

Declared concentrations are converted to integer molecular counts. This
requires rounding. If a positive initial concentration would round to zero, the
discretization preflight rejects the model rather than silently dropping the
species. Accepted runs report the realized discrete state.

### One trajectory is one stochastic realization

With a fixed model and seed, slimmc provides deterministic reproducibility of
the pseudorandom trajectory. That does not turn one trajectory into an
ensemble average. Studies that require uncertainty estimates should use
independent seeds/replicates and analyse between-run variability.

### Rare events

A chemically allowed channel with a very small propensity may fire rarely or
not at all in a finite trajectory. Absence of a realized firing is not by
itself proof that the underlying propensity is zero.

### Distribution smoothing

MWD/CLD smoothing (`hist`, Gaussian, KDE where available) is a representation
of the stored discrete population. Exact moments are calculated from the
population itself. Plotting choices can alter the appearance of the curve and
must not be interpreted as new molecular information.

### Mass model assumptions

`mass_model repeat_units` uses repeat-unit masses. `with_end_groups` additionally
uses declared end-group contributions. Missing or inappropriate end-group
masses can make absolute molar-mass interpretation incorrect even when chain DP
is correct.

### Chain mass spectrum is not m/z

`chain_mass_spectrum()` reports neutral chain masses derived from simulated
chains. It is not an ionization/adduct/isotope model and does not directly
predict experimental mass-spectrometric `m/z` spectra.

### Stored sequence information

`sequence_mode composition` stores composition and supported summary
microstructure information without literal full sequences. Analyses requiring
complete monomer order need `sequence_mode full` and an appropriate
`save_chains` snapshot.

### Numeric ceilings

DP and population representations have explicit implementation limits, exposed
and validated by the engine contracts. For example, `dp_max` cannot exceed the
32-bit DP ceiling documented by the engine. Do not raise such limits in model
files beyond the accepted parser/state contract.

## Engine differences are intentional where documented

Homo and copo share a common model/storage family but are not identical
implementations. Terminal/penultimate pool logic, sequence state and some
mechanism availability are engine-specific. The exact references are the
source of truth when a shared guide and an engine-specific rule differ.

## Before designing a study

Ask four questions:

1. Can the required chemistry be expressed with the implemented linear-chain
   channels?
2. Does the requested observable require full sequences or saved chains?
3. Is `kmc_volume` large enough for the observable to converge?
4. Does the scientific claim require an ensemble rather than one seeded run?

If the answer to the first question is no, changing analysis code cannot make
the missing chemistry valid; the engine model itself is outside scope.

## See also

- [`CONCEPTS.md`](CONCEPTS.md) — Core concepts
- [`THEORY.md`](THEORY.md) — Theory and assumptions
- [`MODEL_SYNTAX.md`](MODEL_SYNTAX.md) — Model syntax
