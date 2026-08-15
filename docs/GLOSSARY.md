# slimmc and pyslimmc glossary

This glossary collects chemical, mathematical, simulation and software terms
used across the slimmc family. Definitions describe how the term is used in
this project and are written primarily for a polymer chemist using the engines
and `pyslimmc`.

## Polymer chemistry

### Active chain / live chain

A polymer molecule carrying a reactive chain end and eligible for further
propagation or another radical reaction. “Live” here means active in the
kinetic model; it does not imply living/controlled polymerization.

### Backbiting

Intramolecular hydrogen abstraction by a growing radical, often producing a
mid-chain radical and potentially short-chain branching. It is not represented
by the current linear-chain slimmc state model.

### Block

A maximal uninterrupted sequence of one monomer type within a copolymer chain.
Block-length statistics describe sequence organization, not separate physical
polymer molecules.

### Branching

Formation of a polymer topology in which one repeat-unit path divides into two
or more paths. Short-/long-chain branching and explicit graph topology are not
implemented in the current engines.

### Chain transfer agent (CTA)

A species that reacts with a growing radical, terminates that chain and may
create a new radical. A CTA can control chain length only when transfer and
subsequent reinitiation kinetics are represented appropriately.

### Combination

A termination mode in which two live radicals form one dead chain. The product
DP and repeat-unit inventory are approximately the sums of the two reactant
chains, subject to the declared end-group mass model.

### Controlled/living radical polymerization

Mechanisms such as ATRP, RAFT, NMP or iodine-transfer polymerization that use
reversible activation/transfer to regulate radical populations and chain
growth. They are not implemented as dedicated slimmc reaction families.

### Copolymer

A polymer containing repeat units from at least two different monomers. The
current copo engine supports binary copolymers and terpolymers with explicit
sequence-aware chain records.

### Crosslinking

Formation of covalent links between polymer chains, eventually producing a
network or gel. Current slimmc chains are linear, so crosslink topology and the
gel fraction cannot be simulated explicitly.

### Cyclization

Intramolecular formation of a ring by reaction of two sites on the same
molecule. It requires molecular connectivity/topology and is outside the
current linear-chain implementation.

### Dead chain

A polymer molecule that no longer carries an active radical in the declared
model. Dead chains may still contain structural information relevant to mass,
end groups, sequence and product distributions.

### Degree of polymerization (DP)

Number of repeat units in one polymer molecule. In a copolymer it is the sum of
the per-monomer repeat-unit counts, not the abundance of that molecular type.

### Depropagation

Removal of the terminal repeat unit from an active chain, returning monomer and
shortening the chain by one. Eligibility depends on DP and terminal/penultimate
identity.

### Disproportionation

A termination mode in which two live radicals produce two dead chains rather
than one combined molecule. It conserves chain count differently from
combination and therefore changes MWD even at similar radical-consumption rate.

### Dyad

Two adjacent repeat units in a copolymer sequence, for example `A|B`. Dyad
fractions quantify nearest-neighbour sequence structure.

### End group

A chemical group at a polymer-chain end originating from initiation, transfer,
termination or another reaction. Its molar-mass contribution affects chain
masses only with the `with_end_groups` mass model.

The conventional labels `H` and `U` are used for the two products of radical
disproportionation: `H` is the hydrogen-terminated product and `U` the
unsaturated product. In the homo engine their built-in mass contributions are
`+1.008` and `-1.008` g/mol relative to the active-end baseline, respectively,
so the pair has zero net end-group mass change. The copo engine uses the same
labels for disproportionation products, but their mass contributions must be
declared explicitly when `with_end_groups` masses are required. `ACTIVE` is the
active radical-end label and normally has a `0.0` g/mol contribution.

### Gel effect / Trommsdorff effect

Autoacceleration associated mainly with diffusion-limited radical termination
as viscosity increases. slimmc can impose conversion-/time-dependent rate
constants phenomenologically, but it does not derive the gel effect from
spatial diffusion.

### Homopolymer

A polymer formed from one monomer type. Its linear chain can be stored compactly
because a full monomer-identity sequence would contain no additional
information.

### Initiation

Creation of the first repeat-unit-containing active chain by reaction of a
radical with monomer. Initiator decomposition and macromolecular initiation
are separate channels when both steps are declared.

### Long-chain branching (LCB)

A branch whose length is comparable to the main chains and strongly influences
rheology and molecular topology. LCB requires a graph-like molecular model and
is not represented by current slimmc chains.

### Monomer

A small molecule consumed to add one repeat unit to a chain. The free-monomer
population and incorporated repeat-unit inventory are stored separately.

### Oligomer

A low-DP polymer molecule below a declared analysis cutoff. The cutoff is an
analysis/model choice; “oligomer” is not a universal DP boundary.

### Polymer pool

A kinetic population used by the engine to decide which chains are eligible for
a reaction channel. A pool answers **where a chain is now in the kinetic
model** (for example an active terminal pool or a dead product pool). It can
change as reactions move a chain between kinetic states.

### Chain origin

A provenance label stored with a chain record that identifies the mechanism
that created that record, such as `init`, `transfer_m`, `term_c`, `term_d_H`,
`term_d_U`, or `term_x`. Origin answers **how the chain record was formed**; it
is not the current kinetic pool and not a log of the most recent reaction.

### Pendant double bond

An unreacted polymerizable double bond attached to an existing chain. Its
reaction can create branching/crosslinking and is not represented by the
current linear-chain engines.

### Penultimate model

A copolymerization model in which propagation depends on both terminal and
penultimate repeat-unit identities. It requires a larger rate tensor and
cannot generally be reduced to two terminal reactivity ratios.

### Propagation

Addition of monomer to an active chain, increasing DP by one and updating the
terminal/penultimate state. In copo, the incoming monomer also updates sequence
and composition statistics.

### Pulsed-laser polymerization (PLP)

Polymerization initiated by timed radical pulses, commonly used for propagation
coefficient measurements. slimmc's exact-time concentration/action boundaries
are particularly suitable for representing idealized PLP pulse schedules.

### Radical pool

A group of active chains sharing kinetic eligibility, such as a terminal or
penultimate identity. Pool membership is a reaction-state classification, not
necessarily a physically separated phase.

### Reactivity ratio

For a binary terminal model, ratios such as
`r_A=k_AA/k_AB` and `r_B=k_BB/k_BA`. They express relative preference of a
terminal radical for like versus unlike monomer.

### Reinitiation

Formation of a new active chain by a radical created in transfer or another
reaction. Keeping reinitiation separate from transfer allows its efficiency and
delay to be modelled explicitly.

### Repeat unit

The monomer-derived structural unit incorporated in a polymer chain. Its molar
mass need not equal the complete chain mass contribution when end groups or
eliminated fragments are considered.

### Sequence

Ordered list of monomer identities along a copolymer chain. Aggregate
composition can remain known when a stored literal sequence is truncated, but
order-dependent analyses then become incomplete.

### Short-chain branching (SCB)

Relatively short branches, often produced by backbiting followed by
propagation. SCB is not represented by the current linear-chain state.

### Terminal model

A copolymerization model in which the propagation rate depends only on the
identity of the terminal repeat unit and incoming monomer. It underlies the
classical binary Mayo–Lewis equation.

### Termination

Reaction consuming active radicals and producing dead chains. slimmc
distinguishes at least combination and disproportionation because their effects
on molecule count and distributions differ.

### Terpolymer

A copolymer containing three monomer types. The current copo engine supports
terminal terpolymer kinetics; binary formulas must not be applied without
appropriate generalization.

### Transfer

Reaction in which a growing radical becomes dead and radical character is
transferred to another species. Transfer changes molecular-weight statistics
even when overall monomer consumption is little affected.

## Composition and microstructure

### Free-monomer fraction ($f_i$)

Mole fraction of monomer `i` among all unreacted monomers at a snapshot. It is
the feed/medium composition used by instantaneous copolymerization theory.

### Instantaneous composition ($F_i^{\mathrm{ins}}$)

Theoretical fraction of newly incorporated repeat units at one state, for
example from Mayo–Lewis. It is not the composition accumulated over a finite
time interval.

### Interval composition ($F_i^{\mathrm{interval}}$)

Actual fraction of net repeat units incorporated between two snapshots. It is
undefined for an interval with no appropriate net incorporation denominator.

### Cumulative composition ($F_i^{\mathrm{cum}}$)

Fraction of monomer `i` among all repeat units incorporated up to a snapshot.
It averages over the entire previous composition drift.

### Composition drift

Change of polymer composition during conversion as free-monomer composition
and/or kinetic parameters evolve. It can be studied by comparing interval
composition with monomer composition and instantaneous theory.

### Mayo–Lewis equation

Classical binary terminal-model expression relating free-monomer fractions and
reactivity ratios to instantaneous copolymer composition. Its assumptions must
be checked before comparison with penultimate, terpolymer or strongly
time-varying systems.

### Transition probability

Probability that a chain in a specified terminal/penultimate state adds a
particular monomer next. Predicted probabilities arise from rates and monomer
availability; observed probabilities arise from realized firing counts.

### Triad

Three adjacent repeat units, for example `A|B|A`. Triads contain more sequence
information than composition or dyads and require sufficiently complete stored
sequences.

## Stochastic simulation

### Chemical master equation (CME)

Equation governing the probability distribution over discrete chemical states
for a specified reaction network. slimmc samples trajectories consistent with
the declared discrete channel system rather than solving the CME distribution
directly.

### Channel

One elementary or macromolecular event type with a rate, eligibility rule,
propensity and state mutation. Similar chemistry in different terminal pools
appears as distinct channels.

### Direct SSA / Gillespie algorithm

Exact stochastic simulation algorithm that samples the waiting time and next
channel from current propensities. “Exact” refers to the discrete declared
model, not to perfect representation of all real chemistry.

### Event / firing

One realized execution of a stochastic channel. Cumulative firing counts are
trajectory observations and fluctuate between seeds.

### Firing share

Fraction of realized firings assigned to a channel over a selected interval or
cumulatively. It should not be confused with an instantaneous propensity share.

### Kinetic Monte Carlo (kMC)

Event-driven stochastic simulation terminology used here largely synonymously
with SSA. Time advances between reaction events rather than in fixed steps.

### Propensity

Current stochastic event intensity `a_j(x)` for channel `j`. It combines the
declared rate constant with counts of eligible reactants and the simulation
volume convention.

### Propensity/rate share

One channel's propensity divided by total propensity at a state. It describes
instantaneous competition, whereas firing share describes realized events.

### Random seed

Integer initializing the pseudo-random generator and identifying one
reproducible trajectory. A fixed seed does not represent statistical
uncertainty; independent seeds are needed for an ensemble.

### Stochastic trajectory

One time-ordered sample path of states and events. It can be physically
plausible yet differ visibly from another trajectory of the same model.

### Ensemble

Collection of independent trajectories, usually with different seeds. Ensemble
statistics estimate expected behaviour and stochastic variability.

### Waiting time `tau`

Random simulated time from the current state to the next stochastic event. It
is sampled from an exponential distribution with rate equal to total
propensity.

### Simulation volume `V_kMC`

Physical volume used to convert molar concentrations into integer molecule
counts and bimolecular propensities. Too small a volume creates severe
discretization and noise; larger volume increases computational population.

### Discretization preflight

Before simulation, checks whether declared positive concentrations map to
meaningful integer populations and reports relative rounding/population risks.
It prevents silent disappearance of a positive initial species.

### Time barrier

Exact simulated-time boundary for a scheduled action. The SSA cannot fire an
event using old parameters beyond that boundary.

### `at`, `every`, `when`

`at` schedules an action once at a time, `every` repeats it periodically and
`when` fires once when a state condition first becomes true. One `when` line
may contain any number of predicates joined by `and`; all are evaluated on the
same state. Separate `when` lines provide OR semantics. They are model language
constructs, not Python callbacks.

### `stop` action

One-shot conditional action that requests clean termination after the current
conditional-action scan. It is valid only with `when`, does not write `save` or
`save_chains` implicitly, and records `termination_reason=stop_condition` when
final validation succeeds.

### Parameter state

Complete kinetic/temperature parameter set valid after a parameter-changing
action. Stored parameter states allow time-dependent kinetics to be
reconstructed during analysis.

## Populations, moments and distributions

### Chain population

Selected set of polymer molecules at one chain snapshot, possibly restricted
by live/dead pool, composition, sequence or end group. It is the common source
for CLD, MWD and chain spectra.

### Compressed chain row

One stored record representing many chemically identical chains. Its
`abundance` is essential; counting the row once would bias every statistic.

### Abundance

Number of physical polymer molecules represented by one compressed row. It is
the basic number weight used in chain-population sums.

### DP counts

Exact aggregation of compressed-chain multiplicity by degree of polymerization.
`dp_counts()` returns the discrete source projection `DP -> N` without
normalization.

### Mass counts

Exact aggregation by neutral chain molar mass. `mass_counts()` returns
`M -> N`. In copolymers, the same DP can contribute to several different
masses.

### Chain-length distribution (CLD)

Normalized discrete distribution on degree of polymerization. `weighting` is
`number`, `mass`, or `z`. Mass CLD weights each DP class by the actual polymer
mass carried by chains in that class.

### Exact mass distribution

`mass_distribution()` returns a normalized discrete distribution on the actual
neutral molar masses present in the selected chain population. `weighting` is
`number`, `mass`, or `z`. This is the source-faithful representation for
discrete/oligomeric populations.

### Molar-mass distribution (MWD)

In pyslimmc, `mwd()` specifically denotes the reconstructed mass-weighted
logarithmic density `dW/dlog10(M)`. It is derived from the exact chain-mass
measure by mcPolymer-style linear interpolation and is normalized by numerical
area in `log10(M)`. It is not the source of exact `Mn`, `Mw`, or `Mz`.

### SEC distribution

Continuous apparent mass distribution produced by Gaussian instrumental
broadening in `log10(M)`, applied directly to the exact mass measure. Its
ordinate is `dW_app/dlog10(M)`. It is independent of the MWD interpolation grid.

### Chain-mass spectrum

Exact neutral-mass stick spectrum obtained from distinct chain masses and
abundances. It is not automatically an experimental MS spectrum because it
contains no ionization, isotope, adduct, fragmentation or charge-state model.

### Neutral mass versus `m/z`

Neutral molar mass describes the uncharged molecule; `m/z` is mass-to-charge
ratio of an ion. Converting one to the other requires an explicit ion/adduct and
charge model, which chain spectra do not assume.

### Number-average molar mass `M_n`

Total polymer mass divided by number of polymer molecules. It is particularly
sensitive to the numerous low-mass chains. The default run-level Slimmc
moments describe the complete stored polymer population represented by the
selected population scope; for the normal `all` scope this includes both live
and dead chains. Use the live/dead population views when a comparison requires
one class only.

### Weight-average molar mass `M_w`

Ratio of the second to first molar-mass moment. It weights heavier molecules
more strongly than `M_n`.

### Z-average molar mass `M_z`

Ratio of the third to second molar-mass moment. It is highly sensitive to the
high-mass tail and requires reliable rare-chain statistics.

### Dispersity `Đ`

Ratio `M_w/M_n` (or `DP_w/DP_n` for chain length). It summarizes breadth but
cannot reveal multimodality or tail shape by itself.

### `DP_n`, `DP_w`, `DP_z`

Number-, weight- and z-style averages calculated using DP rather than molar
mass. For copolymers/end-group-rich oligomers, DP- and mass-weighting are not
interchangeable.

### Sticks

Exact discrete representation with a vertical line at every DP/mass value. It
preserves simulation resolution and is the reference before binning/smoothing.

### Bin

Interval of the x coordinate into which discrete contributions are grouped.
Each bin represents a range, not one exact molecular species.

### Binning / histogram

Grouping exact sticks into bins to reduce noise/detail. Apparent modes and
heights depend on bin width and placement, so bin settings must be reported.

### Bin width

Size of a histogram interval in linear or log10 coordinate. A logarithmic width
of `0.02` is 0.02 decade, not 0.02 g/mol.

### Kernel density estimate (KDE)

Smooth density formed by placing a kernel at every weighted observation. KDE
is a representation of the simulated population and should not be confused
with experimental instrument broadening.

### Bandwidth / `sigma`

Width of the smoothing kernel in the construction coordinate. Too small gives
a noisy curve; too large merges chemically meaningful modes.

### Gaussian broadening

Convolution-like smoothing with Gaussian functions. In the Gaussian
method it is applied to a histogram, while KDE acts directly on weighted exact
observations.

### Grid step

Distance between numerical points at which a smooth KDE curve is evaluated.
It changes curve resolution but not the chosen physical smoothing bandwidth.

### Linear coordinate

Axis/construction coordinate in DP or molar mass itself. Equal spacing
corresponds to equal absolute increments.

### Logarithmic coordinate / `log10(M)`

Coordinate useful for distributions spanning orders of magnitude. Densities per
unit mass and per unit log-mass differ by a Jacobian and are not interchangeable.

### Jacobian

Factor required when transforming a density between coordinates. For
`u=log10(M)`, `dM/du=ln(10)M`.

### Amount

Unnormalized contribution per stick/bin under the selected basis. It retains
absolute stochastic population scale.

### Fraction

Dimensionless normalized contribution. For a discrete distribution its
contributions sum to one under per-series normalization.

### Density

Contribution per unit coordinate. Its area, not usually its point sum, carries
the normalization meaning.

### Absolute normalization

No division by the total of each series. It preserves changes in chain number
or polymer mass between populations/runs.

### Per-series normalization

Every curve is independently normalized. It compares distribution shape while
discarding absolute yield/population differences.

### Base-peak normalization

Scaling in which the highest spectrum line equals 100. It is a presentation
convention and removes absolute abundance information.

### Reference normalization

Scaling multiple series relative to one named reference series. It is useful
for comparing contributions but must not be mistaken for independent unit-area
normalization.

### Overlay

Several curves drawn on the same axes. Overlay supports shape comparison but
can hide smaller curves or make overlapping series look additive when they are
not.

### Stacked distribution

Series drawn cumulatively above one another. It is physically meaningful as a
total only when the series are non-overlapping parts of the same population.

## Output and analysis objects

### Snapshot

Saved state at a specific `snapshot_id`, event count and simulated time. Some
snapshots contain only state; full snapshots additionally contain chains and
other aligned outputs.

### Final versus last

`final` is the authoritative termination state; `last` is the last explicitly
saved member of a series. They may differ, so endpoint analysis should usually
request `final` explicitly.

### Partial output

Run directory missing some optional or not-yet-written files. `pyslimmc` keeps
usable parts accessible while reporting unavailable analyses explicitly.

### Output contract

Documented agreement between engines and readers: filenames/logical keys,
columns, units, identities and completion metadata. It prevents analysis code
from guessing file meaning.

### Schema

Machine-readable description/version of data organization. A schema identifies
how to interpret output; it is different from a chemical mechanism.

### Metadata

Data describing the run, model, producer and output rather than a kinetic
measurement itself. Good metadata makes a trajectory reproducible and
searchable.

### Mass audit

Independent check that chain masses can be interpreted completely and/or agree
with declared chemistry. Missing end groups can make structural MWD unsafe even
when repeat-unit-only MWD is available.

### Invariant

Condition that must always hold for a valid engine state/output, such as
non-negative populations or DP equalling summed composition. An invariant is a
strict correctness requirement, not a statistical expectation.

### Balance ledger

Independent engine counter used to verify that free plus incorporated monomer
matches externally added/removed inventory. Independence gives the validation
check real power.

### Capability

One analysis described by separate flags: implemented in software, applicable
to the chemical model and supported by available data. These must not be
collapsed into one vague available/unavailable flag.

### Validation report

Structured result describing strict validity, completeness, missing/invalid
outputs and warnings. A completed process status does not override failed
scientific/data invariants.

### Diagnostic

Information useful for assessing a run without necessarily making it invalid,
for example small stochastic populations or a statistical SSA sanity check.

## Python analysis terminology

The public Python API intentionally uses ordinary Python/NumPy/Matplotlib
terminology. Generic definitions of concepts such as *method*, *cache*,
*NumPy array*, *DPI*, or Matplotlib backends are outside this glossary.
See [`PYSLIMMC.md`](PYSLIMMC.md) for task-oriented usage and
[`reference/PYSLIMMC_API.md`](reference/PYSLIMMC_API.md) for the exact API.

## See also

- [`THEORY.md`](THEORY.md) for equations and assumptions;
- [`pyslimmc API`](reference/PYSLIMMC_API.md) for the exact analysis contract;
- [`homo reference`](reference/HOMO.md) and
  [`copo reference`](reference/COPO.md) for engine-specific model syntax;
- [`architecture`](development/ARCHITECTURE.md), [`Storage`](reference/STORAGE.md), and
  [`testing`](development/TESTING.md) for implementation contracts.


## Storage, validation, and API

### Slimmc Storage

The canonical directory-based results format written by both engines. Tables are directories, numeric columns are `.npy` files, and metadata and dictionaries are JSON or JSONL.

### Storage format version

Version of the on-disk Slimmc Storage contract. It is independent of the engine and `pyslimmc` versions.

### Engine version

Version of the executable that generated a run.

### pyslimmc version

Version of the Python analysis package reading and analysing Slimmc Storage.

### Validation

Automatic checks of Storage structure and scientific invariants, such as equal column lengths, valid foreign keys, count-unit conversions, composition reconstructing DP, and moments reconstructing from chains.

### Audit

A targeted assessment of whether a selected interpretation is complete and safe. An audit is narrower than validation and may depend on analysis options.

### Mass audit

An independent `pyslimmc` check that required repeat-unit and end-group mass definitions are available for a selected molecular-mass model.

### Completion marker

`RESULTS_COMPLETE`, created only after a completed run has been finalized and validated.

### Partial run

A running, failed, or interrupted run without the completion marker. It may be opened only with explicit `allow_incomplete=True`.

### Snapshot reason

Dictionary value describing why a snapshot was written, for example initial, scheduled, action, or final.

### Fire share

Fraction of realized SSA/KMC events assigned to a channel over an interval or cumulatively.

### Rate share / propensity share

A channel's propensity divided by total propensity at a snapshot. It describes instantaneous competition and is not the same as a fire share.

### Capability

A structured statement separating implementation, chemical applicability, data availability, and final availability of an analysis.

### Read-only NumPy leaf

A numeric endpoint such as `run.t` or `run.conv["A"]` that returns a NumPy array protected against in-place modification.

### Composition sequence mode

`sequence_mode=composition` stores full per-chain composition and engine-aggregated microstructure counters, but not complete symbol order for every chain.

### Full sequence mode

`sequence_mode=full` stores complete linear sequences and supports independent sequence-derived consistency checks.
