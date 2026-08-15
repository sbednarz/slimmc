# Theory of slimmc and pyslimmc

This document states the mathematical and chemical conventions used by the
`slimmc` homo and copo engines and by `pyslimmc`. It is a theory reference, not
a command-language manual. Model syntax is documented in
[`reference/HOMO.md`](reference/HOMO.md) and [`reference/COPO.md`](reference/COPO.md); the
analysis interface is documented in [`PYSLIMMC_API.md`](reference/PYSLIMMC_API.md).
Project terminology is collected in [`GLOSSARY.md`](GLOSSARY.md).

## 1. Scope and assumptions

The engines simulate a spatially homogeneous, well-mixed stochastic reactor.
The state is discrete: small species are molecule counts and polymers are
populations of molecular objects. Temperature and kinetic parameters may
change at scheduled or conditional boundaries, but spatial diffusion and a
continuous viscosity field are not solved.

The current molecular model assumes linear chains. Branching, cyclization and
network topology are outside its state representation. The copo engine stores
sequence information; the homo engine exploits the fact that a single-monomer
linear chain can be represented compactly by its degree of polymerization and
end-group state.

## 2. Notation and units

| Symbol | Meaning |
|---|---|
| `V` | stochastic simulation volume in litres |
| `N_A` | Avogadro constant |
| `n_i` | integer molecule count of species `i` |
| `c_i` | molar concentration of species `i` |
| `k_j` | rate constant of channel `j` |
| `a_j(x)` | propensity of channel `j` in state `x` |
| `a_0` | sum of all propensities |
| `t` | simulated time in seconds |
| `DP` | number of repeat units in one chain |
| `M` | molar mass of one chain |
| `N_q` | multiplicity of a compressed chain row `q` |

Concentration and molecule count are related by

$$
n_i=c_iVN_A,
\qquad
c_i=\frac{n_i}{VN_A}.
$$

The conversion to an integer population necessarily rounds the initial count.
If a positive declared initial concentration would round to zero molecules,
the engine's discretization preflight rejects the model rather than silently
removing the species. For accepted models the written output reports the
realized discrete state, not an imaginary fractional molecule population.

## 3. Direct stochastic simulation algorithm

### 3.1 Chemical master-equation view

Let `x` be the current discrete state. A reaction channel `j` has a state
change vector `nu_j` and a propensity `a_j(x)` such that

$$
a_j(x)\,dt
$$

is the probability, to first order in `dt`, that channel `j` fires in the next
interval `dt`. The probability distribution of the state is governed by the
chemical master equation. The engines generate statistically exact sample
paths of the declared discrete channel system using the direct SSA.

### 3.2 Propensity conventions

For a first-order channel

$$
A\rightarrow\cdots,
\qquad
a=k n_A.
$$

For two different reactants

$$
A+B\rightarrow\cdots,
\qquad
a=\frac{k n_A n_B}{N_A V}.
$$

For two reactants drawn from the same population

$$
2A\rightarrow\cdots,
\qquad
a=\frac{k n_A(n_A-1)}{N_A V}.
$$

In slimmc's declared-rate convention the same-species symmetry factor is
carried by the rate constant; the engine does not insert an additional
`1/2`. Therefore a declared elementary channel `2 A -> B k` reproduces the
macroscopic disappearance law

$$
-\frac{d[A]}{dt}=2k[A]^2.
$$

For bimolecular radical termination this is the usual polymer-kinetics form
`R_t = 2 k_t [P^\bullet]^2`. If a literature source instead defines its
reported constant through `-d[P^\bullet]/dt = k_t[P^\bullet]^2`, its numeric
constant differs by a factor of two from the convention above. Always compare
the defining rate equation, not only the symbol `k_t`.

Polymer channels use the same counting principle. For example, propagation
from a pool containing `n_P` eligible live chains and `n_M` monomer molecules
has propensity

$$
a_p=\frac{k_p n_P n_M}{N_A V}.
$$

Eligibility is part of the channel definition. A depropagation channel, for
example, counts only chains whose length and terminal state permit removal of
the terminal repeat unit.

### 3.3 Waiting time and channel selection

For

$$
a_0(x)=\sum_j a_j(x),
$$

two independent uniform random values `u_1,u_2 in (0,1)` give

$$
\tau=-\frac{\ln u_1}{a_0},
$$

and the selected channel `mu` satisfies

$$
\sum_{j<\mu}a_j < u_2a_0
\leq \sum_{j\leq\mu}a_j.
$$

If `a_0=0`, no stochastic event is possible. Termination metadata distinguishes
this condition from reaching `t_end`, reaching `max_steps`, a memory policy or
a validation failure.

### 3.4 Determinism and statistical reproducibility

A fixed model, engine version and seed should reproduce the same trajectory on
the supported implementation. A seed does not remove stochastic uncertainty:
it merely identifies one sample path. Estimating an expectation or confidence
interval requires independent runs with different seeds.

For a statistic `Y` obtained from `R` independent trajectories, the sample
mean and its estimated standard error are

$$
\bar Y=\frac{1}{R}\sum_{r=1}^{R}Y_r,
\qquad
SE(\bar Y)=\frac{s_Y}{\sqrt R}.
$$

## 4. Simulated-time boundaries and actions

An `at` or `every` action that changes parameters is a time barrier. If the
next stochastic waiting time crosses the boundary, the engine advances to the
boundary, executes the action and recomputes the affected propensities before
sampling the next event. The event is not fired using parameters from the
wrong side of the boundary.

Conditional `when` actions are evaluated from simulated state. They may be
used for phenomenological relations such as a stepwise termination coefficient
dependent on conversion. Such a rule is imposed by the model; it is not an
emergent diffusion calculation.

Concentration dosing changes the discrete count consistently with `V`. Rate
and temperature histories are written so the time-dependent process can be
reconstructed during analysis.

## 5. Homopolymerization mechanisms

### 5.1 Initiation

A conventional initiator decomposition may be represented as an elementary
reaction followed by macromolecular initiation:

$$
I\xrightarrow{k_d}2R,
\qquad
R+M\xrightarrow{k_i}P_1^\bullet.
$$

An efficiency may be represented by the declared elementary-reaction
stoichiometry/efficiency convention documented in the language manual. The
propensity is calculated from the declared `k` and reactant populations; it is
**not** multiplied by `f`. Once that channel fires, the left-hand-side
reactants are consumed and products are formed with probability `f`. Thus `f`
represents the productive fraction of reaction attempts (for example a cage
efficiency), rather than replacing the kinetic constant by `f*k`. The engine
does not silently assume an initiator efficiency that is absent from the model.

### 5.2 Propagation and depropagation

Propagation is

$$
P_n^\bullet+M\xrightarrow{k_p}P_{n+1}^\bullet.
$$

Terminal depropagation is the reverse length-changing event

$$
P_n^\bullet\xrightarrow{k_{dp}}P_{n-1}^\bullet+M,
\qquad n\geq2.
$$

The minimum-DP rule prevents creation of a zero-length polymer object. A model
containing both directions can approach a kinetic balance, but detailed
balance is obtained only when the declared constants and state definition are
thermodynamically consistent.

### 5.3 Termination

Combination produces one dead molecule:

$$
P_n^\bullet+P_m^\bullet
\xrightarrow{k_{tc}}D_{n+m}.
$$

Disproportionation produces two dead molecules:

$$
P_n^\bullet+P_m^\bullet
\xrightarrow{k_{td}}D_n+D_m.
$$

These channels have different effects on molecule counts and the molecular
weight distribution even when they consume radicals at the same total rate.

### 5.4 Transfer and reinitiation

A generic transfer step has the conceptual form

$$
P_n^\bullet+X\xrightarrow{k_{tr}}D_n+R_X.
$$

If `R_X` initiates another chain, reinitiation is represented by a separate
initiation channel. Keeping transfer and reinitiation separate permits delayed
or inefficient reinitiation and makes the corresponding firing counts visible.

Transfer to monomer similarly terminates one growing chain and creates a
monomer-derived radical according to the declared model. End-group identities
affect structural mass only when the `with_end_groups` mass model has complete
end-group definitions.

## 6. Copolymerization and terpolymerization

### 6.1 Terminal model

For monomers `A` and `B`, a terminal model distinguishes the four propagation
constants

$$
k_{AA},\quad k_{AB},\quad k_{BA},\quad k_{BB},
$$

where the first index is the terminal repeat unit and the second is the added
monomer. A chain ending in `A`, for example, has competing propensities

$$
a_{AA}=\frac{k_{AA}n_{P_A}n_A}{N_AV},
\qquad
a_{AB}=\frac{k_{AB}n_{P_A}n_B}{N_AV}.
$$

The binary reactivity ratios are

$$
r_A=\frac{k_{AA}}{k_{AB}},
\qquad
r_B=\frac{k_{BB}}{k_{BA}}.
$$

For free-monomer mole fractions `f_A` and `f_B=1-f_A`, the terminal
Mayo--Lewis instantaneous fraction of incorporated `A` is

$$
F_A^{ins}=
\frac{r_A f_A^2+f_Af_B}
     {r_A f_A^2+2f_Af_B+r_Bf_B^2},
\qquad
F_B^{ins}=1-F_A^{ins}.
$$

This is an instantaneous terminal-model prediction. It is not identical to
the composition accumulated over a finite interval or over the complete run.

### 6.2 Four composition concepts

The analysis deliberately separates:

- `f_i`: fraction of free monomer `i` in the current monomer mixture;
- `F_i^ins`: theoretical instantaneous incorporation fraction;
- `F_i^interval`: actual fraction incorporated between two snapshots;
- `F_i^cum`: cumulative fraction in all repeat units incorporated so far.

With free-monomer amounts `m_i`,

$$
f_i=\frac{m_i}{\sum_j m_j}.
$$

If `Delta N_i` units of monomer `i` were incorporated over an interval,

$$
F_i^{interval}=\frac{\Delta N_i}{\sum_j\Delta N_j}.
$$

For total polymerized units `N_i^poly`,

$$
F_i^{cum}=\frac{N_i^{poly}}{\sum_jN_j^{poly}}.
$$

Quantities with a zero denominator are chemically undefined. `pyslimmc` does
not silently replace every undefined composition by zero.

### 6.3 Explicit penultimate kinetics

The penultimate model distinguishes both the terminal and penultimate units.
A pool `P_AB` represents a chain whose penultimate unit is `A` and terminal
unit is `B`. Adding `C` changes the remembered pair from `AB` to `BC`:

$$
P_{AB}^\bullet+C\xrightarrow{k_{AB,C}}P_{BC}^\bullet.
$$

The rate tensor therefore contains more information than terminal reactivity
ratios. Terminal and penultimate predictions must not be compared as if they
were the same model. `pyslimmc` reports capabilities and raises a domain error
when the required parameter family is incomplete.

### 6.4 Terpolymers

For three monomers, the terminal model has up to nine terminal/addition
combinations. The event competition follows the same propensity law; there is
no binary Mayo--Lewis reduction that represents all three components without
additional assumptions. Composition is consequently handled as named vectors
and ratio arrays rather than hard-coded `A/B` scalars.

### 6.5 Sequence memory and truncation

Sequence-aware chains support composition, terminal state, dyads, triads,
blocks and oligomer analysis. A stored dead-chain sequence may be limited by
`sequence_mode` (`composition` or `full`). The total DP and composition can remain known even
when the literal sequence is incomplete. Analyses requiring the omitted order
must then report incomplete sequence data instead of fabricating motifs.

## 7. Chain mass models

For repeat-unit masses `m_i` and counts `n_i` in a chain, the repeat-unit-only
mass is

$$
M_{repeat}=\sum_i n_i m_i.
$$

For the structural end-group model,

$$
M_{struct}=M_{repeat}+m_{left}+m_{right}+\cdots,
$$

where the additional terms are determined by the reaction history and declared
end groups. Missing required end-group masses make the structural mass
incomplete. The mass audit exposes the missing names and affected rows.

Mass conservation and structural mass completeness are related but distinct:
a reaction may conserve the declared molecular bookkeeping while the absolute
mass remains unavailable because an end group has no assigned molar mass.

## 8. Compressed populations and weighted statistics

An output row may represent `N_q` identical chains. Every population statistic
must use this multiplicity. Treating each row as one molecule biases moments
toward rare unique structures.

For a chain property `x_q`, the number-weighted mean is

$$
\langle x\rangle_n=
\frac{\sum_qN_qx_q}{\sum_qN_q}.
$$

For degree of polymerization,

$$
DP_n=\frac{\sum_qN_qDP_q}{\sum_qN_q},
\qquad
DP_w=\frac{\sum_qN_qDP_q^2}{\sum_qN_qDP_q}.
$$

For molar mass,

$$
M_n=\frac{\sum_qN_qM_q}{\sum_qN_q},
$$

$$
M_w=\frac{\sum_qN_qM_q^2}{\sum_qN_qM_q},
$$

$$
M_z=\frac{\sum_qN_qM_q^3}{\sum_qN_qM_q^2},
\qquad
Đ=\frac{M_w}{M_n}.
$$

Empty populations and zero mass denominators do not have meaningful ordinary
moments. Callers must handle the documented unavailable/undefined result.

## 9. Exact chain projections, CLD, MWD, and SEC

Let compressed chain record `r` have degree of polymerization `D_r`, neutral
molar mass `M_r`, and multiplicity `N_r`.

### 9.1 Exact count projections

The exact DP-count projection is

```text
N_D(D) = sum_{r: D_r=D} N_r,
```

and the exact mass-count projection is

```text
N_M(M) = sum_{r: M_r=M} N_r.
```

These are returned by `dp_counts()` and `mass_counts()` and are not normalized.
Because copolymers or different end groups can give different masses at the
same DP, MWD must be projected independently onto actual `M_r`; it is not a
transformation of aggregated CLD.

### 9.2 Discrete CLD forms

For unique DP classes `D_i` with counts `N_i`, number CLD is

```text
p_i = N_i / sum_j N_j.
```

Mass-weighted CLD uses the actual mass carried by each DP class,

```text
A_i = sum_{r: D_r=D_i} N_r M_r,
w_i = A_i / sum_j A_j.
```

The z-weighted CLD is

```text
z_i = D_i^2 N_i / sum_j D_j^2 N_j.
```

The logarithmic CLD uses the same exact mass fractions as the mass form, but
places them on `log10(D_i)` support.

### 9.3 Discrete MWD forms

For unique neutral masses `M_i` and exact counts `N_i`, number MWD is

```text
p_i = N_i / sum_j N_j,
```

mass MWD is

```text
w_i = M_i N_i / sum_j M_j N_j,
```

and z MWD is

```text
z_i = M_i^2 N_i / sum_j M_j^2 N_j.
```

The logarithmic MWD uses the same exact mass fractions `w_i` on support
`u_i = log10(M_i)`.

### 9.4 Continuous logarithmic density

If a continuous mass density is written as `w(M) = dW/dM` and
`u = log10(M)`, conservation of measure gives

```text
w(M) dM = g(u) du,
g(u) = M ln(10) w(M).
```

The Jacobian belongs to a change of variables for continuous densities. Exact
discrete atom weights are unchanged when their support is re-expressed in
`log10(M)`.

### 9.5 SEC instrumental broadening

pyslimmc models Buback-style SEC broadening as a Gaussian response in
`u = log10(M)` applied directly to exact mass fractions:

```text
S(u) = sum_i w_i /(sigma sqrt(2 pi))
       * exp(-(u-u_i)^2/(2 sigma^2)).
```

Thus `S(u)` is a continuous apparent `dW_app/dlog10(M)` distribution and
integrates analytically to unity. Under linear calibration
`log10(M) = a - b v`, the public parameter `sigma_log10M` corresponds to
`b sigma_v`.

This is an instrumental response model, not generic smoothing. Histogram, KDE,
and arbitrary Gaussian smoothing are not definitions of CLD/MWD in pyslimmc.

### 9.6 Exact moments

Mass moments are evaluated from exact source sums,

```text
Mn = sum N_i M_i / sum N_i,
Mw = sum N_i M_i^2 / sum N_i M_i,
Mz = sum N_i M_i^3 / sum N_i M_i^2,
```

with analogous DP moments. They are invariant with respect to distribution
form and SEC broadening.

## 10. Conversion, rates and firing fractions

For monomer `i`, with initial amount `n_{i,0}` and current free amount `n_i`,

$$
X_i=1-\frac{n_i}{n_{i,0}}.
$$

A total conversion weighted by all initial monomer units is

$$
X_{total}=
\frac{\sum_i(n_{i,0}-n_i)}{\sum_i n_{i,0}},
$$

not the arithmetic mean of the individual `X_i` values.

For cumulative firing count `C_j(t)`, the count over interval `(t_1,t_2]` is

$$
\Delta C_j=C_j(t_2)-C_j(t_1).
$$

The firing share is

$$
s_j=\frac{\Delta C_j}{\sum_k\Delta C_k}.
$$

Firing shares describe realized stochastic events. Propensity or rate shares
describe instantaneous competition. They approach one another statistically
over suitable ensembles/intervals but are not identical observations.

## 11. Microstructure statistics

For fully observed sequences, the fraction of motif `alpha` is its weighted
count divided by the weighted total of motifs of the same order. Dyads count
adjacent pairs; triads count overlapping triples. A sequence of length `DP`
contains `DP-1` dyad positions and `DP-2` triad positions when `DP` is large
enough.

A block is a maximal consecutive run of one monomer identity. Number-average
block length for component `A` is

$$
L_{n,A}=\frac{\text{number of A repeat units}}
              {\text{number of A blocks}}.
$$

Statistics must be weighted by compressed-chain multiplicity. Sequence
truncation changes the available motif denominator and is reported through
capability/completeness metadata.

## 12. Numerical interpretation and limitations

- A single SSA path is not a smooth kinetic law.
- Smaller simulation volumes increase discreteness and typically increase
  relative noise.
- Snapshot spacing limits time resolution of interval-derived quantities.
- A final engine state and the last user-requested snapshot are distinct
  anchors; `pyslimmc` preserves that distinction.
- Smoothing parameters alter a plotted distribution and must be reported.
- Conversion-triggered parameter changes are model assumptions, not automatic
  physical laws.
- Stored sequences may be incomplete even when aggregate composition is exact.
- Direct comparison with SEC/GPC requires accounting for calibration,
  detector response and instrumental broadening; an arbitrary smoothed curve is not itself an
  SEC forward model.

## 13. Validation principles

Theory claims are protected at several levels:

- Nim unit tests check parser, channel and event bookkeeping;
- deterministic engine validation cases check stoichiometry, boundaries,
  moments and output invariants;
- `pyslimmc` fixtures check compressed weights, normalization, transformations,
  composition and error behaviour;
- literature cases, when present, compare scientific results separately from
  the fast default test suite.

The testing layout, Storage contract, implementation boundaries, and source
map are documented in [`development/ARCHITECTURE.md`](development/ARCHITECTURE.md).

## See also

- [`CONCEPTS.md`](CONCEPTS.md) — Conceptual guide
- [`reference/HOMO.md`](reference/HOMO.md) — Homo engine reference
- [`reference/COPO.md`](reference/COPO.md) — Copo engine reference
