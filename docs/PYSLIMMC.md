# pyslimmc user guide

`pyslimmc` is the read-only Python analysis layer for slimmc Storage. It is a
regular Python library, not a GUI. It can be used from Python scripts or the
standard interpreter, IPython, Jupyter Notebook/JupyterLab, Marimo, VS Code
(including `.py` files and `# %%` cells), and other environments that run
Python.

Use this page when you know **what scientific result you want** and need the
shortest path to the relevant object or analysis. For a user-oriented map of
the objects and their relationships, see
[`reference/PYSLIMMC_API_TREE.md`](reference/PYSLIMMC_API_TREE.md). For the
exhaustive public API, see
[`reference/PYSLIMMC_API.md`](reference/PYSLIMMC_API.md); exact generated
signatures are in
[`reference/PYSLIMMC_SIGNATURES.md`](reference/PYSLIMMC_SIGNATURES.md).

The core runtime dependency is NumPy. Matplotlib is optional and is required
only for plotting and graphical reports (`pyslimmc[plot]`).

```python
import pyslimmc as sl

run = sl.open("results/run_000001")
```

## I want to...

| Task | Start with |
|---|---|
| inspect a run | `run.info()` / `run.summary()` |
| plot conversion | `run.t`, `run.conv[...]`, `run.conv.total` |
| get concentrations | `run.conc[...]` |
| get physical moles | `run.moles[...]` |
| get final `M_n`, `M_w`, `M_z`, `Đ` | `run.final.mn`, `.mw`, `.mz`, `.dispersity` |
| get moment series | `run.mn`, `run.mw`, `run.mz`, `run.dispersity` |
| build an MWD | `run.mwd()` |
| build a CLD | `run.cld()` |
| inspect exact DP/mass counts | `run.dp_counts()` / `run.mass_counts()` |
| compare live/dead chains | `run.mwd_series(series=("live", "dead"))` |
| inspect copolymer composition | `run.f`, `run.F` |
| calculate dyads/triads/blocks | `run.microstructure` |
| inspect channel competition | `run.firings.fire_shares()`; `.rate_shares()` currently requires copo `channel_propensities` |
| inspect feed/balances | `run.feeds`, `run.balance`, `run.volume` |
| validate a run | `run.validate(strict=True)` |
| audit chain masses | `run.mass_audit()` |
| verify provenance/reproducibility | `run.reproducibility.verify()` |
| find several runs | `sl.scan(...)` |
| select runs by ID pattern | `runs.match(...)` |
| make a parameter grid | `runs.sweep(...)` |
| access raw Storage | `run.raw` |

## Run identity and versions

```python
sl.__version__                 # pyslimmc
run.engine_version             # producer engine
run.storage_format_version     # slimmc Storage
run.run_id
run.status
```

The component versions are independent. Reproducibility metadata are exposed
without mutating the stored run:

```python
run.reproducibility.input_hash
run.reproducibility.model_hash
run.reproducibility.binary_hash
run.reproducibility.storage_hash
run.reproducibility.git_commit
run.reproducibility.git_dirty
run.reproducibility.verify()
```

See [`SIMULATION_RESULTS.md`](SIMULATION_RESULTS.md) for the distinction between
scientific results, run history, and provenance.

## Snapshot navigation

```python
run.first
run.last
run.final
run.snapshots_with_chains
run.first_with_chains
run.last_with_chains

run.at_time(1.0, method="before")
run.at_event(1000, method="before")
run.at_conversion(0.5, monomer="A", method="nearest")
run.at_temperature(333.15)  # tuple: several snapshots may be equally near
```

Run-wide numeric leaves are read-only NumPy arrays. Snapshot leaves are
scalars.

`final` means a successfully finalized terminal snapshot; `last` means the
last snapshot physically present. Do not substitute one for the other for an
interrupted/failed run.

## State, amounts, and conversion

```python
run.t
run.event
run.sid

run.count["A"]
run.conc["A"]
run.conv["A"]
run.conv.total

snap = run.final
snap.count["A"]
snap.conc["A"]
snap.conv["A"]
```

`run.conv.total` is weighted by initial monomer counts, not the arithmetic mean
of component conversions.

Physical amounts require `param init_volume`:

```python
run.moles["A"]       # physical reactor amount
run.state.moles["A"] # technical amount represented by KMC counts
```

These are deliberately different concepts.

## Semibatch volume, feeds, and balances

```python
run.kmc_volume       # stochastic representation volume, L
run.volume           # physical reactor volume, L

run.c0["A"]
run.count0["A"]
run.moles0["A"]
```

Feed definitions and realized feed events are available directly:

```python
run.feeds.names
run.feeds["F"].concentration
run.feeds["F"].fraction
run.feeds["F"].events
run.feeds["F"].volume_cum
run.feeds["F"].moles_cum
```

Physical balance terms (mol):

```python
run.balance.initial["A"]
run.balance.dosed["A"]
run.balance.total["A"]
run.balance.free["A"]
run.balance.consumed["A"]
run.balance.incorporated["A"]  # monomers
```

`run.balance` and physical `run.moles` require `init_volume`. If `set_c` has
invalidated the physical balance for a species, pyslimmc raises
`AnalysisNotApplicableError` instead of returning a fabricated value.

Chain counts are available at ordinary `save` snapshots:

```python
run.chain_count.live
run.chain_count.dead
run.chain_count.total
```

Chain-resolved concentrations and distributions require `save_chains`.

## Copolymer composition

```python
run.f["A"]
run.f0["A"]
run.F.ins["A"]
run.F.int["A"]
run.F.cum["A"]
```

Interpretation:

- `f` — current free-monomer mole fraction;
- `f0` — initial free-monomer mole fraction;
- `F.ins` — theoretical instantaneous incorporation fraction;
- `F.int` — actual interval polymer composition;
- `F.cum` — cumulative polymer composition.

For higher-level comparisons:

```python
cp = run.copolymerization
cp.reactivity_ratios()
cp.mayo_lewis()
cp.compare_mayo_lewis()
cp.terminal_diagnostics()
cp.penultimate_parameters()
cp.penultimate_composition()
cp.compare_penultimate()
cp.penultimate_diagnostics()
```

Use capability checks when writing code intended for several model classes.

## Chain populations

```python
chains = run.final.chains

chains.all
chains.live
chains.dead
chains.pool("terminal_A")
chains.origin("term_c")
chains.where(dp_min=10, dp_max=100)

chains.dp
chains.molar_mass
chains.abundance
chains.composition.counts["A"]
chains.composition.fractions["A"]
```

One row is one compressed structural record. `chains.n_records` is the number
of stored compressed rows; `chains.n_chains` is the represented physical chain
population. pyslimmc weights compressed rows by their multiplicity.

Alternative chain masses can be reconstructed where supported:

```python
chains.masses(mass_model="repeat_units")
chains.masses(mass_model="with_end_groups")
```

## Exact moments

Simple run-level moments remain available:

```python
run.dpn
run.dpw
run.mn
run.mw
run.mz
run.dispersity       # Mw / Mn
run.dp_dispersity    # DPw / DPn
```

For an explicitly selected population and mass model:

```python
m = run.moments(
    snapshot="final",
    population="dead",
    mass_model="repeat_units",
)

m.dpn
m.dpw
m.dpz
m.mn
m.mw
m.mz
m.dp_dispersity
m.mass_dispersity
```

The fully general path is through a selected `ChainPopulation`:

```python
m = run.final.chains.dead.moments(mass_model="repeat_units")
```

Moments are calculated from exact source-chain information or stored aggregate
moments. They are never recomputed from plotted CLD/MWD/SEC curves. `dpz` is
`NaN` when only an older aggregate moments table is available and the required
third DP moment was not stored.

## Exact DP and mass counts

Exact source projections are separate from normalized distributions:

```python
dp = run.dp_counts(pool="dead")
print(dp.dp)
print(dp.count)

mass = run.mass_counts(pool="dead")
print(mass.mass)
print(mass.count)
```

`DPCounts` represents exact `DP -> chain count`; `MassCounts` represents exact
neutral molar mass `M -> chain count`. They are unnormalized source-of-truth
projections and expose semantic properties (`dp`/`mass`, `count`) rather than
generic `x`/`y` aliases.

For copolymers or end-group-aware masses, several different masses may occur at
the same DP. Therefore general MWD is built independently from actual chain
masses, not by transforming an already aggregated CLD.

## CLD and MWD forms

The distribution API uses one literature-oriented keyword:

```python
form="number" | "mass" | "z" | "log"
```

Defaults are:

```python
cld = run.cld()   # form="number"
mwd = run.mwd()   # form="log"
```

All CLD/MWD results are exact **discrete** normalized distributions. There is no
histogram, KDE, generic Gaussian smoothing, midpoint construction, or generic
`interpolate=` step in the core API.

For CLD:

```python
run.cld(form="number")
run.cld(form="mass")
run.cld(form="z")
run.cld(form="log")
```

For MWD:

```python
run.mwd(form="number")
run.mwd(form="mass")
run.mwd(form="z")
run.mwd(form="log")
```

`form="mass"` CLD accumulates the **actual chain mass** carried by each DP
class. It is not generally equal to weighting by `DP * count` for copolymers.

For exact logarithmic forms, the support is transformed but atomic weights are
unchanged:

```python
m_mass = run.mwd(form="mass")
m_log = run.mwd(form="log")

# same exact mass fractions
# m_log.x == log10(m_mass.mass)
# m_log.y == m_mass.y
```

Thus `mwd(form="log")` is a discrete mass-weighted distribution on
`log10(M)` support. It should not be confused with a finite continuous density.
Changing a Matplotlib axis with `plt.xscale("log")` is only a display choice and
does not change the mathematical distribution.

Distribution moments remain exact source moments for every form:

```python
mwd = run.mwd(form="log")
mwd.mn
mwd.mw
mwd.mz
mwd.dispersity

cld = run.cld(form="mass")
cld.dpn
cld.dpw
cld.dpz
cld.dispersity
```

CLD intentionally does not expose misleading `mn`/`mw`/`mz` aliases.

## Multi-series distributions

Multi-series comparison is a separate composition layer:

```python
g = run.mwd_series(
    series=("live", "dead"),
    form="log",
    normalization="per_series",
)
```

Accepted normalizations are only:

```text
per_series
combined
```

`per_series` normalizes every series independently to 1 and is intended for
shape comparison. `combined` preserves relative physical contributions using
exact source totals appropriate to the selected form and requires pairwise-
disjoint populations; for example `live + dead` is valid, while `all + dead`
is rejected.

Each series keeps its own exact support. pyslimmc does not construct a common
grid, interpolate between supports, or insert synthetic bins.

## SEC broadening

Experimental SEC broadening is intentionally separate from MWD:

```python
sec = run.sec(
    pool="dead",
    sigma_log10M=0.05,
)
```

The model applies a normalized Gaussian instrumental response in
`u = log10(M)` directly to exact source mass fractions:

```text
S(u) = sum_i w_i G_sigma(u - log10(M_i))
```

`sec.y` is therefore a genuine continuous apparent density
`dW_app / dlog10(M)`, unlike the atomic weights returned by discrete
`mwd(form="log")`.

`sigma_log10M` is required because it is an instrumental/model parameter; no
silent default is assumed. `step_log10M=` may be supplied to control numerical
sampling of the output curve. Exact `Mn`, `Mw`, `Mz`, and dispersity attached
to the SEC result still come from the source chain population and do not depend
on broadening or grid spacing.

## Microstructure

Aggregated analyses:

```python
run.microstructure.dyads()
run.microstructure.triads()
run.microstructure.run_lengths("A")
run.microstructure.blockiness()
```

Composition filters/maps work with composition counters in both sequence modes.
Full-sequence transitions, motifs, n-grams and position-resolved analyses
require `sequence_mode full`. See the exact API for the complete method list.

## Channels, firings, and kinetics

```python
run.channels.event_count["prop_AA"]
run.firings.fire_shares()
run.firings.rate_shares()
run.temp
run.k["kp_aa"]
run.actions
```

Fire shares summarize realized events. Rate shares summarize propensity
competition at snapshots.

## Validation, mass audit, and reproducibility

These answer different questions:

```python
validation = run.validate(strict=True)
audit = run.mass_audit()
audit.raise_if_failed()
repro = run.reproducibility.verify()
```

- validation checks Storage/run invariants;
- mass audit checks a chain-mass interpretation;
- reproducibility verification checks hashes/provenance;
- none of these by itself proves stochastic or `V_kMC` convergence.

A practical pre-publication workflow is in [`COOKBOOK.md`](COOKBOOK.md).

## Raw data and diagnostics

```python
run.raw.table("chains")
run.raw.dictionary("monomers")
run.diagnostics.validation
run.diagnostics.memory
run.diagnostics.run_log
run.diagnostics.debug_log
run.summary("summary.json")
```

Public analysis methods raise explicit pyslimmc exceptions for unavailable or
chemically inapplicable results rather than returning misleading zeros.

## Collections of runs

```python
runs = sl.scan("results")

runs.completed
runs.failed
runs.interrupted

run = runs.one(run_id="run_000001")
subset = runs.match("feed_*_T3")
subset = subset.filter(status="completed", var_name="T", var_value=333.15)
```

String/run-ID selection and path selection are intentionally distinct; use
`runs.by_path(...)` for path-based lookup.

For a regular multidimensional study:

```python
sweep = subset.sweep("T", "feed_fraction")
sweep.info()
```

`runs.pack(...)`, `runs.as_table()` and `runs.model_diff()` support interactive
comparison workflows; see the exact API for their complete contracts.

## Reports

```python
report = sl.report("Run report")
report.text(run.summary().to_text())
report.plot(run.mwd())
report.save("report.pdf")

`Report` is optional output formatting; it does not replace the Python analysis API and it is not a GUI.
```

## See also

- common practical tasks: [`COOKBOOK.md`](COOKBOOK.md)
- meaning of KMC/distribution concepts: [`CONCEPTS.md`](CONCEPTS.md)
- API tree: [`reference/PYSLIMMC_API_TREE.md`](reference/PYSLIMMC_API_TREE.md)
- exact API: [`reference/PYSLIMMC_API.md`](reference/PYSLIMMC_API.md)
- exact signatures/defaults: [`reference/PYSLIMMC_SIGNATURES.md`](reference/PYSLIMMC_SIGNATURES.md)
- raw Storage contract: [`reference/STORAGE.md`](reference/STORAGE.md)
