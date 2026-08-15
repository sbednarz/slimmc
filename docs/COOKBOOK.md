# slimmc cookbook

This cookbook contains **short task recipes**, not a second examples
repository. Complete chemistry studies and literature reproductions belong in
the separate examples/literature repository. Recipes stay intentionally short. Most state the goal, the code, and the one
requirement, interpretation point, or caveat that matters for that task rather
than repeating a rigid template.

For the first complete runnable model, start with [`QUICKSTART.md`](QUICKSTART.md).

## 1. Open a run and inspect its endpoint

**Goal.** Load a completed run and inspect its main endpoint quantities.

```python
import pyslimmc as sl

run = sl.open("results/run_000001")
snap = run.final

print(run.status)
print(snap.conv.total)
print(snap.mn, snap.mw, snap.dispersity)
```

**Requirement.** `M_n`, `M_w` and `Đ` require chain data at the selected
snapshot.

**Caveat.** `run.final` is not interchangeable with `run.last` for incomplete
runs.

## 2. Plot conversion versus time

```python
import matplotlib.pyplot as plt

plt.plot(run.t / 60, run.conv.total)
plt.xlabel("Time, min")
plt.ylabel("Total monomer conversion")
plt.tight_layout()
plt.show()
```

For a copolymer:

```python
plt.plot(run.t / 60, run.conv["A"], label="A")
plt.plot(run.t / 60, run.conv["B"], label="B")
plt.plot(run.t / 60, run.conv.total, "--", label="total")
plt.legend()
```

`run.conv.total` is weighted by initial monomer counts.

## 3. Plot a final MWD

The default MWD is the exact discrete mass-weighted logarithmic form:

```python
mwd = run.mwd()
plt.plot(mwd.x, mwd.y)
plt.xlabel("log$_{10}$(M / g mol$^{-1}$)")
plt.ylabel("Polymer mass fraction")
plt.tight_layout()
plt.show()
```

For physical molar-mass support use another form explicitly:

```python
mwd = run.mwd(form="mass")
plt.plot(mwd.mass, mwd.y)
plt.xscale("log")       # display choice only
```

Exact unnormalized mass counts are a separate source projection:

```python
counts = run.mass_counts(pool="dead")
plt.vlines(counts.mass, 0, counts.count)
```

## 4. Plot a CLD

```python
cld = run.cld()  # form="number"
plt.plot(cld.x, cld.y)
plt.xlabel("Degree of polymerization, DP")
plt.ylabel("Chain number fraction")
plt.tight_layout()
plt.show()
```

Use `form="mass"` when the question is how polymer mass is distributed among
DP classes. pyslimmc then accumulates actual chain masses in each DP class.

## 5. Compare live and dead chains

```python
group = run.mwd_series(
    series=("live", "dead"),
    form="log",
    normalization="per_series",
)

for name, dist in group.series.items():
    plt.plot(dist.x, dist.y, label=name)
plt.legend()
plt.show()
```

`per_series` compares shapes. To preserve relative physical contributions use
`normalization="combined"`; the selected series must then be pairwise disjoint.

## 6. SEC broadening

```python
sec = run.sec(pool="dead", sigma_log10M=0.05)
plt.plot(sec.x, sec.y)
plt.xlabel("log$_{10}$(M / g mol$^{-1}$)")
plt.ylabel("$dW_{app}/d\\log_{10}M$")
plt.show()
```

SEC is a continuous instrumental-response model. It is deliberately separate
from the exact discrete MWD; histogram, KDE, and generic Gaussian smoothing are
not MWD methods.

## 7. Inspect copolymer composition drift

```python
plt.plot(run.t, run.f["A"], label="free monomer f(A)")
plt.plot(run.t, run.F.cum["A"], label="cumulative polymer F(A)")
plt.xlabel("Time, s")
plt.ylabel("Mole fraction")
plt.legend()
plt.show()
```

Higher-level terminal-model comparison is available through:

```python
run.copolymerization.compare_mayo_lewis()
```

Use the diagnostic/capability API before applying terminal or penultimate
analyses to a generic collection of runs.

## 8. Dyads, triads, and block statistics

```python
dyads = run.microstructure.dyads()
triads = run.microstructure.triads()
blocks_a = run.microstructure.run_lengths("A")
blockiness = run.microstructure.blockiness()
```

Aggregated statistics are available where the required counters were stored.
Position-resolved motifs, n-grams and full-sequence analyses require
`sequence_mode full`.

## 9. Compare realized events with kinetic competition

```python
fire = run.firings.fire_shares()
rate = run.firings.rate_shares()
```

**Requirement.** `rate_shares()` requires the `channel_propensities` Storage
table, which is currently written by the copo engine. `fire_shares()` can be
available without that table.

**Interpretation.** Fire shares summarize realized stochastic events; rate
shares summarize snapshot propensities. They answer different questions.

## 10. Semibatch: define and schedule a feed

A feed is a named constant-composition mixture:

```text
param init_volume 100 mL

feed F A 1.2
feed F B 0.8
feed F I 0.01
```

Feed concentrations are mol/L. Dose volumes without a suffix are litres; `mL`,
`ml`, and `ML` are accepted explicit millilitre forms.

Single dose:

```text
at 400 feed F 1 mL
```

Finite constant portion-wise feed:

```text
from 600 repeat 60 every 10 feed F 0.01 mL
```

This fires at `600, 610, ..., 1190 s` (subject to `t_end`). More complicated
profiles can be written as several explicit dose blocks; the engine
intentionally does not hide a ramp/expression evaluator inside `feed`.

## 11. Semibatch: two feeds / two pumps

```text
feed monomers A 1.5
feed monomers B 0.5
feed initiator I 0.1

from 0 repeat 60 every 10 feed monomers 0.1 mL
from 5 repeat 60 every 10 feed initiator 0.01 mL
```

`init_volume` is required. Physical reactor volume and KMC representation
volume are tracked separately and consistently.

Inspect the result:

```python
run.volume
run.kmc_volume
run.feeds["monomers"].events
run.feeds["monomers"].volume_cum
run.feeds["monomers"].moles_cum
```

## 12. `add_c` versus `set_c`

A physical inventory adjustment at constant volume:

```text
at 300 add_c A 0.10
at 600 add_c A -0.05
```

A forced concentration override:

```text
at 300 set_c A 0.10
```

`set_c` invalidates the physical balance for that species. pyslimmc raises
`AnalysisNotApplicableError` if you later ask for that balance; this is
intentional.

## 13. Physical semibatch balance

```python
run.balance.initial["A"]
run.balance.dosed["A"]
run.balance.total["A"]
run.balance.free["A"]
run.balance.consumed["A"]
run.balance.incorporated["A"]
```

Amounts are physical mol and require `param init_volume`.

## 14. Convert a small-species concentration to mg/mL

For a known molar mass:

```python
M_AIBN = 164.21  # g/mol

aibn_mg_mL = run.conc["AIBN"] * M_AIBN
final_aibn_mg_mL = run.final.conc["AIBN"] * M_AIBN
```

Because `mol/L * g/mol = g/L = mg/mL`, no additional factor is needed. This
works in batch and semibatch because `run.conc` already uses the current
reactor concentration. A generic `species` declaration does not store molar
mass, so the user supplies it.

## 15. Live/dead polymer concentration in mg/mL

Use a snapshot that contains chains:

```python
import numpy as np

snap = run.last_with_chains
dead = snap.chains.dead
live = snap.chains.live

dead_mg_mL = np.sum(dead.conc * dead.molar_mass)
live_mg_mL = np.sum(live.conc * live.molar_mass)
total_mg_mL = dead_mg_mL + live_mg_mL
```

`conc` already includes compressed-row multiplicity, so the sum is directly
mg/mL. `molar_mass` follows the stored mass-model interpretation.

If needed and supported, reconstruct explicitly:

```python
dead_repeat_units_mg_mL = np.sum(
    dead.conc * dead.masses(mass_model="repeat_units")
)
```

## 16. Select runs by `run_id` pattern

```python
import pyslimmc as sl

runs = sl.scan("results")
study = runs.match("f5a_DBI*_I3")

for one in study:
    print(one.run_id)
```

Matching uses shell-style glob syntax on the complete `run_id` and is
case-sensitive.

## 17. Parameter sweep

```python
runs = sl.scan("results")
sweep = runs.sweep("feed_IA", "temperature")

for one in sweep:
    print(
        one.var["feed_IA"].value,
        one.var["temperature"].value,
        one.mw[-1],
    )

sweep.info()
```

`sweep.info()` reports dimensions, axis values, duplicates, and missing grid
points. Sweep metadata should describe scientific parameter variation, not
random replicate identity.

## 18. Replicate trajectories and stochastic uncertainty

slimmc intentionally treats one model execution as one stochastic trajectory.
For an ensemble, create otherwise identical model files with distinct `seed`
and distinct `output_dir` values, run them, then collect the results.

Example endpoint analysis after the replicate runs exist:

```python
import numpy as np
import pyslimmc as sl

runs = sl.scan("results/replicates").completed
mw = np.asarray([one.final.mw for one in runs], dtype=float)

print("n =", len(mw))
print("mean Mw =", mw.mean())
print("SD Mw =", mw.std(ddof=1))
print("SE Mw =", mw.std(ddof=1) / np.sqrt(len(mw)))
```

Choose the number of replicates and statistical interval according to the
scientific question. slimmc does not claim that a single fixed number of
replicates is universally sufficient.

Do not mix a parameter sweep with a seed ensemble when estimating purely
stochastic variability.

## 19. Check `V_kMC` convergence

A successful run is not automatically a converged finite-population result.
Repeat the same scientific model at progressively larger `kmc_volume` and
compare the quantity you intend to report.

Typical targets include:

```python
run.final.conv.total
run.final.mn
run.final.mw
run.final.dispersity
run.chain_count.total[-1]
```

For distributions compare their scientifically relevant features, not only a
single moment. There is no universal chain-count threshold built into this
recipe.

## 20. Before trusting a result

A useful pre-publication checklist is:

```python
run.validate(strict=True)
run.mass_audit().raise_if_failed()
report = run.reproducibility.verify()
print(report.overall)
```

Then also ask questions that cannot be answered by a structural validator:

1. Did the model use the intended kinetic convention and units?
2. Is the selected mass model appropriate?
3. Is the required snapshot backed by `save_chains`?
4. Is the observable converged with respect to `kmc_volume`?
5. Is one stochastic trajectory sufficient, or are replicate seeds required?
6. Were any `set_c` operations used that invalidate a physical balance?
7. Are terminal/penultimate/microstructure assumptions appropriate for the
   chemistry being interpreted?

Validation, convergence, and scientific model adequacy are separate questions.

## 21. Compare with a classical kinetic limit

Simple limiting cases are useful scientific sanity checks because they test the
meaning of a model independently of Storage/API validation. For example, for a
well-mixed homo FRP model with one initiator-decay channel and one effective
termination coefficient, the stationary-state radical concentration has the
classical scale

$$
[P^\bullet] \sim \sqrt{\frac{f k_d[I]}{k_t}},
$$

when `k_t` is defined using the slimmc disappearance convention
`-d[P^bullet]/dt = 2 k_t [P^bullet]^2`. The corresponding propagation-rate
scale is

$$
R_p \sim k_p[M][P^\bullet].
$$

Use the exact channels and conventions declared in your model when making the
comparison. With multiple termination channels, transfer, depropagation,
non-steady initiation, or strong conversion effects the textbook limit is not
expected to match exactly. Other useful checks include Mayo transfer,
Mayo-Lewis copolymer composition, PLP characteristic chain length, and
`[M]_eq = k_dep/k_p` for a compatible simple prop/deprop model.

The purpose is not to force KMC onto an analytical approximation: it is to ask
whether a deliberately simplified model approaches the known limit when its
assumptions are satisfied.

## 22. Work with an interrupted run

```python
run = sl.open("results/interrupted_run", allow_incomplete=True)
snap = run.last
```

Never relabel `last` as `final`; incomplete Storage is intentionally explicit.

## 23. Build an optional PDF report

```python
r = sl.report("Polymerization report")
r.text(run.summary().to_text())
r.plot(run.mwd(), span="column")
r.save("report.pdf")
```

`Report` is a convenience PDF builder, not a GUI and not the primary analysis interface. Analysis still happens through `Run`, snapshots, tables and plotting objects in Python. For publication figures, plain Matplotlib remains available and gives complete control over axes, labels and layout.

## See also

- [`CONCEPTS.md`](CONCEPTS.md) — `V_kMC`, ensemble, pools, distributions,
  coordinates and normalization
- [`MODEL_SYNTAX.md`](MODEL_SYNTAX.md) — exact user-facing model statements
- [`PYSLIMMC.md`](PYSLIMMC.md) — task-oriented analysis guide
- [`reference/PYSLIMMC_API.md`](reference/PYSLIMMC_API.md) — exhaustive API
