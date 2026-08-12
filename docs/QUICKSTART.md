# Quick start

This page takes a new user from a `.model` file to a checked simulation, a
human-readable summary, and a first Python plot. It intentionally uses only a
small part of the model language. The complete shared syntax is in
[`MODEL_SYNTAX.md`](MODEL_SYNTAX.md); engine-specific details are in the
[reference](reference/).

Slimmc has no GUI. You run models from the command line and use Python with
[`pyslimmc`](PYSLIMMC.md) to inspect, analyse and plot the results. The Python
analysis can be run as ordinary scripts or in the standard Python interpreter,
IPython, Jupyter Notebook/JupyterLab, Marimo, VS Code, and other Python
environments.

## 1. Install slimmc

Native release archives contain `slimmc` and `slimmc-summary` for three targets:

- Linux x86-64 with glibc 2.28 or newer;
- static Linux x86-64 using musl;
- Windows 10/11 x86-64.

Download the archive for your system from the project Releases page, unpack it,
and make the `bin/` directory available on your command path. Source-build and
release-engineering details are documented in
[`development/RELEASES.md`](development/RELEASES.md).

Check the installation:

```bash
slimmc --version
```

## 2. Save a first model

Create `styrene.model`:

```model
# Styrene free-radical polymerization with AIBN at 80 °C.
desc "Styrene FRP with AIBN at 80 C"

param output_dir "results/styrene"
param kmc_volume 1.0e-15
param temperature 353
param t_end 21600
param seed 1000

monomer Sty 1.0 104.15
species AIBN 1.0e-3
species R 0.0

polymer P active
polymer D dead

rate kd arr 2.89e15 130230
rate ki const 331.0
rate kp arr 3.09e7 31700
rate ktc const 1.0e7
rate ktd const 1.0e7

# The final factor in the initiator-decomposition channel is efficiency f.
rxn AIBN -> R + R kd 0.6
macro init R + Sty -> P ki
macro prop P + Sty -> P kp
macro term_c P + P -> D ktc
macro term_d P + P -> D + D ktd

every 360 save
at 21600 save_chains
```

Two Arrhenius pairs in this compact model are traceable literature/data-sheet
values. The styrene propagation Arrhenius pair (`A = 3.09e7
L mol^-1 s^-1`, `Ea = 31.7 kJ/mol`) is the revised IUPAC PLP benchmark from
Beuermann *et al.*, *Polymer Chemistry* **13** (2022) 1891-1900,
[DOI 10.1039/D2PY00147K](https://doi.org/10.1039/D2PY00147K). The AIBN
decomposition pair (`A = 2.89e15 s^-1`, `Ea = 130.23 kJ/mol`) is from the
[Nouryon Perkadox AIBN product data sheet](https://www.nouryon.com/globalassets/inriver/resources/pds-perkadox-aibn-acrylics-production-glo-en.pdf).
`ki`, `ktc`, and `ktd` here are compact illustrative choices rather than a claim of one universal styrene
parameter set; real studies should use conditions-appropriate primary or
critically evaluated data.

The model has five blocks:

1. `param` defines the simulation scale, temperature, time, random seed, and output;
2. `monomer` and `species` define initial concentrations;
3. `polymer` declares live and dead chain pools;
4. `rate` and reaction lines define the kinetic model;
5. `save`/`save_chains` define what is retained for later analysis.

Concentrations are in mol/L, temperature is in K, time is in s, molar masses are
in g/mol, and `kmc_volume` is in L. For the physical and numerical consequences of `kmc_volume`, see
[`CONCEPTS.md`](CONCEPTS.md#simulation-volume-v_kmc) and [`LIMITATIONS.md`](LIMITATIONS.md). Treat it as the
discrete stochastic simulation volume, not the laboratory reactor volume.

## 3. Check before running

```bash
slimmc --check styrene.model
```

`--check` parses and validates the model and runs the discretization preflight,
but does not simulate it. Fix reported errors before starting a long run.

## 4. Run

```bash
slimmc styrene.model
```

The relative `output_dir "results/styrene"` is resolved from the directory
containing the model. Using a run-specific subdirectory avoids collisions when
you add a second model. A successful run writes a Slimmc Storage directory
there.

On a typical current desktop/laptop this deliberately non-trivial example is an
order-of-a-minute run and executes roughly $2\times10^8$ SSA events. Runtime
depends strongly on CPU, build mode, and platform; a quiet terminal between the
initial and final `[run]` messages does not mean the simulation has stopped.

The unified CLI dispatches the model to the homo or copo engine from its
contents; users normally invoke `slimmc`, not the engine executable directly.

## 5. Inspect the result without Python

```bash
slimmc-summary results/styrene
```

The summary reports the engine and versions, run status, validation status,
seed, final time/event, snapshot count, final DP/molar-mass moments, peak memory,
and result size. JSON output is also available:

```bash
slimmc-summary results/styrene --format json
```

This is the fastest way to verify that a run completed and to inspect its basic
polymer statistics without writing analysis code.

## 6. Open the result with pyslimmc

Install the analysis package and plotting extra from the source/package release
used for your project, then:

```python
import matplotlib.pyplot as plt
import pyslimmc as sl

run = sl.open("results/styrene")
print(run.status)
print(run.conv["Sty"][-1])

mwd = run.final.mwd(method="gaussian", coordinate="log10", basis="mass")
plt.plot(mwd.x, mwd.y)
plt.xscale("log")
plt.xlabel("Molar mass, g/mol")
plt.ylabel("dW/dlog10(M)")
plt.tight_layout()
plt.show()
```

The exact `Mn`, `Mw`, `Mz` and dispersity values are calculated from the stored
discrete chain population. A smoothed MWD is a plotting representation of that
population, not the source of the moments.

## See also

- Change or build a model: [`MODEL_SYNTAX.md`](MODEL_SYNTAX.md)
- Understand what a run contains: [`SIMULATION_RESULTS.md`](SIMULATION_RESULTS.md)
- Analyse results in Python: [`PYSLIMMC.md`](PYSLIMMC.md)
- Understand SSA, propensities and polymer statistics: [`THEORY.md`](THEORY.md)
- Check deliberate scope boundaries before designing a study: [`LIMITATIONS.md`](LIMITATIONS.md)
- Look up exact engine syntax: [`reference/HOMO.md`](reference/HOMO.md) and [`reference/COPO.md`](reference/COPO.md)
