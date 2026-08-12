# slimmc

[![CI](https://github.com/sbednarz/slimmc/actions/workflows/ci.yml/badge.svg)](https://github.com/sbednarz/slimmc/actions/workflows/ci.yml)
[![Release](https://github.com/sbednarz/slimmc/actions/workflows/release.yml/badge.svg)](https://github.com/sbednarz/slimmc/actions/workflows/release.yml)
[![License: GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-green.svg)](LICENSE)

**slimmc** is a stochastic simulator (kMC) for radical polymerization that
follows individual polymer chains as they grow, terminate, transfer, or
depolymerize. It supports homo- and copolymerization and can provide
conversion, chain populations, molecular-weight distributions, copolymer
composition, and selected microstructure information. Results can be explored
with `pyslimmc` in Python.

## What slimmc does — and does not do

- simulation of radical homo- and copolymerization under **batch and
  semibatch** conditions;
- prediction of conversion, living/dead chain populations, molecular-weight
  and chain-length distributions, copolymer composition, and selected
  microstructure;
- analysis across the full molecular-size range, from **small molecules
  through oligomers to linear polymers**, using concentrations, chain
  populations, MWD/CLD, and neutral chain-mass spectra;
- explicit tracking of **individual polymer-chain populations** rather than
  only average quantities;
- propagation, termination, transfer, and supported **depropagation**
  mechanisms;
- terminal and explicit penultimate models for copolymer propagation;
- **programmable time- and state-dependent process changes**, including feeds,
  temperature programs, concentration changes, **PLP/UV pulse sequences**,
  inhibition periods, kinetic switches, and conversion-dependent changes of
  rate constants;
- **in-depth Python analysis** with `pyslimmc`, including distributions,
  molecular-weight moments, chain populations, copolymer composition,
  microstructure, validation, mass-balance checks, and reproducibility
  verification;
- preliminary **optimization** of process and model parameters with
  `pyslimmc-opt`;
- **fully traceable and versioned simulation results**: each run records the
  model, software versions, build information, random seed, and provenance
  metadata, so results can be reproduced, compared, and audited later rather
  than existing only as exported tables.

## What slimmc does not cover

- no dedicated mechanistic engine for **controlled/living radical
  polymerization (RDRP)**, including **ATRP, RAFT, NMP, OMRP, RITP, and
  related reversible-deactivation radical polymerization methods**; some
  simplified kinetic effects may be represented with the available generic
  reaction channels, but slimmc does not model their full
  activation/deactivation, reversible-transfer, or dormant/active-state
  mechanisms;
- linear-chain model only: no general **branching, crosslinking, cyclization,
  gelation, or polymer-network/topology** engine;
- no general **chain-length-dependent kinetic laws**;
- no coupled **non-isothermal reactor energy balance**;
- finite stochastic populations introduce statistical noise, so important
  results should be checked for convergence with respect to `kmc_volume` and,
  where appropriate, across independent trajectories.

See [Scope and limitations](docs/LIMITATIONS.md) before designing a study that
depends on any of these features.

## A short history

slimmc began as a manually written Nim homopolymer/PLP simulator and was later
expanded into the current homo/copo family, shared storage format, and Python
analysis stack. Subsequent development, testing, and documentation were
conducted with extensive assistance from **OpenAI ChatGPT and Anthropic
Claude** as engineering tools under human direction and review; they are not
project authors.

The project concept was influenced by the long tradition of polymer kMC
software, in particular the excellent open-source
[mcPolymer](https://www.itc.tu-clausthal.de/en/research/mcpolymer) developed by
Marco Drache and collaborators.

## How a calculation actually goes

You start in a text editor. A slimmc calculation is described in a plain-text
model file, where you write down the chemistry you want to simulate: the
monomer and its concentration, the initiator, the temperature, the reactions
that can occur and their rate constants, and how long the polymerization
should run. It reads much like a reaction scheme written out line by line, and
nothing else is needed to define a calculation.

Then you give that file to `slimmc`. It first checks the model and reports
mistakes without running the simulation, which is worth doing every time — a
misspelled species costs you a second instead of an hour. When the model is
valid, you run it again and the calculation starts.

When the simulation finishes, the results are written to a directory. To see
what the simulation produced, run `slimmc-summary` on that directory: it
reports whether the run completed, the final conversion, molar-mass averages
and dispersity, and the calculation time. That is one command and requires no
programming; for many questions, it may already provide the answer you need.

When you want more than a summary, `pyslimmc` opens the finished calculation in
Python and gives you access to everything it contains — conversion over time,
concentrations, individual chains, molecular-weight and chain-length
distributions, copolymer composition and composition drift, and sequence
statistics such as dyads and triads. You can plot the results, export them, or
use them in your own analysis. `pyslimmc` is read-only: analysing a run cannot
alter the stored simulation result.

And if your question runs the other way round — not "what happens under these
conditions?" but "which conditions give me the polymer I want?" —
`pyslimmc-opt` can search the parameter space by running slimmc repeatedly.

## Documentation

Choose a path rather than reading the documentation front to back:

- **First calculation:** [Quick start](docs/QUICKSTART.md)
- **Write/edit models:** [Model syntax](docs/MODEL_SYNTAX.md)
- **Understand stored results:** [Simulation results](docs/SIMULATION_RESULTS.md)
- **Analyse with Python:** [pyslimmc guide](docs/PYSLIMMC.md)
- **Check scope before modeling:** [Limitations](docs/LIMITATIONS.md)
- **Equations and kinetic conventions:** [Theory](docs/THEORY.md)
- **Exact engine/API contracts:** [documentation index](docs/README.md)
- **Develop slimmc:**
  [development workflow](docs/development/DEVELOPMENT.md); coding agents should
  start with [AGENTS.md](AGENTS.md)

The `docs/reference/` layer is intended to provide complete public syntax/API
coverage, including generated callable-signature inventories for pyslimmc and
pyslimmc-opt.

## Current releases

- **slimmc 5.0.0** — native simulator binaries for Linux and Windows  
  [Release page](https://github.com/sbednarz/slimmc/releases/tag/slimmc-v5.0.0)

- **pyslimmc 4.0.0** — Python analysis package  
  [Release page](https://github.com/sbednarz/slimmc/releases/tag/pyslimmc-v4.0.0)

- **pyslimmc-opt 1.0.0** — optional optimization package  
  [Release page](https://github.com/sbednarz/slimmc/releases/tag/pyslimmc-opt-v1.0.0)


## Installation and builds

Native slimmc release archives are produced for:

- Linux x86-64, glibc 2.28+ (modern systems);
- Linux x86-64, static musl build (older or minimal systems);
- Windows 10/11 x86-64.

The source tree can also build the simulation engines with Nim 2.2.10.

Python analysis is provided by:

- `pyslimmc` — the main analysis package for reading and analysing slimmc
  results;
- `pyslimmc-opt` — an optional package for parameter studies and preliminary
  optimization.

`pyslimmc` requires NumPy. Matplotlib support is optional and is only needed
for plotting and report generation. `pyslimmc-opt` has additional dependencies
for optimization workflows.

Exact release, installation, and build details are in
[Releases](docs/development/RELEASES.md).

## Examples and case studies

Worked examples and scientific case studies are maintained in separate
repositories so they can develop independently from the core simulation code.

- [`slimmc-examples`](https://github.com/sbednarz/slimmc-examples) — tutorials
  and practical examples
- [`slimmc-case-studies`](https://github.com/sbednarz/slimmc-case-studies) —
  literature-based studies, validation, and benchmarks

## Citation and license

Citation metadata are provided in `CITATION.cff`. When reporting or citing a
calculation, record the slimmc/pyslimmc component versions used for it.

slimmc is distributed under the GNU General Public License v3 or later; see
`LICENSE`.

## Acknowledgements

This work has been supported by the Czech–Polish Lead Agency project funded by
the Czech Science Foundation (**25-15669K**) and the National Science Centre,
Poland, WEAVE-UNISONO grant (**2024/06/Y/ST5/00062**).

Full acknowledgements, including the role of AI-assisted development tools,
are in [ACKNOWLEDGEMENTS.md](ACKNOWLEDGEMENTS.md).
