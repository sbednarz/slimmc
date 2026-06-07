# slimmc
A Monte Carlo simulation program of radical polymerization kinetics


## slimmc — model file syntax

`slimmc` runs a kinetic Monte Carlo (kMC) simulation of free‑radical
polymerization. It reads a plain‑text **model file** that defines the kinetic
parameters and initial concentrations, and schedules **breakpoints** — points
in simulated time at which commands run (write output, change a concentration,
print a message).

A model file is just a sequence of lines. Every non‑blank, non‑comment line is
one of three things: a **parameter assignment**, a **`list` directive**, or a
**breakpoint**. Lines are processed independently, so their order is free — with
one requirement: **at least one breakpoint must be defined**, otherwise slimmc
aborts.

---

## Lexical rules

- **Comment** — a line whose **first character** is `#`. It must be in column 0;
  an indented `#` is *not* treated as a comment.
- **Blank lines** are ignored.
- **Numbers** are floating point; scientific notation is accepted
  (`1e-4`, `5.0e-15`, `-2.1`, `3500`).
- **Time** is in seconds; **concentrations** in mol/L.

---

## Kinetic model

The parameters below describe this reaction scheme:

```
I            ->  2 R•        kd , efficiency f      (initiator decomposition)
R• +  M      ->    P•        ki                     (initiation)
P• +  M      ->    P•        kp                     (propagation)
P• +  P•     ->    D         ktc                    (termination by combination)
P• +  P•     ->    D + D     ktd                    (termination by disproportionation)
```

Species symbols used throughout (parameters, `dc`, output files):

| symbol | meaning                          |
|--------|----------------------------------|
| `I`    | initiator                        |
| `Rx`   | primary radical (R•)             |
| `M`    | monomer                          |
| `Px`   | macroradical / living chain (P•) |
| `D`    | dead chain                       |

---

## Parameters

Syntax: `name = value` (whitespace around `=` optional). Unknown names abort with
an error. Any parameter left out keeps its default.

| name   | meaning                                   | unit        | default   |
|--------|-------------------------------------------|-------------|-----------|
| `kd`   | initiator decomposition rate constant     | 1/s         | `0.0`     |
| `f`    | initiator efficiency                      | —           | `0.0`     |
| `ki`   | initiation rate constant                  | L/(mol·s)   | `0.0`     |
| `kp`   | propagation rate constant                 | L/(mol·s)   | `0.0`     |
| `ktc`  | termination by combination                | L/(mol·s)   | `0.0`     |
| `ktd`  | termination by disproportionation         | L/(mol·s)   | `0.0`     |
| `cI0`  | initial initiator concentration           | mol/L       | `0.0`     |
| `cM0`  | initial monomer concentration             | mol/L       | `0.0`     |
| `cRx0` | initial primary‑radical concentration     | mol/L       | `0.0`     |
| `cPx0` | initial macroradical concentration        | mol/L       | `0.0`     |
| `cD0`  | initial dead‑chain concentration          | mol/L       | `0.0`     |
| `MwM`  | monomer molar mass                        | g/mol       | `0.0`     |
| `V_MC` | simulation control volume                 | L           | `1.0e-17` |
| `seed` | PRNG seed (`0` → seeded from clock)        | uint32      | `0`       |

`V_MC` sets the number of molecules in the box: `n = c · V_MC · N_A`. It is
chosen so the box holds a workable number of radicals (≈10 at a steady‑state
radical concentration of ~1e‑7 mol/L corresponds to `V_MC ≈ 1.7e-16 L`).

---

## `list` directives

Print diagnostic information. One directive per line:

| line                 | effect                                                       |
|----------------------|--------------------------------------------------------------|
| `list parameters`    | print the parsed parameter set                               |
| `list breakpoints`   | print the expanded breakpoint schedule                       |
| `list initialstate`  | print initial molecule counts (`nI`, `nRx`, `nM`, `nPx`, `nD`, total `N`) |

---

## Breakpoints

A breakpoint binds a simulated time to a comma‑separated list of commands.

**Single breakpoint**

```
<time> : <commands>
```

**Loop (a series of breakpoints)**

```
<t0> : <dt> : <N> : <commands>
```

generates **N + 1** breakpoints at `t0 + i·dt` for `i = 0, 1, … , N`, each
carrying the same command list. Example: `0:100:10:conc` produces breakpoints at
`0, 100, 200, … , 1000` (eleven points).

If two breakpoints land on the same time, their command lists are merged.

Commands are separated by commas. A comma **inside** a single‑quoted `print`
string is literal, not a separator:

```
3000000:print 'final dump, including a comma', conc, poly
```

---

## Commands

| command          | action                                                                 |
|------------------|------------------------------------------------------------------------|
| `conc`           | write current concentrations of all species to files (see below)       |
| `poly`           | compute and write the chain‑length / molecular‑weight distributions    |
| `dc <S> <v>`     | change concentration of species `S` (`I`/`M`/`Rx`/`Px`/`D`) by `v` mol/L (`v` may be negative); used for semi‑batch feeding |
| `print`          | print a progress line: `t=… (step/total) elapsed‑seconds`              |
| `print progress` | same as bare `print`                                                    |
| `print '<text>'` | print the literal single‑quoted text                                   |

---

## Output files

`conc` **appends** one row per breakpoint to five files — one per species —
`cI.txt`, `cRx.txt`, `cM.txt`, `cPx.txt`, `cD.txt`. Each row is tab‑separated:

```
step    time(s)    concentration(mol/L)    molecule_count
```

`poly` **overwrites** four files per call, with the breakpoint time encoded in
the filename (`P` = living chains, `D` = dead chains):

- `P.<time>.cld.txt`, `D.<time>.cld.txt` — chain‑length distribution; a
  `#nChains: <count>` header line, then `DP <tab> count` for each populated
  degree of polymerization.
- `P.<time>.HlogM.txt`, `D.<time>.HlogM.txt` — the GPC‑style `w(log M)` trace;
  `log10(M) <tab> dw`, where `dw ∝ n · M²` (here `n · DP² · MwM²`) and
  `M = DP · MwM`.

---

## Complete example

```text
# --- initiator decomposition ---
kd   = 1e-4
f    = 0.45
# --- initiation ---
ki   = 0.0002
# --- propagation ---
kp   = 3.5e3
# --- termination ---
ktc  = 1.2e6
ktd  = 2.1e6
# --- initial concentrations (mol/L) ---
cI0  = 0.1
cM0  = 1.2
cRx0 = 0.0
cPx0 = 0.0
cD0  = 0.0
# --- monomer molar mass and control volume ---
MwM  = 101.2
V_MC = 2.0e-17
# --- PRNG seed (0 => seeded from the clock) ---
seed = 12345

# echo the parsed setup and the initial molecule counts
list parameters
list initialstate
list breakpoints

# sample concentrations every 100 s, from 0 to 1000 s
0:100:10:conc

# sample the MWD a few times
200:200:4:poly

# semi-batch: add monomer at t = 500 s
500:dc M 0.5

# final dump at the end of the run
3000000:print 'end of run', conc, poly
```

Run it with:

```sh
slimmc model.model
``` 
