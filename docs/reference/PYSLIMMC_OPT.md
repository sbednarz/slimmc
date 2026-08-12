# pyslimmc-opt 1.0 public API

`pyslimmc_opt` minimizes one finite scalar objective by changing numeric
declarations in a normal Slimmc model. It performs sequential Gaussian-process
Bayesian optimization with Expected Improvement and stores every real Slimmc
trial. Smaller loss is always better.

## Minimal workflow

```python
from pyslimmc_opt import Float, LogFloat, Study


def objective(run, params):
    conversion = float(run.conv.total[-1])
    return (conversion - 0.90) ** 2


study = Study(
    model="base.model",
    parameters=[
        Float("temperature", 330.0, 370.0),
        LogFloat("kp", 1.0e2, 1.0e5),
    ],
    objective=objective,
)

result = study.optimize(
    runs=40,
    initial_random=8,
    optimizer_seed=42,
    slimmc_seed=1000,
)
```

The matching model declaration must be marked with `var`:

```text
var param temperature K
param temperature 350
```

## `Float` and `LogFloat`

```python
Float(name: str, low: float, high: float)
LogFloat(name: str, low: float, high: float)
```

`Float` samples linearly between its bounds. `LogFloat` samples uniformly in
log space and requires `0 < low < high`; use it for positive parameters that
span orders of magnitude, such as many kinetic constants.

| Field | Meaning |
|---|---|
| `name` | Exact target name from one `var` declaration. |
| `low` | Inclusive lower bound. |
| `high` | Inclusive upper bound; must be greater than `low`. |

Automatic replacement supports numeric `var param`, `var species`,
`var monomer`, `var endgroup`, and fixed `var rate` declarations written as
`rate NAME VALUE` or `rate NAME const VALUE`. An Arrhenius declaration has two
independent numbers and is therefore not silently treated as one `Float`; use
`build=` to choose whether the parameter changes `A`, `Ea`, or both.

## `Study`

```python
Study(
    *,
    model: str | Path,
    parameters: Sequence[Float | LogFloat],
    objective: Callable[[Run, Mapping[str, float]], float],
    output_dir: str | Path = "opt_results",
    slimmc: str | Path = "slimmc",
    build: Callable[[str, Mapping[str, float]], str] | None = None,
    trial_timeout: float | None = None,
)
```

| Parameter | Contract |
|---|---|
| `model` | Existing `.model` template; resolved to an absolute path. |
| `parameters` | Non-empty sequence with unique names. Without `build`, every name must have a matching `var`. |
| `objective` | Called as `objective(run, params)` after a completed run; must return one finite float. |
| `output_dir` | Study directory. Existing compatible trials are resumed. |
| `slimmc` | Unified CLI executable name or path. |
| `build` | Optional callback called as `build(rendered_text, params)`; it must return complete model text as `str`. |
| `trial_timeout` | Optional positive wall-time limit in seconds for each Slimmc subprocess; a timeout becomes a `Trial(status="timeout")`. |

Without `build`, every optimization variable must have exactly one matching
model declaration. The study validates this before the first trial, and the
renderer independently verifies that exactly one substitution was performed; a
missing or duplicate target is an error rather than a silently ineffective
optimization variable.

The automatic renderer first substitutes supported `var` declarations, forces
`param output_dir "result"`, and inserts/replaces the trial seed. `build` then
receives that rendered text, so it may edit unsupported structures such as a
feed composition or Arrhenius pair:

```python
def build(model_text, params):
    return model_text.replace(
        "rate kd arr 1.0e15 125000",
        f"rate kd arr {params['A']:.17g} {params['Ea']:.17g}",
    )
```

## `Study.optimize`

```python
Study.optimize(
    *,
    runs: int,
    initial_random: int = 8,
    optimizer_seed: int = 42,
    slimmc_seed: int = 1000,
    candidates: int = 10000,
) -> Result
```

| Parameter | Contract |
|---|---|
| `runs` | Total trial budget, integer `>=1`; includes resumed trials. |
| `initial_random` | Random-design trials before GP/EI proposals, integer `>=1`. It may exceed `runs`. |
| `optimizer_seed` | Non-negative seed for proposals and candidate pools. |
| `slimmc_seed` | Non-negative base seed; trial `n` uses `slimmc_seed+n-1`. |
| `candidates` | Candidate points evaluated per EI step, integer `>=100`. |

Each trial gets `runs/trial_NNNN/`, a rendered model, `stdout.txt`,
`stderr.txt`, and its Storage `result/` when the engine initializes it. Engine,
Storage, objective, non-finite-loss, or callback errors produce a failed
`Trial`; they do not abort the remaining budget. If no trial completes, `report.txt` is still written with the failed/timeout
reasons and `optimize()` then raises `RuntimeError` containing representative
causes.

### Resume contract

Calling `optimize()` again in the same `output_dir` resumes `trials.tsv` and
replays the optimizer RNG so the next proposal matches an uninterrupted study.
The following settings must remain identical: model path, parameter names,
bounds and linear/log scale, `trial_timeout`, `initial_random`,
`optimizer_seed`, `slimmc_seed`, and `candidates`.
Only `runs` may increase. A smaller budget than the number of stored trials is
rejected. Use a new `output_dir` for a changed study.

## `Trial`

```python
Trial(
    number: int,
    params: dict[str, float],
    loss: float | None,
    status: str,
    run_path: Path,
    slimmc_seed: int | None = None,
    error: str = "",
)
```

`status` is `completed`, `failed`, or `timeout`. A completed trial has a finite
`loss` and readable `run_path`; unsuccessful trials have `loss=None` and a
concise reason in `error`. `params` is the exact proposal used to render the
trial.

## `Result`

```python
result.trials       # tuple[Trial, ...], including failures
result.completed    # completed trials only
result.best_trial   # minimum-loss Trial
result.best_params  # copy of best_trial.params
result.best_loss    # float
result.best_run     # Path
```

Verification uses independent KMC seeds at the best point:

```python
summary = result.verify(repeats=5, slimmc_seed=5000)
```

```python
Result.verify(*, repeats: int = 5, slimmc_seed: int = 10000) -> dict[str, object]
```

`repeats` must be at least two. The returned dictionary contains:

| Key | Value |
|---|---|
| `repeats` | requested repeat count |
| `completed` | successful verification count |
| `mean_loss` | arithmetic mean over completed repeats |
| `std_loss` | sample standard deviation (`ddof=1`) |
| `min_loss`, `max_loss` | observed extrema |
| `trials` | tuple of all verification `Trial` objects |

Verification replaces the study's previous `verification/` directory and
writes `verification.tsv`, `report.txt`, and one run directory per repeat. It
raises if all repeats fail.

## Study artifacts

After optimization the study directory contains:

If scikit-learn reports that Gaussian-process kernel parameters reached their
optimization bounds, pyslimmc-opt captures those `ConvergenceWarning` messages
instead of repeatedly writing them to stderr. Unique messages are recorded once
in the `GP DIAGNOSTICS` section of `report.txt`, together with a note that the
objective may be weakly structured or dominated by stochastic noise. This is a
diagnostic, not a failed trial.

- `trials.tsv`: append-only trial history;
- `study.json`: resume-critical settings;
- `report.txt`: search space, settings, failures, and cumulative best loss;
- `best.model`: rendered best model;
- `best_run.txt`: best trial, loss, and run path;
- `runs/`: all trial directories.

The MVP intentionally supports one model, numeric continuous `Float`/`LogFloat`
parameters, one scalar objective, sequential execution, and one process-local
study writer. The random-candidate EI search is intended primarily for
low-dimensional studies (typically about 1--5 continuous variables).
Multi-objective optimization, categorical/integer spaces, parallel writers,
pruning, and distributed databases are outside this API.

`best_trial` is the best **observed stochastic realization**, not an unbiased
estimate of the expected objective at that point. For scientific conclusions,
repeat the selected point with independent Slimmc seeds via `result.verify()`.

## See also

- [`PYSLIMMC_OPT_SIGNATURES.md`](PYSLIMMC_OPT_SIGNATURES.md) — Exact callable signatures
- [`../PYSLIMMC.md`](../PYSLIMMC.md) — pyslimmc guide
- [`../COOKBOOK.md`](../COOKBOOK.md) — Cookbook
