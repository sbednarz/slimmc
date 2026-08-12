# pyslimmc-opt callable signatures

This is the exhaustive callable inventory for the public `pyslimmc_opt` API.
It is generated from the installed source by `scripts/update_api_signatures.py`;
CI rejects a stale inventory. The human-readable reference remains in
[`PYSLIMMC_OPT.md`](PYSLIMMC_OPT.md).

## Top-level callables
- `pyslimmc_opt.Float(name: 'str', low: 'float', high: 'float') -> None`
- `pyslimmc_opt.LogFloat(name: 'str', low: 'float', high: 'float') -> None`
- `pyslimmc_opt.Result(trials: 'tuple[Trial, ...]', _study: "'Study | None'" = None) -> None`
- `pyslimmc_opt.Study(*, model: 'str | Path', parameters: 'Sequence[Float | LogFloat]', objective: 'Objective', output_dir: 'str | Path' = 'opt_results', slimmc: 'str | Path' = 'slimmc', build: 'Build | None' = None, trial_timeout: 'float | None' = None) -> 'None'`
- `pyslimmc_opt.Trial(number: 'int', params: 'dict[str, float]', loss: 'float | None', status: 'str', run_path: 'Path', slimmc_seed: 'int | None' = None, error: 'str' = '') -> None`

## Object callables

### `Result`

Properties: `best_loss`, `best_params`, `best_run`, `best_trial`, `completed`.

- `Result.verify(*, repeats: 'int' = 5, slimmc_seed: 'int' = 10000) -> 'dict[str, object]'`

### `Study`

- `Study.optimize(*, runs: 'int', initial_random: 'int' = 8, optimizer_seed: 'int' = 42, slimmc_seed: 'int' = 1000, candidates: 'int' = 10000) -> 'Result'`


## See also

- [`PYSLIMMC_OPT.md`](PYSLIMMC_OPT.md) — pyslimmc-opt guide and API semantics
- [`../PYSLIMMC.md`](../PYSLIMMC.md) — analysis of Slimmc run results
