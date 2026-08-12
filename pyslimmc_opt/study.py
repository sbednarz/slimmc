from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np
from scipy.special import ndtr
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.exceptions import ConvergenceWarning


@dataclass(frozen=True)
class Float:
    name: str
    low: float
    high: float

    def __post_init__(self) -> None:
        if not self.name or self.low >= self.high:
            raise ValueError("Float requires a name and low < high")


@dataclass(frozen=True)
class LogFloat:
    name: str
    low: float
    high: float

    def __post_init__(self) -> None:
        if not self.name or self.low <= 0 or self.low >= self.high:
            raise ValueError("LogFloat requires a name and 0 < low < high")


@dataclass(frozen=True)
class Trial:
    number: int
    params: dict[str, float]
    loss: float | None
    status: str
    run_path: Path
    slimmc_seed: int | None = None
    error: str = ""


@dataclass(frozen=True)
class Result:
    trials: tuple[Trial, ...]
    _study: "Study | None" = field(default=None, repr=False, compare=False)

    @property
    def completed(self) -> tuple[Trial, ...]:
        return tuple(t for t in self.trials if t.status == "completed" and t.loss is not None)

    @property
    def best_trial(self) -> Trial:
        done = self.completed
        if not done:
            raise RuntimeError("No completed trials")
        return min(done, key=lambda t: float(t.loss))

    @property
    def best_params(self) -> dict[str, float]:
        return dict(self.best_trial.params)

    @property
    def best_loss(self) -> float:
        return float(self.best_trial.loss)

    @property
    def best_run(self) -> Path:
        return self.best_trial.run_path

    def verify(self, *, repeats: int = 5, slimmc_seed: int = 10_000) -> dict[str, object]:
        """Repeat the best point with independent Slimmc seeds."""
        if self._study is None:
            raise RuntimeError("This Result is not attached to a Study")
        return self._study._verify(self.best_params, repeats=repeats, slimmc_seed=slimmc_seed)


Objective = Callable[[object, Mapping[str, float]], float]
Build = Callable[[str, Mapping[str, float]], str]


class Study:
    """Optimize numeric Slimmc ``var`` declarations with GP + Expected Improvement."""

    def __init__(
        self,
        *,
        model: str | Path,
        parameters: Sequence[Float | LogFloat],
        objective: Objective,
        output_dir: str | Path = "opt_results",
        slimmc: str | Path = "slimmc",
        build: Build | None = None,
        trial_timeout: float | None = None,
    ) -> None:
        self.model = Path(model).resolve()
        self.parameters = tuple(parameters)
        self.objective = objective
        self.build = build
        self.output_dir = Path(output_dir).resolve()
        self.trial_timeout = trial_timeout
        self._gp_fit_warnings: list[str] = []
        if trial_timeout is not None and trial_timeout <= 0:
            raise ValueError("trial_timeout must be > 0 or None")
        slimmc_text = str(slimmc)
        slimmc_path = Path(slimmc).expanduser()
        has_path_separator = any(sep and sep in slimmc_text for sep in (os.sep, os.altsep))
        self.slimmc = str(slimmc_path.resolve()) if (slimmc_path.is_absolute() or has_path_separator) else slimmc_text
        if not self.model.is_file():
            raise FileNotFoundError(self.model)
        if not self.parameters:
            raise ValueError("At least one parameter is required")
        names = [p.name for p in self.parameters]
        if len(names) != len(set(names)):
            raise ValueError("Parameter names must be unique")
        self._template = self.model.read_text(encoding="utf-8")
        self._var_kinds = _parse_var_kinds(self._template)
        missing = [name for name in names if name not in self._var_kinds]
        if missing and self.build is None:
            raise ValueError(f"Parameters not declared with var: {', '.join(missing)}")
        unsupported_rates = [
            name for name in names
            if self._var_kinds.get(name) == "rate"
            and not re.search(
                rf"^\s*rate\s+{re.escape(name)}\s+(?:const\s+)?\S+\s*(?:#.*)?$",
                self._template,
                flags=re.MULTILINE,
            )
        ]
        if unsupported_rates and self.build is None:
            raise ValueError(
                "automatic rate optimization supports fixed rates only; "
                f"use build= for: {', '.join(unsupported_rates)}"
            )
        if self.build is None:
            missing_targets = []
            duplicate_targets = []
            for name in names:
                kind = self._var_kinds.get(name)
                if kind is None:
                    continue
                count = _matching_declaration_count(self._template, kind, name)
                if count == 0:
                    missing_targets.append(name)
                elif count > 1:
                    duplicate_targets.append(name)
            if missing_targets:
                raise ValueError(
                    "optimization variables have no matching model declaration: "
                    + ", ".join(missing_targets)
                )
            if duplicate_targets:
                raise ValueError(
                    "optimization variables match multiple model declarations: "
                    + ", ".join(duplicate_targets)
                )

    def optimize(
        self,
        *,
        runs: int,
        initial_random: int = 8,
        optimizer_seed: int = 42,
        slimmc_seed: int = 1_000,
        candidates: int = 10_000,
    ) -> Result:
        if runs < 1 or initial_random < 1 or candidates < 100:
            raise ValueError("runs >= 1, initial_random >= 1, candidates >= 100 required")
        if optimizer_seed < 0 or slimmc_seed < 0:
            raise ValueError("optimizer_seed and slimmc_seed must be >= 0")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "runs").mkdir(exist_ok=True)
        self._write_study_json(runs, initial_random, optimizer_seed, slimmc_seed, candidates)
        trials = self._load_trials(slimmc_seed)
        if len(trials) > runs:
            raise ValueError(
                f"runs={runs} is smaller than the {len(trials)} existing trials; "
                "resume with an equal or larger budget"
            )
        rng = np.random.default_rng(optimizer_seed)
        self._restore_rng(trials, rng, initial_random, candidates)

        while len(trials) < runs:
            number = len(trials) + 1
            completed = [t for t in trials if t.status == "completed" and t.loss is not None]
            if len(completed) < initial_random:
                x = rng.random(len(self.parameters))
            else:
                x = self._suggest(completed, rng, candidates)
            params = self._decode(x)
            trial_seed = slimmc_seed + number - 1
            trial = self._evaluate(number, params, trial_seed)
            trials.append(trial)
            self._append_trial(trial)

        result = Result(tuple(trials), self)
        self._write_report(result, runs, initial_random, optimizer_seed, slimmc_seed, candidates)
        if not result.completed:
            errors = []
            for trial in result.trials:
                if trial.error and trial.error not in errors:
                    errors.append(trial.error)
            detail = "; ".join(errors[:3])
            raise RuntimeError("No completed trials" + (f": {detail}" if detail else ""))
        return result

    def _restore_rng(
        self,
        trials: Sequence[Trial],
        rng: np.random.Generator,
        initial_random: int,
        candidates: int,
    ) -> None:
        """Replay proposal generation so resume continues the exact optimizer stream."""
        prior: list[Trial] = []
        for trial in trials:
            completed = [t for t in prior if t.status == "completed" and t.loss is not None]
            if len(completed) < initial_random:
                rng.random(len(self.parameters))
            else:
                # _suggest() consumes exactly one candidate-pool draw from this RNG;
                # GP fitting itself is deterministic (random_state=0). Replaying the
                # draw restores the exact stream without refitting every historical GP.
                rng.random((candidates, len(self.parameters)))
            prior.append(trial)

    def _write_report(
        self,
        result: Result,
        runs: int,
        initial_random: int,
        optimizer_seed: int,
        slimmc_seed: int,
        candidates: int,
    ) -> None:
        completed = result.completed
        failed = tuple(t for t in result.trials if t.status != "completed")
        best = min(completed, key=lambda t: float(t.loss)) if completed else None

        if best is not None:
            best_trial_dir = best.run_path.parent
            best_model_source = best_trial_dir / self.model.name
            if best_model_source.is_file():
                shutil.copy2(best_model_source, self.output_dir / "best.model")

            (self.output_dir / "best_run.txt").write_text(
                f"trial\t{best.number}\n"
                f"loss\t{best.loss:.17g}\n"
                f"run_path\t{best.run_path}\n",
                encoding="utf-8",
            )

        lines = [
            "PYSLIMMC-OPT REPORT",
            "===================",
            "",
            "MODEL",
            "-----",
            str(self.model),
            "",
            "SEARCH SPACE",
            "------------",
        ]
        for p in self.parameters:
            scale = "log" if isinstance(p, LogFloat) else "linear"
            lines.append(f"{p.name}: {p.low:.12g} .. {p.high:.12g} ({scale})")

        lines.extend([
            "",
            "OPTIMIZATION",
            "------------",
            "method: Gaussian-process Bayesian optimization",
            "kernel: Matern(nu=2.5) + WhiteKernel",
            "acquisition: Expected Improvement",
            f"initial random trials: {initial_random}",
            f"requested trials: {runs}",
            f"completed trials: {len(completed)}",
            f"failed trials: {len(failed)}",
            f"candidate points per Bayesian step: {candidates}",
            f"optimizer seed: {optimizer_seed}",
            f"Slimmc seed base: {slimmc_seed}",
            f"trial timeout: {self.trial_timeout if self.trial_timeout is not None else 'none'}",
            "",
            "BEST OBSERVED TRIAL",
            "-------------------",
        ])
        if best is None:
            lines.append("none (no completed trials)")
        else:
            lines.extend([f"trial: {best.number}", f"loss: {best.loss:.12g}"])
            for p in self.parameters:
                lines.append(f"{p.name}: {best.params[p.name]:.12g}")
            lines.extend([
                f"Slimmc seed: {best.slimmc_seed}",
                f"run: {best.run_path}",
                "note: this is one stochastic observation; verify the selected point with Result.verify().",
            ])
        lines.extend([
            "",
            "SEARCH PROGRESS",
            "---------------",
            "trial\tstatus\tloss\tbest_loss",
        ])

        running_best = math.inf
        for trial in result.trials:
            loss_text = ""
            if trial.status == "completed" and trial.loss is not None:
                running_best = min(running_best, float(trial.loss))
                loss_text = f"{trial.loss:.12g}"
            best_text = "" if not math.isfinite(running_best) else f"{running_best:.12g}"
            lines.append(f"{trial.number}\t{trial.status}\t{loss_text}\t{best_text}")

        if failed:
            lines.extend(["", "FAILED TRIALS", "-------------"])
            for trial in failed:
                lines.append(f"{trial.number}: {trial.error}")

        if self._gp_fit_warnings:
            lines.extend([
                "",
                "GP DIAGNOSTICS",
                "--------------",
                "The Gaussian-process fit reached one or more kernel bounds. "
                "This can indicate weak structure or a noise-dominated objective.",
            ])
            for message in self._gp_fit_warnings:
                lines.append(f"- {message}")

        lines.extend([
            "",
            "FILES",
            "-----",
            "trials: trials.tsv",
            "study settings: study.json",
            "best model: best.model",
            "best run pointer: best_run.txt",
            "all runs: runs/",
            "",
        ])
        (self.output_dir / "report.txt").write_text("\n".join(lines), encoding="utf-8")

    def _suggest(self, completed: Sequence[Trial], rng: np.random.Generator, candidates: int) -> np.ndarray:
        x_train = np.array([self._encode(t.params) for t in completed], dtype=float)
        y_train = np.array([t.loss for t in completed], dtype=float)
        kernel = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(length_scale=np.ones(len(self.parameters)), nu=2.5) + WhiteKernel(1e-6, (1e-10, 1e-1))
        gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=2, random_state=0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", ConvergenceWarning)
            gp.fit(x_train, y_train)
        for item in caught:
            if issubclass(item.category, ConvergenceWarning):
                message = str(item.message)
                if message not in self._gp_fit_warnings:
                    self._gp_fit_warnings.append(message)
        pool = rng.random((candidates, len(self.parameters)))
        mean, std = gp.predict(pool, return_std=True)
        best = float(np.min(y_train))
        improvement = best - mean
        z = np.divide(improvement, std, out=np.zeros_like(improvement), where=std > 1e-15)
        ei = improvement * ndtr(z) + std * np.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi)
        ei[std <= 1e-15] = 0.0
        return pool[int(np.argmax(ei))]

    def _evaluate(self, number: int, params: dict[str, float], slimmc_seed: int) -> Trial:
        trial_dir = self.output_dir / "runs" / f"trial_{number:04d}"
        trial_dir.mkdir(parents=True, exist_ok=False)
        model_path = trial_dir / self.model.name
        run_rel = "result"
        auto_params = {name: value for name, value in params.items() if name in self._var_kinds}
        text = _render_model(self._template, self._var_kinds, auto_params, run_rel, slimmc_seed=slimmc_seed)
        if self.build is not None:
            text = self.build(text, params)
            if not isinstance(text, str):
                raise TypeError("build must return model text as str")
        model_path.write_text(text, encoding="utf-8")
        try:
            proc = subprocess.run(
                [self.slimmc, model_path.name], cwd=trial_dir, capture_output=True,
                text=True, check=False, timeout=self.trial_timeout
            )
            (trial_dir / "stdout.txt").write_text(proc.stdout, encoding="utf-8")
            (trial_dir / "stderr.txt").write_text(proc.stderr, encoding="utf-8")
            if proc.returncode != 0:
                raise RuntimeError(f"slimmc exited with code {proc.returncode}")
            import pyslimmc
            run_path = trial_dir / run_rel
            run = pyslimmc.open(run_path)
            loss = float(self.objective(run, params))
            if not math.isfinite(loss):
                raise ValueError("objective must return a finite number")
            return Trial(number, params, loss, "completed", run_path, slimmc_seed)
        except subprocess.TimeoutExpired as exc:
            return Trial(number, params, None, "timeout", trial_dir / run_rel, slimmc_seed, f"TimeoutExpired: {exc}")
        except Exception as exc:
            return Trial(number, params, None, "failed", trial_dir / run_rel, slimmc_seed, f"{type(exc).__name__}: {exc}")

    def _verify(self, params: Mapping[str, float], *, repeats: int, slimmc_seed: int) -> dict[str, object]:
        if repeats < 2:
            raise ValueError("repeats must be >= 2")
        verify_dir = self.output_dir / "verification"
        if verify_dir.exists():
            shutil.rmtree(verify_dir)
        (verify_dir / "runs").mkdir(parents=True)
        original_output = self.output_dir
        self.output_dir = verify_dir
        trials: list[Trial] = []
        try:
            for i in range(repeats):
                trial = self._evaluate(i + 1, dict(params), slimmc_seed + i)
                trials.append(trial)
        finally:
            self.output_dir = original_output
        completed = [t for t in trials if t.status == "completed" and t.loss is not None]
        if not completed:
            raise RuntimeError("All verification runs failed")
        values = np.array([float(t.loss) for t in completed])
        fields = ["repeat", "status", "slimmc_seed", "loss", "run_path", "error"]
        with (verify_dir / "verification.tsv").open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for t in trials:
                writer.writerow({"repeat": t.number, "status": t.status, "slimmc_seed": t.slimmc_seed, "loss": "" if t.loss is None else repr(t.loss), "run_path": str(t.run_path), "error": t.error})
        summary = {
            "repeats": repeats,
            "completed": len(completed),
            "mean_loss": float(values.mean()),
            "std_loss": float(values.std(ddof=1)) if len(values) > 1 else 0.0,
            "min_loss": float(values.min()),
            "max_loss": float(values.max()),
            "trials": tuple(trials),
        }
        (verify_dir / "report.txt").write_text(
            "PYSLIMMC-OPT VERIFICATION\n"
            "===========================\n\n"
            f"repeats: {repeats}\n"
            f"completed: {len(completed)}\n"
            f"Slimmc seed base: {slimmc_seed}\n"
            f"mean loss: {summary['mean_loss']:.12g}\n"
            f"standard deviation: {summary['std_loss']:.12g}\n"
            f"minimum loss: {summary['min_loss']:.12g}\n"
            f"maximum loss: {summary['max_loss']:.12g}\n",
            encoding="utf-8",
        )
        return summary

    def _decode(self, x: np.ndarray) -> dict[str, float]:
        out: dict[str, float] = {}
        for p, v in zip(self.parameters, x):
            if isinstance(p, LogFloat):
                out[p.name] = float(p.low * (p.high / p.low) ** float(v))
            else:
                out[p.name] = float(p.low + float(v) * (p.high - p.low))
        return out

    def _encode(self, params: Mapping[str, float]) -> list[float]:
        out: list[float] = []
        for p in self.parameters:
            value = float(params[p.name])
            if isinstance(p, LogFloat):
                out.append(math.log(value / p.low) / math.log(p.high / p.low))
            else:
                out.append((value - p.low) / (p.high - p.low))
        return out

    @property
    def _trials_path(self) -> Path:
        return self.output_dir / "trials.tsv"

    def _append_trial(self, trial: Trial) -> None:
        exists = self._trials_path.exists()
        fields = ["trial", "status", "slimmc_seed", *[p.name for p in self.parameters], "loss", "run_path", "error"]
        with self._trials_path.open("a", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            if not exists:
                writer.writeheader()
            row = {"trial": trial.number, "status": trial.status, "slimmc_seed": "" if trial.slimmc_seed is None else trial.slimmc_seed, "loss": "" if trial.loss is None else repr(trial.loss), "run_path": str(trial.run_path), "error": trial.error}
            row.update({k: repr(v) for k, v in trial.params.items()})
            writer.writerow(row)

    def _load_trials(self, slimmc_seed_base: int) -> list[Trial]:
        if not self._trials_path.exists():
            return []
        out: list[Trial] = []
        with self._trials_path.open(encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                params = {p.name: float(row[p.name]) for p in self.parameters}
                loss = float(row["loss"]) if row["loss"] else None
                number = int(row["trial"])
                seed_text = row.get("slimmc_seed", "")
                trial_seed = int(seed_text) if seed_text else slimmc_seed_base + number - 1
                out.append(Trial(number, params, loss, row["status"], Path(row["run_path"]), trial_seed, row.get("error", "")))
        return out

    def _write_study_json(
        self,
        runs: int,
        initial_random: int,
        optimizer_seed: int,
        slimmc_seed: int,
        candidates: int,
    ) -> None:
        path = self.output_dir / "study.json"
        data = {
            "model": str(self.model),
            "parameters": [
                {"name": p.name, "low": p.low, "high": p.high, "scale": "log" if isinstance(p, LogFloat) else "linear"}
                for p in self.parameters
            ],
            "trial_timeout": self.trial_timeout,
            "runs": runs,
            "initial_random": initial_random,
            "optimizer_seed": optimizer_seed,
            "slimmc_seed": slimmc_seed,
            "candidates": candidates,
        }
        if path.exists():
            previous = json.loads(path.read_text(encoding="utf-8"))
            stable = ("model", "parameters", "trial_timeout", "initial_random", "optimizer_seed", "slimmc_seed", "candidates")
            changed = [name for name in stable if previous.get(name) != data[name]]
            if changed:
                raise ValueError(
                    "cannot resume with changed study settings: " + ", ".join(changed)
                )
            if runs > int(previous.get("runs", 0)):
                previous["runs"] = runs
                path.write_text(json.dumps(previous, indent=2) + "\n", encoding="utf-8")
            return
        path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def _parse_var_kinds(text: str) -> dict[str, str]:
    found: dict[str, str] = {}
    for line in text.splitlines():
        m = re.match(r"^\s*var\s+(rate|param|species|monomer|endgroup)\s+([A-Za-z_][A-Za-z0-9_]*)\s+\S+\s*(?:#.*)?$", line)
        if m:
            found[m.group(2)] = m.group(1)
    return found


def _matching_declaration_count(template: str, kind: str, name: str) -> int:
    return sum(_replace_value(line, kind, name, 0.0) != line for line in template.splitlines())


def _render_model(template: str, kinds: Mapping[str, str], params: Mapping[str, float], output_dir: str, *, slimmc_seed: int | None = None) -> str:
    lines = template.splitlines()
    lines = _replace_or_insert_output_dir(lines, output_dir)
    if slimmc_seed is not None:
        lines = _replace_or_insert_seed(lines, slimmc_seed)
    for name, value in params.items():
        kind = kinds[name]
        replaced = [_replace_value(line, kind, name, value) for line in lines]
        count = sum(a != b for a, b in zip(lines, replaced))
        if count != 1:
            raise ValueError(
                f"optimization variable {name!r} expected exactly one matching {kind} declaration; found {count}"
            )
        lines = replaced
    return "\n".join(lines) + "\n"


def _replace_output_dir(line: str, output_dir: str) -> str:
    if re.match(r"^\s*param\s+output_dir\s+", line):
        indent = line[: len(line) - len(line.lstrip())]
        return f'{indent}param output_dir "{output_dir}"'
    return line


def _replace_or_insert_output_dir(lines: list[str], output_dir: str) -> list[str]:
    if any(re.match(r"^\s*param\s+output_dir\s+", line) for line in lines):
        return [_replace_output_dir(line, output_dir) for line in lines]
    insert_at = 0
    for i, line in enumerate(lines):
        if re.match(r"^\s*desc\b", line):
            insert_at = i + 1
            break
    lines.insert(insert_at, f'param output_dir "{output_dir}"')
    return lines


def _replace_value(line: str, kind: str, name: str, value: float) -> str:
    number = format(value, ".17g")
    patterns = {
        "param": rf"^(\s*param\s+{re.escape(name)}\s+)\S+(.*)$",
        "species": rf"^(\s*species\s+{re.escape(name)}\s+)\S+(.*)$",
        "monomer": rf"^(\s*monomer\s+{re.escape(name)}\s+)\S+(.*)$",
        "endgroup": rf"^(\s*endgroup\s+{re.escape(name)}\s+)\S+(.*)$",
        "rate": rf"^(\s*rate\s+{re.escape(name)}\s+(?:const\s+)?)\S+(.*)$",
    }
    m = re.match(patterns[kind], line)
    return f"{m.group(1)}{number}{m.group(2)}" if m else line


def _replace_or_insert_seed(lines: list[str], seed: int) -> list[str]:
    replacement = f"param seed {seed}"
    for i, line in enumerate(lines):
        if re.match(r"^\s*param\s+seed\s+", line):
            indent = line[: len(line) - len(line.lstrip())]
            lines[i] = indent + replacement
            return lines
    insert_at = 0
    for i, line in enumerate(lines):
        if re.match(r"^\s*(?:desc\b|param\s+output_dir\b)", line):
            insert_at = i + 1
    lines.insert(insert_at, replacement)
    return lines
