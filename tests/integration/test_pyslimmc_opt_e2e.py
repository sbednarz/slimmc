from __future__ import annotations

import json
from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def optimization(tmp_path_factory, slimmc_cli):
    import pyslimmc
    from pyslimmc_opt import Float, Study

    root = tmp_path_factory.mktemp("slimmc-opt-integration")
    source = Path(__file__).with_name("models") / "opt.model"
    model = root / "opt.model"
    model.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    study = Study(
        model=model,
        parameters=[Float("t_end", 0.15, 0.30)],
        objective=lambda run, params: -float(run.conv.total[-1]),
        output_dir=root / "study",
        slimmc=slimmc_cli,
    )
    result = study.optimize(runs=3, initial_random=2, optimizer_seed=81, slimmc_seed=7400, candidates=100)
    return study, result, pyslimmc


@pytest.fixture(scope="module")
def optimization_surface(tmp_path_factory, slimmc_cli):
    from pyslimmc_opt import Float, Study

    root = tmp_path_factory.mktemp("slimmc-opt-surface")
    source = Path(__file__).with_name("models") / "opt_surface.model"
    model = root / "opt_surface.model"
    model.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")

    def build(model_text, params):
        return model_text.replace("species R 1.0e-5", f"species R {params['radical']:.17g}")

    study = Study(
        model=model,
        parameters=[
            Float("kp", 1.0e5, 3.0e5),
            Float("M", 0.03, 0.07),
            Float("radical", 5.0e-6, 2.0e-5),
        ],
        objective=lambda run, params: -float(run.conv.total[-1]),
        build=build,
        output_dir=root / "study",
        slimmc=slimmc_cli,
    )
    result = study.optimize(
        runs=1,
        initial_random=1,
        optimizer_seed=82,
        slimmc_seed=7500,
        candidates=100,
    )
    return study, result


def test_opt_real_trials_complete(optimization):
    _, result, _ = optimization
    assert len(result.trials) == 3 and len(result.completed) == 3


def test_opt_trial_seeds(optimization):
    _, result, _ = optimization
    assert [trial.slimmc_seed for trial in result.trials] == [7400, 7401, 7402]


def test_opt_runs_are_real_storage(optimization):
    _, result, pyslimmc = optimization
    for trial in result.trials:
        run = pyslimmc.open(trial.run_path)
        assert run.status == "completed" and run.kinetic_model == "homo"


def test_opt_best_model_is_rerunnable(optimization, slimmc_cli, tmp_path):
    import subprocess
    study, result, _ = optimization
    best = study.output_dir / "best.model"
    text = best.read_text(encoding="utf-8").replace('param output_dir "result"', 'param output_dir "rerun"')
    model = tmp_path / "best.model"
    model.write_text(text, encoding="utf-8")
    proc = subprocess.run([str(slimmc_cli), model.name], cwd=tmp_path, text=True, capture_output=True, check=False)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (tmp_path / "rerun" / "RESULTS_COMPLETE").is_file()


def test_opt_report_artifacts(optimization):
    study, result, _ = optimization
    assert (study.output_dir / "report.txt").is_file()
    assert (study.output_dir / "trials.tsv").is_file()
    assert (study.output_dir / "study.json").is_file()
    assert (study.output_dir / "best_run.txt").is_file()
    assert result.best_run.is_dir()


def test_opt_study_metadata(optimization):
    study, _, _ = optimization
    data = json.loads((study.output_dir / "study.json").read_text(encoding="utf-8"))
    assert data["optimizer_seed"] == 81 and data["slimmc_seed"] == 7400


def test_opt_resume_does_not_duplicate_trials(optimization):
    study, result, _ = optimization
    resumed = study.optimize(runs=3, initial_random=2, optimizer_seed=81, slimmc_seed=7400, candidates=100)
    assert len(resumed.trials) == len(result.trials) == 3
    assert resumed.best_loss == pytest.approx(result.best_loss)


def test_opt_fixed_rate_monomer_and_build_callback(optimization_surface):
    study, result = optimization_surface
    trial = result.best_trial
    rendered = (trial.run_path.parent / study.model.name).read_text(encoding="utf-8")
    assert f"rate kp const {trial.params['kp']:.17g}" in rendered
    assert f"monomer M {trial.params['M']:.17g} 100.0" in rendered
    assert f"species R {trial.params['radical']:.17g}" in rendered


def test_opt_verification_runs_real_independent_trajectories(optimization_surface):
    _, result = optimization_surface
    summary = result.verify(repeats=2, slimmc_seed=7600)
    assert summary["repeats"] == 2 and summary["completed"] == 2
    assert [trial.slimmc_seed for trial in summary["trials"]] == [7600, 7601]


def test_opt_failed_objective_becomes_failed_trial(tmp_path, slimmc_cli):
    from pyslimmc_opt import Float, Study

    source = Path(__file__).with_name("models") / "opt.model"
    model = tmp_path / "opt.model"
    model.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")

    def broken_objective(run, params):
        raise RuntimeError("intentional objective failure")

    study = Study(
        model=model,
        parameters=[Float("t_end", 0.15, 0.20)],
        objective=broken_objective,
        output_dir=tmp_path / "study",
        slimmc=slimmc_cli,
    )
    with pytest.raises(RuntimeError, match="No completed trials"):
        study.optimize(
            runs=1,
            initial_random=1,
            optimizer_seed=83,
            slimmc_seed=7700,
            candidates=100,
        )


def test_opt_resume_rejects_changed_seed(optimization):
    study, _, _ = optimization
    with pytest.raises(ValueError, match="changed study settings"):
        study.optimize(
            runs=4,
            initial_random=2,
            optimizer_seed=999,
            slimmc_seed=7400,
            candidates=100,
        )
