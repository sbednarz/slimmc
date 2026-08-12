from pathlib import Path
import numpy as np
import pytest
from pyslimmc_opt.study import Float, Result, Study, Trial, _parse_var_kinds, _render_model

TEMPLATE = '''desc "x"\nparam output_dir "results/main"\nparam temperature 300\nvar param temperature K\nvar monomer M mol_L\nmonomer M 1.0 100.0\n'''

def test_render_uses_var_contract():
    kinds = _parse_var_kinds(TEMPLATE)
    text = _render_model(TEMPLATE, kinds, {"temperature": 350.0, "M": 0.25}, "result", slimmc_seed=1234)
    assert 'param output_dir "result"' in text
    assert 'param temperature 350' in text
    assert 'monomer M 0.25 100.0' in text
    assert 'var param temperature K' in text
    assert 'param seed 1234' in text


def test_render_replaces_plain_and_const_fixed_rates():
    plain = "var rate kp 1\nrate kp 100\n"
    const = "var rate kp 1\nrate kp const 100\n"
    assert "rate kp 250" in _render_model(plain, {"kp": "rate"}, {"kp": 250}, "result")
    assert "rate kp const 250" in _render_model(const, {"kp": "rate"}, {"kp": 250}, "result")


def test_render_inserts_missing_output_dir():
    rendered = _render_model('desc "x"\nparam t_end 1\n', {}, {}, "result")
    assert rendered.splitlines()[1] == 'param output_dir "result"'

def test_float_validation():
    with pytest.raises(ValueError):
        Float("x", 1, 1)

def test_result_best():
    r = Result((Trial(1, {"x": 1.0}, 2.0, "completed", Path("a")), Trial(2, {"x": 2.0}, 1.0, "completed", Path("b"))))
    assert r.best_params == {"x": 2.0}
    assert r.best_loss == 1.0

def test_build_allows_non_var_parameter(tmp_path):
    model = tmp_path / "model.model"
    model.write_text(TEMPLATE + "feed F M 0.1\n", encoding="utf-8")

    def build(text, params):
        return text.replace("feed F M 0.1", f"feed F M {params['feed_c']}")

    study = Study(
        model=model,
        parameters=[Float("feed_c", 0.01, 0.2)],
        objective=lambda run, params: 0.0,
        build=build,
    )
    assert study.build is build


def test_report_and_best_model_are_written(tmp_path, monkeypatch):
    model = tmp_path / "base.model"
    model.write_text(
        'param output_dir "result"\n'
        'param x 0.5\n'
        'var param x 1\n',
        encoding="utf-8",
    )

    def fake_evaluate(self, number, params, slimmc_seed):
        trial_dir = self.output_dir / "runs" / f"trial_{number:04d}"
        trial_dir.mkdir(parents=True, exist_ok=False)
        rendered = trial_dir / self.model.name
        rendered.write_text(self._template, encoding="utf-8")
        run_path = trial_dir / "result"
        run_path.mkdir()
        return Trial(number, params, params["x"] ** 2, "completed", run_path, slimmc_seed)

    monkeypatch.setattr(Study, "_evaluate", fake_evaluate)
    out = tmp_path / "out"
    study = Study(
        model=model,
        parameters=[Float("x", 0.0, 1.0)],
        objective=lambda run, params: 0.0,
        output_dir=out,
    )
    result = study.optimize(runs=2, initial_random=2, optimizer_seed=1, slimmc_seed=500, candidates=100)

    report = (out / "report.txt").read_text(encoding="utf-8")
    assert "PYSLIMMC-OPT REPORT" in report
    assert "SEARCH PROGRESS" in report
    assert (out / "best.model").is_file()
    assert (out / "best_run.txt").is_file()
    assert [t.slimmc_seed for t in result.trials] == [500, 501]
    assert "optimizer seed: 1" in report
    assert "Slimmc seed base: 500" in report


def test_verify_uses_independent_slimmc_seeds(tmp_path, monkeypatch):
    model = tmp_path / "base.model"
    model.write_text('param output_dir "result"\nparam seed 1\nparam x 0.5\nvar param x 1\n', encoding="utf-8")

    def fake_evaluate(self, number, params, slimmc_seed):
        trial_dir = self.output_dir / "runs" / f"trial_{number:04d}"
        trial_dir.mkdir(parents=True, exist_ok=False)
        (trial_dir / self.model.name).write_text(self._template, encoding="utf-8")
        run_path = trial_dir / "result"
        run_path.mkdir()
        return Trial(number, params, float(slimmc_seed), "completed", run_path, slimmc_seed)

    monkeypatch.setattr(Study, "_evaluate", fake_evaluate)
    study = Study(model=model, parameters=[Float("x", 0.0, 1.0)], objective=lambda run, params: 0.0, output_dir=tmp_path / "out")
    result = study.optimize(runs=1, initial_random=1, optimizer_seed=2, slimmc_seed=100, candidates=100)
    summary = result.verify(repeats=3, slimmc_seed=900)
    assert [t.slimmc_seed for t in summary["trials"]] == [900, 901, 902]
    assert summary["mean_loss"] == pytest.approx(901.0)
    assert (tmp_path / "out" / "verification" / "verification.tsv").is_file()


def test_resume_replays_optimizer_rng(tmp_path, monkeypatch):
    model = tmp_path / "base.model"
    model.write_text('param output_dir "result"\nparam x 0.5\nvar param x 1\n', encoding="utf-8")
    seen = []

    def fake_evaluate(self, number, params, slimmc_seed):
        seen.append((number, params["x"], slimmc_seed))
        trial_dir = self.output_dir / "runs" / f"trial_{number:04d}"
        trial_dir.mkdir(parents=True, exist_ok=False)
        (trial_dir / self.model.name).write_text(self._template, encoding="utf-8")
        run_path = trial_dir / "result"
        run_path.mkdir()
        return Trial(number, params, params["x"], "completed", run_path, slimmc_seed)

    monkeypatch.setattr(Study, "_evaluate", fake_evaluate)
    out = tmp_path / "out"
    study = Study(model=model, parameters=[Float("x", 0.0, 1.0)], objective=lambda run, params: 0.0, output_dir=out)
    study.optimize(runs=2, initial_random=2, optimizer_seed=7, slimmc_seed=100, candidates=100)
    study.optimize(runs=3, initial_random=2, optimizer_seed=7, slimmc_seed=100, candidates=100)
    resumed_third = seen[-1][1]

    out2 = tmp_path / "out2"
    seen.clear()
    study2 = Study(model=model, parameters=[Float("x", 0.0, 1.0)], objective=lambda run, params: 0.0, output_dir=out2)
    study2.optimize(runs=3, initial_random=2, optimizer_seed=7, slimmc_seed=100, candidates=100)
    uninterrupted_third = seen[-1][1]
    assert resumed_third == pytest.approx(uninterrupted_third)


def test_missing_target_declaration_is_rejected(tmp_path):
    model = tmp_path / "model.model"
    model.write_text('var param temperature K\n', encoding='utf-8')
    with pytest.raises(ValueError, match='no matching model declaration'):
        Study(model=model, parameters=[Float('temperature', 300, 400)], objective=lambda run, params: 0.0)


def test_render_rejects_zero_substitutions():
    with pytest.raises(ValueError, match='expected exactly one'):
        _render_model('var param x 1\n', {'x': 'param'}, {'x': 2.0}, 'result')


def test_relative_slimmc_path_is_resolved(tmp_path):
    model = tmp_path / 'model.model'
    model.write_text('param x 0.5\nvar param x 1\n', encoding='utf-8')
    exe = tmp_path / 'slimmc'
    exe.write_text('', encoding='utf-8')
    study = Study(model=model, parameters=[Float('x', 0, 1)], objective=lambda run, params: 0.0, slimmc=str(exe))
    assert Path(study.slimmc).is_absolute()


def test_logfloat_roundtrip(tmp_path):
    from pyslimmc_opt import LogFloat
    model = tmp_path / 'model.model'
    model.write_text('param kp 1000\nvar param kp 1\n', encoding='utf-8')
    study = Study(model=model, parameters=[LogFloat('kp', 1e2, 1e5)], objective=lambda run, params: 0.0)
    decoded = study._decode(np.array([0.5]))
    assert decoded['kp'] == pytest.approx((1e2 * 1e5) ** 0.5)
    assert study._encode(decoded)[0] == pytest.approx(0.5)


def test_all_failed_trials_still_write_report(tmp_path, monkeypatch):
    model = tmp_path / 'model.model'
    model.write_text('param x 0.5\nvar param x 1\n', encoding='utf-8')
    def fake_evaluate(self, number, params, slimmc_seed):
        trial_dir = self.output_dir / 'runs' / f'trial_{number:04d}'
        trial_dir.mkdir(parents=True, exist_ok=False)
        return Trial(number, params, None, 'failed', trial_dir / 'result', slimmc_seed, 'FileNotFoundError: slimmc')
    monkeypatch.setattr(Study, '_evaluate', fake_evaluate)
    out = tmp_path / 'out'
    study = Study(model=model, parameters=[Float('x', 0, 1)], objective=lambda run, params: 0.0, output_dir=out)
    with pytest.raises(RuntimeError, match='FileNotFoundError'):
        study.optimize(runs=2, initial_random=2, candidates=100)
    report = (out / 'report.txt').read_text(encoding='utf-8')
    assert 'no completed trials' in report
    assert 'FileNotFoundError: slimmc' in report


def test_gp_convergence_warnings_are_captured_for_report(tmp_path, monkeypatch):
    from sklearn.exceptions import ConvergenceWarning
    import warnings

    model = tmp_path / 'model.model'
    model.write_text('param x 0.5\nvar param x 1\n', encoding='utf-8')
    study = Study(model=model, parameters=[Float('x', 0, 1)], objective=lambda run, params: 0.0, output_dir=tmp_path / 'out')

    class FakeGP:
        def __init__(self, *args, **kwargs): pass
        def fit(self, x, y):
            warnings.warn('length scale reached lower bound', ConvergenceWarning)
            return self
        def predict(self, pool, return_std=False):
            return np.zeros(len(pool)), np.ones(len(pool))

    monkeypatch.setattr('pyslimmc_opt.study.GaussianProcessRegressor', FakeGP)
    completed = [
        Trial(1, {'x': 0.2}, 1.0, 'completed', Path('a')),
        Trial(2, {'x': 0.8}, 2.0, 'completed', Path('b')),
    ]
    with warnings.catch_warnings(record=True) as outward:
        warnings.simplefilter('always')
        study._suggest(completed, np.random.default_rng(1), 100)
    assert not [w for w in outward if issubclass(w.category, ConvergenceWarning)]
    assert study._gp_fit_warnings == ['length scale reached lower bound']
