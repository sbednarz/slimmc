from dataclasses import dataclass
from pathlib import Path

import pytest

from pyslimmc.run import Variable
from pyslimmc.runs import Runs, SelectionError


@dataclass(frozen=True)
class FakeRun:
    run_id: str
    path: Path
    var: dict[str, Variable]


def fake(run_id: str, ia: float | None, temp: float | None, replicate: int = 0) -> FakeRun:
    variables: dict[str, Variable] = {}
    if ia is not None:
        variables["IA"] = Variable("monomer", "IA", ia, "mol_L")
    if temp is not None:
        variables["temperature"] = Variable("param", "temperature", temp, "K")
    variables["replicate"] = Variable("param", "replicate", replicate, "1")
    return FakeRun(run_id, Path("/tmp") / run_id, variables)


def collection(*items: FakeRun) -> Runs:
    return Runs(Path("/tmp"), {str(item.path): item for item in items})


def test_multidimensional_sweep_returns_runs_in_numeric_lexicographic_order():
    runs = collection(
        fake("z", 0.2, 353.15),
        fake("a", 0.1, 353.15),
        fake("m", 0.2, 343.15),
        fake("b", 0.1, 343.15),
    )
    sweep = runs.sweep("IA", "temperature")
    assert isinstance(sweep, Runs)
    assert sweep.sweep_variables == ("IA", "temperature")
    assert [(r.var["IA"].value, r.var["temperature"].value) for r in sweep] == [
        (0.1, 343.15),
        (0.1, 353.15),
        (0.2, 343.15),
        (0.2, 353.15),
    ]


def test_duplicate_parameter_points_are_preserved_and_reported():
    runs = collection(
        fake("rep_b", 0.1, 343.15, 2),
        fake("rep_a", 0.1, 343.15, 1),
        fake("other", 0.2, 343.15, 1),
    )
    sweep = runs.sweep("IA", "temperature")
    assert [r.run_id for r in sweep] == ["rep_a", "rep_b", "other"]
    summary = sweep._sweep_summary()
    assert summary["duplicate_runs"] == 1
    assert summary["missing_points"] == 0


def test_incomplete_grid_is_reported():
    sweep = collection(
        fake("a", 0.1, 343.15),
        fake("b", 0.1, 353.15),
        fake("c", 0.2, 343.15),
    ).sweep("IA", "temperature")
    summary = sweep._sweep_summary()
    assert summary["complete"] is False
    assert summary["missing_points"] == 1


def test_missing_variable_is_an_error_and_no_runs_are_silently_dropped():
    runs = collection(fake("ok", 0.1, 343.15), fake("bad", 0.2, None))
    with pytest.raises(SelectionError, match="temperature.*missing in 1 run"):
        runs.sweep("IA", "temperature")


def test_sweep_requires_unique_nonempty_variable_names():
    runs = collection(fake("a", 0.1, 343.15))
    with pytest.raises(TypeError):
        runs.sweep()
    with pytest.raises(ValueError):
        runs.sweep("IA", "IA")
    with pytest.raises(TypeError):
        runs.sweep("")


def test_sweep_metadata_and_order_survive_slice_match_filter_and_pack():
    sweep = collection(
        fake("case_02", 0.2, 343.15),
        fake("case_01", 0.1, 343.15),
        fake("case_03", 0.3, 343.15),
    ).sweep("IA", "temperature")
    assert [r.run_id for r in sweep[1:]] == ["case_02", "case_03"]
    assert sweep[1:].sweep_variables == ("IA", "temperature")
    matched = sweep.match("case_0[12]")
    assert [r.run_id for r in matched] == ["case_01", "case_02"]
    filtered = sweep.filter(var_name="temperature", var_value=343.15)
    assert [r.run_id for r in filtered] == ["case_01", "case_02", "case_03"]
    packed = sweep.pack(key="case_*")
    assert list(packed) == ["01", "02", "03"]
