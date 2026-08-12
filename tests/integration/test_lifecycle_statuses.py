from __future__ import annotations

import pytest


def test_failed_run_is_real_incomplete_storage(failed_run):
    assert failed_run.status == "failed"
    assert not failed_run.output_status.complete
    assert failed_run.termination_reason == "runtime_error"


def test_interrupted_run_is_real_incomplete_storage(interrupted_run):
    assert interrupted_run.status == "interrupted"
    assert interrupted_run.termination_reason == "user_interrupt"
    assert not interrupted_run.output_status.complete


@pytest.mark.parametrize("fixture_name", ["failed_path", "interrupted_path"])
def test_incomplete_run_requires_explicit_opt_in(request, fixture_name):
    import pyslimmc

    path = request.getfixturevalue(fixture_name)
    with pytest.raises(pyslimmc.InvalidOutputError):
        pyslimmc.open(path)


def test_real_status_namespaces(
    integration_root,
    homo_path,
    failed_path,
    interrupted_path,
):
    import pyslimmc

    runs = pyslimmc.scan(integration_root)
    assert homo_path in {run.path for run in runs.completed}
    assert failed_path in {run.path for run in runs.failed}
    assert interrupted_path in {run.path for run in runs.interrupted}
    assert runs.one(status="failed").path == failed_path
