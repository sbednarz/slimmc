from __future__ import annotations

import os
import signal
import shutil
import subprocess
import time
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
MODELS = Path(__file__).with_name("models")


def _slimmc() -> Path:
    configured = os.environ.get("SLIMMC_CLI")
    path = Path(configured).expanduser() if configured else ROOT / "bin" / "slimmc"
    path = path.resolve()
    if not path.is_file():
        pytest.fail(f"Slimmc CLI is missing: {path}; run `make build` first")
    return path


@pytest.fixture(scope="session")
def slimmc_cli() -> Path:
    return _slimmc()


def _run_fixture(slimmc: Path, model_name: str, root: Path) -> Path:
    work = root / model_name.removesuffix(".model")
    work.mkdir()
    model = work / model_name
    shutil.copy2(MODELS / model_name, model)
    checked = subprocess.run(
        [str(slimmc), "--check", model.name], cwd=work, text=True,
        capture_output=True, check=False,
    )
    assert checked.returncode == 0, checked.stdout + checked.stderr
    completed = subprocess.run(
        [str(slimmc), model.name], cwd=work, text=True,
        capture_output=True, check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    result = work / "result"
    assert result.is_dir()
    return result


def _run_incomplete_fixture(
    slimmc: Path,
    model_name: str,
    root: Path,
    *,
    interrupt: bool = False,
) -> Path:
    work = root / model_name.removesuffix(".model")
    work.mkdir()
    model = work / model_name
    shutil.copy2(MODELS / model_name, model)
    if interrupt:
        process = subprocess.Popen(
            [str(slimmc), model.name],
            cwd=work,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        metadata = work / "result" / "run_metadata.json"
        deadline = time.monotonic() + 10.0
        while time.monotonic() < deadline and not metadata.is_file() and process.poll() is None:
            time.sleep(0.01)
        assert process.poll() is None, "fixture ended before it could be interrupted"
        process.send_signal(signal.SIGINT)
        stdout, stderr = process.communicate(timeout=20)
        assert process.returncode == 0, stdout + stderr
    else:
        completed = subprocess.run(
            [str(slimmc), model.name], cwd=work, text=True,
            capture_output=True, check=False,
        )
        assert completed.returncode != 0, "failure fixture unexpectedly completed"
    result = work / "result"
    assert (result / "run_metadata.json").is_file()
    assert not (result / "RESULTS_COMPLETE").exists()
    return result


@pytest.fixture(scope="session")
def integration_root(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("slimmc-integration")


@pytest.fixture(scope="session")
def homo_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_fixture(slimmc_cli, "homo.model", integration_root)


@pytest.fixture(scope="session")
def copo_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_fixture(slimmc_cli, "copo.model", integration_root)


@pytest.fixture(scope="session")
def homo_mechanisms_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_fixture(slimmc_cli, "homo_mechanisms.model", integration_root)


@pytest.fixture(scope="session")
def homo_stop_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_fixture(slimmc_cli, "homo_stop.model", integration_root)


@pytest.fixture(scope="session")
def copo_mechanisms_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_fixture(slimmc_cli, "copo_mechanisms.model", integration_root)


@pytest.fixture(scope="session")
def copo_penultimate_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_fixture(slimmc_cli, "copo_penultimate.model", integration_root)


@pytest.fixture(scope="session")
def copo_terpolymer_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_fixture(slimmc_cli, "copo_terpolymer.model", integration_root)


@pytest.fixture(scope="session")
def copo_composition_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_fixture(slimmc_cli, "copo_composition.model", integration_root)


@pytest.fixture(scope="session")
def failed_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_incomplete_fixture(slimmc_cli, "homo_failed.model", integration_root)


@pytest.fixture(scope="session")
def interrupted_path(slimmc_cli: Path, integration_root: Path) -> Path:
    return _run_incomplete_fixture(
        slimmc_cli, "homo_interrupted.model", integration_root, interrupt=True,
    )


@pytest.fixture(scope="session")
def homo_run(homo_path: Path):
    import pyslimmc
    return pyslimmc.open(homo_path)


@pytest.fixture(scope="session")
def copo_run(copo_path: Path):
    import pyslimmc
    return pyslimmc.open(copo_path)


@pytest.fixture(scope="session")
def homo_mechanisms_run(homo_mechanisms_path: Path):
    import pyslimmc
    return pyslimmc.open(homo_mechanisms_path)


@pytest.fixture(scope="session")
def homo_stop_run(homo_stop_path: Path):
    import pyslimmc
    return pyslimmc.open(homo_stop_path)


@pytest.fixture(scope="session")
def copo_mechanisms_run(copo_mechanisms_path: Path):
    import pyslimmc
    return pyslimmc.open(copo_mechanisms_path)


@pytest.fixture(scope="session")
def copo_penultimate_run(copo_penultimate_path: Path):
    import pyslimmc
    return pyslimmc.open(copo_penultimate_path)


@pytest.fixture(scope="session")
def copo_terpolymer_run(copo_terpolymer_path: Path):
    import pyslimmc
    return pyslimmc.open(copo_terpolymer_path)


@pytest.fixture(scope="session")
def copo_composition_run(copo_composition_path: Path):
    import pyslimmc
    return pyslimmc.open(copo_composition_path)


@pytest.fixture(scope="session")
def failed_run(failed_path: Path):
    import pyslimmc
    return pyslimmc.open(failed_path, allow_incomplete=True)


@pytest.fixture(scope="session")
def interrupted_run(interrupted_path: Path):
    import pyslimmc
    return pyslimmc.open(interrupted_path, allow_incomplete=True)
