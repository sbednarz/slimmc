from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


def _call(cli: Path, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run([str(cli), *args], cwd=cwd, text=True, capture_output=True, check=False)


def test_cli_version(slimmc_cli: Path):
    proc = _call(slimmc_cli, "--version")
    version = (Path(__file__).resolve().parents[2] / "VERSION").read_text(encoding="utf-8").strip()
    assert proc.returncode == 0 and version in proc.stdout


@pytest.mark.parametrize("model_name", [
    "homo.model",
    "homo_mechanisms.model",
    "homo_stop.model",
    "homo_failed.model",
    "homo_interrupted.model",
    "copo.model",
    "copo_mechanisms.model",
    "copo_penultimate.model",
    "copo_terpolymer.model",
    "copo_composition.model",
    "opt_surface.model",
])
def test_cli_check_accepts_fixture_models(slimmc_cli: Path, model_name: str):
    model = Path(__file__).with_name("models") / model_name
    proc = _call(slimmc_cli, "--check", str(model))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "[run]" not in proc.stdout
    assert "[done]" not in proc.stdout



def test_homo_run_reports_start_and_done(slimmc_cli: Path, tmp_path: Path):
    source = Path(__file__).with_name("models") / "homo.model"
    model = tmp_path / "homo_status.model"
    model.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    proc = _call(slimmc_cli, model.name, cwd=tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    lines = [line for line in proc.stdout.splitlines() if line.strip()]
    assert lines and lines[0].startswith("[run] event=0 ")
    assert any(line.startswith("[run] ") and "100.0%" in line for line in lines)
    assert lines[-1].startswith("[done] ")
    assert "output=" in lines[-1]

def test_cli_dispatches_homo(homo_run):
    assert homo_run.engine == "slimmc" and homo_run.kinetic_model == "homo"


def test_cli_dispatches_copo(copo_run):
    assert copo_run.engine == "slimmc-copo" and copo_run.kinetic_model == "copo"


def test_cli_rejects_missing_model(slimmc_cli: Path, tmp_path: Path):
    proc = _call(slimmc_cli, "missing.model", cwd=tmp_path)
    assert proc.returncode != 0


def test_cli_rejects_invalid_model(slimmc_cli: Path, tmp_path: Path):
    model = tmp_path / "invalid.model"
    model.write_text('desc "invalid"\nmonomer A 1 100\nmonomer B 1 100\n', encoding="utf-8")
    proc = _call(slimmc_cli, "--check", model.name, cwd=tmp_path)
    assert proc.returncode != 0
