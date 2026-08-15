from __future__ import annotations

import importlib.util
import inspect
import tempfile
from pathlib import Path

import numpy as np
import pytest

import pyslimmc
from pyslimmc.tests.test_l2_5_moments_spectra import build


REMOVED_RUN_NAMES = (
    "history", "oligomers", "engine_mass_audit", "parameter_states",
    "snapshot", "conversion", "concentration", "concentrations",
    "temperature", "T", "runinfo",
)


def _open_fixture():
    td = tempfile.TemporaryDirectory()
    root = Path(td.name) / "run"
    root.mkdir()
    build(root)
    return td, root, pyslimmc.open(root)


def test_public_run_contract_is_storage_only():
    td, root, run = _open_fixture()
    try:
        assert isinstance(run, pyslimmc.Run)
        assert run.metadata.storage == "slimmc-storage"
        for name in REMOVED_RUN_NAMES:
            assert not hasattr(run, name), name
        assert hasattr(run, "at_snapshot")
        assert hasattr(run, "last") and hasattr(run, "final")
        assert "conc" in dir(run) and "conv" in dir(run) and "temp" in dir(run)
        assert inspect.signature(pyslimmc.open).parameters.keys() == {"path", "allow_incomplete"}
        assert importlib.util.find_spec("pyslimmc.homo") is None
        assert importlib.util.find_spec("pyslimmc.copo") is None
    finally:
        td.cleanup()


def test_scan_recognizes_only_storage_runs(tmp_path: Path):
    good = tmp_path / "good"
    good.mkdir()
    build(good)
    old = tmp_path / "old"
    old.mkdir()
    (old / "unrelated.json").write_text("{}")
    runs = pyslimmc.scan(tmp_path)
    assert len(runs) == 1
    assert runs[0].path == good


def test_chain_names_are_canonical():
    td, root, run = _open_fixture()
    try:
        chains = run.final.chains
        assert np.array_equal(chains.count, np.array([3, 2], dtype=np.uint64))
        assert not hasattr(chains, "abundance")
        assert hasattr(chains, "counts")  # monomer-composition mapping, not multiplicity
        assert hasattr(chains, "population_activity")
        assert hasattr(chains, "pool")
        assert not hasattr(chains, "population")
    finally:
        td.cleanup()


def test_mwd_contract_and_single_export(tmp_path: Path):
    td, root, run = _open_fixture()
    try:
        snap = run.final
        for kwargs in ({"method":"hist"}, {"basis":"mass"}, {"coordinate":"log10"},
                       {"output":"density"}, {"bins":10}, {"bin_width":0.1},
                       {"sigma":0.2}, {"grid_step":0.01}, {"normalization":"per_series"}):
            with pytest.raises(TypeError):
                snap.mwd(**kwargs)

        dist = snap.mwd(form="mass")
        assert np.isclose(dist.y.sum(), 1.0)
        assert hasattr(dist, "to_tsv")
        assert not hasattr(dist, "to_csv")
        assert not hasattr(dist, "to_gnuplot")
        out = dist.to_tsv(tmp_path / "mwd.tsv")
        assert out.is_file()
        header = next(line for line in out.read_text().splitlines() if not line.startswith("#"))
        assert "\t" in header
    finally:
        td.cleanup()


def test_root_exports_are_frozen():
    expected = {
        "__version__", "open", "scan", "help", "Run", "Variable", "Variables", "Runs",
        "SelectionError", "Report", "report", "PlotStyle", "available_styles",
        "get_style", "figure_size", "MassAuditResult", "PyslimmcError",
        "FeatureUnavailableError", "ChemicalAnalysisNotApplicableError", "AnalysisNotApplicableError",
        "ChemicalModelIncompatibleError", "DataUnavailableError",
        "IncompleteSequenceDataError", "InvalidOutputError", "ValidationFailedError",
        "NumericalAnalysisError", "DataConsistencyError", "UnknownColumnError",
        "UnknownMonomerError", "UnsupportedChainSchema",
        "SnapshotUnavailableError", "FinalSnapshotUnavailableError",
        "MassModelUnavailableError", "InvalidDistributionConfigurationError",
    }
    assert set(pyslimmc.__all__) == expected
