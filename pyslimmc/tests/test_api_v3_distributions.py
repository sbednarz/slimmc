from pathlib import Path
import tempfile
import numpy as np
import pytest
import pyslimmc
from pyslimmc.tests.test_l2_5_moments_spectra import build


def _run():
    td = tempfile.TemporaryDirectory()
    root = Path(td.name) / "run"
    root.mkdir()
    build(root)
    return td, pyslimmc.open(root)


def test_snapshot_contract_first_last_final():
    td, run = _run()
    try:
        assert run.first.id == run.snapshots.first.id
        assert run.last.id == run.snapshots.last.id
        assert run.final.id == run.snapshots.final.id
        assert not hasattr(run, "initial")
    finally:
        td.cleanup()


def test_distribution_help_before_call(capsys):
    td, run = _run()
    try:
        assert "Molar-mass distribution" in run.mwd.help()
        assert "Chain-length distribution" in run.cld.help()
        assert not hasattr(run, "chain_mass_spectrum")
        out = capsys.readouterr().out
        assert "exact and discrete" in out
    finally:
        td.cleanup()


def test_new_distribution_properties_and_removed_legacy_api():
    td, run = _run()
    try:
        mwd = run.mwd(form="log")
        assert np.allclose(mwd.x, np.log10(mwd.mass))
        assert isinstance(mwd.mn, float)
        assert isinstance(mwd.dispersity, float)
        assert mwd.form == "log"
        assert mwd.metadata["representation"] == "discrete"
        assert not hasattr(mwd, "pdi")
        assert not hasattr(mwd, "method")
        assert not hasattr(mwd, "basis")
        assert not hasattr(mwd, "coordinate")
        assert not hasattr(mwd, "output")

        cld = run.cld(form="number")
        assert np.array_equal(cld.x, cld.dp)
        assert np.isclose(cld.y.sum(), 1.0)
        assert hasattr(cld, "dpn") and hasattr(cld, "dpw") and hasattr(cld, "dpz")
        assert not hasattr(cld, "mn") and not hasattr(cld, "mw") and not hasattr(cld, "mz")

        for kwargs in (
            {"method": "hist"},
            {"basis": "mass"},
            {"coordinate": "log10"},
            {"output": "density"},
            {"bin_width": 0.1},
            {"sigma": 0.05},
            {"grid_step": 0.01},
            {"normalization": "per_series"},
        ):
            with pytest.raises(TypeError):
                run.mwd(**kwargs)
    finally:
        td.cleanup()
