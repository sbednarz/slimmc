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


def test_analysis_help_before_call(capsys):
    td, run = _run()
    try:
        assert "Molar-mass distribution" in run.mwd.help()
        assert "Chain-length distribution" in run.cld.help()
        assert "Neutral chain-mass spectrum" in run.chain_mass_spectrum.help()
        assert "Molar-mass distribution" in capsys.readouterr().out
    finally:
        td.cleanup()


def test_distribution_v3_properties_and_info():
    td, run = _run()
    try:
        mwd = run.mwd(method="hist", coordinate="log10", bin_width=0.1)
        assert np.allclose(mwd.log10_x, np.log10(mwd.x))
        assert isinstance(mwd.mn, float)
        assert isinstance(mwd.dispersity, float)
        assert "Molar mass distribution" in mwd.info_text()
        assert "metadata" in mwd.as_dict()
        assert not hasattr(mwd, "pdi")

        cld = run.cld(method="sticks")
        assert cld.dp_n == cld.mn
        assert cld.dp_w == cld.mw
        assert "Chain length distribution" in cld.info_text()

        spectrum = run.chain_mass_spectrum(normalize="base_peak")
        assert np.array_equal(spectrum.mass, spectrum.x)
        assert spectrum.base_peak_intensity == pytest.approx(1.0)
        assert "not an m/z spectrum" in spectrum.info_text()
        assert not hasattr(run, "chain_spectrum")
    finally:
        td.cleanup()
