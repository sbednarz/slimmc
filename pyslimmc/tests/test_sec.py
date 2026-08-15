from __future__ import annotations

import numpy as np
import pytest

import pyslimmc
from pyslimmc.tests.test_distribution_math import (
    _build_oracle_storage,
    EXPECTED_MN,
    EXPECTED_MW,
    EXPECTED_MZ,
)


def _run(tmp_path):
    root = tmp_path / "run"
    root.mkdir()
    _build_oracle_storage(root)
    return pyslimmc.open(root)


def _area(sec):
    return float(np.trapezoid(sec.y, sec.x))


def test_sec_is_continuous_unit_mass_density(tmp_path):
    run = _run(tmp_path)
    sec = run.sec(pool="dead", sigma_log10M=0.05, mass_model="repeat_units")

    assert sec.x.ndim == 1
    assert sec.y.shape == sec.x.shape
    assert sec.x.size > 100
    assert np.all(np.diff(sec.x) > 0)
    assert np.all(sec.y >= 0)
    assert _area(sec) == pytest.approx(1.0, rel=1e-8, abs=1e-8)
    assert sec.metadata["representation"] == "continuous"
    assert sec.metadata["ordinate"] == "dW_app/dlog10M"


def test_sec_uses_mass_weighted_exact_population(tmp_path):
    run = _run(tmp_path)
    sigma = 0.01
    sec = run.sec(pool="dead", sigma_log10M=sigma, mass_model="repeat_units")

    # For well-separated narrow kernels, integrated peak areas recover the
    # exact mass fractions 200, 600, 250, 1600 / 2650.
    support = np.log10(np.array([100.0, 200.0, 250.0, 400.0]))
    expected = np.array([200.0, 600.0, 250.0, 1600.0]) / 2650.0
    boundaries = (support[:-1] + support[1:]) / 2.0
    edges = np.r_[-np.inf, boundaries, np.inf]
    areas = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (sec.x >= lo) & (sec.x < hi)
        areas.append(float(np.trapezoid(sec.y[mask], sec.x[mask])))
    np.testing.assert_allclose(areas, expected, atol=2e-3, rtol=0)


def test_sec_moments_remain_exact_source_moments(tmp_path):
    run = _run(tmp_path)
    for sigma in (0.01, 0.05, 0.2):
        sec = run.sec(pool="dead", sigma_log10M=sigma, mass_model="repeat_units")
        assert sec.mn == pytest.approx(EXPECTED_MN)
        assert sec.mw == pytest.approx(EXPECTED_MW)
        assert sec.mz == pytest.approx(EXPECTED_MZ)


def test_sec_mass_and_log10_mass_are_consistent(tmp_path):
    run = _run(tmp_path)
    sec = run.sec(pool="dead", sigma_log10M=0.05, mass_model="repeat_units")
    np.testing.assert_allclose(np.log10(sec.mass), sec.x)
    np.testing.assert_array_equal(sec.log10_mass, sec.x)


def test_sec_available_on_run_snapshot_and_population(tmp_path):
    run = _run(tmp_path)
    a = run.sec(pool="dead", sigma_log10M=0.05, mass_model="repeat_units")
    b = run.final.sec(pool="dead", sigma_log10M=0.05, mass_model="repeat_units")
    c = run.final.chains.dead.sec(sigma_log10M=0.05, mass_model="repeat_units")
    np.testing.assert_allclose(a.x, b.x)
    np.testing.assert_allclose(a.y, b.y)
    np.testing.assert_allclose(a.x, c.x)
    np.testing.assert_allclose(a.y, c.y)


def test_sec_requires_explicit_positive_sigma(tmp_path):
    run = _run(tmp_path)
    with pytest.raises(TypeError):
        run.sec(pool="dead")
    for value in (0.0, -0.1, np.nan, np.inf):
        with pytest.raises(ValueError):
            run.sec(pool="dead", sigma_log10M=value)


def test_sec_step_is_numerical_only(tmp_path):
    run = _run(tmp_path)
    fine = run.sec(pool="dead", sigma_log10M=0.05, step_log10M=0.001)
    coarse = run.sec(pool="dead", sigma_log10M=0.05, step_log10M=0.005)
    assert fine.x.size > coarse.x.size
    assert _area(fine) == pytest.approx(1.0, abs=1e-8)
    assert _area(coarse) == pytest.approx(1.0, abs=1e-8)
    assert fine.mn == coarse.mn
    assert fine.mw == coarse.mw
    assert fine.mz == coarse.mz


def test_sec_matches_direct_buback_gaussian_mixture(tmp_path):
    run = _run(tmp_path)
    sigma = 0.05
    sec = run.sec(pool="dead", sigma_log10M=sigma, mass_model="repeat_units")

    mass = np.array([100.0, 200.0, 250.0, 400.0])
    weights = np.array([200.0, 600.0, 250.0, 1600.0]) / 2650.0
    support = np.log10(mass)
    delta = (sec.x[:, None] - support[None, :]) / sigma
    expected = (
        np.exp(-0.5 * delta * delta) / (sigma * np.sqrt(2.0 * np.pi))
    ) @ weights

    np.testing.assert_allclose(sec.y, expected, rtol=5e-14, atol=5e-14)


def test_sigma_log10M_is_buback_b_sigma_v_parameter(tmp_path):
    run = _run(tmp_path)
    sec = run.sec(pool="dead", sigma_log10M=0.05, mass_model="repeat_units")
    assert sec.metadata["sigma_log10M"] == pytest.approx(0.05)
    assert sec.metadata["b_sigma_v_equivalent"] == pytest.approx(0.05)


def test_sec_does_not_hide_coarse_grid_by_renormalizing(tmp_path):
    run = _run(tmp_path)
    with pytest.raises(ValueError, match="grid is too coarse"):
        run.sec(
            pool="dead",
            sigma_log10M=0.01,
            step_log10M=0.2,
            mass_model="repeat_units",
        )
