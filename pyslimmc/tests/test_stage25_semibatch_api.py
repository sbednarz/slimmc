from __future__ import annotations
from pathlib import Path
import numpy as np
import pytest
import pyslimmc
from pyslimmc import AnalysisNotApplicableError

RUN = Path(__file__).resolve().parent / 'data' / 'feedh_results'

def test_volume_and_kmc_volume_are_readonly():
    run = pyslimmc.open(RUN)
    np.testing.assert_allclose(run.volume, [0.101]*4)
    np.testing.assert_allclose(run.kmc_volume, [1.01e-18]*4)
    assert not run.volume.flags.writeable
    assert not run.kmc_volume.flags.writeable

def test_initial_state_and_physical_moles():
    run = pyslimmc.open(RUN)
    assert run.c0['M'] == 0.1
    assert run.moles0['M'] == pytest.approx(0.01)
    np.testing.assert_allclose(run.moles['M'], run.conc['M'] * run.volume)
    assert run.state.moles['M'][0] != pytest.approx(run.moles['M'][0])

def test_feeds():
    run = pyslimmc.open(RUN)
    feed = run.feeds['F']
    assert feed.concentration == {'M': 1.0, 'Q': 0.2}
    assert feed.fraction == {'M': 1.0}
    np.testing.assert_allclose(feed.events.dose, [0.001])
    assert feed.volume_cum == pytest.approx(0.001)
    assert feed.moles_cum['M'] == pytest.approx(0.001)

def test_set_c_invalidates_balance_only_for_target():
    run = pyslimmc.open(RUN)
    with pytest.raises(AnalysisNotApplicableError, match='set_c'):
        _ = run.balance.total['M']
