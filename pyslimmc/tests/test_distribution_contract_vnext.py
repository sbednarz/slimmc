from pathlib import Path
import numpy as np
import pytest
import pyslimmc
from pyslimmc.tests.test_distribution_math import _build_oracle_storage, EXPECTED_CLD_NUMBER, EXPECTED_CLD_MASS, EXPECTED_MASS, EXPECTED_MWD_MASS

@pytest.fixture
def run(tmp_path):
    root=tmp_path/'run'; root.mkdir(); _build_oracle_storage(root); return pyslimmc.open(root)

def test_cld(run):
    a=run.cld(pool='dead', weighting='number', mass_model='repeat_units')
    b=run.cld(pool='dead', weighting='mass', mass_model='repeat_units')
    np.testing.assert_allclose(a.y, EXPECTED_CLD_NUMBER)
    np.testing.assert_allclose(b.y, EXPECTED_CLD_MASS)
    assert a.metadata['representation']=='discrete'

def test_mass_distribution(run):
    d=run.mass_distribution(pool='dead', weighting='mass', mass_model='repeat_units')
    np.testing.assert_allclose(d.x, EXPECTED_MASS)
    np.testing.assert_allclose(d.y, EXPECTED_MWD_MASS)
    assert d.y.sum()==pytest.approx(1)
    assert d.metadata['representation']=='discrete'

def test_mwd_copo_exact_mass(run):
    d=run.mwd(pool='dead', mass_model='repeat_units')
    assert d.metadata['representation']=='density'
    assert d.metadata['source']=='mass_counts'
    assert d.metadata['zero_filled'] is False
    assert np.trapezoid(d.y,d.x)==pytest.approx(1, rel=1e-12)
    assert d.metadata['ordinate']=='dW/dlog10(M)'

def test_mwd_no_form(run):
    with pytest.raises(TypeError):
        run.mwd(pool='dead', form='log')

class _Mom:
    mn=1000.; mw=1200.; mz=1400.; dpn=10.; dpw=12.; dpz=14.

class _HomoRun:
    kinetic_model='homo'
    engine='slimmc-homo'
    dictionaries={'monomers': {0: {'name': 'A', 'molar_mass_increment': 100.0}}}

class FakeHomo:
    snapshot_id=1; t=1.0; run=_HomoRun()
    dp=np.array([10,12], dtype=int)
    count=np.array([100,1], dtype=int)
    counts={'A': np.array([10,12])}
    def masses(self, *, mass_model): return self.dp.astype(float)*100.0
    def moments(self, *, mass_model): return _Mom()

def test_homo_zero_fill_internal():
    from pyslimmc.distributions import build_mwd
    d=build_mwd(FakeHomo(), mass_model='repeat_units')
    assert d.metadata['zero_filled'] is True
    assert d.metadata['source']=='dp_counts'
    assert np.trapezoid(d.y,d.x)==pytest.approx(1)
    # DP 11 is explicitly zero before interpolation; on the regular output grid
    # the density therefore approaches zero at log10(1100).
    i=np.argmin(np.abs(d.x-np.log10(1100.0)))
    assert d.y[i] < 0.1*d.y.max()
