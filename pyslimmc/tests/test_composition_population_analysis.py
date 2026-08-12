from __future__ import annotations
import tempfile
from pathlib import Path
import numpy as np
import pyslimmc
from pyslimmc.tests.test_l2_4_chains import build


def _run(full=False):
    td = tempfile.TemporaryDirectory()
    root = Path(td.name) / 'run'; root.mkdir(); build(root, full=full)
    return td, pyslimmc.open(root)


def test_composition_filters_work_in_composition_and_full_modes():
    for full in (False, True):
        td, run = _run(full)
        try:
            c = run.final.chains
            assert c.where_count('A', min=2).count.tolist() == [2, 1]
            assert c.where_count('B', max=1).count.tolist() == [2, 1]
            assert c.where_fraction('A', min=0.6).count.tolist() == [2, 1]
            assert c.component_count.tolist() == [2, 2, 1]
            assert c.where_component_count(min=2, max=2).count.tolist() == [2, 1]
            assert c.where_components(('A', 'B')).count.tolist() == [2, 1]
            assert c.where_components(('B',), exact=True).count.tolist() == [4]
            assert c.where_components(('B',), exact=False).count.tolist() == [2, 1, 4]
            assert c.where_fraction('A', min=0.6).mwd(method='sticks').mn > 0
        finally:
            td.cleanup()


def test_composition_by_dp_is_multiplicity_weighted_and_readonly():
    td, run = _run()
    try:
        result = run.final.chains.composition_by_dp(bins=[1.5, 2.5, 3.5])
        assert result.chain_count.tolist() == [4.0, 3.0]
        assert result.record_count.tolist() == [1, 2]
        assert result.mean['A'][0] == 0.0
        assert result.mean['A'][1] == 2/3
        assert not result.mean['A'].flags.writeable
        ax = result.plot()
        assert ax.get_xlabel() == 'DP'
    finally:
        td.cleanup()


def test_composition_maps_and_component_classes():
    td, run = _run()
    try:
        chains = run.final.chains
        result = chains.composition_dp_map('A', dp_bins=[1.5,2.5,3.5], fraction_bins=[0,0.5,1.0])
        assert np.sum(result.values) == 7
        assert result.values.shape == (2,2)
        assert result.plot().get_xlabel() == 'DP'
        classes = chains.component_classes()
        assert classes.labels == ('B','AB')
        assert classes.chain_count.tolist() == [4.0,3.0]
        assert np.isclose(np.sum(classes.number_fraction),1.0)
        assert classes.plot().get_xlabel() == 'components'
        xy = chains.composition_map('A','B',bins=[0,0.5,1.0])
        assert np.sum(xy.values) == 7
        mm = chains.composition_mass_map('A', mass_bins=3, fraction_bins=4)
        assert np.sum(mm.values) == 7
    finally:
        td.cleanup()

def test_run_plot_composition_shortcuts():
    td, run = _run()
    try:
        assert run.plot.composition_by_dp(bins=[1.5,2.5,3.5]).get_xlabel() == 'DP'
        assert run.plot.composition_dp_map('A', dp_bins=2, fraction_bins=2).get_xlabel() == 'DP'
        assert run.plot.component_classes().get_xlabel() == 'components'
    finally:
        td.cleanup()

def test_run_plot_general_namespace_shortcuts():
    td, run = _run(full=True)
    try:
        assert run.plot.conversion().get_ylabel() == 'conversion'
        assert run.plot.concentrations(entities=('A','B')).get_ylabel() == 'concentration (mol/L)'
        assert run.plot.counts(entities=('A',)).get_ylabel() == 'count'
        assert run.plot.chain_mass_spectrum().get_xlabel()
        assert run.plot.chain_counts().get_xlabel() == 'DP'
    finally:
        td.cleanup()
