from __future__ import annotations

import numpy as np
import pytest


def test_copo_completed(copo_run):
    assert copo_run.status == "completed" and copo_run.output_status.complete


def test_copo_metadata(copo_run):
    assert copo_run.desc == "integration copo"
    assert copo_run.storage_format_version == "1.2.0"


def test_copo_snapshot_alignment(copo_run):
    n = len(copo_run.t)
    assert n >= 3 and copo_run.event.shape == (n,) and copo_run.sid.shape == (n,)
    assert np.all(np.diff(copo_run.t) >= 0)


def test_copo_monomer_names(copo_run):
    assert copo_run.f.names == ("A", "B")
    assert copo_run.F.cum.names == ("A", "B")


@pytest.mark.parametrize("monomer", ["A", "B"])
def test_copo_conversions(copo_run, monomer: str):
    values = copo_run.conv[monomer]
    assert values.shape == copo_run.t.shape
    assert np.all(np.isfinite(values)) and np.all((values >= 0) & (values <= 1))


def test_copo_total_conversion(copo_run):
    assert np.all((copo_run.conv.total >= 0) & (copo_run.conv.total <= 1))
    assert copo_run.conv.total[-1] > 0


def test_copo_remaining_composition(copo_run):
    assert np.allclose(copo_run.f["A"] + copo_run.f["B"], 1.0)


def test_copo_cumulative_composition(copo_run):
    a, b = copo_run.F.cum["A"], copo_run.F.cum["B"]
    mask = np.isfinite(a) & np.isfinite(b)
    assert mask.any() and np.allclose(a[mask] + b[mask], 1.0)


def test_copo_instantaneous_composition(copo_run):
    a, b = copo_run.F.ins["A"], copo_run.F.ins["B"]
    mask = np.isfinite(a) & np.isfinite(b)
    assert mask.any() and np.allclose(a[mask] + b[mask], 1.0)


def test_copo_feed_balance(copo_run):
    feed = copo_run.feeds["BFEED"]
    assert feed.events.time.size == 1
    assert feed.moles_cum["B"] == pytest.approx(0.0002)
    assert copo_run.volume[-1] == pytest.approx(copo_run.volume[0] + 0.001)


def test_copo_temperature(copo_run):
    assert copo_run.temp[0] == pytest.approx(343.15)
    assert copo_run.temp[-1] == pytest.approx(353.15)


def test_copo_chain_composition(copo_run):
    chains = copo_run.final.chains
    assert chains.total_chains > 0
    assert np.array_equal(chains.composition.counts.total, chains.dp)


def test_copo_full_sequences(copo_run):
    chains = copo_run.final.chains
    assert chains.has_sequences
    assert all(len(sequence) == int(dp) for sequence, dp in zip(chains.sequences, chains.dp))


def test_copo_live_dead_partition(copo_run):
    chains = copo_run.final.chains
    assert chains.live.total_chains + chains.dead.total_chains == chains.total_chains


def test_copo_moments(copo_run):
    finite = np.isfinite(copo_run.mn) & np.isfinite(copo_run.mw) & np.isfinite(copo_run.mz)
    assert finite.any()
    assert np.all(copo_run.mw[finite] >= copo_run.mn[finite])
    assert np.all(copo_run.mz[finite] >= copo_run.mw[finite])


def test_copo_distributions(copo_run):
    assert copo_run.final.cld().x.size > 0
    assert copo_run.final.mwd().x.size > 0
    assert copo_run.final.mass_counts().mass.size > 0


def test_copo_channels(copo_run):
    names = copo_run.event_counts.names
    assert {"prop_PA_A", "prop_PA_B", "prop_PB_A", "prop_PB_B"}.issubset(names)
    assert np.all(copo_run.channels.interval_event_counts() >= 0)


def test_copo_kinetics(copo_run):
    assert copo_run.k["kp_aa"].shape == copo_run.t.shape
    assert np.all(copo_run.k["kp_aa"] > 0)


def test_copo_microstructure(copo_run):
    micro = copo_run.microstructure
    assert micro is not None
    assert copo_run.final.chains.has_sequences


def test_copo_validation(copo_run):
    assert copo_run.validate(strict=True).is_valid


def test_runs_scan_and_status_namespaces(homo_path, copo_path, integration_root):
    import pyslimmc
    runs = pyslimmc.scan(integration_root)
    assert len(runs.completed) >= 2
    assert {"homo", "copo"}.issubset({run.kinetic_model for run in runs.completed})
