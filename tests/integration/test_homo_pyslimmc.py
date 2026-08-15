from __future__ import annotations

import numpy as np
import pytest


def test_homo_completed(homo_run):
    assert homo_run.status == "completed" and homo_run.output_status.complete


def test_homo_metadata(homo_run):
    assert homo_run.desc == "integration homo"
    assert homo_run.storage_format_version == "1.2.0"


def test_homo_snapshot_axes(homo_run):
    assert len(homo_run.t) >= 3
    assert np.all(np.diff(homo_run.t) >= 0)
    assert np.all(np.diff(homo_run.event) >= 0)
    assert np.array_equal(homo_run.sid, np.arange(len(homo_run.sid)))


def test_homo_no_artificial_initial_point(homo_run):
    assert homo_run.t[0] > 0
    assert homo_run.conv.total.shape == homo_run.t.shape


def test_homo_snapshot_selection(homo_run):
    assert homo_run.first.id == 0
    assert homo_run.last.id == homo_run.final.id
    assert homo_run.at_snapshot(homo_run.last.id).id == homo_run.last.id


def test_homo_time_selection(homo_run):
    middle = float(homo_run.t[len(homo_run.t) // 2])
    assert homo_run.at_time(middle, method="nearest").t == pytest.approx(middle)


def test_homo_state_series(homo_run):
    assert "M" in homo_run.conc.names
    assert homo_run.conc["M"].shape == homo_run.t.shape
    assert np.all(homo_run.count["M"] >= 0)


def test_homo_initial_conditions(homo_run):
    assert homo_run.c0["M"] == pytest.approx(0.050)
    assert homo_run.moles0["M"] == pytest.approx(0.0025)


def test_homo_conversion(homo_run):
    x = homo_run.conv["M"]
    assert np.all(np.isfinite(x)) and np.all((x >= 0) & (x <= 1))
    assert np.array_equal(x, homo_run.conv.total)


def test_homo_temperature_action(homo_run):
    assert homo_run.temp[0] == pytest.approx(343.15)
    assert homo_run.temp[-1] == pytest.approx(353.15)
    assert any(action.type == "set_temp" for action in homo_run.actions)


def test_homo_feed_and_volume(homo_run):
    feed = homo_run.feeds["MFEED"]
    assert feed.events.time.size == 1
    assert feed.volume_cum == pytest.approx(0.001)
    assert homo_run.volume[-1] == pytest.approx(homo_run.volume[0] + 0.001)


def test_homo_moments(homo_run):
    finite = np.isfinite(homo_run.mn) & np.isfinite(homo_run.mw)
    assert finite.any() and np.all(homo_run.mw[finite] >= homo_run.mn[finite])
    assert np.all(homo_run.dispersity[finite] >= 1.0)


def test_homo_chain_counts(homo_run):
    assert np.all(homo_run.chain_count.total == homo_run.chain_count.live + homo_run.chain_count.dead)
    assert homo_run.final.chains.total_chains > 0


def test_homo_distributions(homo_run):
    cld = homo_run.final.cld()
    mwd = homo_run.final.mwd()
    assert cld.x.size and mwd.x.size
    assert np.all(cld.y >= 0) and np.all(mwd.y >= 0)


def test_homo_channels(homo_run):
    name = next(name for name in homo_run.event_counts.names if "prop" in name)
    assert np.all(np.diff(homo_run.event_counts[name]) >= 0)


def test_homo_validation_and_mass_audit(homo_run):
    assert homo_run.validate(strict=True).is_valid
    audit = homo_run.mass_audit()
    assert audit.ok and audit.checked_chains > 0
