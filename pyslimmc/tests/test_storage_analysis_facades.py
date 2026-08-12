from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

import pyslimmc
from pyslimmc.core import DataUnavailableError
from pyslimmc.storage_analysis import StorageMicrostructure
from pyslimmc.tests.test_l2_6_channels_kinetics_actions import build as build_channels


def test_storage_firings_facade_from_canonical_tables():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td) / "run"
        root.mkdir()
        build_channels(root)
        run = pyslimmc.open(root)
        f = run.firings
        assert f.channels() == ["prop_AA", "prop_AB"]
        assert f.final_fires() == {"prop_AA": 12, "prop_AB": 8}
        assert f.total_fires() == 20
        assert f.delta_fires("prop_AA") == 5
        assert np.array_equal(f.delta_fires_series("prop_AB"), [3, 5])
        assert f.fire_shares() == pytest.approx({"prop_AA": 0.6, "prop_AB": 0.4})
        assert f.validate() is True
        assert "total fires: 20" in f.info_text()


class _Record:
    def __init__(self, sequence, count=1, rid=0):
        self._sequence = sequence
        self.count = count
        self.chain_record_id = rid

    @property
    def sequence(self):
        if isinstance(self._sequence, Exception):
            raise self._sequence
        return self._sequence


class _Chains(list):
    has_sequences = True


class _Snap:
    def __init__(self, sid, chains):
        self.id = sid
        self.chains = chains


class _FakeRun:
    def __init__(self, *, bad_sequence=False):
        records = _Chains([
            _Record(["A", "A", "B"], count=2, rid=1),
            _Record(ValueError("bad sequence") if bad_sequence else ["A", "B"], count=1, rid=2),
        ])
        self.final = _Snap(2, records)
        self.last = self.final
        self.monomer_names = ("A", "B")
        self.tables = {"microstructure_motifs", "block_statistics"}
        self._tables = {
            "microstructure_motifs": {
                "snapshot_id": np.array([2, 2, 2], dtype=np.uint64),
                "motif_order": np.array([2, 2, 3], dtype=np.uint32),
                "motif_id": np.array([0, 1, 0], dtype=np.uint32),
                "count": np.array([2, 3, 2], dtype=np.uint64),
            },
            "block_statistics": {
                "snapshot_id": np.array([2, 2], dtype=np.uint64),
                "monomer_id": np.array([0, 1], dtype=np.uint32),
                "block_length": np.array([2, 1], dtype=np.uint64),
                "block_count": np.array([2, 3], dtype=np.uint64),
            },
        }

    def table(self, name):
        return self._tables[name]

    def dictionary(self, name):
        if name == "microstructure_dyads":
            return {0: {"name": "A|A"}, 1: {"name": "A|B"}}
        if name == "microstructure_triads":
            return {0: {"name": "A|A|B"}}
        raise KeyError(name)


def test_storage_microstructure_engine_and_sequence_paths_agree():
    m = StorageMicrostructure(_FakeRun())
    dyads = {row["motif"]: row["count"] for row in m.dyads().rows()}
    assert dyads == {"A|A": 2.0, "A|B": 3.0}
    seq_dyads = {row["motif"]: row["count"] for row in m.dyads(source="sequences").rows()}
    assert seq_dyads == dyads
    assert m.transition_fraction() == pytest.approx(3 / 5)
    assert m.homodyad_fraction() == pytest.approx(2 / 5)
    assert m.blockiness() == pytest.approx({"homodyad_fraction": 2 / 5, "transition_fraction": 3 / 5})
    assert m.check_sequence_consistency() == {"dyads_match": True, "triads_match": True}
    rows = m.run_lengths("A").rows()
    assert len(rows) == 1 and rows[0]["run_length"] == 2


def test_storage_microstructure_unreadable_sequence_is_not_silently_dropped():
    m = StorageMicrostructure(_FakeRun(bad_sequence=True))
    with pytest.raises(DataUnavailableError, match="chain record 2"):
        m.dyads(source="sequences")
