from __future__ import annotations
import json
from pathlib import Path
import tempfile
import numpy as np
import pyslimmc


def write_table(root: Path, name: str, columns: dict[str, np.ndarray]) -> None:
    d = root / name
    d.mkdir(parents=True)
    for key, value in columns.items():
        np.save(d / f"{key}.npy", value, allow_pickle=False)


def build(root: Path) -> None:
    meta = {
        "run_id": root.name,
        "storage": "slimmc-storage",
        "storage_format_version": "1.2.0",
        "run_status": "completed",
        "validation_error_count": 0,
        "engine": "slimmc-copo",
        "kinetic_model": "copo",
    }
    (root / "run_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    records = [
        {"record_type":"schema_header","schema_name":"slimmc-storage","schema_version":"1.2.0"},
        {"record_type":"table","name":"snapshots","required":True},
        {"record_type":"table","name":"state","required":True},
        {"record_type":"column","table":"snapshots","name":"time","unit":"s"},
        {"record_type":"column","table":"state","name":"count","unit":"1"},
        {"record_type":"column","table":"state","name":"moles","unit":"mol"},
        {"record_type":"column","table":"state","name":"concentration","unit":"mol/L"},
        {"record_type":"dictionary_entry","dictionary":"snapshot_reasons","id":0,"name":"initial"},
        {"record_type":"dictionary_entry","dictionary":"snapshot_reasons","id":4,"name":"final"},
        {"record_type":"dictionary_entry","dictionary":"state_entities","id":0,"name":"A","kind":"monomer"},
        {"record_type":"dictionary_entry","dictionary":"state_entities","id":1,"name":"live_chains","kind":"aggregate"},
    ]
    (root / "schema.jsonl").write_text("".join(json.dumps(x,separators=(",",":"))+"\n" for x in records))
    write_table(root, "snapshots", {
        "snapshot_id": np.array([0,1,2], dtype="<u8"),
        "time": np.array([0.0,0.5,1.0], dtype="<f8"),
        "kmc_event": np.array([0,5,12], dtype="<u8"),
        "snapshot_reason_id": np.array([0,1,4], dtype="<u4"),
        "is_final": np.array([False,False,True], dtype=np.bool_),
        "has_chains": np.array([False,False,False], dtype=np.bool_),
        "has_sequences": np.array([False,False,False], dtype=np.bool_),
        "kinetic_parameter_set_id": np.array([0,0,0], dtype="<u8"),
    })
    write_table(root, "state", {
        "snapshot_id": np.repeat(np.array([0,1,2], dtype="<u8"), 2),
        "entity_id": np.tile(np.array([0,1], dtype="<u4"), 3),
        "count": np.array([100,0,70,3,40,7], dtype="<u8"),
        "moles": np.array([1.0,0.0,0.7,0.03,0.4,0.07], dtype="<f8"),
        "concentration": np.array([10.0,0.0,7.0,0.3,4.0,0.7], dtype="<f8"),
    })
    (root / "RESULTS_COMPLETE").write_text("slimmc-storage-v1\n")


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        root = Path(td) / "run_000001"; root.mkdir(); build(root)
        run = pyslimmc.open(root)
        assert isinstance(run.t, np.ndarray)
        assert not run.t.flags.writeable
        assert np.array_equal(run.t, [0.0,0.5,1.0])
        assert np.allclose(run.t / 3600.0, np.array([0.0,0.5,1.0]) / 3600.0)
        assert np.array_equal(run.event, [0,5,12])
        assert np.array_equal(run.sid, [0,1,2])
        assert run.conc.names == ("A", "live_chains")
        assert isinstance(run.count["A"], np.ndarray)
        assert not run.count["A"].flags.writeable
        assert np.array_equal(run.count["A"], [100,70,40])
        assert np.allclose(run.state.moles["live_chains"], [0.0,0.03,0.07])
        assert np.allclose(run.conc["A"], [10.0,7.0,4.0])
        assert np.allclose(run.conc["A"] * 1000.0, [10000.0,7000.0,4000.0])
        assert run.state["A"].count[1] == 70
        assert run.state["A"].conc[-1] == 4.0
        snap = run.last
        assert snap.count["A"] == 40
        assert snap.state.moles["live_chains"] == 0.07
        assert snap.conc["A"] == 4.0
        assert snap.state["A"].count == 40
        assert snap.state["A"].conc == 4.0
        assert np.array_equal(run.conc.matrix, [[10.0,0.0],[7.0,0.3],[4.0,0.7]])
        try:
            run.conc["missing"]
        except KeyError:
            pass
        else:
            raise AssertionError("unknown state entity accepted")
    print("pyslimmc L2.2 state series: PASS")



def test_script_contract() -> None:
    main()

if __name__ == "__main__":
    main()
