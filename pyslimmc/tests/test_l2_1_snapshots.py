from __future__ import annotations
import json
from pathlib import Path
import tempfile
import numpy as np
import pyslimmc
from pyslimmc.core import DataUnavailableError


def write_table(root: Path, name: str, columns: dict[str, np.ndarray]) -> None:
    d = root / name
    d.mkdir(parents=True)
    for key, value in columns.items():
        np.save(d / f"{key}.npy", value, allow_pickle=False)


def build(root: Path, status: str = "completed") -> None:
    metadata = {
        "run_id": root.name,
        "storage": "slimmc-storage",
        "storage_format_version": "1.2.0",
        "run_status": status,
        "validation_error_count": 0,
        "engine": "slimmc-copo",
        "kinetic_model": "copo",
    }
    (root / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    records = [
        {"record_type":"schema_header","schema_name":"slimmc-storage","schema_version":"1.2.0"},
        *({"record_type":"table","name":name,"required":True} for name in [
            "snapshots","state","chains","moments","channel_events","kinetic_parameters/values"
        ]),
        {"record_type":"dictionary_entry","dictionary":"snapshot_reasons","id":0,"name":"initial"},
        {"record_type":"dictionary_entry","dictionary":"snapshot_reasons","id":1,"name":"scheduled"},
        {"record_type":"dictionary_entry","dictionary":"snapshot_reasons","id":4,"name":"final"},
    ]
    (root / "schema.jsonl").write_text("".join(json.dumps(x,separators=(",",":"))+"\n" for x in records))
    write_table(root, "snapshots", {
        "snapshot_id": np.array([0,1,2], dtype="<u8"),
        "time": np.array([0.0, 1.0, 3.0], dtype="<f8"),
        "kmc_event": np.array([0,10,30], dtype="<u8"),
        "snapshot_reason_id": np.array([0,1,4], dtype="<u4"),
        "is_final": np.array([False,False,status=="completed"], dtype=np.bool_),
        "has_chains": np.array([False,True,True], dtype=np.bool_),
        "has_sequences": np.array([False,False,True], dtype=np.bool_),
        "kinetic_parameter_set_id": np.array([0,0,1], dtype="<u8"),
    })
    write_table(root, "state", {
        "snapshot_id": np.array([0,1,2], dtype="<u8"),
        "entity_id": np.array([0,0,0], dtype="<u4"),
        "count": np.array([100,80,50], dtype="<u8"),
    })
    write_table(root, "chains", {
        "chain_record_id": np.array([0,1], dtype="<u8"),
        "snapshot_id": np.array([1,2], dtype="<u8"),
        "dp": np.array([5,8], dtype="<u8"),
    })
    write_table(root, "moments", {
        "snapshot_id": np.array([1,2], dtype="<u8"),
        "mn": np.array([500.0,800.0], dtype="<f8"),
    })
    write_table(root, "channel_events", {
        "snapshot_id": np.array([0,1,2], dtype="<u8"),
        "channel_id": np.array([0,0,0], dtype="<u4"),
        "event_count": np.array([0,10,30], dtype="<u8"),
    })
    write_table(root, "kinetic_parameters/values", {
        "kinetic_parameter_set_id": np.array([0,1], dtype="<u8"),
        "kinetic_parameter_id": np.array([0,0], dtype="<u4"),
        "value": np.array([300.0,320.0], dtype="<f8"),
    })
    if status == "completed":
        (root / "RESULTS_COMPLETE").write_text("slimmc-storage-v1\n")


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        root=Path(td)/"run_000001"; root.mkdir(); build(root)
        run=pyslimmc.open(root)
        assert len(run.snapshots)==3
        assert run.first.id==0 and run.last.id==2 and run.final.id==2
        assert run.at_snapshot(1).t==1.0 and run.snapshots[1].event==10
        assert run.at_time(2.5).id==1
        assert run.at_time(2.5, method="after").id==2
        assert run.at_time(2.5, method="nearest").id==2
        assert run.at_event(29).id==1
        assert run.last.reason=="final" and run.last.has_sequences
        assert int(run.last.state.count[0])==50
        assert int(run.last.chains.dp[0])==8
        assert float(run.last.mn)==800.0
        assert int(run.last.channels.event_count[0])==30
        assert float(run.last.kinetics.value[0])==320.0
        try: run.first.chains
        except DataUnavailableError: pass
        else: raise AssertionError("chains unexpectedly available")

        root2=Path(td)/"run_000002"; root2.mkdir(); build(root2, "interrupted")
        interrupted=pyslimmc.open(root2, allow_incomplete=True)
        assert interrupted.last.id==2 and not interrupted.last.is_final
        try: interrupted.final
        except DataUnavailableError: pass
        else: raise AssertionError("interrupted run unexpectedly has final snapshot")
    print("pyslimmc L2.1 snapshots: PASS")


def test_script_contract() -> None:
    main()

if __name__ == "__main__": main()
