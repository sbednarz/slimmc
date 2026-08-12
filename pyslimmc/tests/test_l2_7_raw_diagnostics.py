from __future__ import annotations
import json, tempfile
from pathlib import Path
import numpy as np
import pyslimmc


def wt(root: Path, name: str, cols: dict[str, np.ndarray]) -> None:
    d = root / name
    d.mkdir(parents=True, exist_ok=True)
    for key, value in cols.items():
        np.save(d / f"{key}.npy", value, allow_pickle=False)


def build(root: Path) -> None:
    meta = {
        "run_id": "run", "storage": "slimmc-storage",
        "storage_format_version": "1.2.0", "run_status": "completed",
        "validation_status": "passed", "validation_warning_count": 0,
        "validation_error_count": 0, "engine": "slimmc-copo",
        "kinetic_model": "copo",
    }
    (root / "run_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    schema = [
        {"record_type": "schema_header", "schema_name": "slimmc-storage", "schema_version": "1.2.0"},
        {"record_type": "table", "name": "snapshots", "required": True},
        {"record_type": "table", "name": "state", "required": True},
        {"record_type": "table", "name": "memory", "required": False},
        {"record_type": "column", "table": "snapshots", "name": "time", "unit": "s"},
        {"record_type": "column", "table": "memory", "name": "total_est_B", "unit": "B"},
        {"record_type": "dictionary_entry", "dictionary": "state_entities", "id": 0, "name": "monomer_A", "kind": "monomer"},
    ]
    (root / "schema.jsonl").write_text("".join(json.dumps(x, separators=(",", ":")) + "\n" for x in schema))
    wt(root, "snapshots", {
        "snapshot_id": np.arange(3, dtype="<u8"),
        "time": np.array([0., 1., 2.]),
        "kmc_event": np.array([0, 10, 20], dtype="<u8"),
        "snapshot_reason_id": np.array([0, 1, 4], dtype="<u4"),
        "is_final": np.array([0, 0, 1], dtype=bool),
        "has_chains": np.array([0, 1, 1], dtype=bool),
        "has_sequences": np.zeros(3, dtype=bool),
        "kinetic_parameter_set_id": np.zeros(3, dtype="<u8"),
    })
    wt(root, "state", {
        "snapshot_id": np.arange(3, dtype="<u8"),
        "entity_id": np.zeros(3, dtype="<u4"),
        "count": np.array([100, 80, 60], dtype="<u8"),
        "moles": np.array([100., 80., 60.]),
        "concentration": np.array([100., 80., 60.]),
    })
    wt(root, "memory", {
        "snapshot_id": np.array([1, 2], dtype="<u8"),
        "live_chains": np.array([5, 7], dtype="<u8"),
        "total_est_B": np.array([1024, 2048], dtype="<u8"),
    })
    diagnostics = root / "diagnostics"
    diagnostics.mkdir()
    checks = [
        {"check": "required_files", "status": "pass", "severity": "error"},
        {"check": "optional_note", "status": "warn", "severity": "warning", "message": "example"},
    ]
    (diagnostics / "validation.jsonl").write_text("".join(json.dumps(x) + "\n" for x in checks))
    (diagnostics / "run.log").write_text("start\nfinish\n")
    (root / "RESULTS_COMPLETE").write_text("slimmc-storage-v1\n")


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        root = Path(td) / "run"
        root.mkdir()
        build(root)
        run = pyslimmc.open(root)
        assert run.raw.metadata["run_id"] == "run"
        assert run.raw.table("state") is run.table("state")
        assert run.raw.dictionary("state_entities")[0]["name"] == "monomer_A"
        assert len(run.raw.schema) == 7
        val = run.diagnostics.validation
        assert val.status == "passed" and val.error_count == 0
        assert len(val) == 2 and val["required_files"].passed
        assert val[1].check == "optional_note" and len(val.failed) == 1
        mem = run.diagnostics.memory
        assert mem.columns == ("live_chains", "total_est_B")
        assert np.isnan(mem.total_est_B[0]) and np.array_equal(np.asarray(mem.total_est_B)[1:], [1024., 2048.])
        assert run.column_unit("memory", "total_est_B") == "B"
        assert run.diagnostics.run_log.lines == ("start", "finish")
        assert not run.diagnostics.debug_log.exists
    print("pyslimmc L2.7 raw/diagnostics: PASS")


def test_script_contract() -> None:
    main()

if __name__ == "__main__":
    main()
