from __future__ import annotations
import json
from pathlib import Path
import tempfile
import numpy as np
import pyslimmc

NA = 6.02214076e23


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
    }
    (root / "run_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    (root / "input.model").write_text(
        'desc "regression"\nparam kmc_volume 1e-20\nmonomer IA 1.0 130.10\nmonomer BMA 1.0 142.20\n'
    )
    schema = [{"record_type":"schema_header","schema_name":"slimmc-storage","schema_version":"1.2.0"}]
    for t in ("snapshots", "state", "chains", "chain_composition"):
        schema.append({"record_type":"table","name":t,"required":True})
    for i, m in enumerate(("IA", "BMA")):
        schema.append({"record_type":"dictionary_entry","dictionary":"state_entities","id":i,"name":m,"kind":"monomer"})
        schema.append({"record_type":"dictionary_entry","dictionary":"monomers","id":i,"name":m})
    (root / "schema.jsonl").write_text("".join(json.dumps(x, separators=(",", ":")) + "\n" for x in schema))

    n0 = round(NA * 1e-20)
    times = np.array([900., 1800., 2700.])
    ia = np.array([n0 - 1000, n0 - 2000, n0 - 3000], dtype='<u8')
    bma = np.array([n0 - 2000, n0 - 3500, n0 - 5000], dtype='<u8')
    write_table(root, "snapshots", {
        "snapshot_id": np.arange(3, dtype='<u8'),
        "time": times,
        "kmc_event": np.array([3000, 5500, 8000], dtype='<u8'),
        "snapshot_reason_id": np.ones(3, dtype='<u4'),
        "is_final": np.array([0,0,1], dtype=np.bool_),
        "has_chains": np.array([0,0,1], dtype=np.bool_),
        "has_sequences": np.zeros(3, dtype=np.bool_),
        "kinetic_parameter_set_id": np.zeros(3, dtype='<u8'),
    })
    counts = np.column_stack([ia,bma]).reshape(-1)
    write_table(root, "state", {
        "snapshot_id": np.repeat(np.arange(3, dtype='<u8'), 2),
        "entity_id": np.tile(np.arange(2, dtype='<u4'), 3),
        "count": counts,
        "moles": counts.astype(float) / NA,
        "concentration": counts.astype(float) / NA / 1e-20,
    })
    write_table(root, "chains", {
        "chain_record_id": np.array([0], dtype='<u8'),
        "snapshot_id": np.array([2], dtype='<u8'),
        "count": np.array([1], dtype='<u8'),
    })
    write_table(root, "chain_composition", {
        "chain_record_id": np.array([0,0], dtype='<u8'),
        "monomer_id": np.array([0,1], dtype='<u4'),
        "unit_count": np.array([3000,5000], dtype='<u8'),
    })
    (root / "RESULTS_COMPLETE").write_text("slimmc-storage-v1\n")


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        root = Path(td) / "run"
        root.mkdir()
        build(root)
        run = pyslimmc.open(root)
        assert run.desc == "regression"
        assert np.all(run.t == np.array([900.,1800.,2700.]))
        assert run.conv.total[0] > 0.0
        assert np.allclose(run.F.cum["IA"], [1/3, 2000/5500, 3000/8000])
        assert np.all(np.isfinite(run.F.cum["IA"]))
    print("pyslimmc initial state/F.cum regression: PASS")



def test_script_contract() -> None:
    main()

if __name__ == "__main__":
    main()
