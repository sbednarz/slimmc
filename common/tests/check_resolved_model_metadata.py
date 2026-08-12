#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

REQUIRED_MODEL_KEYS = {
    "schema", "kinetic_model", "desc", "parameters", "memory_policy",
    "monomers", "species", "polymers", "endgroups", "rates",
    "reactions", "actions", "variables", "engine_specific",
}
REQUIRED_PARAMETER_KEYS = {
    "kmc_volume_L", "initial_temperature_K", "t_end_s", "max_events",
    "when_check_events", "seed", "dp_max", "sequence_mode", "mass_model",
}


def run_model(engine: Path, source: Path) -> dict:
    with tempfile.TemporaryDirectory(prefix="slimmc-resolved-") as td:
        root = Path(td)
        model = root / source.name
        model.write_bytes(source.read_bytes())
        completed = subprocess.run(
            [str(engine.resolve()), model.name],
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"engine failed for {source}:\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        run_dir = root / "results" / model.stem
        metadata_path = run_dir / "run_metadata.json"
        if not metadata_path.exists():
            raise AssertionError(f"missing {metadata_path}")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        manifest_path = run_dir / "checksums.sha256"
        assert manifest_path.is_file(), f"missing {manifest_path}"
        manifest = manifest_path.read_text(encoding="utf-8")
        for line in manifest.splitlines():
            digest, relative = line.split("  ", 1)
            payload = (run_dir / relative).read_bytes()
            assert hashlib.sha256(payload).hexdigest() == digest, relative
        expected_storage_hash = hashlib.sha256(
            b"slimmc-storage-hash-v1\n" + manifest.encode("utf-8")
        ).hexdigest()
        storage = metadata["storage_info"]
        assert storage["hash"] == expected_storage_hash
        assert storage["hash_algorithm"] == "sha256"
        assert storage["hash_schema"] == "slimmc-storage-hash-v1"
        assert storage["manifest_file"] == "checksums.sha256"
        assert storage["file_count"] == len(manifest.splitlines())
        return metadata


def check_common(metadata: dict, kinetic_model: str) -> dict:
    model = metadata.get("model")
    assert isinstance(model, dict), "run_metadata.json must contain object field 'model'"
    assert REQUIRED_MODEL_KEYS <= set(model), sorted(REQUIRED_MODEL_KEYS - set(model))
    assert model["schema"] == "slimmc-resolved-model-v1"
    assert model["kinetic_model"] == kinetic_model
    assert REQUIRED_PARAMETER_KEYS <= set(model["parameters"])
    assert isinstance(model["parameters"]["seed"], int)
    assert isinstance(model["parameters"]["dp_max"], int)
    assert isinstance(model["actions"], list)
    assert isinstance(model["variables"], list)
    # Stage 13 defines the resolved structure only. Hash fields are added later.
    assert "hash" not in model
    return model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--homo-engine", required=True, type=Path)
    parser.add_argument("--copo-engine", required=True, type=Path)
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2]
    homo_md = run_model(
        args.homo_engine,
        repo / "homo/tests/validation/phase_c/models/H11_no_t0_save.model",
    )
    homo = check_common(homo_md, "homo")
    periodic = homo["actions"][0]
    assert periodic["trigger"] == {"kind": "periodic", "start_s": 0.0, "step_s": 0.1}
    assert periodic["action"]["kind"] == "save"
    assert any(item["builtin"] for item in homo["endgroups"])

    copo_md = run_model(
        args.copo_engine,
        repo / "copo/tests/validation/phase_b/models/C10_seed_reproducibility.model",
    )
    copo = check_common(copo_md, "copo")
    assert len(copo["monomers"]) == 2
    assert copo["parameters"]["sequence_mode"] == "full"
    assert copo["actions"][0]["trigger"]["start_s"] == 0.0
    assert all("source_line" in reaction for reaction in copo["reactions"])

    print("resolved model metadata: PASS")


if __name__ == "__main__":
    main()
