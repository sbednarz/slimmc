from pathlib import Path

from pyslimmc.run import Run


def test_storage_info_is_exposed_readonly(tmp_path: Path):
    run = Run(
        tmp_path,
        _metadata={
            "run_id": "case_01",
            "storage_info": {
                "name": "slimmc-storage",
                "format_version": "1",
                "complete": True,
                "hash": "abc123",
                "hash_algorithm": "sha256",
                "hash_schema": "slimmc-storage-hash-v1",
                "manifest_file": "checksums.sha256",
                "file_count": 12,
            },
        },
    )
    assert run.storage.hash == "abc123"
    assert run.storage.hash_algorithm == "sha256"
    assert run.storage.hash_schema == "slimmc-storage-hash-v1"
    assert run.storage.manifest_file == "checksums.sha256"


def test_legacy_storage_metadata_fallback(tmp_path: Path):
    run = Run(
        tmp_path,
        _metadata={
            "run_id": "legacy",
            "storage": "slimmc-storage",
            "storage_format_version": "1",
        },
    )
    assert run.storage.name == "slimmc-storage"
    assert run.storage.format_version == "1"
    assert run.storage.hash is None
