import hashlib, json
from pathlib import Path
import pyslimmc as sl


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def make_run(root: Path, value=b"abc"):
    root.mkdir()
    (root/"input.model").write_bytes(b"model")
    (root/"state.npy").write_bytes(value)
    line=f"{sha(root/'state.npy')}  state.npy\n"
    (root/"checksums.sha256").write_text(line)
    storage_hash=hashlib.sha256(b"slimmc-storage-hash-v1\n"+line.encode()).hexdigest()
    meta={"run_id":root.name,"storage":"slimmc-storage","storage_format_version":"1.2.0","run_status":"completed","input":{"file":"input.model","hash":sha(root/'input.model'),"hash_algorithm":"sha256"},
          "model":{"hash":"model123"},"execution":{"binary_hash":"deadbeef","binary_hash_algorithm":"sha256"},
          "storage_info":{"hash":storage_hash,"hash_schema":"slimmc-storage-hash-v1","manifest_file":"checksums.sha256"}}
    (root/"run_metadata.json").write_text(json.dumps(meta))
    (root/"schema.jsonl").write_text("{}\n")
    (root/"RESULTS_COMPLETE").write_text("")
    return sl.open(root)

def test_verify_and_compare(tmp_path):
    a=make_run(tmp_path/'a')
    b=make_run(tmp_path/'b')
    r=a.reproducibility.verify()
    assert r.ok and r.binary_status == "NOT CHECKED"
    assert a.reproducibility.compare(b)["storage"] == "IDENTICAL"
    (b.path/"state.npy").write_bytes(b"changed")
    assert b.reproducibility.verify().overall == "FAILED"
