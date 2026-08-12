from __future__ import annotations
import hashlib, json, os, signal, subprocess, tempfile, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'homo/tests/regression/frp_all_channels_seeded/FRP_ALLCHAN01.model'
BIN = ROOT/'homo/slimmc-stage-h'

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_md(run: Path) -> dict:
    return json.loads((run/'run_metadata.json').read_text(encoding='utf-8'))

with tempfile.TemporaryDirectory() as raw:
    td = Path(raw)
    # completed
    model = td/'completed.model'
    text = BASE.read_text().replace('param t_end 4.0','param t_end 0.02').replace('param max_steps 300000','param max_steps 20000')
    model.write_text(text)
    subprocess.run([str(BIN), str(model)], check=True, capture_output=True, text=True)
    run = td/'results/completed'
    md = read_md(run)
    assert md['run_status'] == 'completed' and md['exit_code'] == 0
    assert md['input_model_sha256'] == digest(run/'input.model')
    assert md['schema_sha256'] == digest(run/'schema.jsonl')
    assert md['platform'] and md['compiler'].startswith('Nim ')
    assert md['build_mode'] in ('debug','release')
    assert (run/'RESULTS_COMPLETE').is_file()
    assert not (run/'.work').exists()
    assert np.count_nonzero(np.load(run/'snapshots/is_final.npy', allow_pickle=False)) == 1

    # failed after initialization: finite input values produce an infinite propensity.
    model = td/'failed.model'
    text = BASE.read_text().replace('rate kd const 1.0e-3','rate kd const 1.0e308')
    model.write_text(text)
    cp = subprocess.run([str(BIN), str(model)], capture_output=True, text=True)
    assert cp.returncode != 0
    run = td/'results/failed'
    md = read_md(run)
    assert md['run_status'] == 'failed' and md['exit_code'] == 1
    assert not (run/'RESULTS_COMPLETE').exists()
    assert (run/'.work').is_dir()
    assert np.count_nonzero(np.load(run/'snapshots/is_final.npy', allow_pickle=False)) == 0
    assert md['input_model_sha256'] == digest(run/'input.model')
    assert md['schema_sha256'] == digest(run/'schema.jsonl')

    # interrupted by Ctrl+C; use a long event-limited run.
    model = td/'interrupted.model'
    text = BASE.read_text().replace('param t_end 4.0','param t_end 1.0e12').replace('param max_steps 300000','param max_steps 2000000000')
    model.write_text(text)
    proc = subprocess.Popen([str(BIN), str(model)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    run = td/'results/interrupted'
    deadline = time.time() + 5
    while time.time() < deadline and not (run/'run_metadata.json').exists() and proc.poll() is None:
        time.sleep(0.01)
    assert proc.poll() is None, 'run ended before interrupt test'
    proc.send_signal(signal.SIGINT)
    out, err = proc.communicate(timeout=15)
    assert proc.returncode == 0, (proc.returncode, out, err)
    md = read_md(run)
    assert md['run_status'] == 'interrupted' and md['exit_code'] == 130
    assert not (run/'RESULTS_COMPLETE').exists()
    assert (run/'.work').is_dir()
    assert np.count_nonzero(np.load(run/'snapshots/is_final.npy', allow_pickle=False)) == 0

print('Homo Slimmc Storage v1 finalization/hash contract: OK')
