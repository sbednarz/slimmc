from __future__ import annotations
import json, subprocess, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / 'homo/tests/regression/frp_all_channels_seeded/FRP_ALLCHAN01.model'
BIN = ROOT / 'homo/slimmc-stage-h'

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    model = td / 'run_000001.model'
    model.write_text(MODEL.read_text(encoding='utf-8'), encoding='utf-8')
    subprocess.run([str(BIN), str(model)], check=True, capture_output=True, text=True)
    run = td / 'results/run_000001'
    records = [json.loads(line) for line in (run/'diagnostics/validation.jsonl').read_text(encoding='utf-8').splitlines()]
    assert records
    assert all(r['status'] == 'pass' for r in records), records
    assert all({'check','status','severity'} <= set(r) for r in records)
    md = json.loads((run/'run_metadata.json').read_text(encoding='utf-8'))
    assert md['validation_status'] == 'passed'
    assert md['validation_error_count'] == 0
    assert md['validation_warning_count'] == 0
    assert md['run_status'] == 'completed'
    assert (run/'RESULTS_COMPLETE').is_file()
print('Homo Slimmc Storage v1 validator contract: OK')
