from __future__ import annotations
import json
from pathlib import Path
import numpy as np


def check_run(root: Path, n_monomers: int) -> None:
    md = json.loads((root / 'run_metadata.json').read_text())
    assert md['volume_mode'] == 'variable'
    assert md['initial_volume_mL'] == 100.0
    assert md['current_volume_mL'] == 103.0
    np.testing.assert_allclose(np.load(root/'snapshots/volume_mL.npy'), [100,101,102,103,103])
    np.testing.assert_allclose(np.load(root/'snapshots/kmc_volume_L.npy'), [1e-18,1.01e-18,1.02e-18,1.03e-18,1.03e-18])
    np.testing.assert_allclose(np.load(root/'feed_events/dose_mL.npy'), [1,1,1])
    np.testing.assert_allclose(np.load(root/'feed_events/volume_before_mL.npy'), [100,101,102])
    np.testing.assert_allclose(np.load(root/'feed_events/volume_after_mL.npy'), [101,102,103])
    assert len(np.load(root/'monomer_balance/snapshot_id.npy')) == 5*n_monomers
    assert not [json.loads(x) for x in (root/'diagnostics/validation.jsonl').read_text().splitlines() if json.loads(x)['status']=='fail']

if __name__ == '__main__':
    base=Path('/mnt/data/stage23_smoke')
    check_run(base/'homo/results/feed',1)
    check_run(base/'copo/results/feed',2)
    print('semibatch storage checks: OK')
