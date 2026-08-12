from __future__ import annotations
import json
import subprocess
import tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / 'homo/tests/regression/frp_all_channels_seeded/FRP_ALLCHAN01.model'
BIN = ROOT / 'homo/slimmc-stage-h'

with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    model = td / 'run_000001.model'
    model.write_text(MODEL.read_text(encoding='utf-8'), encoding='utf-8')
    subprocess.run([str(BIN), str(model)], check=True, capture_output=True, text=True)
    run = td / 'results/run_000001'
    assert (run / 'input.model').read_bytes() == model.read_bytes()
    assert (run / 'RESULTS_COMPLETE').is_file()
    assert not (run / '.work').exists()
    md = json.loads((run / 'run_metadata.json').read_text(encoding='utf-8'))
    assert md['run_status'] == 'completed'
    assert md['run_id'] == 'run_000001'
    assert md['results_format_version'] == '1.0.0'
    assert '\n  "run_id"' in (run / 'run_metadata.json').read_text(encoding='utf-8')

    snap = {p.stem: np.load(p, allow_pickle=False) for p in (run/'snapshots').glob('*.npy')}
    state = {p.stem: np.load(p, allow_pickle=False) for p in (run/'state').glob('*.npy')}
    channels = {p.stem: np.load(p, allow_pickle=False) for p in (run/'channel_events').glob('*.npy')}
    ksets = {p.stem: np.load(p, allow_pickle=False) for p in (run/'kinetic_parameters/sets').glob('*.npy')}
    kvals = {p.stem: np.load(p, allow_pickle=False) for p in (run/'kinetic_parameters/values').glob('*.npy')}
    assert len({len(a) for a in snap.values()}) == 1
    assert len({len(a) for a in state.values()}) == 1
    assert len({len(a) for a in channels.values()}) == 1
    assert len({len(a) for a in ksets.values()}) == 1
    assert len({len(a) for a in kvals.values()}) == 1
    n_snap = len(snap['snapshot_id'])
    assert np.array_equal(snap['snapshot_id'], np.arange(n_snap, dtype=np.uint64))
    assert snap['snapshot_id'].dtype == np.dtype('<u8')
    assert snap['snapshot_reason_id'].dtype == np.dtype('<u4')
    assert snap['time'].dtype == np.dtype('<f8')
    assert snap['is_final'].dtype == np.dtype('bool')
    assert np.count_nonzero(snap['is_final']) == 1
    assert snap['is_final'][-1]
    assert snap['snapshot_id'][0] == 0 and snap['time'][0] == 0 and snap['kmc_event'][0] == 0
    assert np.all(np.diff(snap['time']) >= 0)
    assert np.all(np.diff(snap['kmc_event']) >= 0)

    assert len(state['snapshot_id']) % n_snap == 0
    n_entities = len(state['snapshot_id']) // n_snap
    for i in range(n_snap):
        block = slice(i*n_entities, (i+1)*n_entities)
        assert np.all(state['snapshot_id'][block] == i)
        assert np.array_equal(state['entity_id'][block], np.arange(n_entities, dtype=np.uint32))
    na = md['avogadro_constant_mol_inv']
    vol = md['kmc_volume_L']
    np.testing.assert_allclose(state['moles'], state['count']/na, rtol=1e-12, atol=0)
    np.testing.assert_allclose(state['concentration'], state['moles']/vol, rtol=1e-12, atol=0)
    assert len(channels['snapshot_id']) % n_snap == 0
    n_channels = len(channels['snapshot_id']) // n_snap
    assert n_channels > 0
    assert channels['snapshot_id'].dtype == np.dtype('<u8')
    assert channels['channel_id'].dtype == np.dtype('<u4')
    assert channels['event_count'].dtype == np.dtype('<u8')
    for i in range(n_snap):
        block = slice(i*n_channels, (i+1)*n_channels)
        assert np.all(channels['snapshot_id'][block] == i)
        assert np.array_equal(channels['channel_id'][block], np.arange(n_channels, dtype=np.uint32))
        assert np.all(channels['event_count'][block] ==
                      channels['productive_event_count'][block] +
                      channels['nonproductive_event_count'][block])
        assert int(channels['event_count'][block].sum(dtype=np.uint64)) == int(snap['kmc_event'][i])
    for cid in range(n_channels):
        rows = channels['channel_id'] == cid
        assert np.all(np.diff(channels['event_count'][rows]) >= 0)
        assert np.all(np.diff(channels['productive_event_count'][rows]) >= 0)
        assert np.all(np.diff(channels['nonproductive_event_count'][rows]) >= 0)

    # Complete kinetic-parameter sets, dense by set × parameter id.
    n_sets = len(ksets['kinetic_parameter_set_id'])
    assert n_sets >= 1
    assert np.array_equal(ksets['kinetic_parameter_set_id'], np.arange(n_sets, dtype=np.uint64))
    assert ksets['kinetic_parameter_set_id'].dtype == np.dtype('<u8')
    assert ksets['start_time'].dtype == np.dtype('<f8')
    assert ksets['has_source_action'].dtype == np.dtype('bool')
    assert not ksets['has_source_action'][0]
    assert ksets['source_action_id'][0] == 0
    assert np.all(np.diff(ksets['start_kmc_event']) >= 0)
    assert np.all(np.diff(ksets['start_time']) >= 0)
    assert len(kvals['kinetic_parameter_set_id']) % n_sets == 0
    n_params = len(kvals['kinetic_parameter_set_id']) // n_sets
    assert n_params >= 1
    assert np.all(np.isfinite(kvals['value']))
    for set_id in range(n_sets):
        block = slice(set_id*n_params, (set_id+1)*n_params)
        assert np.all(kvals['kinetic_parameter_set_id'][block] == set_id)
        assert np.array_equal(kvals['kinetic_parameter_id'][block], np.arange(n_params, dtype=np.uint32))
    assert np.all(snap['kinetic_parameter_set_id'] < n_sets)

    actions = {p.stem: np.load(p, allow_pickle=False) for p in (run/'actions').glob('*.npy')}
    assert len({len(a) for a in actions.values()}) == 1
    n_actions = len(actions['action_id'])
    assert n_actions > 0
    assert np.array_equal(actions['action_id'], np.arange(n_actions, dtype=np.uint64))
    assert actions['source_line'].dtype == np.dtype('<u4')
    assert actions['state_changed'].dtype == np.dtype('bool')
    assert np.all(actions['trigger_type_id'] > 0)
    assert np.all(np.isfinite(actions['scheduled_time']))
    assert np.all(actions['has_snapshot'])
    assert np.all(actions['output_written'])
    assert np.all(actions['snapshot_id'] < n_snap)
    assert not np.any(actions['has_kinetic_parameter_set'])
    assert np.all(actions['kinetic_parameter_set_id'] == 0)
    conditions = {p.stem: np.load(p, allow_pickle=False) for p in (run/'action_conditions').glob('*.npy')}
    assert len({len(a) for a in conditions.values()}) == 1
    assert len(conditions['condition_record_id']) == 0
    assert (run/'actions/messages.jsonl').read_text(encoding='utf-8') == ''

    zero = state['count'] == 0
    assert np.all(state['moles'][zero] == 0.0)
    assert np.all(state['concentration'][zero] == 0.0)
print('Homo Slimmc Storage v1 snapshots/state/channel_events/kinetic_parameters/actions contract: OK')
