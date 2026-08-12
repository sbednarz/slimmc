from __future__ import annotations
import json, math, subprocess, tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / 'homo/tests/regression/frp_all_channels_seeded/FRP_ALLCHAN01.model'
BIN = ROOT / 'homo/slimmc-stage-h'

with tempfile.TemporaryDirectory() as raw:
    td = Path(raw)
    model = td / 'run_000001.model'
    model.write_text(MODEL.read_text(encoding='utf-8'), encoding='utf-8')
    subprocess.run([str(BIN), str(model)], check=True, capture_output=True, text=True)
    run = td / 'results/run_000001'
    snap = {p.stem: np.load(p, allow_pickle=False) for p in (run/'snapshots').glob('*.npy')}
    ch = {p.stem: np.load(p, allow_pickle=False) for p in (run/'chains').glob('*.npy')}
    mo = {p.stem: np.load(p, allow_pickle=False) for p in (run/'moments').glob('*.npy')}

    assert len({len(a) for a in mo.values()}) == 1
    n = len(mo['snapshot_id'])
    chain_sids = np.flatnonzero(snap['has_chains']).astype(np.uint64)
    assert n == 6 * len(chain_sids)
    assert mo['snapshot_id'].dtype == np.dtype('<u8')
    assert mo['population_scope_id'].dtype == np.dtype('<u4')
    assert mo['mass_basis_id'].dtype == np.dtype('<u4')
    assert mo['chain_count'].dtype == np.dtype('<u8')
    assert mo['mn'].dtype == np.dtype('<f8')

    expected_keys = [(int(sid), scope, basis)
                     for sid in chain_sids for scope in range(3) for basis in range(2)]
    actual_keys = list(zip(mo['snapshot_id'].tolist(),
                           mo['population_scope_id'].tolist(),
                           mo['mass_basis_id'].tolist()))
    assert actual_keys == expected_keys

    recs = [json.loads(line) for line in (run/'schema.jsonl').read_text(encoding='utf-8').splitlines()]
    end_mass = {r['id']: r.get('molar_mass_contribution')
                for r in recs if r.get('record_type') == 'dictionary_entry'
                and r.get('dictionary') == 'chain_end_types'}
    assert {0,1,2} <= {r['id'] for r in recs if r.get('dictionary') == 'population_scope'}
    assert {0,1} <= {r['id'] for r in recs if r.get('dictionary') == 'mass_bases'}

    for row, (sid, scope, basis) in enumerate(expected_keys):
        mask = ch['snapshot_id'] == sid
        if scope == 1:
            mask &= ch['population_id'] == 0
        elif scope == 2:
            mask &= ch['population_id'] == 1
        counts = ch['count'][mask].astype(np.float64)
        dp = ch['dp'][mask].astype(np.float64)
        masses = ch['molar_mass'][mask].astype(np.float64)
        if basis == 0 and len(masses):
            ends = np.array([
                float(end_mass[int(l)]) + float(end_mass[int(r)])
                for l, r in zip(ch['left_end_id'][mask], ch['right_end_id'][mask])
            ])
            masses = masses - ends

        chain_count = int(np.sum(counts, dtype=np.float64)) if len(counts) else 0
        sum_dp = float(np.sum(counts * dp)) if len(counts) else 0.0
        sum_dp2 = float(np.sum(counts * dp * dp)) if len(counts) else 0.0
        s1 = float(np.sum(counts * masses)) if len(counts) else 0.0
        s2 = float(np.sum(counts * masses * masses)) if len(counts) else 0.0
        s3 = float(np.sum(counts * masses * masses * masses)) if len(counts) else 0.0
        assert int(mo['chain_count'][row]) == chain_count
        np.testing.assert_allclose(mo['sum_dp'][row], sum_dp, rtol=1e-12, atol=0)
        np.testing.assert_allclose(mo['sum_dp2'][row], sum_dp2, rtol=1e-12, atol=0)
        np.testing.assert_allclose(mo['sum_molar_mass'][row], s1, rtol=1e-12, atol=0)
        np.testing.assert_allclose(mo['sum_molar_mass2'][row], s2, rtol=1e-12, atol=0)
        np.testing.assert_allclose(mo['sum_molar_mass3'][row], s3, rtol=1e-12, atol=0)
        if chain_count == 0:
            assert sum_dp == sum_dp2 == s1 == s2 == s3 == 0.0
            for name in ['dp_n','dp_w','mn','mw','mz','dispersity']:
                assert math.isnan(float(mo[name][row]))
        else:
            expected = {
                'dp_n': sum_dp/chain_count,
                'dp_w': sum_dp2/sum_dp,
                'mn': s1/chain_count,
                'mw': s2/s1,
                'mz': s3/s2,
                'dispersity': (s2/s1)/(s1/chain_count),
            }
            for name, value in expected.items():
                np.testing.assert_allclose(mo[name][row], value, rtol=1e-12, atol=0)
                assert np.isfinite(mo[name][row])
        assert not np.any(np.isinf([mo[name][row] for name in ['dp_n','dp_w','mn','mw','mz','dispersity']]))

print('Homo Slimmc Storage v1 moments reconstruction contract: OK')
