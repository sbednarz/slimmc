from __future__ import annotations
import json, subprocess, tempfile
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT/'homo/tests/regression/frp_all_channels_seeded/FRP_ALLCHAN01.model'
BIN = ROOT/'homo/slimmc-stage-h'
with tempfile.TemporaryDirectory() as raw:
    td=Path(raw); model=td/'run_000001.model'
    model.write_text(MODEL.read_text(encoding='utf-8'), encoding='utf-8')
    subprocess.run([str(BIN), str(model)], check=True, capture_output=True, text=True)
    run=td/'results/run_000001'
    snap={p.stem:np.load(p,allow_pickle=False) for p in (run/'snapshots').glob('*.npy')}
    ch={p.stem:np.load(p,allow_pickle=False) for p in (run/'chains').glob('*.npy')}
    seq=np.load(run/'sequences/symbols.npy',allow_pickle=False)
    assert len({len(a) for a in ch.values()})==1
    n=len(ch['chain_record_id']); assert n>0
    assert np.array_equal(ch['chain_record_id'],np.arange(n,dtype=np.uint64))
    assert ch['chain_record_id'].dtype==np.dtype('<u8')
    assert ch['population_id'].dtype==np.dtype('<u4')
    assert ch['molar_mass'].dtype==np.dtype('<f8')
    assert ch['count'].dtype==np.dtype('<u8')
    assert np.all(ch['dp']>=1) and np.all(ch['count']>=1)
    assert np.all(np.isfinite(ch['molar_mass'])) and np.all(ch['molar_mass']>0)
    md=json.loads((run/'run_metadata.json').read_text())
    np.testing.assert_allclose(ch['moles'],ch['count']/md['avogadro_constant_mol_inv'],rtol=1e-12,atol=0)
    np.testing.assert_allclose(ch['concentration'],ch['moles']/md['kmc_volume_L'],rtol=1e-12,atol=0)
    assert np.all(ch['sequence_offset']==0) and np.all(ch['sequence_length']==0)
    assert seq.dtype==np.dtype('<u4') and len(seq)==0
    chain_sids=np.unique(ch['snapshot_id'])
    assert np.all(snap['has_chains'][chain_sids])
    assert set(chain_sids.tolist())==set(np.flatnonzero(snap['has_chains']).tolist())
    # contiguous blocks and deterministic order
    assert np.all(np.diff(ch['snapshot_id'])>=0)
    for sid in chain_sids:
        rows=np.flatnonzero(ch['snapshot_id']==sid)
        assert np.array_equal(rows,np.arange(rows[0],rows[-1]+1))
        keys=list(zip(ch['population_id'][rows],ch['pool_id'][rows],ch['dp'][rows],ch['left_end_id'][rows],ch['right_end_id'][rows],ch['origin_id'][rows]))
        assert keys==sorted(keys)
        assert len(keys)==len(set(keys))
    # dictionary bounds from schema
    recs=[json.loads(x) for x in (run/'schema.jsonl').read_text().splitlines()]
    def ids(name): return {r['id'] for r in recs if r.get('record_type')=='dictionary_entry' and r.get('dictionary')==name}
    assert set(ch['population_id'].tolist())<=ids('chain_populations')
    assert set(ch['pool_id'].tolist())<=ids('chain_pools')
    assert set(ch['origin_id'].tolist())<=ids('chain_origins')
    assert set(ch['left_end_id'].tolist())<=ids('chain_end_types')
    assert set(ch['right_end_id'].tolist())<=ids('chain_end_types')
print('Homo Slimmc Storage v1 chains contract: OK')
