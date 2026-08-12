from __future__ import annotations
import json
from pathlib import Path
import tempfile
import numpy as np
import pyslimmc
from pyslimmc.core import DataUnavailableError


def write_table(root: Path, name: str, columns: dict[str, np.ndarray]) -> None:
    d=root/name; d.mkdir(parents=True)
    for key,value in columns.items(): np.save(d/f'{key}.npy', value, allow_pickle=False)


def build(root: Path, *, homo: bool=False) -> None:
    monomers=['M'] if homo else ['A','B']
    meta={"run_id":root.name,"storage":"slimmc-storage","storage_format_version":"1.2.0","run_status":"completed","validation_error_count":0,"engine":"slimmc-homo" if homo else "slimmc-copo","kinetic_model":"homo" if homo else "copo"}
    (root/'run_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    schema=[{"record_type":"schema_header","schema_name":"slimmc-storage","schema_version":"1.2.0"}]
    for t in ('snapshots','state','chains','chain_composition'): schema.append({"record_type":"table","name":t,"required":True})
    schema += [{"record_type":"dictionary_entry","dictionary":"state_entities","id":i,"name":m,"kind":"monomer"} for i,m in enumerate(monomers)]
    schema += [{"record_type":"dictionary_entry","dictionary":"monomers","id":i,"name":m} for i,m in enumerate(monomers)]
    (root/'schema.jsonl').write_text(''.join(json.dumps(x,separators=(',',':'))+'\n' for x in schema))
    n=3
    write_table(root,'snapshots',{
      'snapshot_id':np.arange(n,dtype='<u8'),'time':np.array([0.,1.,2.],dtype='<f8'),'kmc_event':np.array([0,10,20],dtype='<u8'),
      'snapshot_reason_id':np.array([0,1,4],dtype='<u4'),'is_final':np.array([0,0,1],dtype=np.bool_),
      'has_chains':np.array([1,1,1],dtype=np.bool_),'has_sequences':np.array([0,0,0],dtype=np.bool_),'kinetic_parameter_set_id':np.zeros(n,dtype='<u8')})
    if homo:
      counts=np.array([100,70,40],dtype='<u8'); moles=counts.astype(float)
    else:
      counts=np.array([80,20,50,10,20,5],dtype='<u8'); moles=counts.astype(float)
    ns=len(monomers)
    write_table(root,'state',{'snapshot_id':np.repeat(np.arange(n,dtype='<u8'),ns),'entity_id':np.tile(np.arange(ns,dtype='<u4'),n),'count':counts,'moles':moles,'concentration':moles})
    # one compressed polymer record at each snapshot; cumulative polymerized units
    write_table(root,'chains',{'chain_record_id':np.arange(n,dtype='<u8'),'snapshot_id':np.arange(n,dtype='<u8'),'count':np.ones(n,dtype='<u8')})
    if homo:
      unit=np.array([0,30,60],dtype='<u8'); mids=np.zeros(n,dtype='<u4'); rids=np.arange(n,dtype='<u8')
    else:
      # cumulative: [0,0], [30,10], [60,15]
      rids=np.repeat(np.arange(n,dtype='<u8'),2); mids=np.tile(np.array([0,1],dtype='<u4'),n); unit=np.array([0,0,30,10,60,15],dtype='<u8')
    write_table(root,'chain_composition',{'chain_record_id':rids,'monomer_id':mids,'unit_count':unit})
    (root/'RESULTS_COMPLETE').write_text('slimmc-storage-v1\n')


def main():
  with tempfile.TemporaryDirectory() as td:
    root=Path(td)/'copo'; root.mkdir(); build(root)
    run=pyslimmc.open(root)
    assert np.allclose(run.conv['A'],[0,.375,.75])
    assert np.allclose(run.conv['B'],[0,.5,.75])
    assert np.allclose(run.conv.total,[0,.4,.75])
    assert np.allclose(run.f['A'],[.8,5/6,.8])
    assert np.allclose(run.f0['A'],[.8,.8,.8])
    assert np.allclose(run.F.cum['A'][1:],[.75,.8])
    assert np.isfinite(run.F.cum['A'][0])
    assert np.allclose(run.F.int['A'][1:],[.75,30/35])
    try: run.F.ins
    except DataUnavailableError: pass
    else: raise AssertionError('binary F.ins accepted without kinetic adapter')
  with tempfile.TemporaryDirectory() as td:
    root=Path(td)/'homo'; root.mkdir(); build(root,homo=True)
    run=pyslimmc.open(root)
    assert np.allclose(run.conv['M'],[0,.3,.6])
    assert np.allclose(run.conv.total,[0,.3,.6])
    assert np.allclose(run.f['M'],[1,1,1])
    assert np.allclose(run.F.ins['M'],[1,1,1])
    assert np.allclose(run.F.int['M'][1:],[1,1])
    assert np.allclose(run.F.cum['M'][1:],[1,1])
  print('pyslimmc L2.3 composition: PASS')

def test_script_contract():
    main()

if __name__=='__main__': main()
