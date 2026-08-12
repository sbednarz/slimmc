from __future__ import annotations
import json
from pathlib import Path
import tempfile
import numpy as np
import pyslimmc
from pyslimmc.core import DataUnavailableError


def write_table(root: Path, name: str, columns: dict[str, np.ndarray]) -> None:
    d = root / name; d.mkdir(parents=True)
    for key, value in columns.items(): np.save(d / f"{key}.npy", value, allow_pickle=False)


def build(root: Path, *, full: bool) -> None:
    meta = {"run_id":root.name,"storage":"slimmc-storage","storage_format_version":"1.2.0","run_status":"completed","validation_error_count":0,"engine":"slimmc-copo","kinetic_model":"copo","sequence_mode":"full" if full else "composition"}
    (root/'run_metadata.json').write_text(json.dumps(meta, indent=2)+'\n')
    schema=[{"record_type":"schema_header","schema_name":"slimmc-storage","schema_version":"1.2.0"}]
    for t in ('snapshots','state','chains','chain_composition','sequences'): schema.append({"record_type":"table","name":t,"required":True})
    for i,n in enumerate(('A','B')):
        schema.append({"record_type":"dictionary_entry","dictionary":"monomers","id":i,"name":n})
        schema.append({"record_type":"dictionary_entry","dictionary":"sequence_symbols","id":i,"name":n})
        schema.append({"record_type":"dictionary_entry","dictionary":"state_entities","id":i,"name":n,"kind":"monomer"})
    for d, vals in {
      'chain_populations':('live','dead'), 'chain_pools':('not_applicable','terminal_A','terminal_B'),
      'chain_origins':('unknown','init','term_c'), 'chain_end_types':('not_applicable','unknown','init_end','radical_end','dead_end')}.items():
        for i,n in enumerate(vals): schema.append({"record_type":"dictionary_entry","dictionary":d,"id":i,"name":n})
    (root/'schema.jsonl').write_text(''.join(json.dumps(x,separators=(',',':'))+'\n' for x in schema))
    write_table(root,'snapshots',{
      'snapshot_id':np.array([0,1],dtype='<u8'),'time':np.array([0.,1.],dtype='<f8'),'kmc_event':np.array([0,10],dtype='<u8'),
      'snapshot_reason_id':np.array([0,4],dtype='<u4'),'is_final':np.array([0,1],dtype=np.bool_),
      'has_chains':np.array([0,1],dtype=np.bool_),'has_sequences':np.array([0,full],dtype=np.bool_),'kinetic_parameter_set_id':np.zeros(2,dtype='<u8')})
    write_table(root,'state',{'snapshot_id':np.repeat(np.arange(2,dtype='<u8'),2),'entity_id':np.tile(np.arange(2,dtype='<u4'),2),'count':np.array([10,10,4,6],dtype='<u8'),'moles':np.array([10.,10.,4.,6.]),'concentration':np.array([10.,10.,4.,6.])})
    # three compressed records in final snapshot
    seqs=[[0,0,1],[0,1,0],[1,1]]
    offsets=[]; lengths=[]; symbols=[]
    for seq in seqs:
      offsets.append(len(symbols)); lengths.append(len(seq)); symbols.extend(seq if full else [])
    write_table(root,'chains',{
      'chain_record_id':np.arange(3,dtype='<u8'),'snapshot_id':np.ones(3,dtype='<u8'),
      'population_id':np.array([0,0,1],dtype='<u4'),'pool_id':np.array([1,2,0],dtype='<u4'),'origin_id':np.array([1,1,2],dtype='<u4'),
      'dp':np.array([3,3,2],dtype='<u8'),'molar_mass':np.array([300.,310.,220.]),'count':np.array([2,1,4],dtype='<u8'),'moles':np.array([2.,1.,4.]),'concentration':np.array([2.,1.,4.]),
      'left_end_id':np.array([2,2,2],dtype='<u4'),'right_end_id':np.array([3,3,4],dtype='<u4'),
      'has_first_monomer':np.ones(3,dtype=np.bool_),'first_monomer_id':np.array([0,0,1],dtype='<u4'),
      'has_penultimate_monomer':np.ones(3,dtype=np.bool_),'penultimate_monomer_id':np.array([0,1,1],dtype='<u4'),
      'has_last_monomer':np.ones(3,dtype=np.bool_),'last_monomer_id':np.array([1,0,1],dtype='<u4'),
      'has_sequence':np.full(3,full,dtype=np.bool_),'sequence_offset':np.array(offsets,dtype='<u8') if full else np.zeros(3,dtype='<u8'),'sequence_length':np.array(lengths,dtype='<u8') if full else np.zeros(3,dtype='<u8')})
    write_table(root,'chain_composition',{
      'chain_record_id':np.repeat(np.arange(3,dtype='<u8'),2),'monomer_id':np.tile(np.arange(2,dtype='<u4'),3),
      'unit_count':np.array([2,1,2,1,0,2],dtype='<u8')})
    write_table(root,'sequences',{'symbols':np.array(symbols,dtype='<u4')})
    (root/'RESULTS_COMPLETE').write_text('slimmc-storage-v1\n')


def main():
  for full in (False, True):
    with tempfile.TemporaryDirectory() as td:
      root=Path(td)/'run'; root.mkdir(); build(root,full=full)
      run=pyslimmc.open(root)
      c=run.last.chains
      assert len(c)==3 and len(c.live)==2 and len(c.dead)==1
      assert np.array_equal(c.live.pool('terminal_A').dp,[3])
      assert np.array_equal(c.origin('term_c').count,[4])
      assert len(c.where(dp_min=3,dp_max=3))==2
      assert np.array_equal(c.composition.counts['A'],[2,2,0])
      assert np.array_equal(c.composition.counts.total,c.dp)
      assert np.allclose(c.composition.fractions['A'],[2/3,2/3,0])
      assert c.first_monomer.tolist()==['A','A','B']
      assert c.penultimate_monomer.tolist()==['A','B','B']
      assert c.last_monomer.tolist()==['B','A','B']
      assert c.population_activity_names.tolist()==['live','live','dead']
      assert c.pool_names.tolist()==['terminal_A','terminal_B','not_applicable']
      assert run.chains.last.record(2).composition.counts['B']==2
      assert run.chains.record(1).count==1
      if full:
        assert c.has_sequences
        assert c[0].sequence==('A','A','B')
        assert c.sequences[1]==('A','B','A')
      else:
        assert not c.has_sequences
        try: _=c.sequences
        except DataUnavailableError: pass
        else: raise AssertionError('composition mode exposed full sequences')
  print('pyslimmc L2.4 chains: PASS')

def test_script_contract():
    main()

if __name__=='__main__': main()
