from __future__ import annotations
import json, tempfile
from pathlib import Path
import numpy as np
import pyslimmc


def wt(root: Path, name: str, cols):
    d=root/name; d.mkdir(parents=True)
    for k,v in cols.items(): np.save(d/f'{k}.npy',v,allow_pickle=False)


def build(root: Path):
    meta={'run_id':'run','storage':'slimmc-storage','storage_format_version':'1.2.0','run_status':'completed','validation_error_count':0,'engine':'slimmc-copo','kinetic_model':'copo','sequence_mode':'composition'}
    (root/'run_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    schema=[{'record_type':'schema_header','schema_name':'slimmc-storage','schema_version':'1.2.0'}]
    for t in ('snapshots','state','chains','chain_composition','sequences','moments'): schema.append({'record_type':'table','name':t,'required':True})
    for i,n in enumerate(('A','B')):
        schema += [
            {'record_type':'dictionary_entry','dictionary':'monomers','id':i,'name':n,'molar_mass_increment':100.0+10*i},
            {'record_type':'dictionary_entry','dictionary':'state_entities','id':i,'name':n,'kind':'monomer'}]
    for d,vals in {'chain_populations':('live','dead'),'chain_pools':('not_applicable','terminal_A'),'chain_origins':('unknown','init'),'chain_end_types':('not_applicable','unknown'),'population_scope':('all','live','dead'),'mass_bases':('repeat_units','with_end_groups')}.items():
        for i,n in enumerate(vals): schema.append({'record_type':'dictionary_entry','dictionary':d,'id':i,'name':n})
    (root/'schema.jsonl').write_text(''.join(json.dumps(x,separators=(',',':'))+'\n' for x in schema))
    wt(root,'snapshots',{'snapshot_id':np.array([0,1],'<u8'),'time':np.array([0.,1.]),'kmc_event':np.array([0,10],'<u8'),'snapshot_reason_id':np.array([0,4],'<u4'),'is_final':np.array([0,1],bool),'has_chains':np.array([0,1],bool),'has_sequences':np.array([0,0],bool),'kinetic_parameter_set_id':np.zeros(2,'<u8')})
    wt(root,'state',{'snapshot_id':np.repeat(np.arange(2,dtype='<u8'),2),'entity_id':np.tile(np.arange(2,dtype='<u4'),2),'count':np.array([10,10,4,6],'<u8'),'moles':np.array([10.,10.,4.,6.]),'concentration':np.array([10.,10.,4.,6.])})
    wt(root,'chains',{'chain_record_id':np.arange(2,dtype='<u8'),'snapshot_id':np.ones(2,'<u8'),'population_id':np.array([0,1],'<u4'),'pool_id':np.array([1,0],'<u4'),'origin_id':np.array([1,1],'<u4'),'dp':np.array([2,4],'<u8'),'molar_mass':np.array([210.,430.]),'count':np.array([3,2],'<u8'),'moles':np.array([3.,2.]),'concentration':np.array([3.,2.]),'left_end_id':np.ones(2,'<u4'),'right_end_id':np.ones(2,'<u4'),'has_first_monomer':np.ones(2,bool),'first_monomer_id':np.array([0,0],'<u4'),'has_penultimate_monomer':np.ones(2,bool),'penultimate_monomer_id':np.array([1,1],'<u4'),'has_last_monomer':np.ones(2,bool),'last_monomer_id':np.array([1,1],'<u4'),'has_sequence':np.zeros(2,bool),'sequence_offset':np.zeros(2,'<u8'),'sequence_length':np.zeros(2,'<u8')})
    wt(root,'chain_composition',{'chain_record_id':np.repeat(np.arange(2,dtype='<u8'),2),'monomer_id':np.tile(np.arange(2,dtype='<u4'),2),'unit_count':np.array([1,1,2,2],'<u8')})
    wt(root,'sequences',{'symbols':np.array([],dtype='<u4')})
    # 6 rows: all/live/dead x RU/with ends. Values chosen for easy API checks.
    sid=np.ones(6,'<u8'); pop=np.repeat(np.arange(3,dtype='<u4'),2); basis=np.tile(np.arange(2,dtype='<u4'),3)
    wt(root,'moments',{'snapshot_id':sid,'population_scope_id':pop,'mass_basis_id':basis,'chain_count':np.array([5,5,3,3,2,2],'<u8'),'sum_dp':np.zeros(6),'sum_dp2':np.zeros(6),'dp_n':np.array([2.8,2.8,2.,2.,4.,4.]),'dp_w':np.array([3.142857,3.142857,2.,2.,4.,4.]),'sum_molar_mass':np.zeros(6),'sum_molar_mass2':np.zeros(6),'sum_molar_mass3':np.zeros(6),'mn':np.array([280.,298.,200.,210.,400.,430.]),'mw':np.array([314.2857,329.8658,200.,210.,400.,430.]),'mz':np.array([350.,365.,200.,210.,400.,430.]),'dispersity':np.array([1.12245,1.10693,1.,1.,1.,1.])})
    (root/'RESULTS_COMPLETE').write_text('slimmc-storage-v1\n')


def main():
    with tempfile.TemporaryDirectory() as td:
        root=Path(td)/'run'; root.mkdir(); build(root)
        run=pyslimmc.open(root); snap=run.last
        assert np.isnan(np.asarray(run.mn)[0]) and np.isclose(np.asarray(run.mn)[1],298.)
        assert np.isclose(snap.mn,298.) and np.isclose(snap.moments(population='live', mass_model='repeat_units').dpn,2.)
        assert np.isclose(run.moments(snapshot='final', population='dead', mass_model='with_end_groups').mw,430.)

        dp_counts=snap.dp_counts()
        assert np.array_equal(dp_counts.dp,np.array([2,4])) and dp_counts.total_chains == 5
        mass_counts=snap.mass_counts()
        assert np.array_equal(mass_counts.mass,np.array([210.,430.])) and mass_counts.total_chains == 5

        cld=snap.cld(form='number')
        assert np.array_equal(cld.x,np.array([2.,4.])) and np.isclose(cld.y.sum(),1.)
        mwd=snap.mwd(form='number')
        assert np.array_equal(mwd.x,np.array([210.,430.])) and np.isclose(mwd.y.sum(),1.)
        log_mwd=snap.mwd(form='log')
        assert np.allclose(log_mwd.x,np.log10(mwd.x)) and np.allclose(log_mwd.y,snap.mwd(form='mass').y)
        assert np.array_equal(run.mwd(form='number').x,mwd.x)
        assert not hasattr(run,'chain_mass_spectrum')
    print('pyslimmc L2.5 moments/distributions: PASS')

def test_script_contract():
    main()

if __name__=='__main__': main()
