from __future__ import annotations
import json,tempfile
from pathlib import Path
import numpy as np
import pyslimmc
from pyslimmc.core import DataUnavailableError

def wt(root,name,cols):
 d=root/name; d.mkdir(parents=True)
 for k,v in cols.items(): np.save(d/f'{k}.npy',v,allow_pickle=False)

def build(root):
 meta={'run_id':'run','storage':'slimmc-storage','storage_format_version':'1.2.0','run_status':'completed','validation_error_count':0,'engine':'slimmc-copo','kinetic_model':'copo'}
 (root/'run_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
 schema=[{'record_type':'schema_header','schema_name':'slimmc-storage','schema_version':'1.2.0'}]
 for t in ('snapshots','state','channel_events','actions','action_conditions','kinetic_parameters/sets','kinetic_parameters/values'): schema.append({'record_type':'table','name':t,'required':True})
 ds={
 'monomers':[(0,'A',{}),(1,'B',{})],
 'state_entities':[(0,'monomer_A',{'kind':'monomer'}),(1,'monomer_B',{'kind':'monomer'})],
 'channels':[(0,'prop_AA',{}),(1,'prop_AB',{})],
 'action_types':[(0,'set_k',{}),(1,'print',{})],
 'action_triggers':[(0,'unknown',{}),(1,'at',{}),(2,'every',{}),(3,'when',{})],
 'condition_operators':[(1,'>',{})],
 'condition_observables':[(1,'conversion',{})],
 'kinetic_parameter_definitions':[(0,'temperature_K',{'kind':'temperature','unit':'K'}),(1,'kp_aa',{'kind':'rate_constant'}),(2,'kp_ab',{'kind':'rate_constant'}),(3,'kp_ba',{'kind':'rate_constant'}),(4,'kp_bb',{'kind':'rate_constant'})]}
 for d,vals in ds.items():
  for i,n,e in vals: schema.append({'record_type':'dictionary_entry','dictionary':d,'id':i,'name':n,**e})
 (root/'schema.jsonl').write_text(''.join(json.dumps(x,separators=(',',':'))+'\n' for x in schema))
 wt(root,'snapshots',{'snapshot_id':np.arange(3,dtype='<u8'),'time':np.array([0.,1.,2.]),'kmc_event':np.array([0,10,20],'<u8'),'snapshot_reason_id':np.array([0,1,4],'<u4'),'is_final':np.array([0,0,1],bool),'has_chains':np.zeros(3,bool),'has_sequences':np.zeros(3,bool),'kinetic_parameter_set_id':np.array([0,0,1],'<u8')})
 wt(root,'state',{'snapshot_id':np.repeat(np.arange(3,dtype='<u8'),2),'entity_id':np.tile(np.arange(2,dtype='<u4'),3),'count':np.array([80,20,60,20,40,20],'<u8'),'moles':np.array([80.,20.,60.,20.,40.,20.]),'concentration':np.array([80.,20.,60.,20.,40.,20.])})
 wt(root,'channel_events',{'snapshot_id':np.repeat(np.arange(3,dtype='<u8'),2),'channel_id':np.tile(np.arange(2,dtype='<u4'),3),'event_count':np.array([0,0,7,3,12,8],'<u8'),'productive_event_count':np.array([0,0,7,3,12,8],'<u8'),'nonproductive_event_count':np.zeros(6,'<u8')})
 wt(root,'kinetic_parameters/sets',{'kinetic_parameter_set_id':np.array([0,1],'<u8'),'start_kmc_event':np.array([0,15],'<u8'),'start_time':np.array([0.,1.5]),'has_source_action':np.array([0,1],bool),'source_action_id':np.array([0,0],'<u8')})
 vals=np.array([[300.,2.,1.,1.,4.],[320.,4.,2.,1.,4.]])
 wt(root,'kinetic_parameters/values',{'kinetic_parameter_set_id':np.repeat(np.arange(2,dtype='<u8'),5),'kinetic_parameter_id':np.tile(np.arange(5,dtype='<u4'),2),'value':vals.ravel()})
 wt(root,'actions',{'action_id':np.array([0],'<u8'),'kmc_event':np.array([15],'<u8'),'time':np.array([1.5]),'source_line':np.array([10],'<u4'),'action_type_id':np.array([0],'<u4'),'trigger_type_id':np.array([3],'<u4'),'scheduled_time':np.array([np.nan]),'target_id':np.array([0],'<u4'),'requested_value':np.array([4.]),'before_value':np.array([2.]),'after_value':np.array([4.]),'state_changed':np.array([1],bool),'output_written':np.array([0],bool),'has_snapshot':np.array([0],bool),'snapshot_id':np.array([0],'<u8'),'has_kinetic_parameter_set':np.array([1],bool),'kinetic_parameter_set_id':np.array([1],'<u8')})
 wt(root,'action_conditions',{'condition_record_id':np.array([0],'<u8'),'action_id':np.array([0],'<u8'),'condition_index':np.array([0],'<u4'),'observable_id':np.array([1],'<u4'),'operator_id':np.array([1],'<u4'),'threshold':np.array([.1]),'observed_value':np.array([.2]),'condition_met':np.array([1],bool)})
 (root/'actions'/'messages.jsonl').write_text(json.dumps({'action_id':'0','message':'changed'})+'\n')
 (root/'RESULTS_COMPLETE').write_text('slimmc-storage-v1\n')

def main():
 with tempfile.TemporaryDirectory() as td:
  root=Path(td)/'run'; root.mkdir(); build(root); run=pyslimmc.open(root)
  assert np.array_equal(run.event_counts['prop_AA'],[0,7,12])
  assert run.last.channels.event_count['prop_AB']==8
  assert np.array_equal(run.channels.interval_event_counts(),[[0,0],[7,3],[5,5]])
  assert np.allclose(run.temp,[300,300,320]) and run.last.temp==320
  assert isinstance(run.temp,np.ndarray) and not run.temp.flags.writeable
  assert np.allclose(run.temp-273.15,[26.85,26.85,46.85])
  assert isinstance(run.event_counts["prop_AA"],np.ndarray)
  assert not run.event_counts["prop_AA"].flags.writeable
  assert np.allclose(run.k['kp_aa'],[2,2,4]) and run.last.k['kp_aa']==4
  assert len(run.actions)==1 and run.actions[0].type=='set_k' and run.actions[0].trigger=='when'
  assert run.actions[0].message=='changed' and run.actions[0].conditions[0].observable=='conversion'
  try: run.F.ins['A']
  except DataUnavailableError: pass
  else: raise AssertionError('F.ins accepted without terminal propagation declarations')
 print('pyslimmc L2.6 channels/kinetics/actions: PASS')
def test_script_contract():
    main()

if __name__=='__main__': main()
