from __future__ import annotations
import json, tempfile
from pathlib import Path
import numpy as np
import pyslimmc


def wt(root: Path, name: str, cols: dict[str, np.ndarray]) -> None:
    d=root/name; d.mkdir(parents=True,exist_ok=True)
    for k,v in cols.items(): np.save(d/f"{k}.npy",v,allow_pickle=False)


def build(root: Path, value: float) -> None:
    root.mkdir(parents=True)
    meta={"run_id":root.name,"storage":"slimmc-storage","storage_format_version":"1.2.0",
          "run_status":"completed","validation_status":"passed","validation_error_count":0,
          "validation_warning_count":0,"engine":"slimmc-copo","kinetic_model":"copo",
          "variables":[{"kind":"rate","name":"kp","value":value,"unit":"L_mol_s"},
                       {"kind":"param","name":"temperature","value":333.15,"unit":"K"}]}
    (root/'run_metadata.json').write_text(json.dumps(meta))
    schema=[
      {"record_type":"schema_header","schema_name":"slimmc-storage","schema_version":"1.2.0"},
      {"record_type":"table","name":"snapshots","required":True},
      {"record_type":"table","name":"state","required":True},
      {"record_type":"table","name":"channel_events","required":False},
      {"record_type":"table","name":"moments","required":False},
      {"record_type":"table","name":"kinetic_parameters/values","required":False},
      {"record_type":"table","name":"actions","required":False},
      {"record_type":"table","name":"memory","required":False},
      {"record_type":"dictionary_entry","dictionary":"state_entities","id":0,"name":"monomer_A","kind":"monomer"},
      {"record_type":"dictionary_entry","dictionary":"monomers","id":0,"name":"A"},
      {"record_type":"dictionary_entry","dictionary":"channels","id":0,"name":"prop_AA"},
      {"record_type":"dictionary_entry","dictionary":"population_scope","id":0,"name":"all"},
      {"record_type":"dictionary_entry","dictionary":"mass_bases","id":1,"name":"with_end_groups"},
      {"record_type":"dictionary_entry","dictionary":"kinetic_parameter_definitions","id":0,"name":"temperature","kind":"temperature","unit":"K"},
    ]
    (root/'schema.jsonl').write_text(''.join(json.dumps(x)+'\n' for x in schema))
    wt(root,'snapshots',{"snapshot_id":np.array([0,1],dtype='<u8'),"time":np.array([0.,1.]),
       "kmc_event":np.array([0,2],dtype='<u8'),"is_final":np.array([0,1],dtype=bool),
       "has_chains":np.array([0,0],dtype=bool),"has_sequences":np.array([0,0],dtype=bool),
       "kinetic_parameter_set_id":np.array([0,0],dtype='<u8')})
    wt(root,'state',{"snapshot_id":np.array([0,1],dtype='<u8'),"entity_id":np.array([0,0],dtype='<u4'),
       "count":np.array([100,50],dtype='<u8'),"moles":np.array([1.,.5]),"concentration":np.array([1.,.5])})
    wt(root,'channel_events',{"snapshot_id":np.array([0,1],dtype='<u8'),"channel_id":np.array([0,0],dtype='<u4'),
       "event_count":np.array([0,2],dtype='<u8'),"productive_event_count":np.array([0,2],dtype='<u8'),"nonproductive_event_count":np.array([0,0],dtype='<u8')})
    wt(root,'moments',{"snapshot_id":np.array([1],dtype='<u8'),"population_scope_id":np.array([0],dtype='<u4'),"mass_basis_id":np.array([1],dtype='<u4'),
       "dpn":np.array([10.]),"dpw":np.array([12.]),"mn":np.array([1000.]),"mw":np.array([1200.]),"mz":np.array([1400.]),"dispersity":np.array([1.2])})
    wt(root,'kinetic_parameters/values',{"kinetic_parameter_set_id":np.array([0],dtype='<u8'),"kinetic_parameter_id":np.array([0],dtype='<u4'),"value":np.array([333.15])})
    wt(root,'actions',{"action_id":np.array([],dtype='<u8'),"time":np.array([],dtype='<f8'),"kmc_event":np.array([],dtype='<u8')})
    wt(root,'memory',{"snapshot_id":np.array([0,1],dtype='<u8'),"total_est_B":np.array([100.,120.])})
    (root/'input.model').write_text(f'desc "sweep case"\nvar rate kp L_mol_s\nrate kp {value}\n')
    (root/'diagnostics').mkdir(); (root/'diagnostics'/'validation.jsonl').write_text('')
    (root/'RESULTS_COMPLETE').write_text('ok\n')


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        root=Path(td)/'results'
        build(root/'group'/'run_1',100.0); build(root/'group'/'run_2',200.0)
        runs=pyslimmc.scan(root)
        assert len(runs)==2
        assert runs['run_1'].run_id=='run_1'
        assert runs.var.keys()==('kp','temperature')
        sw=runs.sweep('kp'); assert [r.var['kp'].value for r in sw]==[100.0,200.0]; assert sw.sweep_variables==('kp',); assert sw.info().startswith('Runs')
        text=runs.info(); assert 'group/run_1/' in text and str(root.resolve()) not in text
        run=runs[0]; assert run.var['kp'].value==100.0 and run.var['temperature'].kind=='param'
        assert not hasattr(run, 'var_name') and run.var.info().startswith('Variables')
        info=run.info(); assert 'Common next steps:' in info and 'conversion total: 0.5' in info
        snap=run.last; assert snap.conv.total==0.5 and 'snap.mwd()' in snap.info()
        series=run.conv.total; assert isinstance(series,np.ndarray) and not series.flags.writeable
        assert run.moments.info().startswith('MomentsView')
        assert run.channels.info().startswith('ChannelsView')
        assert run.kinetics.info().startswith('KineticsView')
        assert run.actions.info().startswith('ActionsView')
        assert run.diagnostics.info().startswith('DiagnosticsView')
        assert run.raw.info().startswith('RawView')
    print('pyslimmc API audit / Slimmc Storage scan, sweep and info: PASS')

def test_script_contract():
    main()

if __name__=='__main__': main()
