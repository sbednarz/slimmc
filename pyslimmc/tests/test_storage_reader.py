from __future__ import annotations
import json
from pathlib import Path
import tempfile
import numpy as np
import pyslimmc
from pyslimmc._storage import IncompleteResultsError


def write_npy_table(root: Path, name: str, columns: dict[str, np.ndarray]) -> None:
    d=root/name; d.mkdir(parents=True)
    for key,value in columns.items(): np.save(d/f'{key}.npy', value, allow_pickle=False)

def test_storage_reader_contract():
    with tempfile.TemporaryDirectory() as td:
        root=Path(td)/'run_000001'; root.mkdir()
        metadata={"run_id":"run_000001","storage":"slimmc-storage","storage_format_version":"1.2.0","run_status":"completed","validation_error_count":0,"engine":"slimmc-copo","kinetic_model":"copo"}
        (root/'run_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
        schema=[
          {"record_type":"schema_header","schema_name":"slimmc-storage","schema_version":"1.2.0"},
          {"record_type":"table","name":"snapshots","required":True},
          {"record_type":"table","name":"state","required":True},
          {"record_type":"dictionary_entry","dictionary":"state_entities","id":0,"name":"monomer_A","kind":"monomer"},
        ]
        (root/'schema.jsonl').write_text(''.join(json.dumps(x,separators=(',',':'))+'\n' for x in schema))
        write_npy_table(root,'snapshots',{"snapshot_id":np.array([0],dtype='<u8'),"time":np.array([0.0],dtype='<f8')})
        write_npy_table(root,'state',{"snapshot_id":np.array([0],dtype='<u8'),"entity_id":np.array([0],dtype='<u4'),"count":np.array([4],dtype='<u8')})
        (root/'RESULTS_COMPLETE').write_text('slimmc-storage-v1\n')
        run=pyslimmc.open(root)
        assert run.status=='completed' and run.is_ok and run.kinetic_model=='copo'
        assert run.snapshots.n_rows==1 and int(run.state.count['A'][0])==4
        assert run.dictionary('state_entities')[0]['name']=='monomer_A'
        assert isinstance(run.state.count['A'], np.ndarray) and not run.state.count['A'].flags.writeable

        metadata['run_status']='interrupted'; (root/'run_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n'); (root/'RESULTS_COMPLETE').unlink()
        try: pyslimmc.open(root)
        except IncompleteResultsError: pass
        else: raise AssertionError('incomplete run accepted without opt-in')
        run=pyslimmc.open(root,allow_incomplete=True); assert run.status=='interrupted'
    print('pyslimmc Slimmc Storage L1 reader: PASS')
