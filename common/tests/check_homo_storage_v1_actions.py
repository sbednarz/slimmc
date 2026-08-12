from __future__ import annotations
import json, subprocess, tempfile
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'homo/tests/regression/frp_all_channels_seeded/FRP_ALLCHAN01.model'
BIN = ROOT/'homo/slimmc-stage-h'
with tempfile.TemporaryDirectory() as td:
    td=Path(td); model=td/'run_actions.model'
    text=BASE.read_text().replace('every 1.0 save\nevery 1.0 save_chains\n','')
    text += '\nat 0.20 set_k kp 450\nat 0.30 add_temp 5\nat 0.40 set_c CTA 0.0003\nat 0.50 print "stage f message"\nat 0.60 save\nwhen X M > -0.1 and c M < 1000 print "and fired"\n'
    model.write_text(text)
    subprocess.run([str(BIN),str(model)],check=True,capture_output=True,text=True)
    run=td/'results/run_actions'
    a={p.stem:np.load(p,allow_pickle=False) for p in (run/'actions').glob('*.npy')}
    assert np.array_equal(a['action_id'],np.arange(6,dtype=np.uint64))
    assert np.array_equal(a['action_type_id'],np.array([0,5,8,9,0,2],dtype=np.uint32))
    assert np.array_equal(a['trigger_type_id'],np.array([3,1,1,1,1,1],dtype=np.uint32))
    assert np.allclose(a['scheduled_time'][1:],[.2,.3,.4,.5,.6])
    assert np.isnan(a['scheduled_time'][0])
    assert np.array_equal(a['state_changed'],[False,True,True,True,False,False])
    assert np.array_equal(a['has_kinetic_parameter_set'],[False,True,True,False,False,False])
    assert np.array_equal(a['kinetic_parameter_set_id'][1:3],[1,2])
    assert np.array_equal(a['has_snapshot'],[False,True,True,False,False,True])
    assert np.array_equal(a['output_written'],[False,True,True,False,False,True])
    assert np.isclose(a['requested_value'][1],450)
    assert np.isclose(a['requested_value'][2],5)
    assert np.isclose(a['requested_value'][3],0.0003)
    recs=[json.loads(x) for x in (run/'actions/messages.jsonl').read_text().splitlines()]
    assert recs==[{'action_id':'0','message':'and fired'},{'action_id':'4','message':'stage f message'}]
    c={p.stem:np.load(p,allow_pickle=False) for p in (run/'action_conditions').glob('*.npy')}
    assert np.array_equal(c['condition_record_id'],np.array([0,1],dtype=np.uint64))
    assert np.array_equal(c['action_id'],np.array([0,0],dtype=np.uint64))
    assert np.array_equal(c['condition_index'],np.array([0,1],dtype=np.uint32))
    assert np.array_equal(c['observable_id'],np.array([1,2],dtype=np.uint32))
    assert np.array_equal(c['operator_id'],np.array([1,2],dtype=np.uint32))
    assert np.allclose(c['threshold'],[-0.1,1000.0])
    assert np.all(c['condition_met'])
    assert np.all(c['observed_value'][0:] == c['observed_value'][0:])
    ks=np.load(run/'kinetic_parameters/sets/source_action_id.npy',allow_pickle=False)
    hk=np.load(run/'kinetic_parameters/sets/has_source_action.npy',allow_pickle=False)
    assert np.array_equal(ks[hk],[1,2])
print('Homo Slimmc Storage v1 actions detailed contract: OK')
