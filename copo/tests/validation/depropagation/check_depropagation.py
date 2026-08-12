from __future__ import annotations
import argparse, subprocess, sys, tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
REPO=ROOT.parent
sys.path.insert(0,str(REPO))
import pyslimmc as sl
MODELS=Path(__file__).resolve().parent/'models'

def run_model(engine:Path, model:Path, work:Path):
    work.mkdir(parents=True,exist_ok=True)
    local=work/model.name
    local.write_text(model.read_text())
    cp=subprocess.run([str(engine),str(local)],capture_output=True,text=True)
    assert cp.returncode==0,(model.name,cp.stdout,cp.stderr)
    return sl.open(work/'results'/model.stem)

def chain_units(snapshot):
    c=snapshot.chains.all
    a=int(np.sum(np.asarray(c.count,dtype=np.int64)*np.asarray(c.counts['A'],dtype=np.int64)))
    b=int(np.sum(np.asarray(c.count,dtype=np.int64)*np.asarray(c.counts['B'],dtype=np.int64)))
    return np.array([a,b],dtype=np.int64)

def deprop_matrix(run):
    raw=run.channels.raw
    sid=np.asarray(raw['snapshot_id'],dtype=np.int64)
    cid=np.asarray(raw['channel_id'],dtype=np.int64)
    val=np.asarray(raw['event_count'],dtype=np.int64)
    n=len(run.t); m=int(cid.max())+1
    mat=np.zeros((n,m),dtype=np.int64)
    mat[sid,cid]=val
    # declaration order: 0-1 init, 2-5 prop, 6-9 deprop
    return mat[:,6:10]

def snapshot_index_at(run,t,require_chains=False):
    idx=np.flatnonzero(np.isclose(np.asarray(run.t,float),t,rtol=0,atol=1e-12))
    if require_chains:
        has=np.asarray(run.snapshots.raw['has_chains'],dtype=bool)
        idx=idx[has[idx]]
    assert idx.size>=1,(t,run.t)
    return int(idx[-1])

def snapshot_at(run,t):
    return run.snapshots[snapshot_index_at(run,t,require_chains=True)]

def check_isolated_and_control(engine,work):
    test=run_model(engine,MODELS/'CDEP01_isolated.model',work/'test')
    ctrl=run_model(engine,MODELS/'CDEP02_control.model',work/'control')
    s0=snapshot_at(test,.10); sf=test.last
    u0=chain_units(s0); uf=chain_units(sf)
    dep=deprop_matrix(test)
    i0=snapshot_index_at(test,.10,require_chains=True)
    d=dep[-1]-dep[i0]
    # channels 6,7 release A; 8,9 release B
    released=np.array([d[0]+d[1],d[2]+d[3]],dtype=np.int64)
    np.testing.assert_array_equal(u0-uf,released)
    assert np.all(released>0),released
    # Every transition type AA, BA, AB and BB must occur in this mixed-sequence run.
    assert np.all(d>0),d
    # Composition and terminal bookkeeping remain internally consistent.
    c=sf.chains.all
    np.testing.assert_array_equal(np.asarray(c.dp),np.asarray(c.counts['A'])+np.asarray(c.counts['B']))
    for pool,last_id in [('PA',0),('PB',1)]:
        sub=sf.chains.pool(pool)
        if len(sub.dp):
            assert np.all(np.asarray(sub.raw['last_monomer_id'],dtype=np.int64)==last_id)
    # Matched no-deprop control: no deprop fires and chain-unit inventory does not fall after switch.
    c0=chain_units(snapshot_at(ctrl,.10)); cf=chain_units(ctrl.last)
    assert int(deprop_matrix(ctrl)[-1].sum())==0
    assert np.all(cf>=c0),(c0,cf)

def check_scaling(engine,work):
    rates=[]
    for tag,kd in [('10',10.),('20',20.),('40',40.)]:
        r=run_model(engine,MODELS/f'CDEP03_kd_{tag}.model',work/tag)
        mat=deprop_matrix(r)
        i0=snapshot_index_at(r,.10,require_chains=True)
        n=int((mat[-1]-mat[i0]).sum())
        rates.append(n/.04)
    # Statistical black-box scaling: monotonic and approximately proportional.
    assert rates[0]>0 and rates[0]<rates[1]<rates[2],rates
    assert 1.45 < rates[1]/rates[0] < 2.75,rates
    assert 1.45 < rates[2]/rates[1] < 2.75,rates

def check_competition(engine,work):
    rp=run_model(engine,MODELS/'CDEP04_prop_dominant.model',work/'prop')
    rd=run_model(engine,MODELS/'CDEP04_deprop_dominant.model',work/'deprop')
    p0=chain_units(snapshot_at(rp,.10)).sum(); pf=chain_units(rp.last).sum()
    d0=chain_units(snapshot_at(rd,.10)).sum(); df=chain_units(rd.last).sum()
    assert pf>p0,(p0,pf)
    assert df<d0,(d0,df)
    assert deprop_matrix(rd)[-1].sum()>deprop_matrix(rp)[-1].sum()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--engine',type=Path,required=True); ns=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='slimmc_copo_deprop_') as td:
        w=Path(td)
        check_isolated_and_control(ns.engine.resolve(),w/'D01_D02')
        check_scaling(ns.engine.resolve(),w/'D03')
        check_competition(ns.engine.resolve(),w/'D04')
    print('Copo detailed depropagation validation: PASS')
if __name__=='__main__': main()
