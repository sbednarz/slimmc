from __future__ import annotations
import argparse, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; REPO=ROOT.parent
sys.path.insert(0,str(REPO)); import pyslimmc as sl
MODELS=Path(__file__).resolve().parent/'models'

def run_model(engine,model,work):
    work.mkdir(parents=True,exist_ok=True); local=work/model.name; local.write_text(model.read_text())
    cp=subprocess.run([str(engine),str(local)],capture_output=True,text=True)
    assert cp.returncode==0,(model.name,cp.stdout,cp.stderr)
    return sl.open(work/'results'/model.stem)

def events(run):
    raw=run.channels.raw; sid=np.asarray(raw['snapshot_id'],int); cid=np.asarray(raw['channel_id'],int); val=np.asarray(raw['event_count'],np.int64)
    m=np.zeros((len(run.t),int(cid.max())+1),np.int64); m[sid,cid]=val; return m

def idx_at(run,t):
    ii=np.flatnonzero(np.isclose(np.asarray(run.t,float),t,atol=1e-12,rtol=0)); assert ii.size
    hc=np.asarray(run.snapshots.raw['has_chains'],bool); jj=ii[hc[ii]]; assert jj.size,(t,run.t,hc)
    return int(jj[-1])
def at(run,t): return run.snapshots[idx_at(run,t)]

def summary(s):
    c=s.chains.all; n=np.asarray(c.count,np.int64); dp=np.asarray(c.dp,np.int64)
    A=np.asarray(c.counts['A'],np.int64); B=np.asarray(c.counts['B'],np.int64)
    live=int(np.asarray(c.live.count,np.int64).sum()) if len(c.live.dp) else 0
    dead=int(np.asarray(c.dead.count,np.int64).sum()) if len(c.dead.dp) else 0
    return dict(chains=int(n.sum()), units=int((n*dp).sum()), A=int((n*A).sum()), B=int((n*B).sum()), live=live, dead=dead)

def check(engine,work):
    r=run_model(engine,MODELS/'C01_init.model',work/'C01'); e=events(r)[-1]; s=summary(r.last)
    assert e[:2].sum()>0 and s['chains']==int(e[:2].sum()) and s['units']==s['chains']
    assert s['A']==int(e[0]) and s['B']==int(e[1]) and s['dead']==0

    r=run_model(engine,MODELS/'C02_prop.model',work/'C02'); e=events(r); i=idx_at(r,.04)
    s0,sf=summary(at(r,.04)),summary(r.last); d=e[-1]-e[i]
    assert s0['live']==sf['live'] and sf['dead']==0
    assert sf['A']-s0['A']==int(d[2]+d[4]); assert sf['B']-s0['B']==int(d[3]+d[5])
    assert sf['units']-s0['units']==int(d[2:6].sum())

    r=run_model(engine,MODELS/'C03_term_c.model',work/'C03'); e=events(r); i=idx_at(r,.05)
    s0,sf=summary(at(r,.05)),summary(r.last); n=int((e[-1]-e[i])[6:9].sum()); assert n>0
    assert s0['live']-sf['live']==2*n and sf['dead']-s0['dead']==n and sf['units']==s0['units']

    r=run_model(engine,MODELS/'C04_term_d.model',work/'C04'); e=events(r); i=idx_at(r,.05)
    s0,sf=summary(at(r,.05)),summary(r.last); n=int((e[-1]-e[i])[6:9].sum()); assert n>0
    assert s0['live']-sf['live']==2*n and sf['dead']-s0['dead']==2*n and sf['units']==s0['units']

    r=run_model(engine,MODELS/'C05_transfer_m.model',work/'C05'); e=events(r); i=idx_at(r,.05)
    s0,sf=summary(at(r,.05)),summary(r.last); d=(e[-1]-e[i])[6:10]; n=int(d.sum()); assert n>0
    assert sf['live']==s0['live'] and sf['dead']-s0['dead']==n and sf['chains']-s0['chains']==n
    assert sf['A']-s0['A']==int(d[0]+d[2]); assert sf['B']-s0['B']==int(d[1]+d[3]); assert sf['units']-s0['units']==n
    for x in (s0,sf): assert x['units']==x['A']+x['B'] and min(x.values())>=0

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--engine',type=Path,required=True); ns=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='copo_phase_a_') as td: check(ns.engine.resolve(),Path(td))
    print('Copo phase A black-box chemistry: PASS')
if __name__=='__main__': main()
