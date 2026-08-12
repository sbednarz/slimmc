from __future__ import annotations
import argparse, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]; REPO=ROOT.parent
sys.path.insert(0,str(REPO)); import pyslimmc as sl
MODELS=Path(__file__).resolve().parent/'models'

def run_model(engine,model,work):
 work.mkdir(parents=True,exist_ok=True); local=work/model.name; local.write_text(model.read_text()); cp=subprocess.run([str(engine),str(local)],capture_output=True,text=True); assert cp.returncode==0,(model.name,cp.stdout,cp.stderr); return sl.open(work/'results'/model.stem)
def ev(run):
 raw=run.channels.raw; sid=np.asarray(raw['snapshot_id'],int); cid=np.asarray(raw['channel_id'],int); v=np.asarray(raw['event_count'],np.int64); m=np.zeros((len(run.t),int(cid.max())+1),np.int64); m[sid,cid]=v; return m[-1]
def conc(run,name): return float(np.asarray(run.conc[name],float)[-1])
def chains(run):
 c=run.last.chains.all; live=int(np.asarray(c.live.count,np.int64).sum()) if len(c.live.dp) else 0; dead=int(np.asarray(c.dead.count,np.int64).sum()) if len(c.dead.dp) else 0; return live,dead

def check(engine,w):
 r=run_model(engine,MODELS/'C15_reinit.model',w/'C15'); e=ev(r); assert e[8:10].sum()>0 and e[2:4].sum()>0; live,dead=chains(r); assert live>0 and dead>0 and conc(r,'CTA')<.010
 r=run_model(engine,MODELS/'C16_transfer_h.model',w/'C16'); e=ev(r); n=int(e[6:8].sum()); assert n>0; live,dead=chains(r); assert dead==n and conc(r,'Rcta')>0
 r=run_model(engine,MODELS/'C17_term_x.model',w/'C17'); e=ev(r); n=int(e[6:8].sum()); assert n>0; live,dead=chains(r); assert dead==n and conc(r,'Cap')<.010
 for stem in ('C18_rxn_uni','C18_rxn_bidiff','C18_rxn_bisame'):
  r=run_model(engine,MODELS/f'{stem}.model',w/stem); assert int(ev(r)[0])>0 and conc(r,'X')>0
 vals={}
 for tag in ('0','025','1'):
  r=run_model(engine,MODELS/f'C19_eff_{tag}.model',w/tag); raw=r.channels.raw; vals[tag]=(int(ev(r)[0]),conc(r,'X'),conc(r,'Q'))
 assert vals['0'][1]==0.0 and vals['1'][1]>vals['025'][1]>0.0
 assert vals['0'][2]<.020 and vals['1'][2]<.020
 # Same event cap and rate: product yield follows efficiency while substrate consumption follows firings.
 ratio=vals['025'][1]/vals['1'][1]; assert .15<ratio<.35,vals

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--engine',type=Path,required=True); ns=ap.parse_args()
 with tempfile.TemporaryDirectory(prefix='copo_phase_d_') as td: check(ns.engine.resolve(),Path(td))
 print('Copo phase D black-box chemistry: PASS')
if __name__=='__main__': main()
