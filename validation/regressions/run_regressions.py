#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, shutil, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
CASES={
 'R001_homo_termination_tend': ('t_end',),
 'R002_homo_termination_stop': ('stop_condition',),
 'R003_homo_semibatch_chain_volume': ('passed',),
 'R004_copo_zero_propensity_future_feed': ('feed',),
 'R005_copo_full_feed_transfer_reinit': ('passed',),
 'R006_copo_terpoly_feed_species_padding': ('passed',),
}

def metadata(case: Path) -> dict:
    p=case/'results'/'main'/'run_metadata.json'
    if not p.exists(): raise AssertionError(f'missing {p}')
    return json.loads(p.read_text())

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--slimmc',default='slimmc'); a=ap.parse_args()
    failures=[]
    for name in CASES:
        case=ROOT/name; out=case/'results'
        if out.exists(): shutil.rmtree(out)
        model=case/'model.model'
        cp=subprocess.run([a.slimmc,str(model)],text=True,capture_output=True)
        if cp.returncode:
            failures.append((name,f'run rc={cp.returncode}: {cp.stderr[-500:]}')); continue
        try:
            m=metadata(case)
            if name=='R001_homo_termination_tend': assert m.get('termination_reason')=='t_end',m
            elif name=='R002_homo_termination_stop': assert m.get('termination_reason')=='stop_condition',m
            elif name=='R003_homo_semibatch_chain_volume': assert m.get('validation_status')=='passed',m
            elif name=='R004_copo_zero_propensity_future_feed':
                import numpy as np
                dose=np.load(case/'results'/'main'/'feed_events'/'dose_mL.npy')
                assert dose.size==1 and abs(float(dose[0])-0.01)<1e-12,dose
            else: assert m.get('validation_status')=='passed',m
            print(f'[PASS] {name}')
        except Exception as e: failures.append((name,repr(e)))
    if failures:
        for n,e in failures: print(f'[FAIL] {n}: {e}',file=sys.stderr)
        return 1
    return 0
if __name__=='__main__': raise SystemExit(main())
