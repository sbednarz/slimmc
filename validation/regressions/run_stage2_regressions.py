#!/usr/bin/env python3
from __future__ import annotations
import argparse, os, shutil, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def run(cmd, *, env=None):
    cp=subprocess.run(cmd,text=True,capture_output=True,env=env)
    if cp.returncode:
        raise RuntimeError(f"{' '.join(map(str,cmd))}\n{cp.stdout[-1000:]}\n{cp.stderr[-1000:]}")
    return cp

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--slimmc',required=True); args=ap.parse_args()
    run([sys.executable,str(ROOT/'run_regressions.py'),'--slimmc',args.slimmc])
    for name in ('R007_copo_missing_prop_zero','R008_pyslimmc_Fins_topology'):
        case=ROOT/name; shutil.rmtree(case/'results',ignore_errors=True)
        if name.startswith('R007'):
            check=run([args.slimmc,'--check',str(case/'model.model')])
            assert 'missing macro prop transitions are treated as k=0' in check.stdout
        run([args.slimmc,str(case/'model.model')])
    env=dict(os.environ); env['PYTHONPATH']=str(ROOT.parents[1])
    run([sys.executable,str(ROOT/'R008_pyslimmc_Fins_topology/check.py')],env=env)
    run([sys.executable,str(ROOT/'R009_pyslimmc_helpers/check.py')],env=env)
    run([sys.executable,str(ROOT/'R009_pyslimmc_helpers/check_feed.py')],env=env)
    print('[PASS] stage2 regressions R001-R009')
    return 0
if __name__=='__main__': raise SystemExit(main())
