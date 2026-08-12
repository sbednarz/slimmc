#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, subprocess, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
MANIFEST=ROOT/'engine_channels'/'manifest.tsv'

def entries(group: str):
    with MANIFEST.open(encoding='utf-8') as f:
        for row in csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t', fieldnames=['group','model']):
            if row['group']==group:
                yield (MANIFEST.parent/row['model']).resolve()

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--engine', required=True)
    ap.add_argument('--group', default='channels', choices=('channels','pyslimmc_live'))
    ap.add_argument('--run', action='store_true', help='run models instead of parser --check')
    ns=ap.parse_args()
    engine=Path(ns.engine).resolve()
    selected=list(entries(ns.group))
    if not selected:
        raise SystemExit(f'no validation models for group {ns.group}')
    for model in selected:
        if not model.is_file():
            raise SystemExit(f'missing validation model: {model}')
        do_run = ns.run or ns.group == "pyslimmc_live"
        cmd=[str(engine)] + ([] if do_run else ['--check']) + [str(model)]
        print('[copo-validation]', ' '.join(cmd), flush=True)
        subprocess.run(cmd, check=True, cwd=ROOT.parent.parent)
        if ns.group == "pyslimmc_live":
            result_dir = model.parent / "results" / model.stem
            code = (
                "import sys; from pathlib import Path; "
                f"sys.path.insert(0, {str(ROOT.parent.parent)!r}); "
                "import pyslimmc; "
                f"r=pyslimmc.open({str(result_dir)!r}); "
                "assert len(r.t) > 0; assert r.last is not None"
            )
            subprocess.run([sys.executable, '-c', code], check=True)
    print(f'[copo-validation] {ns.group}: {len(selected)} PASS')
    return 0
if __name__=='__main__':
    sys.exit(main())
