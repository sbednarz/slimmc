from __future__ import annotations
import argparse, subprocess, sys, tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
REPO=ROOT.parent
sys.path.insert(0,str(REPO))
import pyslimmc as sl
from pyslimmc.core import DataUnavailableError
MODELS=Path(__file__).resolve().parent/'models'

def run_model(engine, model, work):
    work.mkdir(parents=True,exist_ok=True)
    local=work/model.name; local.write_text(model.read_text())
    cp=subprocess.run([str(engine),str(local)],cwd=work,capture_output=True,text=True)
    assert cp.returncode==0,(model.name,cp.stdout,cp.stderr)
    return sl.open(work/'results'/model.stem)

def table_counts(tab,key='motif'):
    return {str(r[key]):int(r['count']) for r in tab.rows()}

def block_counts(run):
    return {(str(r['monomer']),int(r['run_length'])):int(r['count']) for r in run.microstructure.run_lengths().rows()}

def sequence_block_counts(run):
    out={}
    for rec in run.final.chains:
        seq=rec.sequence
        if not seq: continue
        n=int(rec.count); start=0
        for i in range(1,len(seq)+1):
            if i==len(seq) or seq[i]!=seq[start]:
                key=(str(seq[start]),i-start); out[key]=out.get(key,0)+n; start=i
    return out

def assert_pool_metadata(run, penultimate):
    c=run.final.chains.all
    pools=np.asarray(c.raw['pool_id'],dtype=np.int64)
    last=np.asarray(c.raw['last_monomer_id'],dtype=np.int64)
    has_last=np.asarray(c.raw['has_last_monomer'],dtype=bool)
    prev=np.asarray(c.raw['penultimate_monomer_id'],dtype=np.int64)
    has_prev=np.asarray(c.raw['has_penultimate_monomer'],dtype=bool)
    dp=np.asarray(c.dp,dtype=np.int64)
    ca=np.asarray(c.counts['A'],dtype=np.int64); cb=np.asarray(c.counts['B'],dtype=np.int64)
    np.testing.assert_array_equal(dp,ca+cb)
    # Dictionary pool ids follow declaration order; dead pool is excluded from terminal assertions.
    for i in range(len(dp)):
        if not has_last[i]: continue
        if penultimate and has_prev[i] and dp[i]>=2 and 2 <= pools[i] < 6:
            assert pools[i] == 2 + int(prev[i])*2+int(last[i]),(i,pools[i],prev[i],last[i])
        elif penultimate and pools[i] < 2:
            assert pools[i] == int(last[i]),(i,pools[i],last[i])
        elif not penultimate and pools[i] < 2:
            assert pools[i] == int(last[i]),(i,pools[i],last[i])

def assert_all_prop_channels_fired(run, expected):
    names=list(run.channels.event_count.names)
    prop=[n for n in names if n.startswith('prop_')]
    assert len(prop)==expected,(prop,expected)
    missing=[n for n in prop if int(run.channels.event_count[n][-1])<=0]
    assert not missing,missing

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--engine',type=Path,required=True); ns=ap.parse_args()
    with tempfile.TemporaryDirectory(prefix='slimmc_copo_tm_') as td:
        w=Path(td); engine=ns.engine.resolve()
        full=run_model(engine,MODELS/'CTM01_penultimate_full.model',w/'full')
        comp=run_model(engine,MODELS/'CTM02_penultimate_composition.model',w/'composition')
        term=run_model(engine,MODELS/'CTM03_terminal_full.model',w/'terminal')

        assert_all_prop_channels_fired(full,12)
        assert_all_prop_channels_fired(term,4)
        assert_pool_metadata(full,True)
        assert_pool_metadata(comp,True)
        assert_pool_metadata(term,False)

        consistency=full.microstructure.check_sequence_consistency()
        assert consistency=={'dyads_match':True,'triads_match':True},consistency
        assert block_counts(full)==sequence_block_counts(full)

        # sequence_mode must not alter chemistry or aggregate microstructure.
        np.testing.assert_allclose(full.t,comp.t,rtol=0,atol=0)
        np.testing.assert_array_equal(full.event,comp.event)
        for name in full.channels.event_count.names:
            np.testing.assert_array_equal(full.channels.event_count[name],comp.channels.event_count[name])
        np.testing.assert_allclose(full.conc['A'],comp.conc['A'],rtol=0,atol=0)
        np.testing.assert_allclose(full.conc['B'],comp.conc['B'],rtol=0,atol=0)
        assert table_counts(full.microstructure.dyads())==table_counts(comp.microstructure.dyads())
        assert table_counts(full.microstructure.triads())==table_counts(comp.microstructure.triads())
        assert block_counts(full)==block_counts(comp)
        assert full.final.chains.has_sequences
        assert not comp.final.chains.has_sequences
        try:
            comp.microstructure.dyads(source='sequences')
        except DataUnavailableError:
            pass
        else:
            raise AssertionError('composition mode unexpectedly exposes literal sequences')

    print('Copo terminal/penultimate and microstructure validation: PASS')
if __name__=='__main__': main()
