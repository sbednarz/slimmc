from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

from .core import (DataUnavailableError, FeatureUnavailableError,
                   ChemicalAnalysisNotApplicableError,
                   ChemicalModelIncompatibleError)
from .run import DataConsistencyError
from .snapshots import NamedValues, _readonly
from .copolymerization import (
    Capability, Capabilities, PairValues, ReactivityRatioSeries,
    MayoLewisSeries, MayoLewisComparison, CompositionResult,
    SnapshotCompositionSeries, IntervalCompositionSeries, CompositionDrift,
    TerminalTransitionDiagnostics, TerminalBlockDiagnostics, TerminalDiagnostics,
    PenultimateParameterSeries, PenultimateCompositionSeries, PenultimateComparison,
    PenultimateDiagnostics, TripleValues, _penultimate_prediction,
)


class StorageMicrostructure:
    """Aggregate dyads, triads and block histograms stored by the copo engine.

    These data are independent of literal sequence storage and therefore remain
    available in ``sequence_mode=composition``.  In ``full`` mode the same
    object can independently reconstruct motifs from saved sequences.
    """
    def __init__(self, run):
        self.run = run

    def _snapshot_id(self, snapshot=None) -> int:
        if snapshot is None or snapshot == "final":
            return int(self.run.final.id)
        if snapshot == "last":
            return int(self.run.last.id)
        if hasattr(snapshot, "id"):
            return int(snapshot.id)
        return int(snapshot)

    def _motif_counts(self, order: int, snapshot=None) -> dict[str, float]:
        from .table import Table
        if "microstructure_motifs" not in self.run.tables:
            raise DataUnavailableError("microstructure_motifs table is unavailable")
        tab = self.run.table("microstructure_motifs")
        sid = self._snapshot_id(snapshot)
        mask = (np.asarray(tab["snapshot_id"], dtype=np.uint64) == sid) & (np.asarray(tab["motif_order"], dtype=np.uint32) == order)
        ids = np.asarray(tab["motif_id"], dtype=np.uint32)[mask]
        counts = np.asarray(tab["count"], dtype=np.uint64)[mask]
        dictionary = self.run.dictionary("microstructure_dyads" if order == 2 else "microstructure_triads")
        out = {}
        for motif_id, count in zip(ids, counts):
            meta = dictionary.get(int(motif_id))
            if meta is None:
                raise DataConsistencyError(f"unknown order-{order} motif_id {int(motif_id)}")
            out[str(meta["name"])] = float(count)
        return out

    @staticmethod
    def _motif_table(counts: Mapping[str, float], name: str):
        from .table import Table
        total = float(sum(counts.values()))
        rows = [(motif, value, value / total if total > 0 else 0.0) for motif, value in sorted(counts.items())]
        return Table(("motif", "count", "fraction"), rows, name=name)

    def _sequence_counts(self, width: int, snapshot=None) -> dict[str, float]:
        chains = self.run.final.chains if snapshot in (None, "final") else self.run.last.chains if snapshot == "last" else self.run.at_snapshot(self._snapshot_id(snapshot)).chains
        if not chains.has_sequences:
            raise DataUnavailableError("complete stored sequences are required for source='sequences'")
        out: dict[str, float] = {}
        for record in chains:
            try:
                seq = record.sequence
            except (DataUnavailableError, DataConsistencyError, KeyError, ValueError) as exc:
                rid = getattr(record, "chain_record_id", "?")
                raise DataUnavailableError(
                    f"cannot reconstruct sequence statistics: chain record {rid} has unreadable sequence data"
                ) from exc
            if seq is None or len(seq) < width:
                continue
            tokens = [str(x) for x in seq]
            for i in range(len(tokens) - width + 1):
                motif = "|".join(tokens[i:i+width])
                out[motif] = out.get(motif, 0.0) + float(record.count)
        return out

    def dyads(self, *, source="engine", snapshot=None):
        if source in {"engine", "table"}:
            counts = self._motif_counts(2, snapshot)
        elif source in {"sequences", "chains"}:
            counts = self._sequence_counts(2, snapshot)
        else:
            raise ValueError("source must be 'engine' or 'sequences'")
        return self._motif_table(counts, f"dyads_{source}")

    def triads(self, *, source="engine", snapshot=None):
        if source in {"engine", "table"}:
            counts = self._motif_counts(3, snapshot)
        elif source in {"sequences", "chains"}:
            counts = self._sequence_counts(3, snapshot)
        else:
            raise ValueError("source must be 'engine' or 'sequences'")
        return self._motif_table(counts, f"triads_{source}")

    def run_lengths(self, monomer: str | None = None, *, snapshot=None):
        from .table import Table
        if "block_statistics" not in self.run.tables:
            raise DataUnavailableError("block_statistics table is unavailable")
        tab = self.run.table("block_statistics")
        sid = self._snapshot_id(snapshot)
        names = tuple(self.run.monomer_names)
        mask = np.asarray(tab["snapshot_id"], dtype=np.uint64) == sid
        mids = np.asarray(tab["monomer_id"], dtype=np.uint32)[mask]
        lengths = np.asarray(tab["block_length"], dtype=np.uint64)[mask]
        counts = np.asarray(tab["block_count"], dtype=np.uint64)[mask]
        rows=[]
        total=float(np.sum(counts))
        for mid, length, count in zip(mids, lengths, counts):
            if int(mid) >= len(names):
                raise DataConsistencyError(f"invalid block monomer_id {int(mid)}")
            name=names[int(mid)]
            if monomer is None or name == monomer:
                rows.append((name, int(length), int(count), float(count)/total if total else 0.0))
        return Table(("monomer", "run_length", "count", "fraction"), rows, name="run_lengths")

    def transition_fraction(self, *, source="engine") -> float:
        rows=self.dyads(source=source).rows(); total=sum(float(r["count"]) for r in rows)
        return (sum(float(r["count"]) for r in rows if str(r["motif"]).split('|')[0] != str(r["motif"]).split('|')[1]) / total) if total else 0.0

    def homodyad_fraction(self, *, source="engine") -> float:
        rows=self.dyads(source=source).rows(); total=sum(float(r["count"]) for r in rows)
        return (sum(float(r["count"]) for r in rows if str(r["motif"]).split('|')[0] == str(r["motif"]).split('|')[1]) / total) if total else 0.0

    def blockiness(self, *, source="engine") -> dict[str, float]:
        return {"homodyad_fraction": self.homodyad_fraction(source=source), "transition_fraction": self.transition_fraction(source=source)}

    def check_sequence_consistency(self, *, snapshot=None) -> dict[str, bool]:
        def counts(table): return {str(r["motif"]): float(r["count"]) for r in table.rows()}
        def equal(a,b): return all(abs(a.get(k,0.0)-b.get(k,0.0)) <= 1e-12 for k in set(a)|set(b))
        return {
            "dyads_match": equal(counts(self.dyads(snapshot=snapshot)), counts(self.dyads(source="sequences", snapshot=snapshot))),
            "triads_match": equal(counts(self.triads(snapshot=snapshot)), counts(self.triads(source="sequences", snapshot=snapshot))),
        }

    def info_text(self):
        return (f"Microstructure\n  dyad motifs: {len(self.dyads())}\n  triad motifs: {len(self.triads())}\n"
                f"  transition fraction: {self.transition_fraction():.6g}\n\nCommon next steps:\n"
                "  run.microstructure.dyads()\n  run.microstructure.triads()\n  run.microstructure.run_lengths()")

    def info(self):
        text=self.info_text(); print(text); return text


class StorageFirings:
    """Firing-count facade reconstructed directly from ``channel_events/``.

    Counts are cumulative on the snapshot axis.  Fire shares describe realised
    SSA events; they are deliberately not labelled as rate/propensity shares.
    """
    def __init__(self, run):
        self.run = run
        self._channels = run.channels

    @property
    def _series(self):
        return self._channels

    def channels(self) -> list[str]:
        return list(self._series.event_count.names)

    def rows(self) -> list[dict[str, Any]]:
        rows=[]
        names=self.channels()
        for i,(sid,t) in enumerate(zip(self.run.sid,self.run.t)):
            row={"snapshot_id":int(sid),"t":float(t)}
            for name in names:
                row["fires_"+name]=int(self._series.event_count[name][i])
            row["total_fires"]=sum(row["fires_"+name] for name in names)
            rows.append(row)
        return rows

    def final_row(self) -> dict[str, Any]:
        rows=self.rows()
        return rows[-1] if rows else {}

    def final_fires(self) -> dict[str,int]:
        row=self.final_row()
        return {name:int(row.get("fires_"+name,0)) for name in self.channels()}

    def total_fires(self) -> int:
        values=self.final_fires()
        total=sum(values.values())
        if len(self.run.event) and total != int(self.run.event[-1]):
            raise DataConsistencyError(
                f"channel event total {total} disagrees with final kmc_event {int(self.run.event[-1])}"
            )
        return total

    def _check(self, channel: str) -> None:
        if channel not in self.channels():
            raise KeyError(f"unknown channel {channel!r}; declared channels: {self.channels()}")

    def channel_fires(self, channel: str | None=None):
        if channel is None:
            return self.final_fires()
        self._check(channel)
        return int(self._series.event_count[channel][-1])

    def delta_fires_series(self, channel: str) -> np.ndarray:
        self._check(channel)
        values=np.asarray(self._series.event_count[channel],dtype=np.int64)
        out=np.diff(values)
        out.flags.writeable=False
        return out

    def delta_fires(self, channel: str | None=None):
        if channel is not None:
            values=self.delta_fires_series(channel)
            return int(values[-1]) if len(values) else 0
        return {name:self.delta_fires(name) for name in self.channels()}

    def fire_shares(self) -> dict[str,float]:
        values=self.final_fires(); total=sum(values.values())
        return {name:(value/total if total else np.nan) for name,value in values.items()}

    def fire_shares_series(self) -> dict[str,np.ndarray]:
        names=self.channels()
        deltas={name:self.delta_fires_series(name).astype(float) for name in names}
        total=np.sum(np.column_stack([deltas[n] for n in names]),axis=1) if names else np.empty(0)
        result={}
        for name in names:
            arr=np.divide(deltas[name],total,out=np.full(total.shape,np.nan),where=total>0)
            arr.flags.writeable=False; result[name]=arr
        return result

    def rate_shares_series(self) -> dict[str, np.ndarray]:
        if "channel_propensities" not in self.run.tables:
            raise DataUnavailableError("channel_propensities table is unavailable")
        tab=self.run.table("channel_propensities")
        sids=np.asarray(tab["snapshot_id"],dtype=np.uint64); cids=np.asarray(tab["channel_id"],dtype=np.uint32)
        prop=np.asarray(tab["propensity"],dtype=float); total=np.asarray(tab["total_propensity"],dtype=float)
        names=self.channels(); out={}
        for cid,name in enumerate(names):
            mask=cids==cid; values=prop[mask]; totals=total[mask]
            if len(values)!=len(self.run.sid) or not np.array_equal(sids[mask],np.asarray(self.run.sid,dtype=np.uint64)):
                raise DataConsistencyError(f"channel_propensities axis mismatch for {name}")
            share=np.divide(values,totals,out=np.full(values.shape,np.nan),where=totals>0); share.flags.writeable=False; out[name]=share
        return out

    def rate_shares(self) -> dict[str,float]:
        series=self.rate_shares_series()
        return {name:(float(values[-1]) if len(values) else np.nan) for name,values in series.items()}

    propensity_shares = rate_shares
    propensity_shares_series = rate_shares_series

    def validate(self):
        counts=np.column_stack([np.asarray(self._series.event_count[n],dtype=np.int64) for n in self.channels()])
        if counts.size and np.any(np.diff(counts,axis=0)<0):
            raise DataConsistencyError("channel event counts are not monotonically non-decreasing")
        totals=np.sum(counts,axis=1) if counts.size else np.zeros(len(self.run.event),dtype=np.int64)
        if not np.array_equal(totals,np.asarray(self.run.event,dtype=np.int64)):
            raise DataConsistencyError("sum of channel event counts disagrees with kmc_event")
        return True

    def info_text(self):
        return (f"Firings\n  channels: {len(self.channels())}\n  total fires: {self.total_fires()}\n\n"
                "Common next steps:\n  run.firings.fire_shares()\n"
                "  run.firings.fire_shares_series()\n  run.firings.delta_fires_series(channel)")

    def info(self):
        text=self.info_text(); print(text); return text


class StorageCopolymerization:
    """Copolymer analysis facade backed only by canonical Storage tables.

    Composition, terminal reactivity ratios and Mayo--Lewis are implemented.
    Explicit penultimate diagnostics remain capability-gated until Storage
    carries the required model/channel semantic dictionary and motif tables.
    """
    def __init__(self, run):
        self._run=run
        self._composition=None
        self._rr=None
        self._ml=None

    def _model(self):
        rates={}; props=[]; operations=[]
        pools={}
        for raw in self._run._input_model_lines():
            line=raw.split('#',1)[0].strip()
            if not line: continue
            p=line.split()
            if len(p)>=3 and p[0] in {'rate','arrhenius'}:
                rates[p[1]]=p[1]
            elif len(p)>=3 and p[0]=='polymer':
                pools[p[1]]=p[2]
            elif len(p)>=2 and p[0]=='macro':
                operations.append(p[1])
                if p[1]=='prop' and '->' in p:
                    arrow=p.index('->')
                    # macro prop PA + B -> PB kp_ab
                    if arrow>=4 and len(p)>arrow+2:
                        props.append(dict(pool=p[2],incoming=p[4],product=p[arrow+1],rate=p[arrow+2]))
        return rates,props,operations,pools

    def _terminal_rate_series(self):
        monomers=tuple(self._run.monomer_names)
        if len(monomers)!=2:
            raise ChemicalAnalysisNotApplicableError("terminal analysis requires exactly two monomers")
        _,props,_,_=self._model()
        if len(props)!=4:
            raise ChemicalModelIncompatibleError(
                f"binary terminal model requires four propagation declarations; found {len(props)}"
            )
        terminal_by_pool={}
        for prop in props:
            terminal_by_pool.setdefault(prop['product'],prop['incoming'])
        mapping={}
        for prop in props:
            terminal=terminal_by_pool.get(prop['pool'])
            incoming=prop['incoming']
            if terminal in monomers and incoming in monomers:
                key=(terminal,incoming)
                if key in mapping: raise ChemicalModelIncompatibleError(f"duplicate terminal pair {key}")
                mapping[key]=prop['rate']
        expected={(a,b) for a in monomers for b in monomers}
        if set(mapping)!=expected:
            raise ChemicalModelIncompatibleError(f"cannot infer terminal propagation pairs; missing {sorted(expected-set(mapping))}")
        values={}
        for pair,name in mapping.items():
            try: values[pair]=_readonly(np.asarray(self._run.k[name],dtype=float))
            except KeyError as exc: raise DataUnavailableError(f"kinetic rate {name!r} unavailable") from exc
        return monomers,PairValues(values)

    @property
    def capabilities(self):
        data=('channel_events' in self._run.tables and 'state' in self._run.tables)
        terminal=True
        try: self._terminal_rate_series()
        except Exception: terminal=False
        try:
            pue=all(x=='explicit' for x in self.penultimate_parameters().classification)
        except Exception:
            pue=False
        return Capabilities({
            'composition':Capability('composition',True,len(self._run.monomer_names)>=2,data,None if data else 'missing_data'),
            'reactivity_ratios':Capability('reactivity_ratios',True,terminal,data,None if terminal and data else 'unsupported_variant'),
            'mayo_lewis':Capability('mayo_lewis',True,terminal,data,None if terminal and data else 'unsupported_variant'),
            'compare_mayo_lewis':Capability('compare_mayo_lewis',True,terminal,data,None if terminal and data else 'unsupported_variant'),
            'terminal_diagnostics':Capability('terminal_diagnostics',True,terminal,data and 'microstructure_motifs' in self._run.tables,None if terminal and data and 'microstructure_motifs' in self._run.tables else 'missing_data'),
            'penultimate_parameters':Capability('penultimate_parameters',True,len(self._run.monomer_names)==2,data,None if len(self._run.monomer_names)==2 and data else 'unsupported_variant'),
            'penultimate_composition':Capability('penultimate_composition',True,pue,data,None if pue and data else 'unsupported_variant'),
            'compare_penultimate':Capability('compare_penultimate',True,pue,data,None if pue and data else 'unsupported_variant'),
            'penultimate_diagnostics':Capability('penultimate_diagnostics',True,pue,data and 'microstructure_motifs' in self._run.tables,None if pue and data and 'microstructure_motifs' in self._run.tables else 'missing_data'),
        })

    def _channel_monomer_ledger(self):
        monomers=tuple(self._run.monomer_names)
        inserted={m:np.zeros(len(self._run.sid),dtype=float) for m in monomers}
        removed={m:np.zeros(len(self._run.sid),dtype=float) for m in monomers}
        # Map channel names from input declarations to incoming monomer.
        mapping={}
        for raw in self._run._input_model_lines():
            p=raw.split('#',1)[0].split()
            if len(p)<2 or p[0]!='macro' or '->' not in p: continue
            op=p[1]; arrow=p.index('->'); rate=p[-1]
            if op in {'init','prop','transfer_m','deprop'}:
                incoming=None
                if '+' in p[:arrow]: incoming=p[p.index('+')+1]
                elif op=='deprop' and arrow+2 < len(p): incoming=p[arrow+2]
                if incoming in monomers:
                    mapping[(op,rate,incoming)]=incoming
        channel_names=self._run.channels.event_count.names
        # Prefer direct rate/channel-name substring matching; writer channel names are model-derived.
        for channel in channel_names:
            lower=channel.lower()
            matched=None; is_remove=False
            for (op,rate,mon),_ in mapping.items():
                if rate.lower() in lower or channel==rate:
                    matched=mon; is_remove=(op=='deprop'); break
            if matched is None:
                # conventional names: prop_AA, init_A, deprop_A
                for mon in monomers:
                    if lower.endswith('_'+mon.lower()) or ('_'+mon.lower()+'_') in lower:
                        if any(k in lower for k in ('prop','init','transfer_m','deprop')):
                            matched=mon; is_remove='deprop' in lower; break
            if matched is not None:
                target=removed if is_remove else inserted
                target[matched]+=np.asarray(self._run.channels.event_count[channel],dtype=float)
        return inserted,removed

    def _build(self):
        if self._composition is not None: return self._composition
        monomers=tuple(self._run.monomer_names)
        counts={m:np.asarray(self._run.count[m],dtype=float) for m in monomers}
        free_matrix=np.column_stack([counts[m] for m in monomers]); free_total=free_matrix.sum(axis=1)
        free_frac=np.divide(free_matrix,free_total[:,None],out=np.full(free_matrix.shape,np.nan),where=free_total[:,None]>0)
        free_defined=free_total>0
        inserted,removed=self._channel_monomer_ledger()
        net={m:inserted[m]-removed[m] for m in monomers}
        net_matrix=np.column_stack([net[m] for m in monomers]); net_total=net_matrix.sum(axis=1)
        cum_frac=np.divide(net_matrix,net_total[:,None],out=np.full(net_matrix.shape,np.nan),where=net_total[:,None]>0)
        cum_defined=net_total>0
        delta={m:np.diff(inserted[m]-removed[m]) for m in monomers}
        delta_matrix=np.column_stack([delta[m] for m in monomers]); delta_total=delta_matrix.sum(axis=1)
        int_frac=np.divide(delta_matrix,delta_total[:,None],out=np.full(delta_matrix.shape,np.nan),where=delta_total[:,None]>0)
        int_defined=delta_total>0
        ids=_readonly(np.asarray(self._run.sid,dtype=np.int64)); times=_readonly(np.asarray(self._run.t,dtype=float)); conv=_readonly(np.asarray(self._run.conv.total,dtype=float))
        free=SnapshotCompositionSeries(ids,times,conv,NamedValues({m:_readonly(free_frac[:,i]) for i,m in enumerate(monomers)}),_readonly(free_defined),NamedValues({m:_readonly(counts[m]) for m in monomers}),'free_monomer_mole_fraction')
        cumulative=SnapshotCompositionSeries(ids,times,conv,NamedValues({m:_readonly(cum_frac[:,i]) for i,m in enumerate(monomers)}),_readonly(cum_defined),NamedValues({m:_readonly(net[m]) for m in monomers}),'cumulative_repeat_unit_fraction')
        incremental=IntervalCompositionSeries(_readonly(ids[:-1]),_readonly(ids[1:]),_readonly((times[:-1]+times[1:])/2),_readonly(np.diff(times)),_readonly((conv[:-1]+conv[1:])/2),NamedValues({m:_readonly(int_frac[:,i]) for i,m in enumerate(monomers)}),_readonly(int_defined),NamedValues({m:_readonly(delta[m]) for m in monomers}))
        self._composition=CompositionResult(free,incremental,cumulative,NamedValues({m:_readonly(inserted[m]) for m in monomers}),NamedValues({m:_readonly(removed[m]) for m in monomers}),NamedValues({m:_readonly(net[m]) for m in monomers}))
        return self._composition

    def composition(self): return self._build()
    def monomer_composition(self): return self._build().free
    def incremental_composition(self): return self._build().incremental
    def cumulative_composition(self): return self._build().cumulative
    def polymer_composition(self): return self._build().cumulative

    def reactivity_ratios(self):
        if self._rr is not None: return self._rr
        monomers,rates=self._terminal_rate_series(); a,b=monomers
        ratios={}; defined={}
        for name,num,den in ((a,(a,a),(a,b)),(b,(b,b),(b,a))):
            n=np.asarray(rates[num]); d=np.asarray(rates[den]); val=np.divide(n,d,out=np.full(n.shape,np.nan),where=d!=0)
            val[(d==0)&(n>0)]=np.inf; ratios[name]=_readonly(val); defined[name]=_readonly(~((d==0)&(n==0)))
        self._rr=ReactivityRatioSeries(monomers,_readonly(self._run.sid),_readonly(self._run.t),rates,NamedValues(ratios),NamedValues(defined)); return self._rr

    def mayo_lewis(self):
        if self._ml is not None: return self._ml
        rr=self.reactivity_ratios(); comp=self._build(); a,b=rr.monomers
        f_a=np.asarray(comp.free.fractions[a]); f_b=np.asarray(comp.free.fractions[b])
        r_a=np.asarray(rr.reactivity_ratios[a]); r_b=np.asarray(rr.reactivity_ratios[b])
        den=r_a*f_a*f_a+2*f_a*f_b+r_b*f_b*f_b
        defined=np.isfinite(den)&(den>0)&comp.free.is_defined
        F_a=np.divide(r_a*f_a*f_a+f_a*f_b,den,out=np.full(den.shape,np.nan),where=defined); F_b=1-F_a; F_b[~defined]=np.nan
        _,_,ops,_=self._model()
        assumptions=MappingProxyType({'binary':True,'terminal':True,'has_depropagation':'deprop' in ops,'has_transfer':any(x.startswith('transfer') for x in ops),'has_termination':any(x.startswith('term') for x in ops)})
        self._ml=MayoLewisSeries(rr.monomers,rr.snapshot_ids,rr.times,comp.free.conversion,rr.rate_constants,rr.reactivity_ratios,comp.free.fractions,NamedValues({a:_readonly(F_a),b:_readonly(F_b)}),_readonly(defined),assumptions); return self._ml

    def compare_mayo_lewis(self, monomer_reference='start', parameter_reference='start'):
        if monomer_reference not in {'start','end','midpoint'} or parameter_reference not in {'start','end','midpoint'}: raise ValueError('references must be start, end, or midpoint')
        ml=self.mayo_lewis(); comp=self._build(); inc=comp.incremental
        def ref(v,mode):
            x=np.asarray(v); return x[:-1] if mode=='start' else x[1:] if mode=='end' else (x[:-1]+x[1:])/2
        f={m:ref(comp.free.fractions[m],monomer_reference) for m in ml.monomers}; r={m:ref(ml.reactivity_ratios[m],parameter_reference) for m in ml.monomers}
        a,b=ml.monomers; den=r[a]*f[a]*f[a]+2*f[a]*f[b]+r[b]*f[b]*f[b]; defined=(den>0)&np.isfinite(den)&inc.is_defined
        pred_a=np.divide(r[a]*f[a]*f[a]+f[a]*f[b],den,out=np.full(den.shape,np.nan),where=defined); pred_b=1-pred_a; pred_b[~defined]=np.nan
        pred=NamedValues({a:_readonly(pred_a),b:_readonly(pred_b)}); diff=NamedValues({m:_readonly(np.asarray(inc.fractions[m])-pred[m]) for m in ml.monomers})
        return MayoLewisComparison(ml.monomers,inc.start_snapshot_ids,inc.end_snapshot_ids,_readonly(self._run.t[:-1]),_readonly(self._run.t[1:]),inc.t_mid,inc.dt,inc.conversion,monomer_reference,parameter_reference,NamedValues({m:_readonly(f[m]) for m in ml.monomers}),pred,inc.fractions,diff,_readonly(defined))

    def composition_drift(self, monomer_reference='start'):
        if monomer_reference not in {'start','end','midpoint'}: raise ValueError('monomer_reference must be start, end, or midpoint')
        result=self._build(); inc=result.incremental; refs={}; diffs={}; defined=inc.is_defined.copy()
        for m in self._run.monomer_names:
            x=np.asarray(result.free.fractions[m]); v=x[:-1] if monomer_reference=='start' else x[1:] if monomer_reference=='end' else (x[:-1]+x[1:])/2
            refs[m]=_readonly(v); diffs[m]=_readonly(np.asarray(inc.fractions[m])-v); defined &= np.isfinite(v)
        return CompositionDrift(inc.start_snapshot_ids,inc.end_snapshot_ids,inc.t_mid,inc.dt,inc.conversion,inc.fractions,NamedValues(refs),NamedValues(diffs),_readonly(defined),monomer_reference)

    def _propagation_records(self):
        records=[]
        for raw in self._run._input_model_lines():
            p=raw.split('#',1)[0].split()
            if len(p)>=8 and p[0]=='macro' and p[1]=='prop' and '->' in p:
                arrow=p.index('->')
                records.append({'pool':p[2], 'incoming':p[4], 'product':p[arrow+1],
                                'rate':p[arrow+2], 'channel':f'prop_{p[2]}_{p[4]}'})
        return records

    def _pool_states(self):
        monomers=tuple(self._run.monomer_names); states={}
        # Initiation products identify terminal-only states.
        for raw in self._run._input_model_lines():
            p=raw.split('#',1)[0].split()
            if len(p)>=8 and p[0]=='macro' and p[1]=='init' and '->' in p:
                arrow=p.index('->'); incoming=p[4]; product=p[arrow+1]
                if incoming in monomers: states[product]=(None,incoming)
        props=self._propagation_records(); changed=True
        while changed:
            changed=False
            for rec in props:
                if rec['pool'] not in states: continue
                _,terminal=states[rec['pool']]
                value=(terminal,rec['incoming'])
                if rec['product'] not in states:
                    states[rec['product']]=value; changed=True
                elif states[rec['product']] != value:
                    if rec['product'] == rec['pool'] and states[rec['product']][0] is None:
                        continue
                    raise ChemicalModelIncompatibleError(f"pool {rec['product']!r} has inconsistent sequence state")
        return states

    def terminal_diagnostics(self, monomer_reference='start', parameter_reference='start'):
        if monomer_reference not in {'start','end','midpoint'} or parameter_reference not in {'start','end','midpoint'}:
            raise ValueError('references must be start, end, or midpoint')
        ml=self.mayo_lewis(); comp=self._build(); a,b=ml.monomers; times=np.asarray(ml.times)
        def ref(v,mode):
            x=np.asarray(v); return x[:-1] if mode=='start' else x[1:] if mode=='end' else (x[:-1]+x[1:])/2
        f={m:ref(comp.free.fractions[m],monomer_reference) for m in (a,b)}
        k={pair:ref(ml.rate_constants[pair],parameter_reference) for pair in ml.rate_constants.keys()}
        records=self._propagation_records(); terminal_by_pool={}
        for rec in records:
            if rec['incoming'] in (a,b): terminal_by_pool.setdefault(rec['product'],rec['incoming'])
        channels={}
        for rec in records:
            terminal=terminal_by_pool.get(rec['pool'])
            if terminal in (a,b) and rec['incoming'] in (a,b):
                channels[(terminal,rec['incoming'])]=rec['channel']
        expected={(x,y) for x in (a,b) for y in (a,b)}
        if set(channels)!=expected:
            raise ChemicalModelIncompatibleError('cannot map all four terminal propagation channels')
        predicted={}; observed={}; differences={}; outgoing={}; defined_rows={}
        for terminal in (a,b):
            den=sum(k[(terminal,j)]*f[j] for j in (a,b))
            fires={j:np.asarray(self._run.firings.delta_fires_series(channels[(terminal,j)]),dtype=float) for j in (a,b)}
            total=fires[a]+fires[b]; outgoing[terminal]=_readonly(total)
            row_defined=(total>0)&np.isfinite(den)&(den>0); defined_rows[terminal]=_readonly(row_defined)
            for j in (a,b):
                pred=np.divide(k[(terminal,j)]*f[j],den,out=np.full(den.shape,np.nan),where=den>0)
                obs=np.divide(fires[j],total,out=np.full(total.shape,np.nan),where=total>0)
                predicted[(terminal,j)]=_readonly(pred); observed[(terminal,j)]=_readonly(obs)
                differences[(terminal,j)]=_readonly(np.where(row_defined,obs-pred,np.nan))
        transitions=TerminalTransitionDiagnostics((a,b),ml.snapshot_ids[:-1],ml.snapshot_ids[1:],_readonly(times[:-1]),_readonly(times[1:]),_readonly((times[:-1]+times[1:])/2),_readonly(np.diff(times)),PairValues(predicted),PairValues(observed),PairValues(differences),NamedValues(outgoing),NamedValues(defined_rows))
        pred_ln={}; pred_lw={}; pred_d={}; obs_ln={}; obs_lw={}; obs_d={}; counts={}; block_defined={}
        block_rows=self._run.microstructure.run_lengths().rows()
        for terminal in (a,b):
            other=b if terminal==a else a
            den=np.asarray(ml.rate_constants[(terminal,terminal)])*np.asarray(ml.monomer_mole_fractions[terminal])+np.asarray(ml.rate_constants[(terminal,other)])*np.asarray(ml.monomer_mole_fractions[other])
            p=np.divide(np.asarray(ml.rate_constants[(terminal,terminal)])*np.asarray(ml.monomer_mole_fractions[terminal]),den,out=np.full(den.shape,np.nan),where=den>0)
            q=1-p; ok=np.isfinite(q)&(q>0)
            pred_ln[terminal]=_readonly(np.divide(1.0,q,out=np.full(q.shape,np.nan),where=ok)); pred_lw[terminal]=_readonly(np.divide(1+p,q,out=np.full(q.shape,np.nan),where=ok)); pred_d[terminal]=_readonly(np.where(ok,1+p,np.nan))
            ln=[]; lw=[]; disp=[]; nblocks=[]; defs=[]
            for sid in ml.snapshot_ids:
                rows=self._run.microstructure.run_lengths(terminal,snapshot=int(sid)).rows()
                n0=sum(int(r['count']) for r in rows); n1=sum(int(r['run_length'])*int(r['count']) for r in rows); n2=sum(int(r['run_length'])**2*int(r['count']) for r in rows)
                nblocks.append(n0); defs.append(n0>0 and n1>0); ln.append(n1/n0 if n0 else np.nan); lw.append(n2/n1 if n1 else np.nan); disp.append((n2/n1)/(n1/n0) if n0 and n1 else np.nan)
            obs_ln[terminal]=_readonly(ln); obs_lw[terminal]=_readonly(lw); obs_d[terminal]=_readonly(disp); counts[terminal]=_readonly(nblocks); block_defined[terminal]=_readonly(defs)
        blocks=TerminalBlockDiagnostics((a,b),ml.snapshot_ids,ml.times,NamedValues(pred_ln),NamedValues(pred_lw),NamedValues(pred_d),NamedValues(obs_ln),NamedValues(obs_lw),NamedValues(obs_d),NamedValues(counts),NamedValues(block_defined))
        return TerminalDiagnostics(transitions,blocks)

    def penultimate_parameters(self):
        monomers=tuple(self._run.monomer_names)
        if len(monomers)!=2: raise ChemicalAnalysisNotApplicableError('penultimate analysis requires exactly two monomers')
        propagation_records=self._propagation_records()
        if len(propagation_records) < 8:
            raise FeatureUnavailableError('explicit binary penultimate analysis requires eight propagation channels')
        index={name:i for i,name in enumerate(monomers)}; states=self._pool_states(); records=[]
        for rec in propagation_records:
            state=states.get(rec['pool'])
            if state and state[0] in index and state[1] in index and rec['incoming'] in index:
                records.append((state[0],state[1],rec['incoming'],rec['rate'],rec['channel']))
        classification='explicit' if len(records)==8 else 'undefined'
        n=len(self._run.sid); tensor=np.full((n,2,2,2),np.nan)
        for previous,terminal,incoming,rate,_ in records:
            tensor[:,index[previous],index[terminal],index[incoming]]=np.asarray(self._run.k[rate],dtype=float)
        defined=np.isfinite(tensor).all(axis=(1,2,3))
        def ratio(num,den):
            x=tensor[:,num[0],num[1],num[2]]; y=tensor[:,den[0],den[1],den[2]]; out=np.full(n,np.nan); np.divide(x,y,out=out,where=y>0); out[(y==0)&(x>0)]=np.inf; return _readonly(out)
        a,b=monomers; ia,ib=index[a],index[b]
        r=NamedValues({a:ratio((ia,ia,ia),(ia,ia,ib)),b:ratio((ib,ib,ib),(ib,ib,ia))})
        rp=NamedValues({a:ratio((ib,ia,ia),(ib,ia,ib)),b:ratio((ia,ib,ib),(ia,ib,ia))})
        sv=NamedValues({a:ratio((ib,ia,ia),(ia,ia,ia)),b:ratio((ia,ib,ib),(ib,ib,ib))})
        tensor.flags.writeable=False
        variable=bool(n>1 and defined.all() and np.any(np.diff(tensor,axis=0)!=0))
        return PenultimateParameterSeries(monomers,_readonly(self._run.sid),_readonly(self._run.t),tensor,tuple(classification for _ in range(n)),variable,r,rp,sv,_readonly(defined))

    def penultimate_composition(self):
        params=self.penultimate_parameters()
        if any(x!='explicit' for x in params.classification): raise ChemicalModelIncompatibleError('eight explicit binary penultimate propagation channels are required')
        comp=self._build(); monomers=params.monomers; n=len(params); fractions=np.full((n,2),np.nan); states=np.full((n,4),np.nan); transitions=np.full((n,4,2),np.nan); defined=np.zeros(n,dtype=bool)
        for i in range(n):
            f=np.asarray([comp.free.fractions[m][i] for m in monomers],dtype=float)
            if np.isfinite(f).all(): fractions[i],states[i],transitions[i],defined[i]=_penultimate_prediction(params.tensor[i],f)
        pairs=tuple((p,t) for p in monomers for t in monomers); triples=tuple((p,t,j) for p,t in pairs for j in monomers)
        return PenultimateCompositionSeries(monomers,params.snapshot_ids,params.times,comp.free.conversion,comp.free.fractions,NamedValues({m:_readonly(fractions[:,i]) for i,m in enumerate(monomers)}),PairValues({pair:_readonly(states[:,i]) for i,pair in enumerate(pairs)}),TripleValues({key:_readonly(transitions[:,i//2,i%2]) for i,key in enumerate(triples)}),_readonly(defined))

    def compare_penultimate(self, monomer_reference='start', parameter_reference='start'):
        if monomer_reference not in {'start','end','midpoint'} or parameter_reference not in {'start','end','midpoint'}: raise ValueError('references must be start, end, or midpoint')
        theory=self.penultimate_composition(); params=self.penultimate_parameters(); comp=self._build(); inc=comp.incremental
        def ref(v,mode):
            x=np.asarray(v); return x[:-1] if mode=='start' else x[1:] if mode=='end' else (x[:-1]+x[1:])/2
        f=np.column_stack([ref(comp.free.fractions[m],monomer_reference) for m in theory.monomers]); k=ref(params.tensor,parameter_reference); pred=np.full((len(inc),2),np.nan); defined=np.zeros(len(inc),dtype=bool)
        for i in range(len(inc)): pred[i],_,_,defined[i]=_penultimate_prediction(k[i],f[i])
        defined &= inc.is_defined & (inc.dt>0); pred[~defined]=np.nan
        vals=NamedValues({m:_readonly(pred[:,i]) for i,m in enumerate(theory.monomers)}); diff=NamedValues({m:_readonly(np.asarray(inc.fractions[m])-vals[m]) for m in theory.monomers}); times=np.asarray(comp.free.times)
        return PenultimateComparison(theory.monomers,inc.start_snapshot_ids,inc.end_snapshot_ids,_readonly(times[:-1]),_readonly(times[1:]),inc.t_mid,inc.dt,inc.conversion,monomer_reference,parameter_reference,vals,inc.fractions,diff,_readonly(defined))

    def penultimate_diagnostics(self, monomer_reference='start', parameter_reference='start'):
        theory=self.penultimate_composition(); params=self.penultimate_parameters(); comparison=self.compare_penultimate(monomer_reference,parameter_reference); monomers=theory.monomers; pairs=tuple((p,t) for p in monomers for t in monomers); triples=tuple((p,t,j) for p,t in pairs for j in monomers); idx={m:i for i,m in enumerate(monomers)}
        def ref(v,mode):
            x=np.asarray(v); return x[:-1] if mode=='start' else x[1:] if mode=='end' else (x[:-1]+x[1:])/2
        f={m:ref(self._build().free.fractions[m],monomer_reference) for m in monomers}; kt=ref(params.tensor,parameter_reference); states=self._pool_states(); recmap={}
        for rec in self._propagation_records():
            state=states.get(rec['pool']); key=(state[0],state[1],rec['incoming']) if state else None
            if key in triples: recmap[key]=rec['channel']
        if set(recmap)!=set(triples): raise ChemicalModelIncompatibleError('cannot map all eight penultimate firing channels')
        predicted={}; observed={}; differences={}; outgoing={}; tdefined={}
        for p,t in pairs:
            ip,it=idx[p],idx[t]; den=sum(kt[:,ip,it,idx[j]]*f[j] for j in monomers); fires={j:np.asarray(self._run.firings.delta_fires_series(recmap[(p,t,j)]),dtype=float) for j in monomers}; total=fires[monomers[0]]+fires[monomers[1]]; outgoing[(p,t)]=_readonly(total); ok=(total>0)&(den>0)&np.isfinite(den); tdefined[(p,t)]=_readonly(ok)
            for j in monomers:
                pr=np.divide(kt[:,ip,it,idx[j]]*f[j],den,out=np.full(den.shape,np.nan),where=den>0); ob=np.divide(fires[j],total,out=np.full(total.shape,np.nan),where=total>0); predicted[(p,t,j)]=_readonly(pr); observed[(p,t,j)]=_readonly(ob); differences[(p,t,j)]=_readonly(np.where(ok,ob-pr,np.nan))
        predicted_states={pair:theory.radical_state_fractions[pair] for pair in pairs}; observed_states={pair:_readonly(np.full(len(theory),np.nan)) for pair in pairs}; state_diff={pair:_readonly(np.full(len(theory),np.nan)) for pair in pairs}
        pred_tri={}; obs_tri={}; tri_diff={}
        for key in triples:
            pr=np.asarray(theory.radical_state_fractions[(key[0],key[1])])*np.asarray(theory.transition_probabilities[key]); obs=[]
            for sid in theory.snapshot_ids:
                rows=self._run.microstructure.triads(snapshot=int(sid)).rows(); d={r['motif']:float(r['count']) for r in rows}; total=sum(d.values()); obs.append(d.get('|'.join(key),0.0)/total if total else np.nan)
            pred_tri[key]=_readonly(pr); obs_tri[key]=_readonly(obs); tri_diff[key]=_readonly(np.asarray(obs)-pr)
        return PenultimateDiagnostics(monomers,comparison.start_snapshot_ids,comparison.end_snapshot_ids,TripleValues(predicted),TripleValues(observed),TripleValues(differences),PairValues(outgoing),PairValues(tdefined),theory.snapshot_ids,PairValues(predicted_states),PairValues(observed_states),PairValues(state_diff),TripleValues(pred_tri),TripleValues(obs_tri),TripleValues(tri_diff),False)
