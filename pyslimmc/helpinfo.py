"""Uniform ``.help()`` / ``.info()`` contract for public pyslimmc objects."""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from .operations import analysis_operation


def _public_names(obj: Any) -> list[str]:
    names=[]
    for name in dir(obj):
        if name.startswith('_') or name in {'help','info'}:
            continue
        try:
            value=getattr(obj,name)
        except Exception:
            continue
        if callable(value):
            continue
        names.append(name)
    return names


def _generic_info(self: Any) -> str:
    cls=type(self).__name__.removeprefix('Storage')
    lines=[cls]
    try:
        if isinstance(self, Mapping):
            keys=list(self.keys())
            lines.append(f"Items: {len(keys)}")
            if keys: lines.append("Keys: " + ", ".join(map(str, keys[:20])) + (" ..." if len(keys)>20 else ""))
        elif isinstance(self, Sequence) and not isinstance(self,(str,bytes,bytearray)):
            lines.append(f"Items: {len(self)}")
    except Exception:
        pass
    names=_public_names(self)
    if names:
        lines.append("Data: " + ", ".join(names[:24]) + (" ..." if len(names)>24 else ""))
    text='\n'.join(lines)
    print(text)
    return text


def _branch_help(title: str, body: str):
    text=f"{title}\n{'-'*len(title)}\n{body.strip()}"
    def help_method(self: Any) -> str:
        print(text)
        return text
    return help_method


def _wrap_operation(cls: type, name: str, text: str) -> None:
    current=cls.__dict__.get(name)
    if current is None or isinstance(current, analysis_operation) or isinstance(current, property):
        return
    if not callable(current):
        return
    setattr(cls, name, analysis_operation(text)(current))


def install() -> None:
    # Imports are intentionally local: installation runs after package modules exist.
    from . import _storage as s
    from . import storage_analysis as a
    from . import copolymerization as c
    from . import runs as r
    from .report import Report
    from .summary import RunSummary
    from .distributions import ChainLengthDistribution, MassDistribution, MolarMassDistribution
    from .counts import DPCounts, MassCounts
    from .run import MassAuditResult

    branch_specs={
      s.StorageRun:("Run", """Obtain a run with sl.open(path), runs.<run_id>, runs.one(run_id=...), runs.match(...)[index], or runs.pack(...)[key]["run"].

Inspect snapshots and axes with run.first, run.last, run.final, run.t, run.event and run.sid. Use run.count, run.moles, run.conc, run.conv, run.mn, run.mw, run.mz and run.dispersity for state and moments.

Variables are available through run.var, run.var["name"].value, run.var.info() and run.var.help(). For distributions call run.mwd.help() or run.cld.help() before the analysis."""),
      s.StorageSnapshots:("Snapshots", "Collection of saved snapshots. Use .first, .last, .final, .at_time(...) or .at_event(...)."),
      s.StorageSnapshot:("Snapshot", "One saved simulation state. Inspect state, chains, moments, kinetics and distributions for this snapshot."),
      s.StorageStateSeries:("State", "Time-dependent state data. Inspect available entities, counts, moles and concentrations."),
      s.SeriesView:("Series", "Named numerical series aligned to saved snapshots. Use keys/indexing and .info() to inspect coverage."),
      s.ConversionSeries:("Conversion", "Monomer and total conversion series aligned to saved snapshots."),
      s.PolymerCompositionSeries:("Polymer composition F", "Polymer composition: instantaneous, interval and cumulative series. Use .info() on each branch."),
      s.StorageChains:("Chains", "Chain population data. Filter by activity/pool/origin, inspect columns, or calculate exact DP/mass counts, MWD and CLD."),
      s.StorageMomentsSeries:("Moments", "Time-dependent number-, weight- and z-average DP/molar-mass moments by population and mass model."),
      s.StorageMomentsSnapshot:("Moments", "Snapshot moments by population and mass model."),
      s.StorageChannelsSeries:("Channels", "Reaction-channel event data. Inspect event counts and shares by channel."),
      s.StorageChannelsSnapshot:("Channels", "Reaction-channel data at one snapshot."),
      a.StorageFirings:("Firings", "Cumulative and interval channel firings, fire shares and rate shares."),
      s.StorageKineticsSeries:("Kinetics", "Kinetic definitions, temperatures and rate constants through time."),
      s.StorageKineticsSnapshot:("Kinetics", "Kinetic values at one snapshot."),
      s.StorageActions:("Actions", "Recorded model actions and their triggers, times, conditions and before/after values."),
      s.StorageAction:("Action", "One executed action record."),
      a.StorageCopolymerization:("Copolymerization", "Copolymer composition, reactivity ratios, Mayo-Lewis, terminal and penultimate analyses. Call an analysis method's .help() before use."),
      a.StorageMicrostructure:("Microstructure", "Dyads, triads, run lengths, transition fractions, blockiness and sequence consistency."),
      s.StorageValidation:("Validation", "Validation checks for this run. Inspect PASS/WARN/FAIL results or call run.validate.help()."),
      s.StorageDiagnostics:("Diagnostics", "Memory, logs and channel trace diagnostics."),
      s.StorageRaw:("Raw", "Low-level schema, tables and dictionaries. Prefer higher-level API unless inspecting storage details."),
      r.Runs:("Runs", "Collection of runs. Filter/select, build a sweep, compare models or export a table."),
      Report:("Report", "Rendered report assembled from one or more data blocks."),
      RunSummary:("Summary", "Compact machine-readable and human-readable summary of one run."),
    }
    for cls,(title,body) in branch_specs.items():
        if 'help' not in cls.__dict__:
            setattr(cls,'help',_branch_help(title,body))
        if 'info' not in cls.__dict__:
            setattr(cls,'info',_generic_info)

    # Every result/data package must expose info().
    result_classes=[ChainLengthDistribution,MassDistribution,MolarMassDistribution,DPCounts,MassCounts,MassAuditResult]
    result_classes += [v for v in vars(c).values() if isinstance(v,type) and v.__module__==c.__name__]
    for cls in result_classes:
        if 'info' not in cls.__dict__:
            setattr(cls,'info',_generic_info)

    operation_groups={
      s.StorageRun:['mwd','cld','mass_distribution','dp_counts','mass_counts','validate','mass_audit','summary'],
      s.StorageSnapshot:['mwd','cld','mass_distribution','dp_counts','mass_counts','validate'],
      s.StorageChains:['mwd','cld','mass_distribution','dp_counts','mass_counts','masses','where'],
      a.StorageMicrostructure:['dyads','triads','run_lengths','transition_fraction','homodyad_fraction','blockiness','check_sequence_consistency'],
      a.StorageFirings:['channels','rows','final_row','final_fires','total_fires','channel_fires','delta_fires','fire_shares','rate_shares','validate'],
      a.StorageCopolymerization:['composition','monomer_composition','incremental_composition','cumulative_composition','polymer_composition','reactivity_ratios','mayo_lewis','compare_mayo_lewis','composition_drift','terminal_diagnostics','penultimate_parameters','penultimate_composition','compare_penultimate','penultimate_diagnostics'],
      s.StorageChannelsSeries:['interval_event_counts','fire_shares'],
      s.StorageKineticsSeries:['by_kind'],
      s.StorageRaw:['table','dictionary'],
      r.Runs:['filter','one','first','sweep','as_table','model_diff'],
    }
    for cls,names in operation_groups.items():
        for name in names:
            _wrap_operation(cls,name,f"{type.__getattribute__(cls,'__name__')}.{name}(...)\n\n{getattr(cls,name).__doc__ or 'See the method signature and documentation for parameters and returned data.'}")
