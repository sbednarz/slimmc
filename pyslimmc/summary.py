from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import json
from typing import Any
import numpy as np


def _last(series: Any) -> float | int | None:
    try:
        a = np.asarray(series)
        if a.size == 0:
            return None
        v = a[-1]
        return v.item() if hasattr(v, 'item') else v
    except Exception:
        return None


def _fmt(v: Any) -> str:
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)

@dataclass(frozen=True)
class RunSummary:
    path: str
    engine: str
    engine_version: str | None
    cli_version: str | None
    status: str
    validation: str
    validation_warnings: int
    validation_errors: int
    seed: int | None
    final_time: float | None
    final_event: int | None
    snapshots: int
    chain_snapshots: int
    final_conversion_total: float | None
    final_conversion: dict[str, float]
    dpn: float | None
    dpw: float | None
    mn: float | None
    mw: float | None
    mz: float | None
    dispersity: float | None
    live_chains: int | None
    dead_chains: int | None
    temperature_initial: float | None
    temperature_final: float | None
    temperature_min: float | None
    temperature_max: float | None
    peak_memory_B: float | None
    results_size_B: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, ensure_ascii=False, sort_keys=False)

    def to_text(self) -> str:
        rows = [
            "Slimmc run summary", "",
            f"Run:                 {self.path}",
            f"Engine:              {self.engine}",
            f"Engine version:      {_fmt(self.engine_version)}",
            f"CLI version:         {_fmt(self.cli_version)}",
            f"Status:              {self.status}",
            f"Validation:          {self.validation}",
            f"Warnings:            {self.validation_warnings}",
            f"Errors:              {self.validation_errors}",
            f"Seed:                {_fmt(self.seed)}", "",
            f"Final time:          {_fmt(self.final_time)}",
            f"Final KMC event:     {_fmt(self.final_event)}",
            f"Snapshots:           {self.snapshots}",
            f"Chain snapshots:     {self.chain_snapshots}", "",
            "Final conversion:",
            f"  total              {_fmt(self.final_conversion_total)}",
        ]
        for name, value in self.final_conversion.items():
            rows.append(f"  {name:<18} {_fmt(value)}")
        rows += ["", "Final moments:",
                 f"  DPN                {_fmt(self.dpn)}",
                 f"  DPW                {_fmt(self.dpw)}",
                 f"  Mn                 {_fmt(self.mn)}",
                 f"  Mw                 {_fmt(self.mw)}",
                 f"  Mz                 {_fmt(self.mz)}",
                 f"  Dispersity         {_fmt(self.dispersity)}", "",
                 "Chains:",
                 f"  live               {_fmt(self.live_chains)}",
                 f"  dead               {_fmt(self.dead_chains)}", "",
                 "Temperature:",
                 f"  initial            {_fmt(self.temperature_initial)}",
                 f"  final              {_fmt(self.temperature_final)}",
                 f"  min                {_fmt(self.temperature_min)}",
                 f"  max                {_fmt(self.temperature_max)}", "",
                 f"Peak memory:         {_fmt(self.peak_memory_B)} B",
                 f"Results size:        {self.results_size_B} B"]
        return "\n".join(rows)

    def write(self, path: str | Path) -> Path:
        p = Path(path)
        if p.suffix.lower() == ".json":
            p.write_text(self.to_json() + "\n", encoding="utf-8")
        else:
            p.write_text(self.to_text() + "\n", encoding="utf-8")
        return p

    def __str__(self) -> str:
        return self.to_text()


def build_summary(run: Any) -> RunSummary:
    md = run.metadata.raw() if hasattr(run, "metadata") else {}
    validation = getattr(getattr(run, "diagnostics", None), "validation", None)
    warnings = int(getattr(validation, "warning_count", 0) or 0)
    errors = int(getattr(validation, "error_count", 0) or 0)
    if getattr(run, "status", None) == "completed":
        validation_status = "PASS" if errors == 0 else "FAIL"
    else:
        validation_status = "PARTIAL" if validation is not None else "-"

    conversion: dict[str, float] = {}
    for name in getattr(run, "monomer_names", ()):
        try:
            v = _last(run.conv[name])
            if v is not None: conversion[name] = float(v)
        except Exception:
            pass

    try:
        chain_sids = np.asarray(run.table("chains")["snapshot_id"], dtype=np.int64)
        chain_snapshots = int(np.unique(chain_sids).size)
    except Exception:
        chain_snapshots = 0

    def moment(name: str) -> float | None:
        try:
            v = _last(getattr(run, name))
            return None if v is None else float(v)
        except Exception:
            return None

    live = dead = None
    try:
        ch = run.last.chains
        live = int(np.sum(np.asarray(ch.live.count, dtype=np.int64)))
        dead = int(np.sum(np.asarray(ch.dead.count, dtype=np.int64)))
    except Exception:
        pass

    ti = tf = tmin = tmax = None
    try:
        arr = np.asarray(run.temp, dtype=float)
        if arr.size:
            ti, tf, tmin, tmax = map(float, (arr[0], arr[-1], np.nanmin(arr), np.nanmax(arr)))
    except Exception:
        pass

    peak_mem = None
    try:
        arr = np.asarray(run.diagnostics.memory.total_est_B, dtype=float)
        if arr.size: peak_mem = float(np.nanmax(arr))
    except Exception:
        pass

    size = 0
    try:
        size = sum(p.stat().st_size for p in Path(run.path).rglob("*") if p.is_file())
    except Exception:
        pass

    seed = md.get("seed")
    if seed is None:
        seed = md.get("rng_seed")
    try: seed = int(seed) if seed is not None else None
    except Exception: seed = None

    conv_total = None
    try:
        v = _last(run.conv.total)
        conv_total = None if v is None else float(v)
    except Exception:
        pass

    return RunSummary(
        path=getattr(run, "_display_path", lambda: Path(run.path).name + "/")(),
        engine=str(getattr(run, "engine", md.get("engine", "-"))),
        engine_version=md.get("engine_version") or md.get("version"),
        cli_version=md.get("cli_version"),
        status=str(getattr(run, "status", md.get("run_status", "unknown"))),
        validation=validation_status,
        validation_warnings=warnings,
        validation_errors=errors,
        seed=seed,
        final_time=(lambda x: None if x is None else float(x))(_last(getattr(run, "t", []))),
        final_event=(lambda x: None if x is None else int(x))(_last(getattr(run, "event", []))),
        snapshots=len(getattr(run, "snapshots", [])),
        chain_snapshots=chain_snapshots,
        final_conversion_total=conv_total,
        final_conversion=conversion,
        dpn=moment("dpn"), dpw=moment("dpw"), mn=moment("mn"), mw=moment("mw"), mz=moment("mz"), dispersity=moment("dispersity"),
        live_chains=live, dead_chains=dead,
        temperature_initial=ti, temperature_final=tf, temperature_min=tmin, temperature_max=tmax,
        peak_memory_B=peak_mem, results_size_B=size,
    )
