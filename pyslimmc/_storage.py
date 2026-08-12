"""Minimal, engine-neutral reader for Slimmc Storage 1.2.0."""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterator, Mapping, Any
import math

import numpy as np

from .run import Run, MassAuditResult, MassEntry, Variables
from .chains import ChainPopulation
from .core import (
    InvalidOutputError, DataUnavailableError, AnalysisNotApplicableError, OutputStatus,
    ValidationReport, ValidationFailedError,
    FinalSnapshotUnavailableError, SnapshotUnavailableError,
)
from .operations import analysis_operation, MWD_HELP, CLD_HELP, SPECTRUM_HELP


AVOGADRO = 6.02214076e23


class IncompleteResultsError(InvalidOutputError):
    """Raised when a non-completed Slimmc Storage run is opened without opt-in."""


@dataclass(frozen=True)
class StorageTable(Mapping[str, np.ndarray]):
    path: Path
    _columns: dict[str, np.ndarray]

    @classmethod
    def open(cls, path: Path, *, mmap_mode: str | None = "r") -> "StorageTable":
        cols: dict[str, np.ndarray] = {}
        if path.is_dir():
            for file in sorted(path.glob("*.npy")):
                cols[file.stem] = np.load(file, allow_pickle=False, mmap_mode=mmap_mode)
        lengths = {len(v) for v in cols.values()}
        if len(lengths) > 1:
            raise InvalidOutputError(f"Unequal column lengths in table {path}: {sorted(lengths)}")
        return cls(path, cols)

    def __getitem__(self, key: str) -> np.ndarray:
        return self._columns[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._columns)

    def __len__(self) -> int:
        return len(self._columns)

    @property
    def n_rows(self) -> int:
        return len(next(iter(self._columns.values()))) if self._columns else 0

    @property
    def columns(self) -> tuple[str, ...]:
        return tuple(self._columns)

    def filtered(self, mask: np.ndarray) -> "StorageTable":
        mask = np.asarray(mask, dtype=bool)
        if self._columns and len(mask) != self.n_rows:
            raise ValueError("filter mask length does not match table row count")
        return StorageTable(self.path, {name: values[mask] for name, values in self._columns.items()})

    def _ipython_key_completions_(self):
        return list(self._columns)

    def __dir__(self):
        return sorted(set(super().__dir__()) | set(self._columns))

    def __getattr__(self, name: str) -> np.ndarray:
        try:
            return self._columns[name]
        except KeyError as exc:
            raise AttributeError(name) from exc


def _readonly_array(values: Any, *, dtype=None) -> np.ndarray:
    """Return a NumPy array suitable for public numeric L2 endpoints.

    The returned array is always read-only. Memory-mapped source columns remain
    zero-copy where possible; computed arrays are protected against accidental
    mutation before exposure through the public API.
    """
    array = np.asarray(values, dtype=dtype)
    if array.flags.writeable:
        array.setflags(write=False)
    return array


@dataclass(frozen=True)
class SeriesView:
    """Lightweight NumPy-compatible time series without a pandas dependency."""

    values: np.ndarray
    time: np.ndarray
    snapshot_id: np.ndarray
    name: str
    unit: str | None = None

    def __array__(self, dtype=None):
        return np.asarray(self.values, dtype=dtype)

    def __len__(self) -> int:
        return len(self.values)

    def __iter__(self):
        return iter(self.values)

    def __getitem__(self, key):
        return self.values[key]

    @property
    def shape(self):
        return self.values.shape

    @property
    def dtype(self):
        return self.values.dtype

    def to_numpy(self, *, copy: bool = False) -> np.ndarray:
        return np.array(self.values, copy=True) if copy else np.asarray(self.values)

    def info_text(self) -> str:
        values = np.asarray(self.values)
        finite = values[np.isfinite(values)] if np.issubdtype(values.dtype, np.number) else np.asarray([])
        lines = [
            "SeriesView",
            f"  source: {self.name}",
            f"  length: {len(self)}",
            f"  unit: {self.unit or '-'}",
        ]
        if len(values):
            lines.append(f"  initial: {values[0]}")
            lines.append(f"  final: {values[-1]}")
        if finite.size:
            lines.append(f"  min: {np.min(finite)}")
            lines.append(f"  max: {np.max(finite)}")
        lines += [
            "", "Common next steps:",
            "  series.values", "  series.time", "  series.snapshot_id",
            "  np.asarray(series)",
        ]
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text

    def __repr__(self) -> str:
        unit = f", unit={self.unit!r}" if self.unit else ""
        return f"SeriesView(name={self.name!r}, shape={self.shape}{unit})"


class _NamedSeries(Mapping[str, np.ndarray]):
    def __init__(self, values: dict[str, np.ndarray]):
        self._values = values
    def __getitem__(self, key: str) -> np.ndarray:
        return self._values[key]
    def __iter__(self):
        return iter(self._values)
    def __len__(self):
        return len(self._values)
    @property
    def names(self) -> tuple[str, ...]:
        return tuple(self._values)
    @property
    def matrix(self) -> np.ndarray:
        if not self._values:
            return np.empty((0, 0), dtype=float)
        return np.column_stack([np.asarray(self._values[name], dtype=float) for name in self._values])
    def keys(self):
        return self.names
    def _ipython_key_completions_(self):
        return list(self.names)


class ConversionSeries(_NamedSeries):
    def __init__(self, values: dict[str, np.ndarray], total: np.ndarray):
        super().__init__(values)
        self.total = total


class PolymerCompositionSeries:
    def __init__(self, run: "StorageRun"):
        self._run = run
    @property
    def instantaneous(self) -> _NamedSeries:
        monomers = self._run.monomer_names
        if len(monomers) == 1:
            values = np.ones(len(self._run.snapshots), dtype=float)
            values[np.asarray(self._run.conc[monomers[0]], dtype=float) <= 0] = np.nan
            return _NamedSeries({monomers[0]: self._run._series(values, f"F_ins_{monomers[0]}", None)})
        if len(monomers) != 2:
            raise DataUnavailableError("F.ins currently supports homo and binary terminal Mayo-Lewis models")
        a,b=monomers; rates=self._run.k
        # Infer terminal propagation pairs from the declared macro topology,
        # never from user-chosen kinetic-parameter names.
        props=[]
        for raw in self._run._input_model_lines():
            line=raw.split("#",1)[0].strip()
            if not line: continue
            parts=line.split()
            if len(parts)>=8 and parts[0]=="macro" and parts[1]=="prop" and "->" in parts:
                arrow=parts.index("->")
                if arrow>=4 and len(parts)>arrow+2:
                    props.append(dict(pool=parts[2], incoming=parts[4], product=parts[arrow+1], rate=parts[arrow+2]))
        terminal_by_pool={}
        for prop in props:
            terminal_by_pool.setdefault(prop["product"], prop["incoming"])
        rate_by_pair={}
        for prop in props:
            terminal=terminal_by_pool.get(prop["pool"])
            incoming=prop["incoming"]
            if terminal in monomers and incoming in monomers:
                rate_by_pair[(terminal,incoming)]=prop["rate"]
        def find(x,y):
            name=rate_by_pair.get((x,y))
            if name is None:
                raise DataUnavailableError(f"F.ins requires a terminal propagation declaration for {x}->{y}")
            try:
                return np.asarray(rates[name],dtype=float)
            except KeyError as exc:
                raise DataUnavailableError(f"F.ins kinetic rate {name!r} is unavailable") from exc
        kaa,kab,kba,kbb=find(a,a),find(a,b),find(b,a),find(b,b)
        fa=np.asarray(self._run.f[a],dtype=float); fb=np.asarray(self._run.f[b],dtype=float)
        ra=np.divide(kaa,kab,out=np.full_like(kaa,np.nan),where=kab!=0); rb=np.divide(kbb,kba,out=np.full_like(kbb,np.nan),where=kba!=0)
        den=ra*fa*fa+2*fa*fb+rb*fb*fb
        Fa=np.divide(ra*fa*fa+fa*fb,den,out=np.full_like(fa,np.nan),where=den!=0)
        return _NamedSeries({a:self._run._series(Fa,f"F_ins_{a}"),b:self._run._series(1-Fa,f"F_ins_{b}")})
    @property
    def interval(self) -> _NamedSeries:
        return self._run._polymer_composition(kind="interval")
    @property
    def cumulative(self) -> _NamedSeries:
        return self._run._polymer_composition(kind="cumulative")
    @property
    def ins(self) -> _NamedSeries:
        return self.instantaneous
    @property
    def int(self) -> _NamedSeries:
        return self.interval
    @property
    def cum(self) -> _NamedSeries:
        return self.cumulative


class _StateFieldSeries(Mapping[str, np.ndarray]):
    def __init__(self, state: "StorageStateSeries", column: str):
        self._state = state
        self._column = column

    def __iter__(self):
        return iter(self._state.names)

    def __len__(self) -> int:
        return len(self._state.names)

    def __getitem__(self, name: str | int) -> np.ndarray:
        entity_id = int(name) if isinstance(name, (int, np.integer)) else self._state.entity_id(name)
        display_name = self._state.names[entity_id] if entity_id < len(self._state.names) else str(entity_id)
        values = self._state.matrix(self._column)[:, entity_id]
        return _readonly_array(values)

    @property
    def names(self) -> tuple[str, ...]:
        return self._state.names

    @property
    def matrix(self) -> np.ndarray:
        return self._state.matrix(self._column)

    def keys(self):
        return self.names


    def _ipython_key_completions_(self):
        return list(self.keys())

class _StateFieldSnapshot(Mapping[str, Any]):
    def __init__(self, state: "StorageStateSnapshot", column: str):
        self._state = state
        self._column = column

    def __iter__(self):
        return iter(self._state.names)

    def __len__(self) -> int:
        return len(self._state.names)

    def __getitem__(self, name: str | int):
        entity_id = int(name) if isinstance(name, (int, np.integer)) else self._state.entity_id(name)
        value = self._state.row_values(self._column)[entity_id]
        return value.item() if hasattr(value, "item") else value

    @property
    def names(self) -> tuple[str, ...]:
        return self._state.names

    @property
    def values(self) -> np.ndarray:
        return self._state.row_values(self._column)

    def keys(self):
        return self.names

    def _ipython_key_completions_(self):
        return list(self.keys())


@dataclass(frozen=True)
class StorageStateSeries:
    run: "StorageRun"

    @property
    def raw(self) -> StorageTable:
        return self.run.table("state")

    @property
    def names(self) -> tuple[str, ...]:
        entries = self.run.dictionaries.get("state_entities", {})
        if entries:
            return tuple(str(entries[i].get("name", i)) for i in sorted(entries))
        raw = self.raw
        if "entity_id" not in raw or raw.n_rows == 0:
            return ()
        n_entities = int(np.max(raw["entity_id"])) + 1
        return tuple(str(i) for i in range(n_entities))

    def entity_id(self, name: str) -> int:
        aliases={str(name),f"monomer_{name}",f"species_{name}"}
        for i, entry in self.run.dictionaries.get("state_entities", {}).items():
            if str(entry.get("name")) in aliases or entry.get("label") == name:
                return int(i)
        available = ", ".join(self.names)
        raise KeyError(f"Unknown state entity {name!r}. Available: {available}")

    def matrix(self, column: str) -> np.ndarray:
        if column not in self.raw:
            raise KeyError(column)
        n_snapshots = len(self.run.snapshots)
        n_entities = len(self.names)
        values = self.raw[column]
        if len(values) != n_snapshots * n_entities:
            raise InvalidOutputError(
                f"state/{column}.npy has {len(values)} rows; expected "
                f"{n_snapshots} snapshots × {n_entities} entities"
            )
        return values.reshape(n_snapshots, n_entities)

    @property
    def counts(self) -> _StateFieldSeries:
        return _StateFieldSeries(self, "count")

    @property
    def moles(self) -> _StateFieldSeries:
        return _StateFieldSeries(self, "moles")

    @property
    def concentrations(self) -> _StateFieldSeries:
        return _StateFieldSeries(self, "concentration")

    @property
    def count(self) -> _StateFieldSeries:
        return self.counts

    @property
    def conc(self) -> _StateFieldSeries:
        return self.concentrations

    def __getitem__(self, name: str):
        return _StateEntitySeries(self, name)


@dataclass(frozen=True)
class _StateEntitySeries:
    state: StorageStateSeries
    name: str

    @property
    def count(self) -> np.ndarray:
        return self.state.counts[self.name]

    @property
    def moles(self) -> np.ndarray:
        return self.state.moles[self.name]

    @property
    def concentration(self) -> np.ndarray:
        return self.state.concentrations[self.name]

    @property
    def conc(self) -> np.ndarray:
        return self.concentration


@dataclass(frozen=True)
class StorageStateSnapshot:
    snapshot: "StorageSnapshot"

    @property
    def run(self) -> "StorageRun":
        return self.snapshot.run

    @property
    def raw(self) -> StorageTable:
        return self.snapshot._rows("state")

    @property
    def names(self) -> tuple[str, ...]:
        return self.run.state.names

    def entity_id(self, name: str) -> int:
        return self.run.state.entity_id(name)

    def row_values(self, column: str) -> np.ndarray:
        ids = np.asarray(self.raw["entity_id"], dtype=np.int64)
        values = np.asarray(self.raw[column])
        expected = np.arange(len(self.names), dtype=np.int64)
        if not np.array_equal(ids, expected):
            raise InvalidOutputError(f"state rows for snapshot {self.snapshot.id} are not dense and ordered")
        return values

    @property
    def counts(self) -> _StateFieldSnapshot:
        return _StateFieldSnapshot(self, "count")

    @property
    def moles(self) -> _StateFieldSnapshot:
        return _StateFieldSnapshot(self, "moles")

    @property
    def concentrations(self) -> _StateFieldSnapshot:
        return _StateFieldSnapshot(self, "concentration")

    @property
    def count(self) -> _StateFieldSnapshot:
        return self.counts

    @property
    def conc(self) -> _StateFieldSnapshot:
        return self.concentrations

    def __getitem__(self, name: str):
        return _StateEntitySnapshot(self, name)


@dataclass(frozen=True)
class _StateEntitySnapshot:
    state: StorageStateSnapshot
    name: str

    @property
    def count(self):
        return self.state.counts[self.name]

    @property
    def moles(self):
        return self.state.moles[self.name]

    @property
    def concentration(self):
        return self.state.concentrations[self.name]

    @property
    def conc(self):
        return self.concentration




class _ChainCompositionCounts(Mapping[str, np.ndarray]):
    def __init__(self, composition: "StorageChainComposition"):
        self._composition = composition
    def __getitem__(self, name: str) -> np.ndarray:
        return self._composition.matrix[:, self._composition.monomer_id(name)]
    def __iter__(self):
        return iter(self._composition.names)
    def __len__(self):
        return len(self._composition.names)
    @property
    def total(self) -> np.ndarray:
        return np.sum(self._composition.matrix, axis=1, dtype=np.uint64)
    @property
    def matrix(self) -> np.ndarray:
        return self._composition.matrix
    @property
    def names(self) -> tuple[str, ...]:
        return self._composition.names


    def _ipython_key_completions_(self):
        return list(self.keys())

class _ChainCompositionFractions(_ChainCompositionCounts):
    def __getitem__(self, name: str) -> np.ndarray:
        counts = np.asarray(super().__getitem__(name), dtype=float)
        total = np.asarray(self.total, dtype=float)
        return np.divide(counts, total, out=np.full(counts.shape, np.nan), where=total != 0)
    @property
    def matrix(self) -> np.ndarray:
        counts = np.asarray(self._composition.matrix, dtype=float)
        total = np.sum(counts, axis=1)
        return np.divide(counts, total[:, None], out=np.full(counts.shape, np.nan), where=total[:, None] != 0)


@dataclass(frozen=True)
class StorageChainComposition:
    chains: "StorageChains"

    @property
    def names(self) -> tuple[str, ...]:
        return self.chains.run.monomer_names

    def monomer_id(self, name: str) -> int:
        for i, candidate in enumerate(self.names):
            if candidate == name:
                return i
        raise KeyError(f"Unknown monomer {name!r}. Available: {', '.join(self.names)}")

    @property
    def matrix(self) -> np.ndarray:
        n = len(self.chains)
        m = len(self.names)
        out = np.zeros((n, m), dtype=np.uint64)
        if n == 0:
            return out
        if "chain_composition" not in self.chains.run.tables:
            if m == 1 and "dp" in self.chains.raw:
                out[:, 0] = np.asarray(self.chains.dp, dtype=np.uint64)
                return out
            raise DataUnavailableError("chain_composition table is unavailable")
        comp = self.chains.run.table("chain_composition")
        rows = {int(rid): i for i, rid in enumerate(np.asarray(self.chains.chain_record_id, dtype=np.uint64))}
        for rid, mid, units in zip(comp["chain_record_id"], comp["monomer_id"], comp["unit_count"]):
            pos = rows.get(int(rid))
            if pos is not None:
                mid_i = int(mid)
                if mid_i < 0 or mid_i >= m:
                    raise InvalidOutputError("chain_composition contains an out-of-range monomer_id")
                out[pos, mid_i] = int(units)
        return out

    @property
    def counts(self) -> _ChainCompositionCounts:
        return _ChainCompositionCounts(self)

    @property
    def fractions(self) -> _ChainCompositionFractions:
        return _ChainCompositionFractions(self)


@dataclass(frozen=True)
class StorageChainRecord:
    chains: "StorageChains"
    index: int

    def _value(self, name: str):
        value = self.chains.raw[name][self.index]
        return value.item() if hasattr(value, "item") else value

    @property
    def chain_record_id(self) -> int: return int(self._value("chain_record_id"))
    @property
    def snapshot_id(self) -> int: return int(self._value("snapshot_id"))
    @property
    def dp(self) -> int: return int(self._value("dp"))
    @property
    def molar_mass(self) -> float: return float(self._value("molar_mass"))
    @property
    def count(self) -> int: return int(self._value("count"))
    @property
    def moles(self) -> float: return float(self._value("moles"))
    @property
    def conc(self) -> float: return float(self._value("concentration"))
    @property
    def population(self) -> str: return self.chains._decoded_at("population_id", "chain_populations", self.index)
    @property
    def pool(self) -> str: return self.chains._decoded_at("pool_id", "chain_pools", self.index)
    @property
    def origin(self) -> str: return self.chains._decoded_at("origin_id", "chain_origins", self.index)
    @property
    def left_end(self) -> str: return self.chains._decoded_at("left_end_id", "chain_end_types", self.index)
    @property
    def right_end(self) -> str: return self.chains._decoded_at("right_end_id", "chain_end_types", self.index)
    @property
    def composition(self):
        matrix = self.chains.composition.matrix[self.index]
        names = self.chains.composition.names
        counts = {name: int(matrix[i]) for i, name in enumerate(names)}
        total = int(np.sum(matrix))
        fractions = {name: (counts[name] / total if total else np.nan) for name in names}
        return _ChainRecordComposition(counts, fractions)
    @property
    def first_monomer(self): return self.chains.first_monomer[self.index]
    @property
    def penultimate_monomer(self): return self.chains.penultimate_monomer[self.index]
    @property
    def last_monomer(self): return self.chains.last_monomer[self.index]
    @property
    def sequence(self) -> tuple[str, ...]:
        return self.chains._sequence_at(self.index)


@dataclass(frozen=True)
class _ChainRecordComposition:
    counts: dict[str, int]
    fractions: dict[str, float]


class StorageChainOrigin:
    def __init__(self, chains: "StorageChains"):
        self.chains = chains

    def __call__(self, name: str) -> "StorageChains":
        return self.chains._filter(self.chains.raw["origin_id"] == self.chains._id_for_name("chain_origins", name))

    def summary(self):
        from .table import Table
        names = self.chains.origin_names
        groups = {}
        for name, count in zip(names, self.chains.count):
            rec = groups.setdefault(str(name), [0, 0])
            rec[0] += 1; rec[1] += int(count)
        total = sum(v[1] for v in groups.values())
        rows = [(name, nr, nc, nc/total if total else float("nan")) for name, (nr, nc) in sorted(groups.items())]
        return Table(["origin", "n_records", "n_chains", "fraction"], rows, name="origin_summary")


class StorageChains(ChainPopulation):
    """High-level view of compressed chain records."""
    def __init__(self, run: "StorageRun", raw: StorageTable | None = None, *,
                 _analysis_root: "StorageChains | None" = None,
                 _root_indices: np.ndarray | None = None):
        self.run = run
        self.raw = run.table("chains") if raw is None else raw
        self._decoded_cache = {}
        # Analysis caches belong to the unfiltered population. Filtered views keep
        # row indices into that population, so expensive per-record statistics can
        # be computed once and sliced without rescanning stored sequences.
        self._analysis_root = self if _analysis_root is None else _analysis_root
        if _root_indices is None:
            _root_indices = np.arange(self.raw.n_rows, dtype=np.int64)
        self._root_indices = np.asarray(_root_indices, dtype=np.int64)
        self._root_indices.setflags(write=False)

    def __len__(self) -> int: return self.raw.n_rows
    @property
    def n_records(self) -> int:
        """Number of compressed chain records in this view."""
        return len(self)
    @property
    def n_chains(self) -> int:
        """Physical chain population represented by the compressed records."""
        return int(np.asarray(self.count, dtype=np.int64).sum()) if len(self) else 0
    def __iter__(self):
        for i in range(len(self)):
            yield StorageChainRecord(self, i)
    def __getitem__(self, index: int) -> StorageChainRecord:
        if index < 0: index += len(self)
        if index < 0 or index >= len(self): raise IndexError(index)
        return StorageChainRecord(self, index)
    def __getattr__(self, name: str):
        if name in self.raw: return self.raw[name]
        raise AttributeError(name)

    def _dict_name(self, dictionary: str, value: int) -> str:
        entry = self.run.dictionaries.get(dictionary, {}).get(int(value))
        return str(entry.get("name", value)) if entry else str(value)
    def _decoded(self, column: str, dictionary: str) -> np.ndarray:
        key=(column,dictionary)
        cached=self._decoded_cache.get(key)
        if cached is not None: return cached
        ids=np.asarray(self.raw[column],dtype=np.int64)
        entries=self.run.dictionaries.get(dictionary,{})
        if ids.size == 0:
            out=np.empty(0,dtype=object)
        else:
            max_id=max(int(ids.max()), max(entries, default=-1))
            lookup=np.empty(max_id+1,dtype=object)
            for i in range(max_id+1):
                entry=entries.get(i)
                lookup[i]=str(entry.get("name",i)) if entry else str(i)
            out=lookup[ids]
        out.setflags(write=False)
        self._decoded_cache[key]=out
        return out
    def _decoded_at(self, column: str, dictionary: str, index: int) -> str:
        return self._dict_name(dictionary, int(self.raw[column][index]))
    def _id_for_name(self, dictionary: str, name: str) -> int:
        for i, entry in self.run.dictionaries.get(dictionary, {}).items():
            if entry.get("name") == name or entry.get("label") == name:
                return int(i)
        raise KeyError(f"Unknown {dictionary} value {name!r}")
    def _filter(self, mask: np.ndarray) -> "StorageChains":
        mask = np.asarray(mask, dtype=bool)
        return StorageChains(
            self.run,
            self.raw.filtered(mask),
            _analysis_root=self._analysis_root,
            _root_indices=self._root_indices[mask],
        )

    def _with_mask(self, mask: np.ndarray) -> "StorageChains":
        return self._filter(mask)

    def _raw_arrays(self) -> Mapping[str, np.ndarray]:
        arrays = dict(self.raw._columns)
        if "molar_mass" in arrays:
            arrays.setdefault("mass", arrays["molar_mass"])
        arrays.setdefault("left_end", self.left_end)
        arrays.setdefault("right_end", self.right_end)
        return arrays

    def _count_arrays(self) -> Mapping[str, np.ndarray]:
        matrix = self.composition.matrix
        return {name: matrix[:, i] for i, name in enumerate(self.composition.names)}

    @property
    def snapshot_id(self) -> int:
        if len(self) == 0:
            return int(self.run.last.id)
        values = np.unique(np.asarray(self.raw["snapshot_id"], dtype=np.uint64))
        if len(values) != 1:
            raise ValueError("distribution methods require chains from exactly one snapshot")
        return int(values[0])

    @property
    def t(self) -> float:
        return float(self.run.snapshots[self.snapshot_id].t)

    @property
    def kmc_event(self) -> int:
        return int(self.run.snapshots[self.snapshot_id].event)

    def masses(self, *, mass_model: str = "with_end_groups") -> np.ndarray:
        if mass_model not in {"with_end_groups", "repeat_units"}:
            raise ValueError("mass_model must be 'repeat_units' or 'with_end_groups'")
        if mass_model == "with_end_groups":
            return np.asarray(self.molar_mass, dtype=float)
        # Reconstruct repeat-unit-only mass from composition dictionary metadata.
        masses = np.zeros(len(self), dtype=float)
        matrix = self.composition.matrix.astype(float)
        entries = self.run.dictionaries.get("monomers", {})
        increments=[]
        for i, name in enumerate(self.composition.names):
            entry=entries.get(i, {})
            value=entry.get("molar_mass_increment", entry.get("molar_mass", None))
            if value is None:
                # Fallback: full stored mass remains scientifically safer than guessing.
                return np.asarray(self.molar_mass, dtype=float)
            increments.append(float(value))
        if increments:
            masses = matrix @ np.asarray(increments, dtype=float)
        return masses

    def select(self, *, pool: str) -> "StorageChains":
        if pool == "all": return self
        if pool in {"live", "active"}: return self.live
        if pool == "dead": return self.dead
        return self.pool(pool)

    @property
    def all(self) -> "StorageChains": return self
    @property
    def live(self) -> "StorageChains": return self.population_activity("live")
    @property
    def dead(self) -> "StorageChains": return self.population_activity("dead")
    def population_activity(self, name: str) -> "StorageChains":
        return self._filter(self.raw["population_id"] == self._id_for_name("chain_populations", name))
    def pool(self, name: str) -> "StorageChains":
        return self._filter(self.raw["pool_id"] == self._id_for_name("chain_pools", name))
    @property
    def origin(self) -> StorageChainOrigin:
        return StorageChainOrigin(self)
    def where(self, *, dp_min: int | None = None, dp_max: int | None = None) -> "StorageChains":
        mask = np.ones(len(self), dtype=bool)
        if dp_min is not None: mask &= np.asarray(self.dp) >= dp_min
        if dp_max is not None: mask &= np.asarray(self.dp) <= dp_max
        return self._filter(mask)
    def where_count(self, monomer: str, *, min: int | None = None, max: int | None = None) -> "StorageChains":
        from .composition_analysis import _validate_bounds
        _validate_bounds(min, max, name="where_count")
        values = np.asarray(self.composition.counts[monomer])
        mask = np.ones(len(self), dtype=bool)
        if min is not None: mask &= values >= int(min)
        if max is not None: mask &= values <= int(max)
        return self._filter(mask)

    def where_fraction(self, monomer: str, *, min: float | None = None, max: float | None = None) -> "StorageChains":
        from .composition_analysis import _validate_bounds
        _validate_bounds(min, max, name="where_fraction")
        values = np.asarray(self.composition.fractions[monomer], dtype=float)
        mask = np.isfinite(values)
        if min is not None: mask &= values >= float(min)
        if max is not None: mask &= values <= float(max)
        return self._filter(mask)

    @property
    def component_count(self) -> np.ndarray:
        out = np.count_nonzero(np.asarray(self.composition.matrix) > 0, axis=1).astype(np.uint16)
        out.setflags(write=False)
        return out

    def where_component_count(self, *, min: int | None = None, max: int | None = None) -> "StorageChains":
        from .composition_analysis import _validate_bounds
        _validate_bounds(min, max, name="where_component_count")
        values = self.component_count
        mask = np.ones(len(self), dtype=bool)
        if min is not None: mask &= values >= int(min)
        if max is not None: mask &= values <= int(max)
        return self._filter(mask)

    def where_components(self, components, *, exact: bool = True) -> "StorageChains":
        requested = tuple(components)
        if not requested:
            raise ValueError("components must not be empty")
        known = self.composition.names
        unknown = [name for name in requested if name not in known]
        if unknown:
            raise KeyError(f"Unknown monomer {unknown[0]!r}. Available: {', '.join(known)}")
        present = np.asarray(self.composition.matrix) > 0
        mask = np.ones(len(self), dtype=bool)
        for name in requested:
            mask &= present[:, known.index(name)]
        if exact:
            mask &= np.count_nonzero(present, axis=1) == len(set(requested))
        return self._filter(mask)

    def composition_by_dp(self, *, bins=None):
        from .composition_analysis import composition_by_dp
        return composition_by_dp(self, bins=bins)

    def composition_dp_map(self, monomer: str, *, dp_bins=None, fraction_bins=None):
        from .composition_analysis import composition_dp_map
        return composition_dp_map(self, monomer, dp_bins=dp_bins, fraction_bins=fraction_bins)

    def composition_mass_map(self, monomer: str, *, mass_model="with_end_groups", mass_bins=None, fraction_bins=None):
        from .composition_analysis import composition_mass_map
        return composition_mass_map(self, monomer, mass_model=mass_model, mass_bins=mass_bins, fraction_bins=fraction_bins)

    def composition_map(self, x: str, y: str, *, bins=None):
        from .composition_analysis import composition_map
        return composition_map(self, x, y, bins=bins)

    def component_classes(self):
        from .composition_analysis import component_classes
        return component_classes(self)

    def sequence_stats(self, *, progress=None):
        from .full_analysis import sequence_stats
        return sequence_stats(self, progress=progress)

    def where_transition_count(self, *, min=None, max=None, progress=None):
        from .composition_analysis import _validate_bounds
        _validate_bounds(min, max, name="where_transition_count")
        values = self.sequence_stats(progress=progress).transition_count
        mask = np.ones(len(self), dtype=bool)
        if min is not None: mask &= values >= min
        if max is not None: mask &= values <= max
        return self._filter(mask)

    def where_transition_fraction(self, *, min=None, max=None, progress=None):
        from .composition_analysis import _validate_bounds
        _validate_bounds(min, max, name="where_transition_fraction")
        values = self.sequence_stats(progress=progress).transition_fraction
        mask = np.ones(len(self), dtype=bool)
        if min is not None: mask &= values >= min
        if max is not None: mask &= values <= max
        return self._filter(mask)

    def where_block_count(self, monomer, *, min=None, max=None, progress=None):
        from .composition_analysis import _validate_bounds
        _validate_bounds(min, max, name="where_block_count")
        values = self.sequence_stats(progress=progress).block_count[monomer]
        mask = np.ones(len(self), dtype=bool)
        if min is not None: mask &= values >= min
        if max is not None: mask &= values <= max
        return self._filter(mask)

    def where_max_block(self, monomer, *, min=None, max=None, progress=None):
        from .composition_analysis import _validate_bounds
        _validate_bounds(min, max, name="where_max_block")
        values = self.sequence_stats(progress=progress).max_block_length[monomer]
        mask = np.ones(len(self), dtype=bool)
        if min is not None: mask &= values >= min
        if max is not None: mask &= values <= max
        return self._filter(mask)

    def _motif_tokens(self, motif):
        if isinstance(motif, str):
            if "|" in motif:
                motif = tuple(part for part in motif.split("|") if part)
            elif all(len(name) == 1 for name in self.composition.names):
                motif = tuple(motif)
            else:
                raise ValueError("string motifs for multi-character monomer names must use '|' separators")
        else:
            motif = tuple(motif)
        if not motif:
            raise ValueError("motif must not be empty")
        unknown = [name for name in motif if name not in self.composition.names]
        if unknown:
            raise KeyError(f"Unknown monomer {unknown[0]!r}")
        return motif

    def contains_motif(self, motif, *, min_occurrences=1):
        motif = self._motif_tokens(motif)
        mask=[]
        for seq in self.sequences:
            n=sum(tuple(seq[i:i+len(motif)]) == motif for i in range(max(0,len(seq)-len(motif)+1)))
            mask.append(n >= min_occurrences)
        return self._filter(np.asarray(mask,dtype=bool))

    def starts_with(self, motif):
        motif = self._motif_tokens(motif)
        return self._filter(np.asarray([tuple(seq[:len(motif)]) == motif for seq in self.sequences],dtype=bool))

    def ends_with(self, motif):
        motif = self._motif_tokens(motif)
        return self._filter(np.asarray([tuple(seq[-len(motif):]) == motif for seq in self.sequences],dtype=bool))

    def block_lengths(self, monomer=None, *, progress=None):
        from .full_analysis import block_lengths
        return block_lengths(self, monomer, progress=progress)

    def dyads_by_dp(self, *, bins=16):
        from .full_analysis import dyads_by_dp
        return dyads_by_dp(self, bins=bins)

    def triads_by_composition(self, monomer: str, *, bins=12):
        from .full_analysis import triads_by_composition
        return triads_by_composition(self, monomer, bins=bins)

    def block_count(self, monomer=None, *, progress=None):
        stats = self.sequence_stats(progress=progress)
        if monomer is not None:
            return stats.block_count[monomer]
        arrays = tuple(stats.block_count.values())
        result = np.sum(np.stack(arrays, axis=0), axis=0) if arrays else np.zeros(len(self), dtype=np.int64)
        return _readonly_array(result)

    def junction_positions(self, left: str, right: str) -> tuple[tuple[int, ...], ...]:
        if left not in self.composition.names or right not in self.composition.names:
            raise KeyError(f"Unknown monomer pair {left!r}, {right!r}")
        return tuple(tuple(i + 1 for i, (a, b) in enumerate(zip(seq[:-1], seq[1:])) if a == left and b == right) for seq in self.sequences)

    def junction_position(self, left: str, right: str) -> np.ndarray:
        values = [positions[0] if len(positions) == 1 else np.nan for positions in self.junction_positions(left, right)]
        return _readonly_array(np.asarray(values, dtype=float))

    def transition_matrix(self, *, normalize=None, progress=None):
        from .full_analysis import transition_matrix
        return transition_matrix(self, normalize=normalize, progress=progress)

    def microstructure_by_dp(self, statistic, *, monomer=None, bins=None, progress=None):
        from .full_analysis import microstructure_by_dp
        return microstructure_by_dp(self, statistic, monomer=monomer, bins=bins, progress=progress)

    def motif_counts(self, motif, *, progress=None):
        from .full_analysis import motif_counts
        return motif_counts(self, motif, progress=progress)

    def ngrams(self, n=4, *, min_count=1, progress=None):
        from .full_analysis import ngrams
        return ngrams(self, n=n, min_count=min_count, progress=progress)

    def position_profile(self, *, bins=20, progress=None):
        from .full_analysis import position_profile
        return position_profile(self, bins=bins, progress=progress)

    def microstructure_map(self, statistic, *, monomer=None, dp_bins=None, value_bins=None, progress=None):
        from .full_analysis import microstructure_map
        return microstructure_map(self, statistic, monomer=monomer, dp_bins=dp_bins, value_bins=value_bins, progress=progress)

    def at_snapshot(self, snapshot_id: int) -> "StorageChains":
        return self._filter(self.raw["snapshot_id"] == snapshot_id)
    @property
    def last(self) -> "StorageChains":
        return self.at_snapshot(self.run.last.id)
    def record(self, chain_record_id: int) -> StorageChainRecord:
        pos = np.flatnonzero(np.asarray(self.chain_record_id) == chain_record_id)
        if len(pos) == 0: raise KeyError(f"Unknown chain_record_id {chain_record_id}")
        return StorageChainRecord(self, int(pos[0]))

    @property
    def chain_record_id(self) -> np.ndarray: return _readonly_array(self.raw["chain_record_id"])
    @property
    def dp(self) -> np.ndarray: return _readonly_array(self.raw["dp"])
    @property
    def molar_mass(self) -> np.ndarray: return _readonly_array(self.raw["molar_mass"])
    @property
    def count(self) -> np.ndarray: return _readonly_array(self.raw["count"])
    @property
    def moles(self) -> np.ndarray: return _readonly_array(self.raw["moles"])
    @property
    def conc(self) -> np.ndarray: return _readonly_array(self.raw["concentration"])
    @property
    def population_activity_names(self) -> np.ndarray: return self._decoded("population_id", "chain_populations")
    @property
    def pool_names(self) -> np.ndarray: return self._decoded("pool_id", "chain_pools")
    @property
    def origin_names(self) -> np.ndarray: return self._decoded("origin_id", "chain_origins")
    @property
    def left_end(self) -> np.ndarray: return self._decoded("left_end_id", "chain_end_types")
    @property
    def right_end(self) -> np.ndarray: return self._decoded("right_end_id", "chain_end_types")
    @property
    def composition(self) -> StorageChainComposition: return StorageChainComposition(self)

    def info_text(self) -> str:
        lines = [
            "ChainsView",
            f"  run: {self.run._display_path()}",
            f"  snapshot_id: {self.snapshot_id if len(self) else '-'}",
            f"  records: {len(self)}",
            f"  total count: {int(np.sum(self.count)) if len(self) else 0}",
        ]
        if len(self):
            lines += [
                f"  DP range: {int(np.min(self.dp))}-{int(np.max(self.dp))}",
                f"  mean DP: {float(np.average(self.dp, weights=self.count)):.6g}",
                f"  sequence data: {'full' if self.has_sequences else 'composition'}",
            ]
        lines += [
            "", "Available fields:",
            "  chains.dp", "  chains.molar_mass", "  chains.count",
            "  chains.composition", "  chains.left_end", "  chains.right_end",
            "", "Common next steps:",
            "  chains.live", "  chains.dead", "  chains.where(dp_min=10)",
            "  chains.mwd()", "  chains.cld()",
        ]
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text

    @analysis_operation(CLD_HELP)
    def cld(self, *, mass_model: str | None = None, series=None, **kwargs):
        return ChainPopulation.cld.__get__(self, type(self))(mass_model=mass_model, series=series, **kwargs)

    @analysis_operation(MWD_HELP)
    def mwd(self, *, mass_model: str | None = None, series=None, **kwargs):
        if str(kwargs.get("method", "gaussian")).lower() == "sticks" and "sigma" not in kwargs:
            kwargs["sigma"] = None
        return ChainPopulation.mwd.__get__(self, type(self))(mass_model=mass_model, series=series, **kwargs)

    @analysis_operation(SPECTRUM_HELP)
    def chain_mass_spectrum(self, *, mass_model: str | None = None, series=None, normalize: str = "count", **kwargs):
        return ChainPopulation.chain_mass_spectrum.__get__(self, type(self))(mass_model=mass_model, series=series, normalize=normalize, **kwargs)

    def _monomer_field(self, prefix: str) -> np.ndarray:
        has_col = f"has_{prefix}"
        id_col = f"{prefix}_id"
        if id_col not in self.raw: return np.full(len(self), None, dtype=object)
        out=[]
        for i, value in enumerate(self.raw[id_col]):
            if has_col in self.raw and not bool(self.raw[has_col][i]): out.append(None)
            else: out.append(self._dict_name("monomers", int(value)))
        return np.asarray(out, dtype=object)
    @property
    def first_monomer(self) -> np.ndarray: return self._monomer_field("first_monomer")
    @property
    def penultimate_monomer(self) -> np.ndarray: return self._monomer_field("penultimate_monomer")
    @property
    def last_monomer(self) -> np.ndarray: return self._monomer_field("last_monomer")
    @property
    def has_sequences(self) -> bool:
        return bool(len(self)) and "has_sequence" in self.raw and bool(np.all(self.raw["has_sequence"]))
    def _sequence_at(self, index: int) -> tuple[str, ...]:
        if "has_sequence" not in self.raw or not bool(self.raw["has_sequence"][index]):
            raise DataUnavailableError("Full sequence is unavailable for this chain record")
        symbols = self.run.table("sequences")["symbols"]
        offset = int(self.raw["sequence_offset"][index]); length = int(self.raw["sequence_length"][index])
        ids = symbols[offset:offset+length]
        dictionary = "sequence_symbols" if "sequence_symbols" in self.run.dictionaries else "monomers"
        return tuple(self._dict_name(dictionary, int(v)) for v in ids)
    @property
    def sequences(self) -> tuple[tuple[str, ...], ...]:
        if len(self) == 0: return ()
        if "has_sequence" not in self.raw or not bool(np.all(self.raw["has_sequence"])):
            raise DataUnavailableError("Full sequences are unavailable for one or more selected chain records")
        return tuple(self._sequence_at(i) for i in range(len(self)))




class _MomentLeaf:
    def __init__(self, owner, population: str, mass_basis: str):
        self.owner=owner; self.population=population; self.mass_basis=mass_basis
    def _series(self, name: str): return self.owner._get(self.population, self.mass_basis, name)
    @property
    def chain_count(self): return self._series("chain_count")
    @property
    def dp_n(self): return self._series("dp_n")
    @property
    def dp_w(self): return self._series("dp_w")
    @property
    def mn(self): return self._series("mn")
    @property
    def mw(self): return self._series("mw")
    @property
    def mz(self): return self._series("mz")
    @property
    def dispersity(self): return self._series("dispersity")

class _MassBasisBranch:
    def __init__(self, owner, population: str): self.owner=owner; self.population=population
    @property
    def repeat_units(self): return _MomentLeaf(self.owner, self.population, "repeat_units")
    @property
    def with_end_groups(self): return _MomentLeaf(self.owner, self.population, "with_end_groups")

class StorageMomentsSeries:
    def __init__(self, run): self.run=run
    @property
    def raw(self): return self.run.table("moments")
    def _id(self, dictionary, name):
        for i,e in self.run.dictionaries.get(dictionary, {}).items():
            if e.get("name")==name: return int(i)
        defaults={"population_scope":{"all":0,"live":1,"dead":2},"mass_bases":{"repeat_units":0,"with_end_groups":1}}
        return defaults[dictionary][name]
    def _get(self,population,mass_basis,column):
        p=self._id("population_scope",population); b=self._id("mass_bases",mass_basis)
        out=np.full(len(self.run.snapshots), np.nan, dtype=float)
        tab=self.raw; sid_to_pos={int(s):i for i,s in enumerate(self.run.sid)}
        mask=(np.asarray(tab["population_scope_id"])==p)&(np.asarray(tab["mass_basis_id"])==b)
        for sid,val in zip(np.asarray(tab["snapshot_id"])[mask], np.asarray(tab[column])[mask]):
            if int(sid) in sid_to_pos: out[sid_to_pos[int(sid)]]=val
        unit={"mn":"g/mol","mw":"g/mol","mz":"g/mol"}.get(column)
        return self.run._series(out, column, unit)
    @property
    def all(self): return _MassBasisBranch(self,"all")
    @property
    def live(self): return _MassBasisBranch(self,"live")
    @property
    def dead(self): return _MassBasisBranch(self,"dead")
    def select(self, *, population_scope="all", mass_basis="with_end_groups"):
        return _MomentLeaf(self,population_scope,mass_basis)
    @property
    def default(self): return self.all.with_end_groups

    def info_text(self) -> str:
        names = ("dpn", "dpw", "mn", "mw", "mz", "dispersity")
        lines = ["MomentsView", "  scope: series over run", "  default: all / with_end_groups",
                 "  populations: all, live, dead", "  mass bases: repeat_units, with_end_groups",
                 "  fields: " + ", ".join(names), "", "Common next steps:",
                 "  run.moments.all.with_end_groups.mw", "  run.moments.live.repeat_units.dpn",
                 "  run.moments.select(population=\"dead\", mass_basis=\"with_end_groups\")"]
        return "\n".join(lines)

    def info(self) -> str:
        text=self.info_text(); print(text); return text


class StorageMomentsSnapshot:
    def __init__(self, snapshot): self.snapshot=snapshot
    @property
    def raw(self):
        tab=self.snapshot.run.table("moments")
        return tab.filtered(np.asarray(tab["snapshot_id"]) == self.snapshot.id)
    def __getattr__(self, name):
        if name in self.raw:
            return self.raw[name]
        raise AttributeError(name)
    def _get(self,population,mass_basis,column):
        series=StorageMomentsSeries(self.snapshot.run)._get(population,mass_basis,column)
        pos=np.flatnonzero(np.asarray(self.snapshot.run.sid)==self.snapshot.id)
        if len(pos)==0: raise DataUnavailableError("snapshot moments are unavailable")
        return np.asarray(series)[int(pos[0])].item()
    @property
    def all(self): return _MassBasisBranch(self,"all")
    @property
    def live(self): return _MassBasisBranch(self,"live")
    @property
    def dead(self): return _MassBasisBranch(self,"dead")
    def select(self, *, population_scope="all", mass_basis="with_end_groups"):
        return _MomentLeaf(self,population_scope,mass_basis)
    @property
    def default(self): return self.all.with_end_groups




class _DictionaryMappedSeries(Mapping[str, np.ndarray]):
    def __init__(self, run: "StorageRun", dictionary: str, matrix: np.ndarray, *, name_prefix: str, unit: str | None = None):
        self.run=run; self.dictionary=dictionary; self._matrix=np.asarray(matrix); self.name_prefix=name_prefix; self.unit=unit
        entries=run.dictionaries.get(dictionary,{})
        self._names=(tuple(str(entries[i].get("name",i)) for i in sorted(entries))
                     if entries else tuple(str(i) for i in range(self._matrix.shape[1])))
    def __iter__(self): return iter(self._names)
    def __len__(self): return len(self._names)
    def __getitem__(self,key):
        if isinstance(key,(int,np.integer)): idx=int(key); name=self._names[idx]
        else:
            name=str(key)
            try: idx=self._names.index(name)
            except ValueError as exc: raise KeyError(name) from exc
        return self.run._series(self._matrix[:,idx],f"{self.name_prefix}_{name}",self.unit)
    @property
    def names(self): return self._names
    @property
    def matrix(self): return self._matrix
    def keys(self): return self._names
    def _ipython_key_completions_(self): return list(self._names)

class StorageChannelsSeries:
    def __init__(self,run):
        self.run=run
        self._matrix_cache={}
    @property
    def raw(self): return self.run.table("channel_events")
    def _matrix(self,column):
        cached=self._matrix_cache.get(column)
        if cached is not None: return cached
        n=len(self.run.snapshots); entries=self.run.dictionaries.get("channels",{})
        nchan=(max(entries)+1) if entries else (int(np.max(self.raw["channel_id"]))+1 if self.raw.n_rows else 0)
        out=np.zeros((n,nchan),dtype=np.uint64); sidpos={int(s):i for i,s in enumerate(self.run.sid)}
        for sid,cid,val in zip(self.raw["snapshot_id"],self.raw["channel_id"],self.raw[column]): out[sidpos[int(sid)],int(cid)]=val
        out.setflags(write=False)
        self._matrix_cache[column]=out
        return out
    @property
    def event_count(self): return _DictionaryMappedSeries(self.run,"channels",self._matrix("event_count"),name_prefix="event_count")
    @property
    def productive(self): return _DictionaryMappedSeries(self.run,"channels",self._matrix("productive_event_count"),name_prefix="productive_event_count")
    @property
    def nonproductive(self): return _DictionaryMappedSeries(self.run,"channels",self._matrix("nonproductive_event_count"),name_prefix="nonproductive_event_count")
    @property
    def productive_event_count(self): return self.productive
    @property
    def nonproductive_event_count(self): return self.nonproductive
    def interval_event_counts(self):
        m=self._matrix("event_count").astype(np.int64); return np.vstack([m[0],np.diff(m,axis=0)]) if len(m) else m
    def fire_shares(self):
        d=self.interval_event_counts().astype(float); den=d.sum(axis=1,keepdims=True)
        return np.divide(d,den,out=np.full(d.shape,np.nan),where=den!=0)

    def info_text(self) -> str:
        names = self.event_count.names
        shown = ", ".join(names[:8]) + (" ..." if len(names) > 8 else "")
        lines=["ChannelsView", "  scope: series over run", f"  channels: {len(names)}", f"  names: {shown or '-'}",
               "", "Common next steps:", "  run.channels.event_count[\"channel_name\"]",
               "  run.channels.productive[\"channel_name\"]", "  run.channels.interval_event_counts()",
               "  run.channels.fire_shares()"]
        return "\n".join(lines)

    def info(self) -> str:
        text=self.info_text(); print(text); return text

class _SnapshotNamedValues(Mapping):
    def __init__(self,names,values): self.names=tuple(names); self.values=np.asarray(values)
    def __iter__(self): return iter(self.names)
    def __len__(self): return len(self.names)
    def __getitem__(self,key):
        if isinstance(key,(int,np.integer)): return self.values[int(key)].item()
        try: i=self.names.index(str(key))
        except ValueError as exc: raise KeyError(key) from exc
        return self.values[i].item()
    def keys(self): return self.names
    def _ipython_key_completions_(self): return list(self.names)

class _SnapshotConversion(_SnapshotNamedValues):
    def __init__(self, names, values, total):
        super().__init__(names, values)
        self.total = float(total)


class _SnapshotPolymerComposition:
    def __init__(self, snapshot):
        self.snapshot = snapshot
    def _at(self, source):
        names = source.names
        return _SnapshotNamedValues(names, [np.asarray(source[name])[self.snapshot.index] for name in names])
    @property
    def instantaneous(self): return self._at(self.snapshot.run.F.ins)
    @property
    def interval(self): return self._at(self.snapshot.run.F.int)
    @property
    def cumulative(self): return self._at(self.snapshot.run.F.cum)
    @property
    def ins(self): return self.instantaneous
    @property
    def int(self): return self.interval
    @property
    def cum(self): return self.cumulative



class StorageChannelsSnapshot:
    def __init__(self,snapshot): self.snapshot=snapshot
    def _mapping(self,column):
        series=getattr(StorageChannelsSeries(self.snapshot.run),column)
        pos=self.snapshot.index
        return _SnapshotNamedValues(series.names,[np.asarray(series[name])[pos] for name in series.names])
    @property
    def event_count(self): return self._mapping("event_count")
    @property
    def productive(self): return self._mapping("productive")
    @property
    def nonproductive(self): return self._mapping("nonproductive")

class StorageKineticsSeries:
    def __init__(self,run): self.run=run
    @property
    def definitions(self): return self.run.dictionaries.get("kinetic_parameter_definitions",{})
    @property
    def names(self): return tuple(str(self.definitions[i].get("name",i)) for i in sorted(self.definitions))
    def _parameter_matrix(self):
        tab=self.run.table("kinetic_parameters/values"); nsets=(int(np.max(tab["kinetic_parameter_set_id"]))+1 if tab.n_rows else 0); npar=len(self.names)
        m=np.full((nsets,npar),np.nan)
        for sid,pid,val in zip(tab["kinetic_parameter_set_id"],tab["kinetic_parameter_id"],tab["value"]): m[int(sid),int(pid)]=val
        setids=np.asarray(self.run.snapshots.kinetic_parameter_set_id,dtype=int)
        return m[setids]
    def __iter__(self): return iter(self.names)
    def __len__(self): return len(self.names)
    def keys(self): return self.names
    def _ipython_key_completions_(self): return list(self.names)
    def __getitem__(self,name):
        try: idx=self.names.index(name)
        except ValueError as exc: raise KeyError(name) from exc
        unit=self.definitions[idx].get("unit"); unit=None if unit in (None,"","1","model_defined","engine_native") else str(unit)
        return self.run._series(self._parameter_matrix()[:,idx],name,unit)
    @property
    def temperature(self):
        for i,e in self.definitions.items():
            if e.get("kind")=="temperature": return self[str(e.get("name",i))]
        raise DataUnavailableError("temperature parameter is unavailable")
    def by_kind(self,kind):
        vals={str(e.get("name",i)):self[str(e.get("name",i))] for i,e in sorted(self.definitions.items()) if e.get("kind")==kind}
        return _NamedSeries(vals)
    @property
    def rate_constants(self): return self.by_kind("rate_constant")
    @property
    def arrhenius_A(self): return self.by_kind("arrhenius_A")
    @property
    def arrhenius_Ea(self): return self.by_kind("arrhenius_Ea")
    @property
    def efficiency(self): return self.by_kind("efficiency")

    def info_text(self) -> str:
        names=self.names; shown=", ".join(names[:10]) + (" ..." if len(names)>10 else "")
        kinds=sorted({str(e.get("kind","unknown")) for e in self.definitions.values()})
        lines=["KineticsView", "  scope: series over run", f"  parameters: {len(names)}", f"  kinds: {', '.join(kinds) or '-'}", f"  names: {shown or '-'}",
               "", "Common next steps:", "  run.kinetics.temperature", "  run.k[\"parameter_name\"]",
               "  run.kinetics.rate_constants", "  run.kinetics.arrhenius_A", "  run.kinetics.arrhenius_Ea"]
        return "\n".join(lines)

    def info(self) -> str:
        text=self.info_text(); print(text); return text


class StorageKineticsSnapshot:
    def __init__(self,snapshot): self.snapshot=snapshot
    @property
    def raw(self):
        return self.snapshot._rows("kinetic_parameters/values","kinetic_parameter_set_id",self.snapshot.kinetic_parameter_set_id)
    def __getattr__(self,name):
        if name in self.raw: return self.raw[name]
        raise AttributeError(name)
    def __getitem__(self,name): return np.asarray(StorageKineticsSeries(self.snapshot.run)[name])[self.snapshot.index].item()
    @property
    def temperature(self): return np.asarray(StorageKineticsSeries(self.snapshot.run).temperature)[self.snapshot.index].item()
    @property
    def k(self):
        rs=StorageKineticsSeries(self.snapshot.run).rate_constants
        return {name:np.asarray(rs[name])[self.snapshot.index].item() for name in rs.names}

@dataclass(frozen=True)
class StorageCondition:
    action: "StorageAction"; row: int
    def _v(self,n): return self.action.run.table("action_conditions")[n][self.row].item()
    @property
    def observable(self):
        i=int(self._v("observable_id")); ds=("action_observables","condition_observables")
        for d in ds:
            e=self.action.run.dictionaries.get(d,{}).get(i)
            if e:return str(e.get("name",i))
        return str(i)
    @property
    def operator(self):
        i=int(self._v("operator_id"))
        for d in ("action_operators","condition_operators"):
            e=self.action.run.dictionaries.get(d,{}).get(i)
            if e:return str(e.get("name",i))
        return str(i)
    @property
    def threshold(self): return float(self._v("threshold"))
    @property
    def observed_value(self): return float(self._v("observed_value"))
    @property
    def met(self): return bool(self._v("condition_met"))

@dataclass(frozen=True)
class StorageAction:
    run: "StorageRun"; row: int
    def _v(self,n): return self.run.table("actions")[n][self.row].item()
    @property
    def id(self): return int(self._v("action_id"))
    @property
    def t(self): return float(self._v("time"))
    @property
    def event(self): return int(self._v("kmc_event"))
    @property
    def source_line(self): return int(self._v("source_line"))
    def _dict_name(self,col,*dicts):
        i=int(self._v(col))
        for d in dicts:
            e=self.run.dictionaries.get(d,{}).get(i)
            if e:return str(e.get("name",i))
        return str(i)
    @property
    def type(self): return self._dict_name("action_type_id","action_types")
    @property
    def trigger(self): return self._dict_name("trigger_type_id","action_trigger_types","action_triggers")
    @property
    def before_value(self): return float(self._v("before_value"))
    @property
    def after_value(self): return float(self._v("after_value"))
    @property
    def state_changed(self): return bool(self._v("state_changed"))
    @property
    def output_written(self): return bool(self._v("output_written"))
    @property
    def conditions(self):
        if "action_conditions" not in self.run.tables:return ()
        t=self.run.table("action_conditions"); rows=np.flatnonzero(np.asarray(t["action_id"])==self.id)
        return tuple(StorageCondition(self,int(r)) for r in rows)
    @property
    def message(self): return self.run._action_messages().get(self.id)
    @property
    def snapshot(self):
        if "has_snapshot" not in self.run.table("actions") or not bool(self._v("has_snapshot")): return None
        return self.run.at_snapshot(int(self._v("snapshot_id")))
    @property
    def kinetic_parameter_set_id(self):
        if "has_kinetic_parameter_set" not in self.run.table("actions") or not bool(self._v("has_kinetic_parameter_set")): return None
        return int(self._v("kinetic_parameter_set_id"))

class StorageActions:
    def __init__(self,run): self.run=run
    def __len__(self): return self.run.table("actions").n_rows
    def __iter__(self): return (StorageAction(self.run,i) for i in range(len(self)))
    def __getitem__(self,i): return StorageAction(self.run,int(i))
    @property
    def raw(self): return self.run.table("actions")

    def info_text(self) -> str:
        kinds=[]; triggers=[]
        for a in self:
            try: kinds.append(str(a.type))
            except Exception: pass
            try: triggers.append(str(a.trigger))
            except Exception: pass
        lines=["ActionsView", f"  actions: {len(self)}", f"  types: {', '.join(sorted(set(kinds))) or '-'}", f"  triggers: {', '.join(sorted(set(triggers))) or '-'}",
               "", "Common next steps:", "  run.actions[0]", "  list(run.actions)", "  action.conditions", "  action.snapshot"]
        return "\n".join(lines)

    def info(self) -> str:
        text=self.info_text(); print(text); return text


@dataclass(frozen=True)
class StorageSnapshot:
    """One point in the run chronology."""

    run: "StorageRun"
    index: int

    def _value(self, name: str):
        return self.run.table("snapshots")[name][self.index].item()

    @property
    def id(self) -> int:
        return int(self._value("snapshot_id"))

    @property
    def t(self) -> float:
        return float(self._value("time"))

    @property
    def time(self) -> float:
        return self.t

    @property
    def event(self) -> int:
        return int(self._value("kmc_event")) if "kmc_event" in self.run.table("snapshots") else 0

    @property
    def kmc_event(self) -> int:
        return self.event

    @property
    def reason_id(self) -> int:
        return int(self._value("snapshot_reason_id")) if "snapshot_reason_id" in self.run.table("snapshots") else 0

    @property
    def reason(self) -> str:
        entry = self.run.dictionaries.get("snapshot_reasons", {}).get(self.reason_id)
        return str(entry.get("name")) if entry else str(self.reason_id)

    @property
    def is_final(self) -> bool:
        return bool(self._value("is_final")) if "is_final" in self.run.table("snapshots") else False

    @property
    def has_chains(self) -> bool:
        return bool(self._value("has_chains")) if "has_chains" in self.run.table("snapshots") else False

    @property
    def has_sequences(self) -> bool:
        return bool(self._value("has_sequences")) if "has_sequences" in self.run.table("snapshots") else False

    @property
    def kinetic_parameter_set_id(self) -> int:
        return int(self._value("kinetic_parameter_set_id")) if "kinetic_parameter_set_id" in self.run.table("snapshots") else 0

    def _rows(self, table_name: str, key: str = "snapshot_id", value: int | None = None) -> StorageTable:
        table = self.run.table(table_name)
        if key not in table:
            return StorageTable(table.path, {})
        target = self.id if value is None else value
        return table.filtered(table[key] == target)

    @property
    def state(self) -> StorageStateSnapshot:
        return StorageStateSnapshot(self)

    @property
    def count(self) -> _StateFieldSnapshot:
        return self.state.counts

    @property
    def kmc_volume(self) -> float:
        return float(self._value("kmc_volume_L"))

    @property
    def volume(self) -> float:
        value = float(self._value("volume_mL"))
        if value <= 0:
            raise AnalysisNotApplicableError("physical volume is unavailable because init_volume was not defined")
        return value / 1000.0

    @property
    def moles(self) -> Mapping[str, float]:
        return {name: float(self.conc[name]) * self.volume for name in self.state.names}

    @property
    def conc(self) -> _StateFieldSnapshot:
        return self.state.concentrations

    @property
    def output_status(self) -> OutputStatus:
        available = tuple(sorted(self.tables))
        required = tuple(sorted(
            str(record["name"]) for record in self.schema_records
            if record.get("record_type") == "table" and bool(record.get("required", False))
        ))
        missing = tuple(name for name in required if name not in self.tables)
        invalid: list[str] = []
        if self.status == "completed" and not (self.path / "RESULTS_COMPLETE").is_file():
            invalid.append("RESULTS_COMPLETE")
        has_final = False
        try:
            has_final = bool(np.count_nonzero(np.asarray(self.snapshots.raw["is_final"], dtype=bool)) == 1)
        except Exception as exc:
            invalid.append(f"snapshots: {exc}")
        return OutputStatus(
            available=available,
            missing=missing,
            invalid=tuple(invalid),
            run_completed=self.status == "completed",
            has_final_snapshot=has_final,
        )

    def validate(self, *, strict: bool = False) -> ValidationReport:
        stored = self.diagnostics.validation
        failed = tuple(item.check for item in stored.failed)
        invalid = list(self.output_status.invalid)
        if stored.error_count:
            invalid.append(f"stored validation reports {stored.error_count} error(s)")
        warnings: list[str] = []
        if stored.warning_count:
            warnings.append(f"stored validation reports {stored.warning_count} warning(s)")
        if not stored.records and stored.status not in {"passed", "pass"}:
            warnings.append("stored validation was not run or has no check records")
        report = ValidationReport(
            is_valid=not invalid and stored.error_count == 0,
            is_complete=self.output_status.complete and stored.status in {"passed", "pass"},
            missing_outputs=self.output_status.missing,
            invalid_outputs=tuple(invalid),
            warnings=tuple(warnings),
            failed_checks=failed,
            details=tuple(dict(record) for record in stored.records),
        )
        if strict and (not report.is_valid or not report.is_complete):
            details = report.invalid_outputs or report.missing_outputs or report.warnings
            raise ValidationFailedError("; ".join(details) or "strict validation failed")
        return report

    @property
    def validation(self) -> ValidationReport:
        return self.validate()

    def refresh(self):
        fresh = open_storage(
            self.path,
            allow_incomplete=self.status != "completed",
            mmap_mode="r",
        )
        self._metadata = fresh._metadata
        self.run_id = fresh.run_id
        self._prefix = fresh._prefix
        self.schema_records = fresh.schema_records
        self.dictionaries = fresh.dictionaries
        self._tables = fresh._tables
        self.__dict__.setdefault("_cache", {}).clear()
        return self

    @property
    def monomers(self) -> Mapping[str, Mapping[str, Any]]:
        return {
            str(entry.get("name", idx)).removeprefix("monomer_"): dict(entry)
            for idx, entry in sorted(self.dictionaries.get("monomers", {}).items())
        }

    @property
    def endgroups(self) -> Mapping[str, Mapping[str, Any]]:
        return {
            str(entry.get("name", idx)): dict(entry)
            for idx, entry in sorted(self.dictionaries.get("chain_end_types", {}).items())
        }

    @property
    def monomer_names(self) -> tuple[str, ...]:
        def public_name(value: object) -> str:
            name = str(value)
            return name[len("monomer_"):] if name.startswith("monomer_") else name

        entries = self.dictionaries.get("monomers", {})
        if entries:
            return tuple(public_name(entries[i].get("name", i)) for i in sorted(entries))
        state_entries = self.dictionaries.get("state_entities", {})
        names = [public_name(entry.get("name", i)) for i, entry in sorted(state_entries.items())
                 if entry.get("kind") == "monomer"]
        return tuple(names)

    def _series(self, values: np.ndarray, name: str, unit: str | None = None) -> np.ndarray:
        # name/unit remain available from the owning facade and schema; numeric
        # leaves deliberately terminate in a plain read-only ndarray.
        return _readonly_array(values)

    @property
    def _conversion(self) -> ConversionSeries:
        monomers = self.monomer_names
        if not monomers:
            raise DataUnavailableError("No monomer entities are declared in schema.jsonl")
        by: dict[str, np.ndarray] = {}
        initial = []
        current = []
        for name in monomers:
            values = np.asarray(self.count[name], dtype=float)
            n0 = float(values[0])
            conv = np.divide(n0 - values, n0, out=np.full(values.shape, np.nan), where=n0 != 0)
            by[name] = self._series(conv, f"conversion_{name}", None)
            initial.append(n0)
            current.append(values)
        initial_arr = np.asarray(initial, dtype=float)
        current_arr = np.column_stack(current)
        denom = float(np.sum(initial_arr))
        total = np.divide(denom - np.sum(current_arr, axis=1), denom,
                          out=np.full(len(self.snapshots), np.nan), where=denom != 0)
        return ConversionSeries(by, self._series(total, "conversion_total", None))

    @property
    def conv(self) -> ConversionSeries:
        return self._conversion

    @property
    def free_monomer_composition(self) -> _NamedSeries:
        monomers = self.monomer_names
        matrix = np.column_stack([np.asarray(self.conc[name], dtype=float) for name in monomers])
        total = np.sum(matrix, axis=1)
        frac = np.divide(matrix, total[:, None], out=np.full(matrix.shape, np.nan), where=total[:, None] != 0)
        return _NamedSeries({name: self._series(frac[:, i], f"f_{name}", None) for i, name in enumerate(monomers)})

    @property
    def f(self) -> _NamedSeries:
        return self.free_monomer_composition

    @property
    def initial_monomer_composition(self) -> _NamedSeries:
        monomers = self.monomer_names
        initial = np.asarray([float(self.conc[name][0]) for name in monomers])
        total = float(np.sum(initial))
        frac = np.divide(initial, total, out=np.full(initial.shape, np.nan), where=total != 0)
        return _NamedSeries({name: self._series(np.full(len(self.snapshots), frac[i]), f"f0_{name}", None)
                             for i, name in enumerate(monomers)})

    @property
    def f0(self) -> _NamedSeries:
        return self.initial_monomer_composition

    def _chain_unit_totals(self) -> tuple[np.ndarray, np.ndarray]:
        monomers = self.monomer_names
        n_snap = len(self.snapshots)
        totals = np.full((n_snap, len(monomers)), np.nan, dtype=float)
        present = np.zeros(n_snap, dtype=bool)
        if "chains" not in self.tables or "chain_composition" not in self.tables:
            return totals, present
        chains = self.table("chains")
        comp = self.table("chain_composition")
        if chains.n_rows == 0:
            return totals, present
        rec_to_snap = np.asarray(chains["snapshot_id"], dtype=np.int64)
        abundance = np.asarray(chains["count"], dtype=float)
        rec_ids = np.asarray(comp["chain_record_id"], dtype=np.int64)
        mon_ids = np.asarray(comp["monomer_id"], dtype=np.int64)
        units = np.asarray(comp["unit_count"], dtype=float)
        sid_to_pos = {int(s): i for i, s in enumerate(np.asarray(self.snapshots.ids))}
        totals[:] = 0.0
        for rid, mid, unit_count in zip(rec_ids, mon_ids, units):
            if rid < 0 or rid >= len(rec_to_snap) or mid < 0 or mid >= len(monomers):
                raise InvalidOutputError("chain_composition contains an out-of-range identifier")
            pos = sid_to_pos.get(int(rec_to_snap[rid]))
            if pos is None:
                raise InvalidOutputError("chains references an unknown snapshot_id")
            totals[pos, mid] += abundance[rid] * unit_count
            present[pos] = True
        totals[~present, :] = np.nan
        return totals, present

    def _polymer_composition(self, *, kind: str) -> _NamedSeries:
        monomers = self.monomer_names
        if len(monomers) == 1:
            name = monomers[0]
            return _NamedSeries({name: self._series(np.ones(len(self.snapshots), dtype=float), f"F_{kind}_{name}", None)})
        cumulative_units, present = self._chain_unit_totals()
        if kind == "cumulative":
            sums = np.nansum(cumulative_units, axis=1)
            frac = np.divide(cumulative_units, sums[:, None], out=np.full(cumulative_units.shape, np.nan),
                             where=sums[:, None] != 0)
            defined = np.flatnonzero(np.all(np.isfinite(frac), axis=1))
            if defined.size:
                first = int(defined[0])
                frac[:first] = frac[first]
                previous = frac[first].copy()
                for row in range(first + 1, len(frac)):
                    if np.all(np.isfinite(frac[row])):
                        previous = frac[row].copy()
                    else:
                        frac[row] = previous
        elif kind == "interval":
            delta = np.full_like(cumulative_units, np.nan)
            previous = None
            for i in range(len(cumulative_units)):
                if not present[i]:
                    continue
                if previous is None:
                    delta[i] = cumulative_units[i]
                else:
                    delta[i] = cumulative_units[i] - cumulative_units[previous]
                previous = i
            sums = np.nansum(delta, axis=1)
            frac = np.divide(delta, sums[:, None], out=np.full(delta.shape, np.nan), where=sums[:, None] > 0)
        else:
            raise ValueError("kind must be 'interval' or 'cumulative'")
        return _NamedSeries({name: self._series(frac[:, i], f"F_{kind}_{name}", None)
                             for i, name in enumerate(monomers)})

    @property
    def polymer_composition(self) -> PolymerCompositionSeries:
        return PolymerCompositionSeries(self)

    @property
    def F(self) -> PolymerCompositionSeries:
        return self.polymer_composition

    @property
    def chains(self) -> StorageChains:
        if not self.has_chains or "chains" not in self.run.tables:
            raise DataUnavailableError(f"Snapshot {self.id} has no stored chains")
        return StorageChains(self.run, self._rows("chains"))

    @property
    def moments(self) -> StorageMomentsSnapshot:
        if "moments" not in self.run.tables:
            raise DataUnavailableError("moments table is unavailable")
        return StorageMomentsSnapshot(self)
    @property
    def dpn(self): return self.moments.default.dp_n
    @property
    def dpw(self): return self.moments.default.dp_w
    @property
    def mn(self): return self.moments.default.mn
    @property
    def mw(self): return self.moments.default.mw
    @property
    def mz(self): return self.moments.default.mz
    @property
    def dispersity(self): return self.moments.default.dispersity
    @analysis_operation(CLD_HELP)
    def cld(self, **kwargs): return self.chains.cld(**kwargs)
    @analysis_operation(MWD_HELP)
    def mwd(self, **kwargs): return self.chains.mwd(**kwargs)
    @analysis_operation(SPECTRUM_HELP)
    def chain_mass_spectrum(self, **kwargs): return self.chains.chain_mass_spectrum(**kwargs)

    @property
    def channels(self) -> StorageChannelsSnapshot:
        return StorageChannelsSnapshot(self)

    @property
    def kinetics(self) -> StorageKineticsSnapshot:
        if "kinetic_parameters/values" not in self.run.tables:
            raise DataUnavailableError("kinetic parameter values are unavailable")
        return StorageKineticsSnapshot(self)

    @property
    def temp(self): return self.kinetics.temperature
    @property
    def k(self): return self.kinetics.k

    @property
    def output_status(self) -> OutputStatus:
        available = tuple(sorted(self.tables))
        required = tuple(sorted(
            str(record["name"]) for record in self.schema_records
            if record.get("record_type") == "table" and bool(record.get("required", False))
        ))
        missing = tuple(name for name in required if name not in self.tables)
        invalid: list[str] = []
        if self.status == "completed" and not (self.path / "RESULTS_COMPLETE").is_file():
            invalid.append("RESULTS_COMPLETE")
        has_final = False
        try:
            has_final = bool(np.count_nonzero(np.asarray(self.snapshots.raw["is_final"], dtype=bool)) == 1)
        except Exception as exc:
            invalid.append(f"snapshots: {exc}")
        return OutputStatus(
            available=available,
            missing=missing,
            invalid=tuple(invalid),
            run_completed=self.status == "completed",
            has_final_snapshot=has_final,
        )

    def validate(self, *, strict: bool = False) -> ValidationReport:
        stored = self.diagnostics.validation
        failed = tuple(item.check for item in stored.failed)
        invalid = list(self.output_status.invalid)
        if stored.error_count:
            invalid.append(f"stored validation reports {stored.error_count} error(s)")
        warnings: list[str] = []
        if stored.warning_count:
            warnings.append(f"stored validation reports {stored.warning_count} warning(s)")
        if not stored.records and stored.status not in {"passed", "pass"}:
            warnings.append("stored validation was not run or has no check records")
        report = ValidationReport(
            is_valid=not invalid and stored.error_count == 0,
            is_complete=self.output_status.complete and stored.status in {"passed", "pass"},
            missing_outputs=self.output_status.missing,
            invalid_outputs=tuple(invalid),
            warnings=tuple(warnings),
            failed_checks=failed,
            details=tuple(dict(record) for record in stored.records),
        )
        if strict and (not report.is_valid or not report.is_complete):
            details = report.invalid_outputs or report.missing_outputs or report.warnings
            raise ValidationFailedError("; ".join(details) or "strict validation failed")
        return report

    @property
    def validation(self) -> ValidationReport:
        return self.validate()

    def refresh(self):
        fresh = open_storage(
            self.path,
            allow_incomplete=self.status != "completed",
            mmap_mode="r",
        )
        self._metadata = fresh._metadata
        self.run_id = fresh.run_id
        self._prefix = fresh._prefix
        self.schema_records = fresh.schema_records
        self.dictionaries = fresh.dictionaries
        self._tables = fresh._tables
        self.__dict__.setdefault("_cache", {}).clear()
        return self

    @property
    def monomers(self) -> Mapping[str, Mapping[str, Any]]:
        return {
            str(entry.get("name", idx)).removeprefix("monomer_"): dict(entry)
            for idx, entry in sorted(self.dictionaries.get("monomers", {}).items())
        }

    @property
    def endgroups(self) -> Mapping[str, Mapping[str, Any]]:
        return {
            str(entry.get("name", idx)): dict(entry)
            for idx, entry in sorted(self.dictionaries.get("chain_end_types", {}).items())
        }

    @property
    def monomer_names(self) -> tuple[str, ...]:
        return self.run.monomer_names

    @property
    def _conversion(self):
        source = self.run.conv
        names = source.names
        return _SnapshotConversion(
            names,
            [np.asarray(source[name])[self.index] for name in names],
            np.asarray(source.total)[self.index],
        )

    @property
    def conv(self): return self._conversion

    @property
    def f(self):
        source = self.run.f
        return _SnapshotNamedValues(source.names, [np.asarray(source[name])[self.index] for name in source.names])

    @property
    def f0(self):
        source = self.run.f0
        return _SnapshotNamedValues(source.names, [np.asarray(source[name])[self.index] for name in source.names])

    @property
    def F(self): return _SnapshotPolymerComposition(self)

    def info_text(self) -> str:
        lines = [
            "Snapshot",
            f"  run: {self.run._display_path()}",
            f"  snapshot_id: {self.id}",
            f"  time: {self.t:.6g}",
            f"  kmc_event: {self.event}",
            f"  final: {str(self.is_final).lower()}",
            f"  chain data: {'available' if self.has_chains else 'unavailable'}",
        ]
        try: lines.append(f"  conversion total: {self.conv.total:.6g}")
        except Exception: pass
        for label, getter in (("Mn", lambda: self.mn), ("Mw", lambda: self.mw), ("dispersity", lambda: self.dispersity)):
            try: lines.append(f"  {label}: {float(getter()):.6g}")
            except Exception: pass
        lines += [
            "", "Main objects:",
            "  snap.state", "  snap.chains", "  snap.moments", "  snap.channels", "  snap.kinetics",
            "", "Common next steps:",
        ]
        if self.monomer_names:
            lines.append(f'  snap.conc["{self.monomer_names[0]}"]')
        lines += ["  snap.chains.live", "  snap.chains.where(dp_min=10)", "  snap.mwd()", "  snap.cld()"]
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text


@dataclass(frozen=True)
class StorageSnapshots:
    run: "StorageRun"

    @property
    def raw(self) -> StorageTable:
        return self.run.table("snapshots")

    def __len__(self) -> int:
        return self.raw.n_rows

    def __iter__(self):
        for i in range(len(self)):
            yield StorageSnapshot(self.run, i)

    def __getitem__(self, snapshot_id: int) -> StorageSnapshot:
        ids = self.raw["snapshot_id"]
        pos = np.flatnonzero(ids == snapshot_id)
        if len(pos) == 0:
            raise KeyError(f"Unknown snapshot_id {snapshot_id}")
        return StorageSnapshot(self.run, int(pos[0]))

    def __getattr__(self, name: str):
        return getattr(self.raw, name)

    @property
    def ids(self) -> np.ndarray:
        return self.raw["snapshot_id"]

    @property
    def time(self) -> np.ndarray:
        return self.raw["time"]

    @property
    def kmc_event(self) -> np.ndarray:
        return self.raw["kmc_event"]

    @property
    def first(self) -> StorageSnapshot:
        if not len(self):
            raise SnapshotUnavailableError("Run has no snapshots")
        return StorageSnapshot(self.run, 0)

    @property
    def last(self) -> StorageSnapshot:
        if not len(self):
            raise SnapshotUnavailableError("Run has no snapshots")
        return StorageSnapshot(self.run, len(self) - 1)

    @property
    def final(self) -> StorageSnapshot:
        if "is_final" not in self.raw:
            raise FinalSnapshotUnavailableError("Snapshot finality is unavailable")
        pos = np.flatnonzero(self.raw["is_final"])
        if len(pos) == 0:
            raise FinalSnapshotUnavailableError("Run has no final snapshot")
        if len(pos) > 1:
            raise InvalidOutputError("Run contains more than one final snapshot")
        return StorageSnapshot(self.run, int(pos[0]))

    def _at_values(self, values: np.ndarray, value: float | int, method: str, label: str) -> StorageSnapshot:
        values = np.asarray(values)
        if len(values) == 0:
            raise SnapshotUnavailableError("Run has no snapshots")
        if method == "before":
            pos = np.flatnonzero(values <= value)
            if len(pos) == 0:
                raise DataUnavailableError(f"No snapshot at or before {label}={value}")
            idx = int(pos[-1])
        elif method == "after":
            pos = np.flatnonzero(values >= value)
            if len(pos) == 0:
                raise DataUnavailableError(f"No snapshot at or after {label}={value}")
            idx = int(pos[0])
        elif method == "nearest":
            idx = int(np.argmin(np.abs(values.astype(float) - float(value))))
        else:
            raise ValueError("method must be 'before', 'after', or 'nearest'")
        return StorageSnapshot(self.run, idx)

    def _at(self, column: str, value: float | int, method: str) -> StorageSnapshot:
        return self._at_values(np.asarray(self.raw[column]), value, method, column)

    def at_time(self, time: float, *, method: str = "before") -> StorageSnapshot:
        return self._at("time", time, method)

    def at_event(self, event: int, *, method: str = "before") -> StorageSnapshot:
        return self._at("kmc_event", event, method)




@dataclass(frozen=True)
class ValidationCheck:
    check: str
    status: str
    severity: str
    details: Mapping[str, Any]

    @property
    def passed(self) -> bool:
        return self.status == "pass"


class StorageValidation:
    def __init__(self, run: "StorageRun"):
        self.run = run
        self.path = run.path / "diagnostics" / "validation.jsonl"
        self._records = _read_jsonl(self.path) if self.path.is_file() and self.path.stat().st_size else []

    def __len__(self) -> int:
        return len(self._records)

    def __iter__(self):
        for record in self._records:
            yield ValidationCheck(
                check=str(record.get("check", record.get("name", ""))),
                status=str(record.get("status", "")),
                severity=str(record.get("severity", "")),
                details=record,
            )

    def __getitem__(self, key: int | str) -> ValidationCheck:
        if isinstance(key, int):
            return tuple(self)[key]
        for item in self:
            if item.check == key:
                return item
        raise KeyError(key)

    @property
    def status(self) -> str:
        return str(self.run._metadata.get("validation_status", "not_run"))

    @property
    def warning_count(self) -> int:
        return int(self.run._metadata.get("validation_warning_count", 0))

    @property
    def error_count(self) -> int:
        return int(self.run._metadata.get("validation_error_count", 0))

    @property
    def passed(self) -> bool:
        return self.status == "passed" and self.error_count == 0

    @property
    def failed(self) -> tuple[ValidationCheck, ...]:
        return tuple(item for item in self if item.status != "pass")

    @property
    def records(self) -> tuple[Mapping[str, Any], ...]:
        return tuple(self._records)


@dataclass(frozen=True)
class TextLog:
    path: Path

    @property
    def exists(self) -> bool:
        return self.path.is_file()

    @property
    def text(self) -> str:
        return self.path.read_text(encoding="utf-8") if self.exists else ""

    @property
    def lines(self) -> tuple[str, ...]:
        return tuple(self.text.splitlines())

    def __str__(self) -> str:
        return self.text


class StorageMemory:
    def __init__(self, run: "StorageRun"):
        self.run = run
        if "memory" not in run.tables:
            raise DataUnavailableError("memory diagnostics table is unavailable")
        self.raw = run.table("memory")

    @property
    def columns(self) -> tuple[str, ...]:
        return tuple(name for name in self.raw.columns if name != "snapshot_id")

    def __getitem__(self, name: str) -> np.ndarray:
        if name not in self.raw:
            raise KeyError(name)
        if "snapshot_id" not in self.raw:
            raise InvalidOutputError("memory table has no snapshot_id column")
        values = np.full(len(self.run.snapshots), np.nan, dtype=float)
        positions = {int(sid): i for i, sid in enumerate(np.asarray(self.run.snapshots.ids))}
        for sid, value in zip(self.raw["snapshot_id"], self.raw[name]):
            pos = positions.get(int(sid))
            if pos is None:
                raise InvalidOutputError("memory table references an unknown snapshot_id")
            values[pos] = value
        return self.run._series(values, name, self.run.column_unit("memory", name))

    def __getattr__(self, name: str) -> np.ndarray:
        if name in self.raw:
            return self[name]
        raise AttributeError(name)




class StorageChannelTrace:
    """Event-level SSA/KMC trace stored in diagnostics/channel_trace/*.npy."""
    def __init__(self, run: "StorageRun"):
        self.run = run

    @property
    def raw(self) -> StorageTable:
        if "diagnostics/channel_trace" not in self.run.tables:
            raise DataUnavailableError("channel trace is unavailable; run the engine with --trace-channels N")
        return self.run.table("diagnostics/channel_trace")

    def __len__(self) -> int:
        return self.raw.n_rows

    @property
    def kmc_event(self) -> np.ndarray: return _readonly_array(self.raw["kmc_event"])
    @property
    def t(self) -> np.ndarray: return _readonly_array(self.raw["time"])
    @property
    def dt(self) -> np.ndarray: return _readonly_array(self.raw["dt"])
    @property
    def channel_id(self) -> np.ndarray: return _readonly_array(self.raw["channel_id"])
    @property
    def rate(self) -> np.ndarray: return _readonly_array(self.raw["rate"])
    @property
    def propensity(self) -> np.ndarray: return _readonly_array(self.raw["propensity"])
    @property
    def total_propensity(self) -> np.ndarray: return _readonly_array(self.raw["total_propensity"])

    @property
    def channel(self) -> np.ndarray:
        entries = self.run.dictionaries.get("channels", {})
        values = [str(entries.get(int(i), {}).get("name", int(i))) for i in self.channel_id]
        return _readonly_array(np.asarray(values, dtype=object))

    @property
    def enabled(self) -> bool:
        return bool(self.run._metadata.get("channel_trace_enabled", "diagnostics/channel_trace" in self.run.tables))

    @property
    def complete(self) -> bool:
        return bool(self.run._metadata.get("channel_trace_complete", False))

    @property
    def truncated(self) -> bool:
        return bool(self.run._metadata.get("channel_trace_truncated", not self.complete))

    @property
    def limit(self) -> int | None:
        value = self.run._metadata.get("channel_trace_limit")
        return None if value is None else int(value)

    def by_channel(self, name: str) -> StorageTable:
        mask = np.asarray(self.channel, dtype=object) == str(name)
        if not np.any(mask):
            available = sorted(set(map(str, self.channel.tolist())))
            raise KeyError(f"unknown traced channel {name!r}; available: {available}")
        return self.raw.filtered(mask)

    def channel_counts(self) -> dict[str, int]:
        names, counts = np.unique(np.asarray(self.channel, dtype=object), return_counts=True)
        return {str(name): int(count) for name, count in zip(names, counts)}

    def info_text(self) -> str:
        return "\n".join([
            "ChannelTraceView",
            f"  enabled: {self.enabled}",
            f"  rows: {len(self) if 'diagnostics/channel_trace' in self.run.tables else 0}",
            f"  complete: {self.complete}",
            f"  truncated: {self.truncated}",
            "",
            "Common next steps:",
            "  trace.kmc_event",
            "  trace.channel",
            "  trace.by_channel(\"prop_A_A\")",
            "  trace.channel_counts()",
        ])

    def info(self) -> str:
        text = self.info_text(); print(text); return text


class StorageDiagnostics:
    def __init__(self, run: "StorageRun"):
        self.run = run

    @property
    def validation(self) -> StorageValidation:
        return StorageValidation(self.run)

    @property
    def memory(self) -> StorageMemory:
        return StorageMemory(self.run)

    @property
    def run_log(self) -> TextLog:
        return TextLog(self.run.path / "diagnostics" / "run.log")

    @property
    def debug_log(self) -> TextLog:
        return TextLog(self.run.path / "diagnostics" / "debug.log")

    @property
    def channel_trace(self) -> StorageChannelTrace:
        return StorageChannelTrace(self.run)

    def info_text(self) -> str:
        v=self.validation
        lines=["DiagnosticsView", f"  validation: {v.status}", f"  checks: {len(v)}", f"  warnings: {v.warning_count}", f"  errors: {v.error_count}",
               f"  memory records: {self.memory.raw.n_rows if 'memory' in self.run.tables else 0}", f"  run.log: {'available' if self.run_log.exists else 'missing'}", f"  debug.log: {'available' if self.debug_log.exists else 'missing'}",
               "", "Common next steps:", "  run.diagnostics.validation", "  run.diagnostics.memory", "  run.diagnostics.channel_trace", "  run.diagnostics.run_log.text"]
        return "\n".join(lines)

    def info(self) -> str:
        text=self.info_text(); print(text); return text


class StorageRaw:
    def __init__(self, run: "StorageRun"):
        self.run = run

    @property
    def metadata(self) -> Mapping[str, Any]:
        return self.run._metadata

    @property
    def schema(self) -> tuple[Mapping[str, Any], ...]:
        return tuple(self.run.schema_records)

    @property
    def tables(self) -> Mapping[str, StorageTable]:
        return self.run.tables

    @property
    def dictionaries(self) -> Mapping[str, Mapping[int, Mapping[str, Any]]]:
        return self.run.dictionaries

    def table(self, name: str) -> StorageTable:
        return self.run.table(name)

    def dictionary(self, name: str) -> Mapping[int, Mapping[str, Any]]:
        return self.run.dictionary(name)


    def info_text(self) -> str:
        tables=tuple(sorted(self.tables)); dictionaries=tuple(sorted(self.dictionaries))
        lines=["RawView", f"  tables: {len(tables)}", "  " + ", ".join(tables), f"  dictionaries: {len(dictionaries)}", "  " + ", ".join(dictionaries),
               "", "Common next steps:", "  run.raw.metadata", "  run.raw.schema", "  run.raw.table(\"chains\")", "  run.raw.dictionary(\"monomers\")"]
        return "\n".join(lines)

    def info(self) -> str:
        text=self.info_text(); print(text); return text



class StorageFeedEvents:
    def __len__(self) -> int:
        return self.raw.n_rows

    @property
    def n_events(self) -> int:
        return len(self)

    def __init__(self, run: "StorageRun", feed_id: int | None = None):
        self.run = run
        self.feed_id = None if feed_id is None else int(feed_id)

    @property
    def raw(self) -> StorageTable:
        table = self.run.table("feed_events")
        if self.feed_id is None:
            return table
        return table.filtered(np.asarray(table["feed_id"], dtype=np.int64) == self.feed_id)

    @property
    def time(self) -> np.ndarray:
        actions = self.run.table("actions")
        ids = np.asarray(self.raw["action_id"], dtype=np.int64)
        return _readonly_array(np.asarray(actions["time"], dtype=float)[ids])

    @property
    def dose(self) -> np.ndarray:
        return _readonly_array(np.asarray(self.raw["dose_mL"], dtype=float) / 1000.0)

    @property
    def dose_mL(self) -> np.ndarray:
        return _readonly_array(self.raw["dose_mL"])

    @property
    def cumulative_volume(self) -> np.ndarray:
        return _readonly_array(np.cumsum(np.asarray(self.dose, dtype=float)))

    @property
    def cumulative_amount(self) -> _NamedSeries:
        if self.feed_id is None:
            feed_ids = np.asarray(self.raw["feed_id"], dtype=np.int64)
            values: dict[str, np.ndarray] = {}
            for name in self.run.state.names:
                increments = np.zeros(len(feed_ids), dtype=float)
                for i, fid in enumerate(feed_ids):
                    feed = next((f for f in self.run.feeds.values() if f.id == int(fid)), None)
                    if feed is not None:
                        increments[i] = float(feed.concentration.get(name, 0.0)) * float(self.dose[i])
                values[name] = _readonly_array(np.cumsum(increments))
            return _NamedSeries(values)
        feed = next(f for f in self.run.feeds.values() if f.id == self.feed_id)
        return _NamedSeries({
            name: _readonly_array(np.cumsum(np.asarray(self.dose, dtype=float) * float(concentration)))
            for name, concentration in feed.concentration.items()
        })

    @property
    def volume_before(self) -> np.ndarray:
        return _readonly_array(np.asarray(self.raw["volume_before_mL"], dtype=float) / 1000.0)

    @property
    def volume_after(self) -> np.ndarray:
        return _readonly_array(np.asarray(self.raw["volume_after_mL"], dtype=float) / 1000.0)

    @property
    def kmc_volume_before(self) -> np.ndarray:
        return _readonly_array(self.raw["kmc_volume_before_L"])

    @property
    def kmc_volume_after(self) -> np.ndarray:
        return _readonly_array(self.raw["kmc_volume_after_L"])


class StorageFeed:
    def __init__(self, run: "StorageRun", feed_id: int, record: Mapping[str, Any]):
        self.run = run
        self.id = int(feed_id)
        self.name = str(record.get("name", feed_id))
        self._record = dict(record)

    @property
    def concentration(self) -> Mapping[str, float]:
        return {str(x["name"]): float(x["concentration_mol_L"]) for x in self._record.get("components", [])}

    @property
    def fraction(self) -> Mapping[str, float]:
        monomers = set(self.run.monomer_names)
        vals = {k: v for k, v in self.concentration.items() if k in monomers}
        total = sum(vals.values())
        return {k: (v / total if total else float("nan")) for k, v in vals.items()}

    @property
    def events(self) -> StorageFeedEvents:
        return StorageFeedEvents(self.run, self.id)

    @property
    def volume_cum(self) -> float:
        return float(np.sum(self.events.dose))

    @property
    def moles_cum(self) -> Mapping[str, float]:
        return {name: concentration * self.volume_cum for name, concentration in self.concentration.items()}


class StorageFeeds(Mapping[str, StorageFeed]):
    def __init__(self, run: "StorageRun"):
        self.run = run
        records = run._metadata.get("model", {}).get("feeds", []) or []
        self._values = {str(rec.get("name", i)): StorageFeed(run, i, rec) for i, rec in enumerate(records)}
    def __getitem__(self, name: str) -> StorageFeed: return self._values[name]
    def __iter__(self): return iter(self._values)
    def __len__(self): return len(self._values)
    @property
    def names(self) -> tuple[str, ...]: return tuple(self._values)
    def _ipython_key_completions_(self): return list(self._values)


class StorageChainCountSeries:
    def __init__(self, run: "StorageRun"):
        self.run = run

    @property
    def live(self) -> np.ndarray:
        return _readonly_array(self.run.snapshots.raw["chain_count_live"])

    @property
    def dead(self) -> np.ndarray:
        return _readonly_array(self.run.snapshots.raw["chain_count_dead"])

    @property
    def total(self) -> np.ndarray:
        return _readonly_array(self.run.snapshots.raw["chain_count_total"])


class StorageBalance:
    """Read-only physical material balance in mol for monomers and free species."""
    def __init__(self, run: "StorageRun"):
        self.run = run

    @property
    def names(self) -> tuple[str, ...]:
        entries = self.run.dictionaries.get("balance_entities", {})
        return tuple(str(entries[i].get("name", i)) for i in sorted(entries))

    def _matrix(self, column: str) -> _NamedSeries:
        raw = self.run._metadata.get("initial_volume_mL")
        if raw is None or float(raw) <= 0:
            raise AnalysisNotApplicableError("physical material balance requires init_volume")
        if "species_balance" not in self.run.tables:
            # Backward-compatible read of Stage 23/25 monomer-only balance.
            if "monomer_balance" not in self.run.tables:
                raise AnalysisNotApplicableError("material balance is not available in this storage run")
            table = self.run.table("monomer_balance")
            names = self.run.monomer_names
            legacy = {
                "initial_moles": "initial_moles",
                "total_moles": "introduced_moles",
                "free_moles": "free_moles",
                "consumed_moles": None,
                "dosed_moles": None,
            }
            source = legacy[column]
            n_snap = len(self.run.snapshots)
            if source is None:
                initial = np.asarray(table["initial_moles"], dtype=float)
                total = np.asarray(table["introduced_moles"], dtype=float)
                free = np.asarray(table["free_moles"], dtype=float)
                values = total - free if column == "consumed_moles" else total - initial
            else:
                values = np.asarray(table[source], dtype=float)
        else:
            table = self.run.table("species_balance")
            names = self.names
            n_snap = len(self.run.snapshots)
            values = np.asarray(table[column], dtype=float)
        if len(values) != n_snap * len(names):
            raise InvalidOutputError(f"species_balance/{column} has invalid row count")
        matrix = values.reshape(n_snap, len(names))
        return _NamedSeries({name: _readonly_array(matrix[:, i]) for i, name in enumerate(names)})

    def _assert_valid(self, name: str) -> None:
        invalid = self.run._balance_invalidations()
        if name in invalid:
            t, line = invalid[name]
            raise AnalysisNotApplicableError(
                f"material balance for {name!r} is invalid because set_c was executed at t={t:g} s"
                + (f" (model line {line})" if line else "")
            )

    class _Checked(Mapping[str, np.ndarray]):
        def __init__(self, owner: "StorageBalance", values: _NamedSeries): self.owner, self.values = owner, values
        def __getitem__(self, name): self.owner._assert_valid(name); return self.values[name]
        def __iter__(self): return iter(self.values)
        def __len__(self): return len(self.values)
        @property
        def names(self): return self.values.names
        def _ipython_key_completions_(self): return list(self.values.names)

    def _checked(self, column: str): return self._Checked(self, self._matrix(column))
    @property
    def initial(self): return self._checked("initial_moles")
    @property
    def dosed(self): return self._checked("dosed_moles")
    @property
    def total(self): return self._checked("total_moles")
    @property
    def free(self): return self._checked("free_moles")
    @property
    def consumed(self): return self._checked("consumed_moles")
    @property
    def incorporated(self):
        if "monomer_balance" not in self.run.tables:
            raise AnalysisNotApplicableError("incorporated balance is available only for monomers")
        table = self.run.table("monomer_balance")
        names = self.run.monomer_names
        values = np.asarray(table["incorporated_moles"], dtype=float)
        n_snap = len(self.run.snapshots)
        if len(values) != n_snap * len(names):
            raise InvalidOutputError("monomer_balance/incorporated_moles has invalid row count")
        matrix = values.reshape(n_snap, len(names))
        vals = _NamedSeries({name: _readonly_array(matrix[:, i]) for i, name in enumerate(names)})
        return self._Checked(self, vals)


@dataclass
class StorageRun(Run):
    """Read-only, memory-mapped view of one canonical Slimmc Storage 1.2.0 run."""

    schema_records: list[dict] = None  # type: ignore[assignment]
    dictionaries: dict[str, dict[int, dict]] = None  # type: ignore[assignment]
    _tables: dict[str, StorageTable] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        super().__post_init__()
        self.schema_records = [] if self.schema_records is None else self.schema_records
        self.dictionaries = {} if self.dictionaries is None else self.dictionaries
        self._tables = {} if self._tables is None else self._tables

    def _input_model_lines(self) -> list[str]:
        path = self.path / "input.model"
        return path.read_text(encoding="utf-8").splitlines() if path.is_file() else []

    @property
    def desc(self) -> str | None:
        value = self._metadata.get("description")
        if value:
            return str(value)
        for line in self._input_model_lines():
            text = line.strip()
            if text.startswith("desc "):
                return text[5:].strip().strip('"')
        return None

    @property
    def var(self) -> Variables:
        records = self._metadata.get("variables", [])
        if records is None:
            records = []
        if not isinstance(records, list):
            raise DataConsistencyError("run_metadata.json field 'variables' must be an array")
        return Variables(records)

    def _display_path(self) -> str:
        rel = getattr(self, "relative_dir", ".")
        if rel and rel != ".":
            return rel.rstrip("/") + "/"
        return self.path.name.rstrip("/") + "/"

    def info_text(self) -> str:
        lines = [
            "Run",
            f"  path: {self._display_path()}",
            f"  pyslimmc: {__import__('pyslimmc').__version__}",
            f"  engine: {self.engine} {self.engine_version}".rstrip(),
            f"  storage: {self.storage_format} {self.storage_format_version}".rstrip(),
            f"  status: {self.status}",
            f"  validation: {'PASS' if self.is_ok else ('FAIL' if self.status == 'completed' else 'partial')}",
            f"  monomers: {', '.join(self.monomer_names) or '-'}",
            f"  snapshots: {len(self.snapshots)}",
            f"  variables: {len(self.var)}",
        ]
        try:
            lines.append(f"  chain snapshots: {len(set(np.asarray(self.table('chains')['snapshot_id'], dtype=int)))}")
        except Exception:
            lines.append("  chain snapshots: 0")
        lines.append("")
        lines.append("Final state:")
        try: lines.append(f"  time: {float(self.t[-1]):.6g} s")
        except Exception: pass
        try: lines.append(f"  conversion total: {float(self.conv.total[-1]):.6g}")
        except Exception: pass
        for label, getter in (("Mn", lambda: self.mn[-1]), ("Mw", lambda: self.mw[-1]), ("dispersity", lambda: self.dispersity[-1])):
            try: lines.append(f"  {label}: {float(getter()):.6g}")
            except Exception: pass
        if self.var:
            lines += ["", "Variables:"]
            for variable in self.var.values():
                lines.append(f"  {variable.name}: {variable.value:.8g} {variable.unit} ({variable.kind})")
        lines += [
            "", "Main objects:",
            "  run.state", "  run.chains", "  run.moments", "  run.channels",
            "  run.kinetics", "  run.actions", "  run.diagnostics",
            "", "Common next steps:",
            "  run.last", "  run.final", "  run.conv.total",
        ]
        if self.monomer_names:
            lines.append(f'  run.conc["{self.monomer_names[0]}"]')
            if len(self.monomer_names) > 1:
                lines.append(f'  run.F.cum["{self.monomer_names[0]}"]')
        lines += ["  run.last.chains.dead", "  run.last.mwd()"]
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text

    def summary(self, path: str | Path | None = None):
        """Return a compact scientific/technical summary of this run.

        When *path* is supplied, write TXT or JSON according to the suffix.
        No validation is rerun; the stored validation result is reported.
        """
        from .summary import build_summary
        result = build_summary(self)
        if path is not None:
            result.write(path)
        return result

    @property
    def raw(self) -> StorageRaw:
        return StorageRaw(self)

    @property
    def diagnostics(self) -> StorageDiagnostics:
        return StorageDiagnostics(self)

    @property
    def status(self) -> str:
        return str(self._metadata.get("run_status", ""))

    @property
    def version(self) -> str:
        """Producer engine version (not the pyslimmc or Storage version)."""
        return self.engine_version

    @property
    def engine_version(self) -> str:
        return str(self._metadata.get("engine_version", ""))

    @property
    def cli_version(self) -> str:
        return str(self._metadata.get("cli_version", ""))

    @property
    def storage_format(self) -> str:
        return str(self._metadata.get("storage", "slimmc-storage"))

    @property
    def storage_format_version(self) -> str:
        return str(self._metadata.get("storage_format_version", ""))

    @property
    def schema(self) -> str:
        return self.storage_format

    @property
    def is_complete(self) -> bool:
        return self.status == "completed" and (self.path / "RESULTS_COMPLETE").is_file()

    @property
    def is_ok(self) -> bool:
        return self.is_complete and self._metadata.get("validation_error_count", 0) == 0

    @property
    def engine(self) -> str:
        return str(self._metadata.get("engine", ""))

    @property
    def kinetic_model(self) -> str:
        return str(self._metadata.get("kinetic_model", ""))

    @property
    def tables(self) -> Mapping[str, StorageTable]:
        return self._tables

    def table(self, name: str) -> StorageTable:
        try:
            return self._tables[name]
        except KeyError as exc:
            raise KeyError(f"Unknown Slimmc Storage table {name!r}. Available: {', '.join(sorted(self._tables))}") from exc

    @property
    def snapshots(self) -> StorageSnapshots:
        return StorageSnapshots(self)

    @property
    def first(self) -> StorageSnapshot:
        return self.snapshots.first

    @property
    def last(self) -> StorageSnapshot:
        return self.snapshots.last

    @property
    def final(self) -> StorageSnapshot:
        return self.snapshots.final

    @property
    def snapshots_with_chains(self) -> tuple[StorageSnapshot, ...]:
        """All snapshots whose chain table was explicitly saved."""
        return tuple(snap for snap in self.snapshots if snap.has_chains)

    @property
    def first_with_chains(self) -> StorageSnapshot:
        snapshots = self.snapshots_with_chains
        if not snapshots:
            raise DataUnavailableError("Run has no snapshot with chains")
        return snapshots[0]

    @property
    def last_with_chains(self) -> StorageSnapshot:
        snapshots = self.snapshots_with_chains
        if not snapshots:
            raise DataUnavailableError("Run has no snapshot with chains")
        return snapshots[-1]

    def at_snapshot(self, snapshot_id: int) -> StorageSnapshot:
        return self.snapshots[snapshot_id]

    def at_time(self, time: float, *, method: str = "before") -> StorageSnapshot:
        return self.snapshots.at_time(time, method=method)

    def at_event(self, event: int, *, method: str = "before") -> StorageSnapshot:
        return self.snapshots.at_event(event, method=method)

    def at_conversion(self, conversion: float, *, monomer: str | None = None, method: str = "before") -> StorageSnapshot:
        values = np.asarray(self.conv.total if monomer is None else self.conv[monomer], dtype=float)
        return self.snapshots._at_values(values, conversion, method, "conversion")

    def at_temperature(self, temperature: float) -> tuple[StorageSnapshot, ...]:
        values = np.asarray(self.temp, dtype=float)
        if len(values) == 0:
            raise SnapshotUnavailableError("Run has no snapshots")
        distance = np.abs(values - float(temperature))
        minimum = np.nanmin(distance)
        return tuple(StorageSnapshot(self, int(i)) for i in np.flatnonzero(np.isclose(distance, minimum, rtol=0.0, atol=0.0)))

    @property
    def state(self) -> StorageStateSeries:
        return StorageStateSeries(self)

    @property
    def t(self) -> np.ndarray:
        return _readonly_array(self.snapshots.time)

    @property
    def event(self) -> np.ndarray:
        return _readonly_array(self.snapshots.kmc_event)

    @property
    def sid(self) -> np.ndarray:
        return _readonly_array(self.snapshots.ids)

    @property
    def kmc_volume(self) -> np.ndarray:
        return _readonly_array(self.snapshots.raw["kmc_volume_L"])

    @property
    def volume(self) -> np.ndarray:
        if "volume_mL" not in self.snapshots.raw:
            raise AnalysisNotApplicableError("physical volume is unavailable because init_volume was not defined")
        values = np.asarray(self.snapshots.raw["volume_mL"], dtype=float)
        if not np.any(values > 0):
            raise AnalysisNotApplicableError("physical volume is unavailable because init_volume was not defined")
        return _readonly_array(values / 1000.0)

    def _initial_concentrations(self) -> dict[str, float]:
        model = self._metadata.get("model", {}) or {}
        result = {}
        for key in ("monomers", "species"):
            for rec in model.get(key, []) or []:
                result[str(rec["name"])] = float(rec.get("initial_concentration_mol_L", 0.0))
        return result

    @property
    def c0(self) -> Mapping[str, float]:
        return self._initial_concentrations()

    @property
    def count0(self) -> Mapping[str, int]:
        v0=float(self._metadata.get("initial_kmc_volume_L", self._metadata.get("kmc_volume_L", 0.0)))
        return {name:int(round(c*AVOGADRO*v0)) for name,c in self.c0.items()}

    @property
    def moles0(self) -> Mapping[str, float]:
        raw = self._metadata.get("initial_volume_mL")
        if raw is None or float(raw) <= 0:
            raise AnalysisNotApplicableError("physical initial moles require init_volume")
        v0 = float(raw) / 1000.0
        return {name: c * v0 for name, c in self.c0.items()}

    @property
    def feeds(self) -> StorageFeeds:
        return StorageFeeds(self)

    @property
    def feed_events(self) -> StorageFeedEvents:
        return StorageFeedEvents(self)

    @property
    def chain_count(self) -> StorageChainCountSeries:
        return StorageChainCountSeries(self)

    def _balance_invalidations(self) -> dict[str, tuple[float, int]]:
        if "actions" not in self.tables: return {}
        action_types=self.dictionaries.get("action_types", {})
        set_ids={int(i) for i,r in action_types.items() if r.get("name")=="set_c"}
        targets=self.dictionaries.get("action_targets", {})
        a=self.table("actions"); out={}
        for i,typ in enumerate(np.asarray(a["action_type_id"],dtype=int)):
            if typ not in set_ids: continue
            target=targets.get(int(a["target_id"][i]),{})
            name=str(target.get("name", ""))
            if name: out.setdefault(name,(float(a["time"][i]),int(a["source_line"][i])))
        return out

    @property
    def balance(self) -> StorageBalance:
        return StorageBalance(self)

    @property
    def count(self) -> _StateFieldSeries:
        return self.state.counts

    @property
    def moles(self) -> _NamedSeries:
        volume = np.asarray(self.volume, dtype=float)
        return _NamedSeries({name: _readonly_array(np.asarray(self.conc[name], dtype=float) * volume) for name in self.state.names})

    @property
    def conc(self) -> _StateFieldSeries:
        return self.state.concentrations

    @property
    def output_status(self) -> OutputStatus:
        available = tuple(sorted(self.tables))
        required = tuple(sorted(
            str(record["name"]) for record in self.schema_records
            if record.get("record_type") == "table" and bool(record.get("required", False))
        ))
        missing = tuple(name for name in required if name not in self.tables)
        invalid: list[str] = []
        if self.status == "completed" and not (self.path / "RESULTS_COMPLETE").is_file():
            invalid.append("RESULTS_COMPLETE")
        has_final = False
        try:
            has_final = bool(np.count_nonzero(np.asarray(self.snapshots.raw["is_final"], dtype=bool)) == 1)
        except Exception as exc:
            invalid.append(f"snapshots: {exc}")
        return OutputStatus(
            available=available,
            missing=missing,
            invalid=tuple(invalid),
            run_completed=self.status == "completed",
            has_final_snapshot=has_final,
        )

    def validate(self, *, strict: bool = False) -> ValidationReport:
        stored = self.diagnostics.validation
        failed = tuple(item.check for item in stored.failed)
        invalid = list(self.output_status.invalid)
        if stored.error_count:
            invalid.append(f"stored validation reports {stored.error_count} error(s)")
        warnings: list[str] = []
        if stored.warning_count:
            warnings.append(f"stored validation reports {stored.warning_count} warning(s)")
        if not stored.records and stored.status not in {"passed", "pass"}:
            warnings.append("stored validation was not run or has no check records")
        report = ValidationReport(
            is_valid=not invalid and stored.error_count == 0,
            is_complete=self.output_status.complete and stored.status in {"passed", "pass"},
            missing_outputs=self.output_status.missing,
            invalid_outputs=tuple(invalid),
            warnings=tuple(warnings),
            failed_checks=failed,
            details=tuple(dict(record) for record in stored.records),
        )
        if strict and (not report.is_valid or not report.is_complete):
            details = report.invalid_outputs or report.missing_outputs or report.warnings
            raise ValidationFailedError("; ".join(details) or "strict validation failed")
        return report

    @property
    def validation(self) -> ValidationReport:
        return self.validate()

    def refresh(self):
        fresh = open_storage(
            self.path,
            allow_incomplete=self.status != "completed",
            mmap_mode="r",
        )
        self._metadata = fresh._metadata
        self.run_id = fresh.run_id
        self._prefix = fresh._prefix
        self.schema_records = fresh.schema_records
        self.dictionaries = fresh.dictionaries
        self._tables = fresh._tables
        self.__dict__.setdefault("_cache", {}).clear()
        return self

    @property
    def monomers(self) -> Mapping[str, Mapping[str, Any]]:
        return {
            str(entry.get("name", idx)).removeprefix("monomer_"): dict(entry)
            for idx, entry in sorted(self.dictionaries.get("monomers", {}).items())
        }

    @property
    def endgroups(self) -> Mapping[str, Mapping[str, Any]]:
        return {
            str(entry.get("name", idx)): dict(entry)
            for idx, entry in sorted(self.dictionaries.get("chain_end_types", {}).items())
        }

    @property
    def monomer_names(self) -> tuple[str, ...]:
        def public_name(value: object) -> str:
            name = str(value)
            return name[len("monomer_"):] if name.startswith("monomer_") else name

        entries = self.dictionaries.get("monomers", {})
        if entries:
            return tuple(public_name(entries[i].get("name", i)) for i in sorted(entries))
        state_entries = self.dictionaries.get("state_entities", {})
        names = [public_name(entry.get("name", i)) for i, entry in sorted(state_entries.items())
                 if entry.get("kind") == "monomer"]
        return tuple(names)

    def _series(self, values: np.ndarray, name: str, unit: str | None = None) -> np.ndarray:
        # name/unit remain available from the owning facade and schema; numeric
        # leaves deliberately terminate in a plain read-only ndarray.
        return _readonly_array(values)

    def _initial_monomer_counts(self) -> dict[str, float] | None:
        """Return model-declared initial monomer counts, when reconstructible.

        Canonical model files declare ``param kmc_volume`` in litres and each
        ``monomer`` concentration in mol/L.  Using these declarations avoids
        treating the first stored snapshot as the initial state when no save
        was requested at t=0.
        """
        volume = None
        concentrations: dict[str, float] = {}
        for raw in self._input_model_lines():
            parts = raw.split("#", 1)[0].split()
            if not parts:
                continue
            if len(parts) >= 3 and parts[0] == "param" and parts[1] == "kmc_volume":
                try:
                    volume = float(parts[2])
                except ValueError:
                    return None
            elif len(parts) >= 3 and parts[0] == "monomer":
                try:
                    concentrations[parts[1]] = float(parts[2])
                except ValueError:
                    return None
        if volume is None or not math.isfinite(volume) or volume <= 0:
            return None
        if any(name not in concentrations for name in self.monomer_names):
            return None
        return {
            name: float(round(concentrations[name] * AVOGADRO * volume))
            for name in self.monomer_names
        }

    @property
    def _conversion(self) -> ConversionSeries:
        monomers = self.monomer_names
        if not monomers:
            raise DataUnavailableError("No monomer entities are declared in schema.jsonl")
        if "monomer_balance" in self.tables:
            table = self.table("monomer_balance")
            values = np.asarray(table["conversion"], dtype=float)
            expected = len(self.snapshots) * len(monomers)
            if len(values) == expected:
                matrix = values.reshape(len(self.snapshots), len(monomers))
                by = {name: self._series(matrix[:, i], f"conversion_{name}", None) for i, name in enumerate(monomers)}
                introduced = np.asarray(table["introduced_moles"], dtype=float).reshape(len(self.snapshots), len(monomers))
                free = np.asarray(table["free_moles"], dtype=float).reshape(len(self.snapshots), len(monomers))
                denom = np.sum(introduced, axis=1)
                total = np.divide(np.sum(introduced-free, axis=1), denom, out=np.full(len(denom), np.nan), where=denom != 0)
                return ConversionSeries(by, self._series(total, "conversion_total", None))
        by: dict[str, np.ndarray] = {}
        initial = []
        current = []
        declared = self._initial_monomer_counts()
        for name in monomers:
            values = np.asarray(self.count[name], dtype=float)
            n0 = float(declared[name]) if declared is not None else float(values[0])
            conv = np.divide(n0 - values, n0, out=np.full(values.shape, np.nan), where=n0 != 0)
            by[name] = self._series(conv, f"conversion_{name}", None)
            initial.append(n0)
            current.append(values)
        initial_arr = np.asarray(initial, dtype=float)
        current_arr = np.column_stack(current)
        denom = float(np.sum(initial_arr))
        total = np.divide(denom - np.sum(current_arr, axis=1), denom,
                          out=np.full(len(self.snapshots), np.nan), where=denom != 0)
        return ConversionSeries(by, self._series(total, "conversion_total", None))

    @property
    def conv(self) -> ConversionSeries:
        return self._conversion

    @property
    def free_monomer_composition(self) -> _NamedSeries:
        monomers = self.monomer_names
        matrix = np.column_stack([np.asarray(self.conc[name], dtype=float) for name in monomers])
        total = np.sum(matrix, axis=1)
        frac = np.divide(matrix, total[:, None], out=np.full(matrix.shape, np.nan), where=total[:, None] != 0)
        return _NamedSeries({name: self._series(frac[:, i], f"f_{name}", None) for i, name in enumerate(monomers)})

    @property
    def f(self) -> _NamedSeries:
        return self.free_monomer_composition

    @property
    def initial_monomer_composition(self) -> _NamedSeries:
        monomers = self.monomer_names
        declared = self._initial_monomer_counts()
        initial = np.asarray(
            [float(declared[name]) for name in monomers]
            if declared is not None
            else [float(self.conc[name][0]) for name in monomers]
        )
        total = float(np.sum(initial))
        frac = np.divide(initial, total, out=np.full(initial.shape, np.nan), where=total != 0)
        return _NamedSeries({name: self._series(np.full(len(self.snapshots), frac[i]), f"f0_{name}", None)
                             for i, name in enumerate(monomers)})

    @property
    def f0(self) -> _NamedSeries:
        return self.initial_monomer_composition

    def _chain_unit_totals(self) -> tuple[np.ndarray, np.ndarray]:
        monomers = self.monomer_names
        n_snap = len(self.snapshots)
        totals = np.full((n_snap, len(monomers)), np.nan, dtype=float)
        present = np.zeros(n_snap, dtype=bool)
        if "chains" not in self.tables or "chain_composition" not in self.tables:
            return totals, present
        chains = self.table("chains")
        comp = self.table("chain_composition")
        if chains.n_rows == 0:
            return totals, present
        rec_to_snap = np.asarray(chains["snapshot_id"], dtype=np.int64)
        abundance = np.asarray(chains["count"], dtype=float)
        rec_ids = np.asarray(comp["chain_record_id"], dtype=np.int64)
        mon_ids = np.asarray(comp["monomer_id"], dtype=np.int64)
        units = np.asarray(comp["unit_count"], dtype=float)
        sid_to_pos = {int(s): i for i, s in enumerate(np.asarray(self.snapshots.ids))}
        totals[:] = 0.0
        for rid, mid, unit_count in zip(rec_ids, mon_ids, units):
            if rid < 0 or rid >= len(rec_to_snap) or mid < 0 or mid >= len(monomers):
                raise InvalidOutputError("chain_composition contains an out-of-range identifier")
            pos = sid_to_pos.get(int(rec_to_snap[rid]))
            if pos is None:
                raise InvalidOutputError("chains references an unknown snapshot_id")
            totals[pos, mid] += abundance[rid] * unit_count
            present[pos] = True
        totals[~present, :] = np.nan
        return totals, present

    def _polymer_composition(self, *, kind: str) -> _NamedSeries:
        monomers = self.monomer_names
        if len(monomers) == 1:
            name = monomers[0]
            return _NamedSeries({name: self._series(np.ones(len(self.snapshots), dtype=float), f"F_{kind}_{name}", None)})
        declared = self._initial_monomer_counts()
        if declared is not None:
            initial = np.asarray([declared[name] for name in monomers], dtype=float)
            current = np.column_stack([np.asarray(self.count[name], dtype=float) for name in monomers])
            cumulative_units = initial[None, :] - current
            # Tiny negative values can only arise from representation/rounding.
            cumulative_units[np.abs(cumulative_units) < 0.5] = 0.0
            present = np.ones(len(self.snapshots), dtype=bool)
        else:
            cumulative_units, present = self._chain_unit_totals()
        if kind == "cumulative":
            sums = np.nansum(cumulative_units, axis=1)
            frac = np.divide(cumulative_units, sums[:, None], out=np.full(cumulative_units.shape, np.nan),
                             where=sums[:, None] > 0)
            defined = np.flatnonzero(np.all(np.isfinite(frac), axis=1))
            if defined.size:
                first = int(defined[0])
                frac[:first] = frac[first]
                previous = frac[first].copy()
                for row in range(first + 1, len(frac)):
                    if np.all(np.isfinite(frac[row])):
                        previous = frac[row].copy()
                    else:
                        frac[row] = previous
        elif kind == "interval":
            delta = np.full_like(cumulative_units, np.nan)
            previous = None
            for i in range(len(cumulative_units)):
                if not present[i]:
                    continue
                if previous is None:
                    delta[i] = cumulative_units[i]
                else:
                    delta[i] = cumulative_units[i] - cumulative_units[previous]
                previous = i
            sums = np.nansum(delta, axis=1)
            frac = np.divide(delta, sums[:, None], out=np.full(delta.shape, np.nan), where=sums[:, None] > 0)
        else:
            raise ValueError("kind must be 'interval' or 'cumulative'")
        return _NamedSeries({name: self._series(frac[:, i], f"F_{kind}_{name}", None)
                             for i, name in enumerate(monomers)})

    @property
    def polymer_composition(self) -> PolymerCompositionSeries:
        return PolymerCompositionSeries(self)

    @property
    def F(self) -> PolymerCompositionSeries:
        return self.polymer_composition

    @property
    def chains(self) -> StorageChains:
        return StorageChains(self)

    @property
    def plot(self):
        from .composition_analysis import RunPlotNamespace
        cached = self.__dict__.get("_plot_namespace")
        if cached is None:
            cached = RunPlotNamespace(self)
            self.__dict__["_plot_namespace"] = cached
        return cached

    @property
    def moments(self) -> StorageMomentsSeries:
        if "moments" not in self.tables:
            raise DataUnavailableError("moments table is unavailable")
        return StorageMomentsSeries(self)
    @property
    def dpn(self): return self.moments.default.dp_n
    @property
    def dpw(self): return self.moments.default.dp_w
    @property
    def mn(self): return self.moments.default.mn
    @property
    def mw(self): return self.moments.default.mw
    @property
    def mz(self): return self.moments.default.mz
    @property
    def dispersity(self): return self.moments.default.dispersity
    @analysis_operation(CLD_HELP)
    def cld(self, *, snapshot="final", **kwargs): return self._resolve_chain_snapshot(snapshot).cld(**kwargs)
    @analysis_operation(MWD_HELP)
    def mwd(self, *, snapshot="final", **kwargs): return self._resolve_chain_snapshot(snapshot).mwd(**kwargs)
    @analysis_operation(SPECTRUM_HELP)
    def chain_mass_spectrum(self, *, snapshot="final", **kwargs): return self._resolve_chain_snapshot(snapshot).chain_mass_spectrum(**kwargs)

    def _resolve_chain_snapshot(self, snapshot):
        if snapshot == "final":
            return self.final
        if snapshot == "last":
            return self.last
        if isinstance(snapshot, (int, np.integer)):
            return self.at_snapshot(int(snapshot))
        if hasattr(snapshot, "chains"):
            return snapshot
        raise ValueError("snapshot must be 'final', 'last', a snapshot_id, or a Snapshot")

    def chain_counts(self, *, snapshot="final", pool="all", grouping: str = "dp"):
        if grouping != "dp":
            raise ValueError("Slimmc Storage chain_counts currently supports grouping='dp' only")
        from .chain_counts import ChainCounts, ChainCountsGroup
        snap = self._resolve_chain_snapshot(snapshot)
        if isinstance(pool, (tuple, list)):
            return ChainCountsGroup({
                str(name): ChainCounts.from_population(snap.chains.select(pool=str(name)), pool=str(name))
                for name in pool
            })
        selected = snap.chains.select(pool=str(pool))
        return ChainCounts.from_population(selected, pool=str(pool))

    def mass_audit(self, *, tolerance: float = 1.0e-9, snapshot="final",
                   mass_model: str | None = None) -> MassAuditResult:
        if mass_model is None:
            model = self._metadata.get("model", {})
            params = model.get("parameters", {}) if isinstance(model, dict) else {}
            mass_model = str(params.get("mass_model", "repeat_units"))
        if mass_model not in {"with_end_groups", "repeat_units"}:
            raise ValueError("mass_model must be 'repeat_units' or 'with_end_groups'")
        pop = self._resolve_chain_snapshot(snapshot).chains
        monomer_entries: list[MassEntry] = []
        monomer_masses: dict[str, float | None] = {}
        for name, meta in self.monomers.items():
            value = meta.get("molar_mass_increment", meta.get("molar_mass"))
            declared = value is not None and np.isfinite(float(value))
            mw = float(value) if declared else None
            monomer_masses[name] = mw
            monomer_entries.append(MassEntry(name, mw, declared, "model" if declared else "missing"))

        end_entries: list[MassEntry] = []
        end_masses: dict[str, float | None] = {}
        for name, meta in self.endgroups.items():
            known = meta.get("has_known_molar_mass_contribution")
            value = meta.get("molar_mass_contribution", meta.get("molar_mass"))
            declared = value is not None and (known is not False) and np.isfinite(float(value))
            mw = float(value) if declared else None
            end_masses[name] = mw
            end_entries.append(MassEntry(name, mw, declared, "model" if declared else "missing"))

        missing_monomers = tuple(e.name for e in monomer_entries if not e.declared)
        used_ends = set(map(str, np.asarray(pop.left_end, dtype=object))) | set(map(str, np.asarray(pop.right_end, dtype=object)))
        ignored = {"not_applicable", "unknown", "0", ""}
        missing_endgroups = tuple(sorted(
            name for name in used_ends if name not in ignored and end_masses.get(name) is None
        )) if mass_model == "with_end_groups" else ()

        expected = np.zeros(len(pop), dtype=float)
        comp = np.asarray(pop.composition.matrix, dtype=float)
        for i, name in enumerate(pop.composition.names):
            mw = monomer_masses.get(name)
            if mw is not None:
                expected += comp[:, i] * mw
        if mass_model == "with_end_groups":
            expected += np.asarray([end_masses.get(str(name)) or 0.0 for name in pop.left_end])
            expected += np.asarray([end_masses.get(str(name)) or 0.0 for name in pop.right_end])
        actual = np.asarray(pop.molar_mass, dtype=float)
        delta = actual - expected
        row_ok = np.abs(delta) <= tolerance
        from .table import Table
        details = Table(
            ("chain_record_id", "dp", "count", "molar_mass", "expected_mass", "delta_mass", "ok"),
            zip(np.asarray(pop.chain_record_id).tolist(), np.asarray(pop.dp).tolist(),
                np.asarray(pop.count).tolist(), actual.tolist(), expected.tolist(),
                delta.tolist(), row_ok.tolist()),
            name="mass_audit",
        )
        numeric_ok = bool(np.all(row_ok)) if len(row_ok) else True
        ok = not missing_monomers and not missing_endgroups and numeric_ok
        warnings: list[str] = []
        if not numeric_ok:
            warnings.append("stored and independently reconstructed chain masses differ")
        return MassAuditResult(
            ok=ok, mass_model=mass_model, checked_records=len(pop),
            checked_chains=int(np.sum(np.asarray(pop.count), dtype=np.int64)),
            missing_monomers=missing_monomers, missing_endgroups=missing_endgroups,
            implicit_zero_monomers=missing_monomers,
            implicit_zero_endgroups=missing_endgroups, warnings=tuple(warnings),
            details=details, entries=tuple(monomer_entries + end_entries),
        )

    def _action_messages(self) -> dict[int,str]:
        path=self.path/"actions"/"messages.jsonl"
        if not path.is_file(): return {}
        out={}
        for obj in _read_jsonl(path): out[int(obj["action_id"])]=str(obj.get("message",""))
        return out

    @property
    def actions(self) -> StorageActions:
        return StorageActions(self)

    @property
    def channels(self) -> StorageChannelsSeries:
        return StorageChannelsSeries(self)
    @property
    def channel_events(self) -> StorageTable:
        return self.table("channel_events")
    @property
    def event_counts(self): return self.channels.event_count

    @property
    def firings(self):
        from .storage_analysis import StorageFirings
        return StorageFirings(self)

    @property
    def microstructure(self):
        if self.kinetic_model not in {"copo", "copolymer", "terpolymer"} and "copo" not in self.engine:
            from .core import ChemicalAnalysisNotApplicableError
            raise ChemicalAnalysisNotApplicableError("microstructure analysis is available only for copolymer runs")
        from .storage_analysis import StorageMicrostructure
        return StorageMicrostructure(self)

    @property
    def copolymerization(self):
        if self.kinetic_model not in {"copo", "copolymer", "terpolymer"} and "copo" not in self.engine:
            from .core import ChemicalAnalysisNotApplicableError
            raise ChemicalAnalysisNotApplicableError("copolymerization analysis is available only for copolymer runs")
        from .storage_analysis import StorageCopolymerization
        return StorageCopolymerization(self)

    @property
    def kinetics(self) -> StorageKineticsSeries:
        if "kinetic_parameters/values" not in self.tables: raise DataUnavailableError("kinetic parameters are unavailable")
        return StorageKineticsSeries(self)
    @property
    def temp(self): return self.kinetics.temperature
    @property
    def k(self): return self.kinetics.rate_constants

    def column_unit(self, table: str, column: str) -> str | None:
        for record in self.schema_records:
            if (record.get("record_type") == "column" and
                    record.get("table") == table and record.get("name") == column):
                unit = record.get("unit")
                return None if unit in (None, "", "1") else str(unit)
        return None

    def dictionary(self, name: str) -> dict[int, dict]:
        try:
            return self.dictionaries[name]
        except KeyError as exc:
            raise KeyError(f"Unknown schema dictionary {name!r}") from exc


def _read_jsonl(path: Path) -> list[dict]:
    records: list[dict] = []
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            raise InvalidOutputError(f"Blank JSONL line in {path} at line {lineno}")
        try:
            obj = json.loads(line)
        except json.JSONDecodeError as exc:
            raise InvalidOutputError(f"Invalid JSONL in {path} at line {lineno}: {exc}") from exc
        if not isinstance(obj, dict):
            raise InvalidOutputError(f"JSONL record in {path} at line {lineno} is not an object")
        records.append(obj)
    return records


def open_storage(path: str | Path, *, allow_incomplete: bool = False,
                    mmap_mode: str | None = "r") -> StorageRun:
    root = Path(path)
    metadata_path = root / "run_metadata.json"
    schema_path = root / "schema.jsonl"
    if not metadata_path.is_file() or not schema_path.is_file():
        raise InvalidOutputError(f"Not a Slimmc Storage 1.2.0 run directory: {root}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    storage_name = metadata.get("storage")
    if storage_name != "slimmc-storage":
        raise InvalidOutputError(f"Unsupported Slimmc Storage identifier in {metadata_path}")
    status = metadata.get("run_status")
    complete_marker = (root / "RESULTS_COMPLETE").is_file()
    if status == "completed" and not complete_marker:
        raise InvalidOutputError("Run claims completed but RESULTS_COMPLETE is missing")
    if status != "completed" and not allow_incomplete:
        raise IncompleteResultsError(
            f"Run status is {status!r}; pass allow_incomplete=True for diagnostic access"
        )

    records = _read_jsonl(schema_path)
    dictionaries: dict[str, dict[int, dict]] = {}
    table_names: list[str] = []
    table_required: dict[str, bool] = {}
    column_records: dict[str, list[dict]] = {}
    for record in records:
        if record.get("record_type") == "table":
            name = str(record["name"])
            table_names.append(name)
            table_required[name] = bool(record.get("required", False))
        elif record.get("record_type") == "column":
            column_records.setdefault(str(record["table"]), []).append(record)
        elif record.get("record_type") == "dictionary_entry":
            dictionaries.setdefault(str(record["dictionary"]), {})[int(record["id"])] = record

    tables: dict[str, StorageTable] = {}
    for name in table_names:
        directory = root / name
        if not directory.is_dir():
            if status == "completed" and table_required.get(name, False):
                raise InvalidOutputError(f"Required Slimmc Storage table directory is missing: {directory}")
            continue
        for record in column_records.get(name, ()):
            file_name = str(record.get("file", f"{record['name']}.npy"))
            column_path = directory / file_name
            if status == "completed" and bool(record.get("required", False)) and not column_path.is_file():
                raise InvalidOutputError(f"Required Slimmc Storage column is missing: {column_path}")
        try:
            table = StorageTable.open(directory, mmap_mode=mmap_mode)
        except (OSError, ValueError, EOFError) as exc:
            raise InvalidOutputError(f"Cannot read Slimmc Storage table {directory}: {exc}") from exc
        for record in column_records.get(name, ()):
            column_name = str(record["name"])
            if column_name not in table:
                continue
            expected = str(record.get("dtype", ""))
            actual = np.asarray(table[column_name]).dtype.name
            aliases = {"bool": "bool", "float64": "float64", "uint64": "uint64", "uint32": "uint32", "int64": "int64", "int32": "int32"}
            if expected in aliases and actual != aliases[expected]:
                raise InvalidOutputError(
                    f"Unexpected dtype for {name}/{column_name}: expected {expected}, got {actual}"
                )
        tables[name] = table

    return StorageRun(
        path=root,
        _metadata=metadata,
        run_id=str(metadata.get("run_id", root.name)),
        _prefix=str(metadata.get("run_id", root.name)),
        schema_records=records,
        dictionaries=dictionaries,
        _tables=tables,
    )
