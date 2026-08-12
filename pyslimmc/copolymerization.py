"""Data-first copolymer composition analysis.

All quantities are reconstructed from engine-written state and cumulative
channel-firing ledgers.  In particular, interval composition uses insertion
events, never a difference of net polymerized inventory.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .run import DataConsistencyError
from .snapshots import NamedValues, _readonly
from .core import ChemicalModelIncompatibleError, ChemicalAnalysisNotApplicableError


class PairValues:
    """Read-only mapping whose keys are ``(terminal, incoming)`` pairs."""

    def __init__(self, values):
        self._values = MappingProxyType(dict(values))

    def __getitem__(self, pair):
        return self._values[pair]

    def keys(self):
        return tuple(self._values)

    def items(self):
        return tuple(self._values.items())

    @property
    def array(self):
        if not self._values:
            return _readonly(np.empty((0, 0)))
        return _readonly(np.column_stack([self._values[key] for key in self._values]))


@dataclass(frozen=True)
class Capability:
    name: str
    implemented: bool
    applicable: bool
    data_available: bool
    reason: str | None = None

    @property
    def available(self) -> bool:
        return self.implemented and self.applicable and self.data_available


class Capabilities:
    def __init__(self, values):
        self._values = MappingProxyType(dict(values))

    def __getitem__(self, name):
        return self._values[name]

    def __getattr__(self, name):
        try:
            return self._values[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def keys(self):
        return tuple(self._values)

    def items(self):
        return tuple(self._values.items())


@dataclass(frozen=True)
class ReactivityRatioSeries:
    monomers: tuple[str, str]
    snapshot_ids: np.ndarray
    times: np.ndarray
    rate_constants: PairValues
    reactivity_ratios: NamedValues
    is_defined: NamedValues

    @property
    def array(self):
        return self.reactivity_ratios.array

    def __len__(self):
        return len(self.snapshot_ids)

    @property
    def final(self):
        return self.at_index(-1)

    def at_index(self, index: int):
        return _Row(self, _index(len(self), index))

    def at_snapshot(self, snapshot_id: int):
        found = np.flatnonzero(self.snapshot_ids == int(snapshot_id))
        if not len(found):
            raise KeyError(f"unknown snapshot_id {snapshot_id}; available: {self.snapshot_ids.tolist()}")
        return _Row(self, int(found[0]))


@dataclass(frozen=True)
class MayoLewisSeries:
    monomers: tuple[str, str]
    snapshot_ids: np.ndarray
    times: np.ndarray
    conversion: np.ndarray
    rate_constants: PairValues
    reactivity_ratios: NamedValues
    monomer_mole_fractions: NamedValues
    instantaneous_repeat_unit_fractions: NamedValues
    is_defined: np.ndarray
    assumptions: Mapping[str, bool]

    @property
    def fraction_array(self):
        return self.instantaneous_repeat_unit_fractions.array

    def __len__(self):
        return len(self.snapshot_ids)

    @property
    def final(self):
        return self.at_index(-1)

    def at_index(self, index: int):
        return _Row(self, _index(len(self), index))

    def at_snapshot(self, snapshot_id: int):
        found = np.flatnonzero(self.snapshot_ids == int(snapshot_id))
        if not len(found):
            raise KeyError(f"unknown snapshot_id {snapshot_id}; available: {self.snapshot_ids.tolist()}")
        return _Row(self, int(found[0]))

    def to_tsv(self, path):
        a, b = self.monomers
        columns = ["snapshot_id", "t", "is_defined"]
        columns += [f"f_{a}", f"f_{b}", f"r_{a}", f"r_{b}", f"F_instantaneous_{a}", f"F_instantaneous_{b}"]
        rows = [[self.snapshot_ids[i], self.times[i], int(self.is_defined[i]),
                 self.monomer_mole_fractions[a][i], self.monomer_mole_fractions[b][i],
                 self.reactivity_ratios[a][i], self.reactivity_ratios[b][i],
                 self.instantaneous_repeat_unit_fractions[a][i],
                 self.instantaneous_repeat_unit_fractions[b][i]] for i in range(len(self))]
        return _write_tsv(path, columns, rows, "mayo_lewis")

    def plot(self, path=None, *, x="conversion", style="screen", ax=None, span=None):
        if x == "time":
            axis, xlabel = self.times, "time"
        elif x == "conversion":
            axis, xlabel = self.conversion, "overall conversion"
        else:
            raise ValueError("x must be 'conversion' or 'time'")
        return _plot_lines(axis, self.instantaneous_repeat_unit_fractions, xlabel,
                           "instantaneous repeat-unit fraction", path, style, ax=ax, span=span)


@dataclass(frozen=True)
class MayoLewisComparison:
    monomers: tuple[str, str]
    start_snapshot_ids: np.ndarray
    end_snapshot_ids: np.ndarray
    t_start: np.ndarray
    t_end: np.ndarray
    t_mid: np.ndarray
    dt: np.ndarray
    conversion: np.ndarray
    monomer_reference: str
    parameter_reference: str
    monomer_mole_fractions: NamedValues
    mayo_lewis_fractions: NamedValues
    incremental_repeat_unit_fractions: NamedValues
    composition_difference: NamedValues
    is_defined: np.ndarray

    def __len__(self):
        return len(self.end_snapshot_ids)

    def at_index(self, index):
        return _Row(self, _index(len(self), index))

    def ending_at_snapshot(self, snapshot_id):
        found = np.flatnonzero(self.end_snapshot_ids == int(snapshot_id))
        if not len(found):
            raise KeyError(f"no interval ending at snapshot_id {snapshot_id}; available: {self.end_snapshot_ids.tolist()}")
        return _Row(self, int(found[0]))

    def to_tsv(self, path):
        columns = ["start_snapshot_id", "end_snapshot_id", "t_start", "t_end", "t_mid", "dt",
                   "conversion", "is_defined"]
        for name in self.monomers:
            columns += [f"f_{name}", f"F_mayo_lewis_{name}", f"F_incremental_{name}",
                        f"composition_difference_{name}"]
        rows = []
        for i in range(len(self)):
            row = [self.start_snapshot_ids[i], self.end_snapshot_ids[i], self.t_start[i], self.t_end[i],
                   self.t_mid[i], self.dt[i], self.conversion[i], int(self.is_defined[i])]
            for name in self.monomers:
                row += [self.monomer_mole_fractions[name][i], self.mayo_lewis_fractions[name][i],
                        self.incremental_repeat_unit_fractions[name][i], self.composition_difference[name][i]]
            rows.append(row)
        return _write_tsv(path, columns, rows,
                          f"compare_mayo_lewis; monomer_reference={self.monomer_reference}; parameter_reference={self.parameter_reference}")

    def plot(self, path=None, *, style="screen", ax=None, span=None):
        return _plot_lines(self.conversion, self.composition_difference, "overall conversion",
                           "incremental minus Mayo-Lewis", path, style, ax=ax, span=span)


@dataclass(frozen=True)
class TerminalTransitionDiagnostics:
    monomers: tuple[str, str]
    start_snapshot_ids: np.ndarray
    end_snapshot_ids: np.ndarray
    t_start: np.ndarray
    t_end: np.ndarray
    t_mid: np.ndarray
    dt: np.ndarray
    predicted_probabilities: PairValues
    observed_probabilities: PairValues
    transition_difference: PairValues
    outgoing_fires: NamedValues
    is_defined: NamedValues


@dataclass(frozen=True)
class TerminalBlockDiagnostics:
    monomers: tuple[str, str]
    snapshot_ids: np.ndarray
    times: np.ndarray
    predicted_Ln: NamedValues
    predicted_Lw: NamedValues
    predicted_dispersity: NamedValues
    observed_Ln: NamedValues
    observed_Lw: NamedValues
    observed_dispersity: NamedValues
    block_count: NamedValues
    is_defined: NamedValues


@dataclass(frozen=True)
class TerminalDiagnostics:
    transitions: TerminalTransitionDiagnostics
    blocks: TerminalBlockDiagnostics


@dataclass(frozen=True)
class PenultimateParameterSeries:
    monomers: tuple[str, str]
    snapshot_ids: np.ndarray
    times: np.ndarray
    tensor: np.ndarray
    classification: tuple[str, ...]
    variable: bool
    r: NamedValues
    r_prime: NamedValues
    s: NamedValues
    is_defined: np.ndarray

    def __len__(self):
        return len(self.snapshot_ids)

    @property
    def final(self):
        return self.at_index(-1)

    def at_index(self, index):
        return _Row(self, _index(len(self), index))

    def at_snapshot(self, snapshot_id):
        found = np.flatnonzero(self.snapshot_ids == int(snapshot_id))
        if not len(found):
            raise KeyError(f"unknown snapshot_id {snapshot_id}; available: {self.snapshot_ids.tolist()}")
        return _Row(self, int(found[0]))


@dataclass(frozen=True)
class PenultimateCompositionSeries:
    monomers: tuple[str, str]
    snapshot_ids: np.ndarray
    times: np.ndarray
    conversion: np.ndarray
    monomer_mole_fractions: NamedValues
    instantaneous_repeat_unit_fractions: NamedValues
    radical_state_fractions: PairValues
    transition_probabilities: object
    is_defined: np.ndarray

    @property
    def fraction_array(self):
        return self.instantaneous_repeat_unit_fractions.array

    def __len__(self):
        return len(self.snapshot_ids)

    @property
    def final(self):
        return self.at_index(-1)

    def at_index(self, index):
        return _Row(self, _index(len(self), index))


@dataclass(frozen=True)
class PenultimateComparison:
    monomers: tuple[str, str]
    start_snapshot_ids: np.ndarray
    end_snapshot_ids: np.ndarray
    t_start: np.ndarray
    t_end: np.ndarray
    t_mid: np.ndarray
    dt: np.ndarray
    conversion: np.ndarray
    monomer_reference: str
    parameter_reference: str
    penultimate_fractions: NamedValues
    incremental_repeat_unit_fractions: NamedValues
    composition_difference: NamedValues
    is_defined: np.ndarray


class TripleValues:
    """Read-only mapping keyed by ``(previous, terminal, incoming)``."""
    def __init__(self, values):
        self._values = MappingProxyType(dict(values))

    def __getitem__(self, key):
        return self._values[key]

    def keys(self):
        return tuple(self._values)

    def items(self):
        return tuple(self._values.items())


@dataclass(frozen=True)
class PenultimateDiagnostics:
    monomers: tuple[str, str]
    start_snapshot_ids: np.ndarray
    end_snapshot_ids: np.ndarray
    predicted_transitions: TripleValues
    observed_transitions: TripleValues
    transition_difference: TripleValues
    outgoing_fires: PairValues
    transition_is_defined: PairValues
    snapshot_ids: np.ndarray
    predicted_radical_states: PairValues
    observed_radical_states: PairValues
    radical_state_difference: PairValues
    predicted_triads: TripleValues
    observed_triads: TripleValues
    triad_difference: TripleValues
    radical_states_available: bool


def _fractions(values: Mapping[str, np.ndarray]) -> tuple[NamedValues, np.ndarray]:
    names = tuple(values)
    if not names:
        return NamedValues({}), _readonly([], dtype=bool)
    matrix = np.column_stack([np.asarray(values[name], dtype=float) for name in names])
    total = matrix.sum(axis=1)
    defined = total > 0
    out = np.full(matrix.shape, np.nan, dtype=float)
    np.divide(matrix, total[:, None], out=out, where=defined[:, None])
    return NamedValues({name: _readonly(out[:, i]) for i, name in enumerate(names)}), _readonly(defined)


def _penultimate_prediction(k: np.ndarray, f: np.ndarray):
    """Return (F, stationary states, transition rows, defined) for binary PUE."""
    states = ((0, 0), (0, 1), (1, 0), (1, 1))
    transition = np.full((4, 2), np.nan)
    matrix = np.zeros((4, 4), dtype=float)
    for row, (previous, terminal) in enumerate(states):
        propensity = np.asarray([k[previous, terminal, incoming] * f[incoming] for incoming in range(2)])
        total = propensity.sum()
        if not np.isfinite(total) or total <= 0:
            return np.full(2, np.nan), np.full(4, np.nan), transition, False
        transition[row] = propensity / total
        for incoming in range(2):
            matrix[row, states.index((terminal, incoming))] += transition[row, incoming]
    # A unique stationary distribution is required. Reducible/ambiguous chains
    # are chemically underdetermined without an initial radical-state mixture.
    if np.linalg.matrix_rank(matrix.T - np.eye(4), tol=1e-10) != 3:
        return np.full(2, np.nan), np.full(4, np.nan), transition, False
    augmented = np.vstack([matrix.T - np.eye(4), np.ones(4)])
    target = np.asarray([0.0, 0.0, 0.0, 0.0, 1.0])
    stationary, *_ = np.linalg.lstsq(augmented, target, rcond=None)
    stationary[np.abs(stationary) < 1e-13] = 0.0
    if np.any(stationary < -1e-10) or not np.isclose(stationary.sum(), 1.0, atol=1e-9):
        return np.full(2, np.nan), np.full(4, np.nan), transition, False
    stationary = np.maximum(stationary, 0.0)
    stationary /= stationary.sum()
    fractions = stationary @ transition
    return fractions, stationary, transition, True


def _index(length: int, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError("index must be an integer")
    i = int(value)
    if i < 0:
        i += length
    if i < 0 or i >= length:
        raise IndexError(value)
    return i


class _Row:
    def __init__(self, series, index: int):
        self._series = series
        self._index = index

    def __getattr__(self, name):
        value = getattr(self._series, name)
        if isinstance(value, NamedValues):
            return MappingProxyType({key: value[key][self._index].item() for key in value.keys()})
        if isinstance(value, PairValues):
            return MappingProxyType({key: value[key][self._index].item() for key in value.keys()})
        if isinstance(value, np.ndarray):
            item = value[self._index]
            return item.item() if isinstance(item, np.generic) else item
        raise AttributeError(name)


@dataclass(frozen=True)
class SnapshotCompositionSeries:
    snapshot_ids: np.ndarray
    times: np.ndarray
    conversion: np.ndarray
    fractions: NamedValues
    is_defined: np.ndarray
    counts: NamedValues
    fraction_name: str

    @property
    def mole_fractions(self) -> NamedValues:
        return self.fractions

    @property
    def repeat_unit_fractions(self) -> NamedValues:
        return self.fractions

    @property
    def fraction_array(self) -> np.ndarray:
        return self.fractions.array

    def __len__(self):
        return len(self.snapshot_ids)

    def at_index(self, index: int) -> _Row:
        return _Row(self, _index(len(self), index))

    def at_snapshot(self, snapshot_id: int) -> _Row:
        if isinstance(snapshot_id, bool) or not isinstance(snapshot_id, (int, np.integer)):
            raise TypeError("snapshot_id must be an integer")
        found = np.flatnonzero(self.snapshot_ids == int(snapshot_id))
        if not len(found):
            raise KeyError(f"unknown snapshot_id {snapshot_id}; available: {self.snapshot_ids.tolist()}")
        return _Row(self, int(found[0]))

    def to_tsv(self, path) -> Path:
        columns = ["snapshot_id", "t", "conversion", "is_defined"]
        columns += [f"count_{name}" for name in self.counts.keys()]
        columns += [f"fraction_{name}" for name in self.fractions.keys()]
        rows = []
        for i in range(len(self)):
            rows.append([self.snapshot_ids[i], self.times[i], self.conversion[i], int(self.is_defined[i])]
                        + [self.counts[name][i] for name in self.counts.keys()]
                        + [self.fractions[name][i] for name in self.fractions.keys()])
        return _write_tsv(path, columns, rows, self.fraction_name)

    def plot(self, path=None, *, x: str = "conversion", style: str = "screen",
             ax=None, span: str | None = None):
        x_values, xlabel = (self.conversion, "overall conversion") if x == "conversion" else (self.times, "time")
        if x not in {"conversion", "time"}:
            raise ValueError("x must be 'conversion' or 'time'")
        return _plot_lines(x_values, self.fractions, xlabel, self.fraction_name, path, style,
                           ax=ax, span=span)


@dataclass(frozen=True)
class IntervalCompositionSeries:
    start_snapshot_ids: np.ndarray
    end_snapshot_ids: np.ndarray
    t_mid: np.ndarray
    dt: np.ndarray
    conversion: np.ndarray
    fractions: NamedValues
    is_defined: np.ndarray
    inserted_counts: NamedValues

    @property
    def repeat_unit_fractions(self) -> NamedValues:
        return self.fractions

    @property
    def fraction_array(self) -> np.ndarray:
        return self.fractions.array

    def __len__(self):
        return len(self.end_snapshot_ids)

    def at_index(self, index: int) -> _Row:
        return _Row(self, _index(len(self), index))

    def ending_at_snapshot(self, snapshot_id: int) -> _Row:
        if isinstance(snapshot_id, bool) or not isinstance(snapshot_id, (int, np.integer)):
            raise TypeError("snapshot_id must be an integer")
        found = np.flatnonzero(self.end_snapshot_ids == int(snapshot_id))
        if not len(found):
            raise KeyError(f"no interval ending at snapshot_id {snapshot_id}; available: {self.end_snapshot_ids.tolist()}")
        return _Row(self, int(found[0]))

    def to_tsv(self, path) -> Path:
        columns = ["start_snapshot_id", "end_snapshot_id", "t_mid", "dt", "conversion", "is_defined"]
        columns += [f"inserted_{name}" for name in self.inserted_counts.keys()]
        columns += [f"fraction_{name}" for name in self.fractions.keys()]
        rows = []
        for i in range(len(self)):
            rows.append([self.start_snapshot_ids[i], self.end_snapshot_ids[i], self.t_mid[i], self.dt[i],
                         self.conversion[i], int(self.is_defined[i])]
                        + [self.inserted_counts[name][i] for name in self.inserted_counts.keys()]
                        + [self.fractions[name][i] for name in self.fractions.keys()])
        return _write_tsv(path, columns, rows, "incremental_repeat_unit_fraction")

    def plot(self, path=None, *, x: str = "conversion", style: str = "screen",
             ax=None, span: str | None = None):
        x_values, xlabel = (self.conversion, "overall conversion") if x == "conversion" else (self.t_mid, "interval midpoint time")
        if x not in {"conversion", "time"}:
            raise ValueError("x must be 'conversion' or 'time'")
        return _plot_lines(x_values, self.fractions, xlabel, "incremental repeat-unit fraction", path, style,
                           ax=ax, span=span)


@dataclass(frozen=True)
class CompositionResult:
    free: SnapshotCompositionSeries
    incremental: IntervalCompositionSeries
    cumulative: SnapshotCompositionSeries
    inserted_counts: NamedValues
    removed_counts: NamedValues
    net_counts: NamedValues


@dataclass(frozen=True)
class CompositionDrift:
    start_snapshot_ids: np.ndarray
    end_snapshot_ids: np.ndarray
    t_mid: np.ndarray
    dt: np.ndarray
    conversion: np.ndarray
    incremental: NamedValues
    monomer_reference_fractions: NamedValues
    composition_difference: NamedValues
    is_defined: np.ndarray
    monomer_reference: str

    def __len__(self):
        return len(self.end_snapshot_ids)

    def at_index(self, index: int) -> _Row:
        return _Row(self, _index(len(self), index))

    def ending_at_snapshot(self, snapshot_id: int) -> _Row:
        if isinstance(snapshot_id, bool) or not isinstance(snapshot_id, (int, np.integer)):
            raise TypeError("snapshot_id must be an integer")
        found = np.flatnonzero(self.end_snapshot_ids == int(snapshot_id))
        if not len(found):
            raise KeyError(f"no interval ending at snapshot_id {snapshot_id}; available: {self.end_snapshot_ids.tolist()}")
        return _Row(self, int(found[0]))

    def to_tsv(self, path) -> Path:
        names = self.composition_difference.keys()
        columns = ["start_snapshot_id", "end_snapshot_id", "t_mid", "dt", "conversion", "is_defined"]
        for name in names:
            columns += [f"incremental_{name}", f"reference_{name}", f"composition_difference_{name}"]
        rows = []
        for i in range(len(self)):
            row = [self.start_snapshot_ids[i], self.end_snapshot_ids[i], self.t_mid[i], self.dt[i],
                   self.conversion[i], int(self.is_defined[i])]
            for name in names:
                row += [self.incremental[name][i], self.monomer_reference_fractions[name][i],
                        self.composition_difference[name][i]]
            rows.append(row)
        return _write_tsv(path, columns, rows, f"composition_drift; monomer_reference={self.monomer_reference}")

    def plot(self, path=None, *, style: str = "screen", ax=None, span: str | None = None):
        return _plot_lines(self.conversion, self.composition_difference, "overall conversion",
                           "composition difference", path, style, ax=ax, span=span)




def _write_tsv(path, columns, rows, object_name: str) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# object: {object_name}",
             f"# exported_at: {datetime.now(timezone.utc).isoformat()}",
             "\t".join(columns)]
    for row in rows:
        lines.append("\t".join(str(value) for value in row))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def _plot_lines(x, values: NamedValues, xlabel: str, ylabel: str, path, style: str,
                *, ax=None, span: str | None = None):
    from .plotting import apply_axes_style, create_axes, require_owned_geometry, style_kwargs
    require_owned_geometry(ax, span)
    if ax is None:
        figure, ax = create_axes(style, span=span)
    else:
        figure = ax.figure
    for i, name in enumerate(values.keys()):
        ax.plot(x, values[name], label=name, **style_kwargs(style, index=i))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if "difference" in ylabel:
        finite = np.concatenate([
            np.asarray(values[name], dtype=float)[np.isfinite(values[name])]
            for name in values.keys()
        ]) if values.keys() else np.asarray([], dtype=float)
        extent = max(0.05, float(np.max(np.abs(finite))) * 1.08) if finite.size else 0.05
        ax.set_ylim(-extent, extent)
    else:
        ax.set_ylim(0.0, 1.05)
    ax.legend()
    apply_axes_style(ax, style)
    if path is not None:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(target)
    return figure
