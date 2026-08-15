from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from .mass_model import record_masses, resolve_mass_model

_FORMS = {"number", "mass", "z", "log"}


def _readonly(values, *, dtype=float):
    a = np.asarray(values, dtype=dtype)
    a.flags.writeable = False
    return a


def _validate_form(form: str) -> str:
    form = str(form).lower()
    if form not in _FORMS:
        raise ValueError("form must be 'number', 'mass', 'z', or 'log'")
    return form


@dataclass(frozen=True)
class _DiscreteDistribution:
    _x: np.ndarray
    _y: np.ndarray
    form: str
    total_chains: int
    snapshot_id: int | None
    t: float | None
    metadata: dict[str, Any]

    def __post_init__(self):
        self._x.flags.writeable = False
        self._y.flags.writeable = False

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def meta(self):
        return self.metadata

    @property
    def is_empty(self):
        return self._x.size == 0

    def to_tsv(self, path: str | Path):
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="") as handle:
            handle.write(f"# form: {self.form}\n")
            handle.write("x\ty\n")
            for x, y in zip(self.x, self.y):
                handle.write(f"{float(x):.17g}\t{float(y):.17g}\n")
        return target

    def plot(self, *, ax=None, **plot_kwargs):
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:  # pragma: no cover
            raise ImportError("plot() requires optional dependency matplotlib") from exc
        if ax is None:
            _, ax = plt.subplots()
        ax.plot(self.x, self.y, **plot_kwargs)
        kind = self.metadata.get("kind")
        if kind == "CLD":
            ax.set_xlabel("log10(DP)" if self.form == "log" else "Degree of polymerization, DP")
            labels = {
                "number": "Chain number fraction",
                "mass": "Polymer mass fraction",
                "z": "z-weighted fraction",
                "log": "Polymer mass fraction",
            }
        else:
            ax.set_xlabel("log10(M / g mol$^{-1}$)" if self.form == "log" else "Molar mass, g mol$^{-1}$")
            labels = {
                "number": "Chain number fraction",
                "mass": "Polymer mass fraction",
                "z": "z-weighted fraction",
                "log": "Polymer mass fraction",
            }
        ax.set_ylabel(labels[self.form])
        return ax

    def info(self):
        return (
            f"{type(self).__name__}\n"
            f"  form: {self.form}\n"
            f"  representation: discrete\n"
            f"  total_chains: {self.total_chains}"
        )

    def help(self):
        return self.info()


@dataclass(frozen=True)
class ChainLengthDistribution(_DiscreteDistribution):
    _dp: np.ndarray
    dpn: float
    dpw: float
    dpz: float

    @property
    def dp(self):
        return self._dp

    @property
    def dispersity(self):
        return self.dpw / self.dpn if self.dpn else float("nan")


@dataclass(frozen=True)
class MolarMassDistribution(_DiscreteDistribution):
    _mass: np.ndarray
    mass_model: str
    mn: float
    mw: float
    mz: float

    @property
    def mass(self):
        return self._mass

    @property
    def dispersity(self):
        return self.mw / self.mn if self.mn else float("nan")


def build_cld(population, *, form="number", mass_model: str | None = None):
    form = _validate_form(form)
    dp_record = np.asarray(population.dp, dtype=float)
    count_record = np.asarray(population.count, dtype=float)
    if np.any(dp_record <= 0):
        raise ValueError("CLD requires strictly positive degree of polymerization")

    unique_dp, inverse = np.unique(dp_record, return_inverse=True)
    grouped_count = np.zeros(unique_dp.size, dtype=float)
    np.add.at(grouped_count, inverse, count_record)

    resolved_model = resolve_mass_model(population, mass_model)
    if form in {"mass", "log"}:
        mass_record, resolved_model = record_masses(population, resolved_model)
        grouped_mass = np.zeros(unique_dp.size, dtype=float)
        np.add.at(grouped_mass, inverse, count_record * mass_record)
        denom = float(np.sum(grouped_mass))
        y = grouped_mass / denom if denom else grouped_mass
    elif form == "number":
        denom = float(np.sum(grouped_count))
        y = grouped_count / denom if denom else grouped_count
    else:  # z
        raw = unique_dp * unique_dp * grouped_count
        denom = float(np.sum(raw))
        y = raw / denom if denom else raw

    x = np.log10(unique_dp) if form == "log" else unique_dp.copy()
    source_moments = population.moments(mass_model=resolved_model)
    dpn, dpw, dpz = source_moments.dpn, source_moments.dpw, source_moments.dpz
    meta = {
        "kind": "CLD",
        "form": form,
        "representation": "discrete",
        "mass_model": resolved_model,
    }
    return ChainLengthDistribution(
        _readonly(x), _readonly(y), form, int(np.sum(count_record)),
        int(population.snapshot_id), population.t, meta,
        _readonly(unique_dp), dpn, dpw, dpz,
    )


def build_mwd(population, *, form="log", mass_model: str | None = None):
    form = _validate_form(form)
    mass_record, model = record_masses(population, mass_model)
    count_record = np.asarray(population.count, dtype=float)

    unique_mass, inverse = np.unique(mass_record, return_inverse=True)
    grouped_count = np.zeros(unique_mass.size, dtype=float)
    np.add.at(grouped_count, inverse, count_record)

    if form == "number":
        raw = grouped_count
    elif form in {"mass", "log"}:
        raw = unique_mass * grouped_count
    else:  # z
        raw = unique_mass * unique_mass * grouped_count
    denom = float(np.sum(raw))
    y = raw / denom if denom else raw
    x = np.log10(unique_mass) if form == "log" else unique_mass.copy()

    source_moments = population.moments(mass_model=model)
    mn, mw, mz = source_moments.mn, source_moments.mw, source_moments.mz
    meta = {
        "kind": "MWD",
        "form": form,
        "representation": "discrete",
        "mass_model": model,
    }
    return MolarMassDistribution(
        _readonly(x), _readonly(y), form, int(np.sum(count_record)),
        int(population.snapshot_id), population.t, meta,
        _readonly(unique_mass), model, mn, mw, mz,
    )


@dataclass(frozen=True)
class DistributionGroup:
    """Named exact discrete distributions sharing one normalization contract.

    Each member keeps its own exact support.  ``per_series`` normalizes every
    member independently.  ``combined`` preserves the members' exact physical
    contributions and is only valid for pairwise-disjoint source populations.
    """
    series: Mapping[str, _DiscreteDistribution]
    form: str
    normalization: str
    kind: str
    series_disjoint: bool
    metadata: dict[str, Any]

    def __post_init__(self):
        object.__setattr__(self, "series", MappingProxyType(dict(self.series)))

    @property
    def series_names(self) -> tuple[str, ...]:
        return tuple(self.series)

    def __getitem__(self, name: str):
        return self.series[name]

    @property
    def x(self):
        raise ValueError("DistributionGroup has multiple supports; use group.series['name'].x")

    @property
    def y(self):
        raise ValueError("DistributionGroup has multiple series; use group.series['name'].y")

    @property
    def is_empty(self) -> bool:
        return all(value.is_empty for value in self.series.values())

    def plot(self, *, ax=None, **plot_kwargs):
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:  # pragma: no cover
            raise ImportError("plot() requires optional dependency matplotlib") from exc
        if ax is None:
            _, ax = plt.subplots()
        for name, value in self.series.items():
            kwargs = dict(plot_kwargs)
            kwargs.setdefault("label", name)
            value.plot(ax=ax, **kwargs)
        return ax

    def info(self):
        return (
            f"{type(self).__name__}\n"
            f"  kind: {self.kind}\n"
            f"  form: {self.form}\n"
            f"  normalization: {self.normalization}\n"
            f"  series: {', '.join(self.series_names)}\n"
            f"  pairwise_disjoint: {self.series_disjoint}"
        )

    def help(self):
        return self.info()


def _select_series_population(population, selector):
    from .chains import ChainPopulation
    if isinstance(selector, ChainPopulation):
        return selector
    if not isinstance(selector, str):
        raise TypeError("each series value must be a ChainPopulation or pool selector string")
    if selector == "all":
        value = population.all
        return value() if callable(value) else value
    if selector in {"live", "active"}:
        value = population.live
        return value() if callable(value) else value
    if selector == "dead":
        value = population.dead
        return value() if callable(value) else value
    return population.pool(selector)


def _resolve_series(population, series) -> dict[str, Any]:
    if isinstance(series, (tuple, list)):
        series = {str(name): str(name) for name in series}
    if not isinstance(series, Mapping) or not series:
        raise TypeError("series must be a non-empty mapping or tuple/list of pool selectors")
    resolved = {}
    for name, selector in series.items():
        if not isinstance(name, str) or not name:
            raise ValueError("series names must be non-empty strings")
        selected = _select_series_population(population, selector)
        if int(selected.snapshot_id) != int(population.snapshot_id):
            raise ValueError("all series must come from the same snapshot_id")
        if selected.t != population.t:
            raise ValueError("all series must come from the same saved time")
        base_run = getattr(population, "run", None)
        selected_run = getattr(selected, "run", None)
        if base_run is not None and selected_run is not None and selected_run is not base_run:
            raise ValueError("all series must come from the same run")
        resolved[name] = selected
    return resolved


def _population_record_ids(population):
    root = getattr(population, "_analysis_root", None)
    indices = getattr(population, "_root_indices", None)
    if root is not None and indices is not None:
        return (id(root), frozenset(int(v) for v in np.asarray(indices, dtype=np.int64)))
    return None


def _row_signatures(population) -> set[tuple[Any, ...]]:
    raw = population._raw_arrays()
    # Prefer persistent record IDs when the storage exposes them.
    if "chain_record_id" in raw:
        return {("chain_record_id", int(v)) for v in np.asarray(raw["chain_record_id"])}
    names = tuple(sorted(name for name in raw if np.asarray(raw[name]).ndim == 1))
    signatures = set()
    for i in range(len(population)):
        values = []
        for name in names:
            value = raw[name][i]
            values.append(value.item() if isinstance(value, np.generic) else value)
        signatures.add(tuple(values))
    return signatures


def _are_disjoint(populations: Mapping[str, Any]) -> bool:
    seen_indices: set[tuple[int, int]] = set()
    can_use_indices = True
    for population in populations.values():
        identity = _population_record_ids(population)
        if identity is None:
            can_use_indices = False
            break
        root_id, indices = identity
        current = {(root_id, index) for index in indices}
        if seen_indices.intersection(current):
            return False
        seen_indices.update(current)
    if can_use_indices:
        return True

    seen = set()
    for population in populations.values():
        current = _row_signatures(population)
        if seen.intersection(current):
            return False
        seen.update(current)
    return True


def _source_total(population, *, kind: str, form: str, mass_model: str | None):
    count = np.asarray(population.count, dtype=float)
    if kind == "CLD":
        if form == "number":
            return float(np.sum(count))
        if form == "z":
            dp = np.asarray(population.dp, dtype=float)
            return float(np.dot(count, dp * dp))
        mass, _ = record_masses(population, mass_model)
        return float(np.dot(count, mass))
    if kind == "MWD":
        if form == "number":
            return float(np.sum(count))
        mass, _ = record_masses(population, mass_model)
        if form in {"mass", "log"}:
            return float(np.dot(count, mass))
        return float(np.dot(count, mass * mass))
    raise ValueError("kind must be CLD or MWD")


def _build_series(population, *, series, kind: str, form: str, normalization: str,
                  mass_model: str | None):
    form = _validate_form(form)
    normalization = str(normalization).lower()
    if normalization not in {"per_series", "combined"}:
        raise ValueError("normalization must be 'per_series' or 'combined'")
    populations = _resolve_series(population, series)
    disjoint = _are_disjoint(populations)
    if normalization == "combined" and not disjoint:
        raise ValueError("normalization='combined' requires pairwise-disjoint series")

    builder = build_cld if kind == "CLD" else build_mwd
    raw = {
        name: builder(selected, form=form, mass_model=mass_model)
        for name, selected in populations.items()
    }
    totals = {
        name: _source_total(selected, kind=kind, form=form, mass_model=mass_model)
        for name, selected in populations.items()
    }
    combined_total = float(sum(totals.values()))
    built = {}
    for name, value in raw.items():
        scale = 1.0
        if normalization == "combined":
            scale = totals[name] / combined_total if combined_total else 0.0
        meta = dict(value.metadata)
        meta.update(
            normalization=normalization,
            series_name=name,
            source_total=totals[name],
            combined_source_total=combined_total if normalization == "combined" else None,
        )
        built[name] = replace(value, _y=_readonly(value.y * scale), metadata=meta)

    metadata = {
        "kind": kind,
        "form": form,
        "representation": "discrete",
        "normalization": normalization,
        "series_names": tuple(built),
        "series_disjoint": disjoint,
        "source_totals": totals,
    }
    return DistributionGroup(built, form, normalization, kind, disjoint, metadata)


def build_cld_series(population, *, series, form="number", normalization="per_series",
                     mass_model: str | None = None):
    return _build_series(
        population, series=series, kind="CLD", form=form,
        normalization=normalization, mass_model=mass_model,
    )


def build_mwd_series(population, *, series, form="log", normalization="per_series",
                     mass_model: str | None = None):
    return _build_series(
        population, series=series, kind="MWD", form=form,
        normalization=normalization, mass_model=mass_model,
    )
