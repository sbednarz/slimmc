from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from .mass_model import record_masses, resolve_mass_model

_WEIGHTINGS = {"number", "mass", "z"}
_DEFAULT_LOG10_STEP = 0.005  # mcPolymer grid spacing


def _readonly(values, *, dtype=float):
    a = np.asarray(values, dtype=dtype)
    a.flags.writeable = False
    return a


def _validate_weighting(weighting: str) -> str:
    weighting = str(weighting).lower()
    if weighting not in _WEIGHTINGS:
        raise ValueError("weighting must be 'number', 'mass', or 'z'")
    return weighting


@dataclass(frozen=True)
class _ExactDistribution:
    _x: np.ndarray
    _y: np.ndarray
    weighting: str
    total_chains: int
    snapshot_id: int | None
    t: float | None
    metadata: dict[str, Any]

    def __post_init__(self):
        self._x.flags.writeable = False
        self._y.flags.writeable = False

    @property
    def x(self): return self._x

    @property
    def y(self): return self._y

    @property
    def meta(self): return self.metadata

    @property
    def representation(self): return "discrete"

    @property
    def is_empty(self): return self._x.size == 0

    def to_tsv(self, path: str | Path):
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="") as handle:
            handle.write(f"# weighting: {self.weighting}\n")
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
        if self.metadata.get("kind") == "CLD":
            ax.set_xlabel("Degree of polymerization, DP")
        else:
            ax.set_xlabel("Molar mass, g mol$^{-1}$")
        labels = {
            "number": "Chain number fraction",
            "mass": "Polymer mass fraction",
            "z": "z-weighted fraction",
        }
        ax.set_ylabel(labels[self.weighting])
        return ax

    def info(self):
        return (
            f"{type(self).__name__}\n"
            f"  weighting: {self.weighting}\n"
            f"  representation: discrete\n"
            f"  total_chains: {self.total_chains}"
        )

    def help(self): return self.info()


@dataclass(frozen=True)
class ChainLengthDistribution(_ExactDistribution):
    _dp: np.ndarray
    dpn: float
    dpw: float
    dpz: float

    @property
    def dp(self): return self._dp

    @property
    def dispersity(self):
        return self.dpw / self.dpn if self.dpn else float("nan")


@dataclass(frozen=True)
class MassDistribution(_ExactDistribution):
    _mass: np.ndarray
    mass_model: str
    mn: float
    mw: float
    mz: float

    @property
    def mass(self): return self._mass

    @property
    def dispersity(self):
        return self.mw / self.mn if self.mn else float("nan")


@dataclass(frozen=True)
class MolarMassDistribution:
    """Reconstructed differential molar-mass distribution.

    The v1 contract is the mass-weighted density with respect to log10(M):
    dW/dlog10(M), normalized to unit area.  It is derived from the exact kMC
    population and is not the exact discrete mass distribution.
    """

    _x: np.ndarray
    _y: np.ndarray
    total_chains: int
    snapshot_id: int | None
    t: float | None
    mass_model: str
    mn: float
    mw: float
    mz: float
    metadata: dict[str, Any]

    def __post_init__(self):
        self._x.flags.writeable = False
        self._y.flags.writeable = False

    @property
    def x(self): return self._x

    @property
    def y(self): return self._y

    @property
    def log10_mass(self): return self._x

    @property
    def log10_x(self): return self._x

    @property
    def mass(self): return np.power(10.0, self._x)

    @property
    def weighting(self): return "mass"

    @property
    def coordinate(self): return "log10"

    @property
    def representation(self): return "density"

    @property
    def dispersity(self):
        return self.mw / self.mn if self.mn else float("nan")

    @property
    def meta(self): return self.metadata

    @property
    def is_empty(self): return self._x.size == 0

    def to_tsv(self, path: str | Path):
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="") as handle:
            handle.write("# ordinate: dW/dlog10(M)\n")
            handle.write("log10_mass\tdW_dlog10M\n")
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
        ax.set_xlabel("log10(M / g mol$^{-1}$)")
        ax.set_ylabel("dW/dlog10(M)")
        return ax

    def info(self):
        return (
            f"{type(self).__name__}\n"
            f"  weighting: mass\n"
            f"  coordinate: log10(M)\n"
            f"  representation: density\n"
            f"  reconstruction: linear interpolation\n"
            f"  total_chains: {self.total_chains}"
        )

    def help(self): return self.info()


def build_cld(population, *, weighting="number", mass_model: str | None = None):
    weighting = _validate_weighting(weighting)
    dp_record = np.asarray(population.dp, dtype=float)
    count_record = np.asarray(population.count, dtype=float)
    if np.any(dp_record <= 0):
        raise ValueError("CLD requires strictly positive degree of polymerization")

    unique_dp, inverse = np.unique(dp_record, return_inverse=True)
    grouped_count = np.zeros(unique_dp.size, dtype=float)
    np.add.at(grouped_count, inverse, count_record)

    resolved_model = resolve_mass_model(population, mass_model)
    if weighting == "mass":
        mass_record, resolved_model = record_masses(population, resolved_model)
        grouped_mass = np.zeros(unique_dp.size, dtype=float)
        np.add.at(grouped_mass, inverse, count_record * mass_record)
        raw = grouped_mass
    elif weighting == "number":
        raw = grouped_count
    else:
        raw = unique_dp * unique_dp * grouped_count

    denom = float(np.sum(raw))
    y = raw / denom if denom else raw
    source_moments = population.moments(mass_model=resolved_model)
    meta = {
        "kind": "CLD",
        "weighting": weighting,
        "coordinate": "DP",
        "representation": "discrete",
        "normalization": "sum",
        "units": "dimensionless",
        "mass_model": resolved_model,
    }
    return ChainLengthDistribution(
        _readonly(unique_dp), _readonly(y), weighting, int(np.sum(count_record)),
        int(population.snapshot_id), population.t, meta,
        _readonly(unique_dp), source_moments.dpn, source_moments.dpw, source_moments.dpz,
    )


def build_mass_distribution(population, *, weighting="mass", mass_model: str | None = None):
    weighting = _validate_weighting(weighting)
    mass_record, model = record_masses(population, mass_model)
    count_record = np.asarray(population.count, dtype=float)
    unique_mass, inverse = np.unique(mass_record, return_inverse=True)
    grouped_count = np.zeros(unique_mass.size, dtype=float)
    np.add.at(grouped_count, inverse, count_record)

    if weighting == "number":
        raw = grouped_count
    elif weighting == "mass":
        raw = unique_mass * grouped_count
    else:
        raw = unique_mass * unique_mass * grouped_count
    denom = float(np.sum(raw))
    y = raw / denom if denom else raw

    source_moments = population.moments(mass_model=model)
    meta = {
        "kind": "mass_distribution",
        "weighting": weighting,
        "coordinate": "M",
        "representation": "discrete",
        "normalization": "sum",
        "units": "dimensionless",
        "mass_model": model,
    }
    return MassDistribution(
        _readonly(unique_mass), _readonly(y), weighting, int(np.sum(count_record)),
        int(population.snapshot_id), population.t, meta,
        _readonly(unique_mass), model, source_moments.mn, source_moments.mw, source_moments.mz,
    )


def _is_homopolymer_population(population) -> bool:
    """Return True only when the run metadata explicitly identifies homo kinetics."""
    run = getattr(population, "run", None)
    if run is None:
        return False
    kinetic_model = str(getattr(run, "kinetic_model", "") or "").lower()
    engine = str(getattr(run, "engine", "") or "").lower()
    return kinetic_model in {"homo", "homopolymer"} or "homo" in engine


def _homo_zero_filled_support(
    population, mass_record: np.ndarray, count_record: np.ndarray, *, mass_model: str
):
    """Return the natural zero-filled homo DP lattice when its mass law is explicit.

    The homo/copo decision comes from run metadata, not from the observed chain
    composition.  Missing DP masses are calculated from the declared mass model.
    If a single DP maps to multiple masses (for example because end groups vary),
    this function deliberately falls back to the exact-mass path.
    """
    if not _is_homopolymer_population(population):
        return None

    dp = np.asarray(population.dp, dtype=np.int64)
    unique_dp = np.unique(dp)
    if unique_dp.size < 2:
        return None

    masses_by_dp = []
    counts_by_dp = []
    for k in unique_dp:
        mask = dp == k
        values = np.asarray(mass_record[mask], dtype=float)
        if not np.allclose(values, values[0], rtol=1e-12, atol=1e-10):
            return None
        masses_by_dp.append(float(values[0]))
        counts_by_dp.append(float(np.sum(count_record[mask])))

    run = population.run
    monomer_entries = getattr(run, "dictionaries", {}).get("monomers", {})
    if len(monomer_entries) != 1:
        return None
    monomer = next(iter(monomer_entries.values()))
    repeat_mass = monomer.get("molar_mass_increment", monomer.get("molar_mass"))
    if repeat_mass is None:
        return None
    repeat_mass = float(repeat_mass)

    if mass_model == "repeat_units":
        offset = 0.0
    else:
        observed_offset = np.asarray(masses_by_dp) - unique_dp.astype(float) * repeat_mass
        if not np.allclose(observed_offset, observed_offset[0], rtol=1e-12, atol=1e-10):
            return None
        offset = float(observed_offset[0])

    # Confirm that the explicit chemical mass law reproduces occupied states.
    expected = unique_dp.astype(float) * repeat_mass + offset
    if not np.allclose(expected, masses_by_dp, rtol=1e-10, atol=1e-8):
        return None

    full_dp = np.arange(int(unique_dp.min()), int(unique_dp.max()) + 1, dtype=np.int64)
    full_mass = full_dp.astype(float) * repeat_mass + offset
    if np.any(full_mass <= 0):
        return None
    full_count = np.zeros(full_dp.size, dtype=float)
    full_count[unique_dp - full_dp[0]] = np.asarray(counts_by_dp, dtype=float)
    return full_mass, full_count

def _mc_polymer_density(source_mass: np.ndarray, source_count: np.ndarray, *, step: float):
    source_mass = np.asarray(source_mass, dtype=float)
    source_count = np.asarray(source_count, dtype=float)
    if source_mass.size < 2:
        raise ValueError(
            "MWD density requires at least two source support points; "
            "use mass_distribution() for a discrete/oligomeric population"
        )
    order = np.argsort(source_mass)
    source_mass = source_mass[order]
    source_count = source_count[order]
    u = np.log10(source_mass)
    h = source_count * source_mass * source_mass

    # mcPolymer uses a global regular log10(M) grid with spacing 0.005.
    first = int(np.ceil(u[0] / step))
    last = int(np.floor(u[-1] / step))
    if last - first < 1:
        raise ValueError(
            "MWD source spans less than two reconstruction grid points; "
            "use mass_distribution() for the exact discrete distribution"
        )
    grid = np.arange(first, last + 1, dtype=float) * step
    raw = np.interp(grid, u, h, left=0.0, right=0.0)
    area = float(np.trapezoid(raw, grid))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("MWD reconstruction has zero or invalid area")
    return grid, raw / area


def build_mwd(population, *, mass_model: str | None = None):
    """Build dW/dlog10(M) by mcPolymer-style linear interpolation.

    For a homopolymer with a single-valued affine M(DP), missing integer DP
    states are explicitly zero-filled before interpolation.  General/copolymer
    populations use the occupied exact-mass support without zero filling.
    """
    mass_record, model = record_masses(population, mass_model)
    count_record = np.asarray(population.count, dtype=float)

    homo = _homo_zero_filled_support(
        population, mass_record, count_record, mass_model=model
    )
    if homo is not None:
        source_mass, source_count = homo
        zero_filled = True
        source = "dp_counts"
    else:
        unique_mass, inverse = np.unique(mass_record, return_inverse=True)
        grouped_count = np.zeros(unique_mass.size, dtype=float)
        np.add.at(grouped_count, inverse, count_record)
        source_mass, source_count = unique_mass, grouped_count
        zero_filled = False
        source = "mass_counts"

    x, y = _mc_polymer_density(source_mass, source_count, step=_DEFAULT_LOG10_STEP)
    source_moments = population.moments(mass_model=model)
    meta = {
        "kind": "MWD",
        "weighting": "mass",
        "coordinate": "log10(M)",
        "representation": "density",
        "normalization": "integral",
        "ordinate": "dW/dlog10(M)",
        "units": "decade^-1",
        "reconstruction": "mcPolymer-style linear interpolation",
        "grid_step_log10M": _DEFAULT_LOG10_STEP,
        "source": source,
        "zero_filled": zero_filled,
        "mass_model": model,
    }
    return MolarMassDistribution(
        _readonly(x), _readonly(y), int(np.sum(count_record)),
        int(population.snapshot_id), population.t, model,
        source_moments.mn, source_moments.mw, source_moments.mz, meta,
    )


@dataclass(frozen=True)
class DistributionGroup:
    series: Mapping[str, Any]
    normalization: str
    kind: str
    series_disjoint: bool
    metadata: dict[str, Any]

    def __post_init__(self):
        object.__setattr__(self, "series", MappingProxyType(dict(self.series)))

    @property
    def series_names(self): return tuple(self.series)

    def __getitem__(self, name: str): return self.series[name]

    @property
    def x(self): raise ValueError("DistributionGroup has multiple supports; use group.series['name'].x")

    @property
    def y(self): raise ValueError("DistributionGroup has multiple series; use group.series['name'].y")

    @property
    def is_empty(self): return all(value.is_empty for value in self.series.values())


def _select_series_population(population, selector):
    from .chains import ChainPopulation
    if isinstance(selector, ChainPopulation): return selector
    if not isinstance(selector, str):
        raise TypeError("each series value must be a ChainPopulation or pool selector string")
    if selector == "all":
        value = population.all; return value() if callable(value) else value
    if selector in {"live", "active"}:
        value = population.live; return value() if callable(value) else value
    if selector == "dead":
        value = population.dead; return value() if callable(value) else value
    return population.pool(selector)


def _resolve_series(population, series):
    if isinstance(series, (tuple, list)):
        series = {str(name): str(name) for name in series}
    if not isinstance(series, Mapping) or not series:
        raise TypeError("series must be a non-empty mapping or tuple/list of pool selectors")
    resolved = {}
    for name, selector in series.items():
        selected = _select_series_population(population, selector)
        resolved[str(name)] = selected
    return resolved


def _row_signatures(population):
    raw = population._raw_arrays()
    if "chain_record_id" in raw:
        return {("chain_record_id", int(v)) for v in np.asarray(raw["chain_record_id"])}
    names = tuple(sorted(name for name in raw if np.asarray(raw[name]).ndim == 1))
    return {tuple((raw[name][i].item() if isinstance(raw[name][i], np.generic) else raw[name][i]) for name in names)
            for i in range(len(population))}


def _are_disjoint(populations):
    seen = set()
    for population in populations.values():
        current = _row_signatures(population)
        if seen.intersection(current): return False
        seen.update(current)
    return True


def _build_exact_series(population, *, series, kind, weighting, normalization, mass_model):
    weighting = _validate_weighting(weighting)
    normalization = str(normalization).lower()
    if normalization not in {"per_series", "combined"}:
        raise ValueError("normalization must be 'per_series' or 'combined'")
    populations = _resolve_series(population, series)
    disjoint = _are_disjoint(populations)
    if normalization == "combined" and not disjoint:
        raise ValueError("normalization='combined' requires pairwise-disjoint series")
    builder = build_cld if kind == "CLD" else build_mass_distribution
    built0 = {name: builder(p, weighting=weighting, mass_model=mass_model) for name, p in populations.items()}
    # derive exact source totals from each normalized distribution's denominator
    totals = {}
    for name, p in populations.items():
        count = np.asarray(p.count, dtype=float)
        if kind == "CLD":
            if weighting == "number": totals[name] = float(count.sum())
            elif weighting == "z":
                dp = np.asarray(p.dp, dtype=float); totals[name] = float(np.dot(count, dp * dp))
            else:
                mass, _ = record_masses(p, mass_model); totals[name] = float(np.dot(count, mass))
        else:
            mass, _ = record_masses(p, mass_model)
            if weighting == "number": totals[name] = float(count.sum())
            elif weighting == "mass": totals[name] = float(np.dot(count, mass))
            else: totals[name] = float(np.dot(count, mass * mass))
    total = float(sum(totals.values()))
    built = {}
    for name, value in built0.items():
        scale = totals[name] / total if normalization == "combined" and total else (0.0 if normalization == "combined" else 1.0)
        meta = dict(value.metadata); meta.update(normalization=normalization, series_name=name)
        built[name] = replace(value, _y=_readonly(value.y * scale), metadata=meta)
    return DistributionGroup(built, normalization, kind, disjoint, {
        "kind": kind, "weighting": weighting, "representation": "discrete",
        "normalization": normalization, "series_names": tuple(built), "series_disjoint": disjoint,
    })


def build_cld_series(population, *, series, weighting="number", normalization="per_series", mass_model=None):
    return _build_exact_series(population, series=series, kind="CLD", weighting=weighting,
                               normalization=normalization, mass_model=mass_model)



def build_mwd_series(population, *, series, normalization="per_series", mass_model=None):
    populations = _resolve_series(population, series)
    disjoint = _are_disjoint(populations)
    if normalization not in {"per_series", "combined"}:
        raise ValueError("normalization must be 'per_series' or 'combined'")
    if normalization == "combined" and not disjoint:
        raise ValueError("normalization='combined' requires pairwise-disjoint series")
    built = {name: build_mwd(p, mass_model=mass_model) for name, p in populations.items()}
    if normalization == "combined":
        totals = {}
        for name, p in populations.items():
            mass, _ = record_masses(p, mass_model)
            totals[name] = float(np.dot(np.asarray(p.count, dtype=float), mass))
        combined = float(sum(totals.values()))
        for name, value in list(built.items()):
            scale = totals[name] / combined if combined else 0.0
            meta = dict(value.metadata); meta.update(normalization="combined", series_name=name)
            built[name] = replace(value, _y=_readonly(value.y * scale), metadata=meta)
    return DistributionGroup(built, normalization, "MWD", disjoint, {
        "kind": "MWD", "weighting": "mass", "representation": "density",
        "normalization": normalization, "series_names": tuple(built), "series_disjoint": disjoint,
    })
