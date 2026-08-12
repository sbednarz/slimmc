from __future__ import annotations

"""Common CLD/MWD construction on :class:`ChainPopulation`.

The engine adapters only resolve a snapshot and (for MWD) its per-chain
masses.  All binning, smoothing, normalization, metadata, export and plotting
is implemented here once for homo- and copolymer results.
"""

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any
from collections.abc import Iterator, Mapping
import warnings

import numpy as np

from .chains import ChainPopulation
from .distribution import Distribution
from .table import Table


def _trapezoid(y, x) -> float:
    """Integrate on NumPy 1.24+ and NumPy 2.x without version branching."""
    integrate = getattr(np, "trapezoid", None)
    if integrate is None:
        integrate = getattr(np, "trapz", None)
    if integrate is None:
        # Reachable once a NumPy release removes both `trapezoid` (unlikely)
        # and the `trapz` alias (already true on recent NumPy 2.x) --
        # previously this looked up `np.trapz` unconditionally and let a
        # raw `AttributeError` escape straight from NumPy's module
        # `__getattr__`, which read like an internal bug rather than an
        # environment/compatibility problem.
        from .core import NumericalAnalysisError
        raise NumericalAnalysisError(
            "no compatible trapezoidal integration function found on this "
            "NumPy install (neither np.trapezoid nor np.trapz is available)"
        )
    with warnings.catch_warnings():
        # NumPy versions carrying both names may deprecate the alias.
        warnings.simplefilter("ignore", DeprecationWarning)
        return float(integrate(y, x))


def _readonly(values) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    result.flags.writeable = False
    return result


# Basis -> single-letter symbol used in the CLD/MWD descriptor (see
# _build_descriptor()). "z" gets its own symbol, Z -- it must never collapse
# onto "W" (mass), which was a real bug: basis="mass" and basis="mass" produced
# the identical descriptor string even though the underlying data differs
# (see _build_descriptor's docstring).
_BASIS_SYMBOL = {"number": "N", "mass": "W"}


def _build_descriptor(*, kind: str, basis: str, coordinate: str, output: str) -> str:
    """Build the canonical CLD/MWD descriptor string.

    This is the single source of truth for the string stored in
    ``metadata["descriptor"]`` -- which doubles as the plot/gnuplot y-axis
    label and the text table/CSV export comment. The notation contract is
    documented here and in docs/PYSLIMMC.md.

    Notation, by ``output``:

    - ``"density"``  -> ``d{SYMBOL}/d{AXIS}``            e.g. ``dW/dlog10(M)``
    - ``"fraction"`` -> ``{SYMBOL}_fraction({AXIS})``     e.g. ``W_fraction(log10(M))``
    - ``"amount"``   -> ``{SYMBOL}_abs({AXIS})``          e.g. ``W_abs(log10(M))``

    ``SYMBOL`` is basis-aware (N/W/Z, never conflating mass and z-weighted
    data). ``AXIS`` is ``DP``/``M`` for ``coordinate="linear"`` and exactly
    ``log10(DP)``/``log10(M)`` for ``coordinate="log10"`` -- always the full,
    unambiguous ``log10`` (never a bare ``"log"``, which could otherwise be
    misread as natural log). Only ``"density"`` uses differential (``d.../d...``)
    notation, matching Chen et al., Ind. Eng. Chem. Res. 2025, 64, 3695-3703,
    the mass-weighted density transformation for coordinate="log10"
    -- ``"fraction"``/``"amount"`` are deliberately NOT written as differentials
    since they are not densities (no division by a coordinate width).

    ``method`` (sticks/hist/kde/gaussian) is intentionally NOT part of this
    string -- it is a numerical-approximation detail, not part of the
    physical quantity's definition (``dW/dlog10(M)`` means the same thing
    whether it was estimated via a histogram or a KDE), and remains
    available separately as ``metadata["method"]``.
    """
    symbol = _BASIS_SYMBOL[basis]
    axis = "DP" if kind == "CLD" else "M"
    axis_expr = f"log10({axis})" if coordinate == "log10" else axis
    if output == "density":
        return f"d{symbol}/d{axis_expr}"
    if output == "fraction":
        return f"{symbol}_fraction({axis_expr})"
    if output == "amount":
        return f"{symbol}_abs({axis_expr})"
    raise ValueError(f"unknown output {output!r}")  # pragma: no cover -- _validate_common already restricts this


def select_pool(population: ChainPopulation, pool: str) -> ChainPopulation:
    if not isinstance(pool, str):
        raise TypeError("pool must be a string")
    if pool == "all":
        value = population.all
        return value() if callable(value) else value
    if pool in {"live", "active"}:
        value = population.live
        return value() if callable(value) else value
    if pool == "dead":
        value = population.dead
        return value() if callable(value) else value
    return population.select(pool=pool)


def _snapshot_meta(population: ChainPopulation, *, pool: str) -> dict[str, Any]:
    raw = population._raw_arrays()
    event = getattr(population, "kmc_event", None)
    if event is None and "kmc_event" in raw and len(raw["kmc_event"]):
        event = int(raw["kmc_event"][0])
    return {
        "snapshot_id": int(population.snapshot_id),
        "t": population.t,
        "kmc_event": event,
        "pool": pool,
        "selection": f"{population.total_chains} chains / {population.compressed_rows} compressed rows",
    }


def _validate_common(*, basis: str, method: str, coordinate: str,
                     output: str, normalization: str) -> tuple[str, str, str, str, str]:
    basis = str(basis).lower()
    if basis not in {"number", "mass"}:
        raise ValueError("basis must be 'number' or 'mass'")
    method = {"histogram": "hist"}.get(str(method).lower(), str(method).lower())
    if method not in {"sticks", "hist", "kde", "gaussian"}:
        raise ValueError("method must be 'sticks', 'hist', 'gaussian', or 'kde'")
    coordinate = {"dp": "linear", "mass": "linear"}.get(str(coordinate).lower(), str(coordinate).lower())
    if coordinate not in {"linear", "log10"}:
        raise ValueError("coordinate must be 'linear' or 'log10'")
    output = str(output).lower()
    if output not in {"amount", "fraction", "density"}:
        raise ValueError("output must be amount, fraction, or density")
    normalization = {"sum": "per_series", "fraction": "per_series"}.get(
        str(normalization).lower(), str(normalization).lower()
    )
    if normalization not in {"absolute", "per_series"}:
        raise ValueError("single-series normalization must be absolute or per_series")
    if output == "amount" and normalization != "absolute":
        raise ValueError("output='amount' requires normalization='absolute'")
    if output == "fraction" and normalization == "absolute":
        raise ValueError("output='fraction' requires a normalized denominator")
    if method == "kde" and output != "density":
        raise ValueError("method='kde' always requires output='density'")
    return basis, method, coordinate, output, normalization


def _representation(base: Distribution, *, method: str, coordinate: str,
                    bins: int | None, bin_width: float | None,
                    sigma: float | None, grid_step: float | None) -> Distribution:
    if not base.x:
        return Distribution(
            base.name, base.x_name, (), (), (), (), (), basis=base.basis,
            method=method, source_stats=base.source_stats,
            meta=dict(base.meta, method=method, coordinate=coordinate),
        )
    if method == "sticks":
        if bins is not None or bin_width is not None or sigma is not None:
            raise ValueError("method='sticks' does not accept bins, bin_width, or sigma")
        return base
    if bins is not None and bin_width is not None:
        raise ValueError("bins and bin_width are mutually exclusive")
    if method == "hist":
        if sigma is not None:
            raise ValueError("method='hist' does not accept sigma")
        if coordinate == "log10":
            return Distribution.from_histogram_log(
                base, bins=bins, bin_width=0.02 if bins is None and bin_width is None else bin_width
            )
        return Distribution.from_histogram(
            base.name, base.x_name, zip(base.x, base.count, base.x),
            bins=bins, bin_width=bin_width, basis=base.basis, meta=dict(base.meta),
        )
    if method == "gaussian":
        if coordinate == "log10":
            hist = Distribution.from_histogram_log(
                base, bins=bins, bin_width=0.01 if bins is None and bin_width is None else bin_width
            )
        else:
            hist = Distribution.from_histogram(
                base.name, base.x_name, zip(base.x, base.count, base.x),
                bins=bins, bin_width=bin_width, basis=base.basis, meta=dict(base.meta),
            )
        default_sigma = 0.04 if coordinate == "log10" else (
            float(hist.meta.get("bin_width", 1.0))
        )
        return hist.gaussian(sigma=float(default_sigma if sigma is None else sigma))
    if bins is not None or bin_width is not None:
        raise ValueError("method='kde' does not accept bins or bin_width")
    return base.kde(
        bandwidth=float(0.04 if sigma is None else sigma),
        transform="log10" if coordinate == "log10" else None,
        grid_step=grid_step,
    )


@dataclass(frozen=True)
class _PopulationDistribution:
    _distribution: Distribution
    _x: np.ndarray
    _y: np.ndarray
    metadata: dict[str, Any]

    def __post_init__(self) -> None:
        self._x.flags.writeable = False
        self._y.flags.writeable = False

    @property
    def x(self) -> np.ndarray:
        """Physical DP or molar mass, even when coordinate='log10'."""
        return self._x

    @property
    def log10_x(self) -> np.ndarray:
        if np.any(self._x <= 0):
            raise ValueError("log10_x requires strictly positive x values")
        return _readonly(np.log10(self._x))

    @property
    def y(self) -> np.ndarray:
        return self._y

    @property
    def meta(self) -> dict[str, Any]:
        return self.metadata

    @property
    def snapshot_id(self) -> int | None:
        value = self.metadata.get("snapshot_id")
        return None if value is None else int(value)

    @property
    def time(self) -> float | None:
        return self.metadata.get("t")

    @property
    def t(self) -> float | None:
        return self.time

    @property
    def kmc_event(self) -> int | None:
        return self.metadata.get("kmc_event")

    @property
    def pool(self) -> str:
        return str(self.metadata["pool"])

    @property
    def basis(self) -> str:
        return str(self.metadata["basis"])

    @property
    def method(self) -> str:
        return str(self.metadata["method"])

    @property
    def coordinate(self) -> str:
        return str(self.metadata["coordinate"])

    @property
    def normalization(self) -> str:
        return str(self.metadata["normalization"])

    @property
    def output(self) -> str:
        return str(self.metadata["output"])

    @property
    def is_empty(self) -> bool:
        return self._distribution.n <= 0

    @property
    def n(self) -> float:
        return self._distribution.n

    @property
    def total_weight(self) -> float:
        return self._distribution.total_weight

    @property
    def mn(self) -> float:
        return self._distribution.mn()

    @property
    def mw(self) -> float:
        return self._distribution.mw()

    @property
    def mz(self) -> float:
        return self._distribution.mz()

    @property
    def dispersity(self) -> float:
        return self._distribution.pdi()

    @property
    def data(self) -> Table:
        x_name = "dp" if self.metadata["kind"] == "CLD" else "mass"
        log_name = "log10_dp" if x_name == "dp" else "log10_mass"
        rows = []
        for x, y in zip(self.x, self.y):
            rows.append((float(x), math.log10(float(x)) if x > 0 else float("nan"), float(y)))
        return Table((x_name, log_name, "value"), rows, name=self.metadata["kind"])

    def as_table(self) -> Table:
        """Legacy full table; canonical x/y table is available as ``.data``."""
        return self._distribution.as_table()

    def as_dict(self) -> dict[str, Any]:
        result = self._distribution.to_dict()
        result.update({"metadata": dict(self.metadata), "x": self.x.tolist(), "y": self.y.tolist()})
        return result

    def _metadata_lines(self) -> list[str]:
        keys = ("kind", "descriptor", "snapshot_id", "t", "kmc_event", "pool",
                "basis", "method", "coordinate", "output", "normalization",
                "sigma", "mass_model", "normalize", "neutral_mass", "mz",
                "ionization_model")
        lines = ["# pyslimmc_export: 1", "# object: distribution"]
        for key in keys:
            value = self.metadata.get(key)
            if value is not None:
                lines.append(f"# {key}: {value}")
        lines.append(f"# integral: {self.metadata.get('integral', 0.0):.17g}")
        return lines

    def _write(self, path: str | Path, *, delimiter: str, metadata: str = "comments") -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        table = self.data
        with target.open("w", encoding="utf-8", newline="") as handle:
            if metadata == "comments":
                handle.write("\n".join(self._metadata_lines()) + "\n")
            handle.write(delimiter.join(table.columns) + "\n")
            for row in table:
                handle.write(delimiter.join(str(value) for value in row) + "\n")
        target.with_suffix(".meta.json").write_text(
            json.dumps(self.metadata, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8"
        )
        return target

    def to_tsv(self, path: str | Path, *, metadata: str = "comments", layout: str = "wide") -> Path:
        if layout not in {"wide", "long"}:
            raise ValueError("layout must be wide or long")
        return self._write(path, delimiter="\t", metadata=metadata)



    def plot(self, ax=None, *, path: str | Path | None = None, dpi: int = 300,
             xscale: str | None = None, yscale: str = "linear",
             normalize: bool | None = None, title: str | None = None,
             style: str = "screen", span: str | None = None, **kwargs):
        from .plotting import apply_axes_style, create_axes, require_owned_geometry, style_kwargs
        require_owned_geometry(ax, span)
        if ax is None:
            _, ax = create_axes(style, span=span)
        # When coordinate="log10", we log-transform x ourselves and hand
        # matplotlib a plain linear axis -- tick locations/labels then come
        # out as ordinary decade numbers (1, 2, 3, ...), not matplotlib's
        # own log-scale formatting (10^1, 10^2, ...) applied to physical x
        # at draw time. xscale="log" would double-log already-log10 data,
        # so it is rejected explicitly rather than silently misinterpreted.
        if self.coordinate == "log10":
            x_plot = np.log10(self.x)
            if xscale is None:
                xscale = "linear"
            elif xscale != "linear":
                raise ValueError(
                    "xscale must be 'linear' (or omitted) when coordinate='log10' -- "
                    "x is already log10-transformed here, so a matplotlib log-scale "
                    "axis on top of it would double-log. Build with coordinate='linear' "
                    "instead if you want physical x on a log-scaled matplotlib axis."
                )
        else:
            x_plot = self.x
            if xscale is None:
                xscale = "linear"
        y = self.y
        if normalize is True:
            area = _trapezoid(y, np.log10(self.x) if self.coordinate == "log10" else self.x) \
                if self.method == "kde" and len(y) > 1 else float(np.sum(y))
            y = y / area if area > 0 else y
        elif normalize is False:
            y = _readonly(self._distribution.y)
        defaults = style_kwargs(style)
        defaults.update(kwargs)
        if self.method == "sticks":
            ax.vlines(x_plot, 0, y, **defaults)
        else:
            ax.plot(x_plot, y, **defaults)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        if self.metadata["kind"] == "CLD":
            ax.set_xlabel("log10(DP)" if self.coordinate == "log10" else "DP")
        else:
            ax.set_xlabel("log10(M)" if self.coordinate == "log10" else "M / g mol$^{-1}$")
        ax.set_ylabel(self.metadata["descriptor"])
        if title is not None:
            ax.set_title(title)
        apply_axes_style(ax, style)
        if path is not None:
            ax.figure.savefig(path, dpi=dpi)
        return ax

    def info_text(self) -> str:
        kind = self.metadata.get("kind")
        title = {
            "MWD": "Molar mass distribution",
            "CLD": "Chain length distribution",
            "chain_mass_spectrum": "Neutral chain mass spectrum",
        }.get(kind, str(kind or "Distribution"))
        final = self.metadata.get("is_final")
        lines = [title, "-" * len(title), "Snapshot:",
                 f"  id:              {self.snapshot_id}",
                 f"  time:            {self.time if self.time is not None else 'not available'}",
                 f"  final:           {'yes' if final else 'no' if final is not None else 'not available'}",
                 "", "Population:",
                 f"  pool:            {self.pool}",
                 f"  physical chains: {self.n:g}",
                 "", "Representation:",
                 f"  method:          {self.method}",
                 f"  basis:           {self.basis}",
                 f"  coordinate:      {self.coordinate}",
                 f"  output:          {self.output}",
                 f"  normalization:   {self.normalization}"]
        for key, label, unit in (("bin_width", "bin width", " decade" if self.coordinate == "log10" else ""),
                                 ("sigma", "sigma", " decade" if self.coordinate == "log10" else ""),
                                 ("grid_step", "grid step", " decade" if self.coordinate == "log10" else "")):
            value = self.metadata.get(key)
            if value is not None:
                lines.append(f"  {label + ':':17} {value:g}{unit}")
        lines += ["", "Axes:"]
        if kind == "CLD":
            lines += ["  x:               degree of polymerization",
                      "  log10_x:         log10(x)"]
        elif kind == "chain_mass_spectrum":
            lines += ["  mass:            neutral chain mass, g/mol",
                      "  intensity:       relative or absolute chain count"]
        else:
            lines += ["  x:               physical molar mass, g/mol",
                      "  log10_x:         log10(x)"]
        lines.append(f"  y:               {self.metadata.get('descriptor', 'distribution value')}")
        if kind == "CLD":
            lines += ["", "Moments:", f"  DPn:             {self.dp_n:.6g}",
                      f"  DPw:             {self.dp_w:.6g}", f"  DPz:             {self.dp_z:.6g}",
                      f"  dispersity:      {self.dispersity:.6g}"]
        elif kind == "MWD":
            lines += ["", "Moments:", f"  Mn:              {self.mn:.6g} g/mol",
                      f"  Mw:              {self.mw:.6g} g/mol", f"  Mz:              {self.mz:.6g} g/mol",
                      f"  dispersity:      {self.dispersity:.6g}"]
        else:
            lines += ["", "Peak:", f"  unique masses:   {len(self.mass)}",
                      f"  base peak mass:  {self.base_peak_mass}",
                      f"  base peak value: {self.base_peak_intensity:.6g}", "",
                      "Note: neutral masses only; this is not an m/z spectrum."]
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text

    def __repr__(self) -> str:
        kind = self.metadata.get("kind", type(self).__name__)
        if kind == "CLD":
            extra = f", DPn={self.dp_n:.6g}, DPw={self.dp_w:.6g}"
        elif kind == "MWD":
            extra = f", Mn={self.mn:.6g}, Mw={self.mw:.6g}"
        else:
            extra = f", peaks={len(self.x)}"
        return (f"{type(self).__name__}(method={self.method!r}, basis={self.basis!r}, "
                f"coordinate={self.coordinate!r}, points={len(self.x)}{extra})")

    def __getattr__(self, name: str):
        if name in {"pdi", "dpn", "dpw", "dpz", "to_dict"}:
            raise AttributeError(name)
        return getattr(self._distribution, name)


class ChainLengthDistribution(_PopulationDistribution):
    @property
    def dp_n(self) -> float:
        return self.mn

    @property
    def dp_w(self) -> float:
        return self.mw

    @property
    def dp_z(self) -> float:
        return self.mz


class MolarMassDistribution(_PopulationDistribution):
    pass


class ChainMassSpectrum(_PopulationDistribution):
    """Exact neutral-chain mass spectrum; not an m/z or ionization model."""

    @property
    def mass(self) -> np.ndarray:
        return self.x

    @classmethod
    def from_table(cls, table: Table, *, normalize: str = "base_peak",
                   snapshot_id: int | None = None) -> "ChainMassSpectrum":
        return _chain_mass_spectrum_from_table(table, normalize=normalize, snapshot_id=snapshot_id)

    @property
    def intensity(self) -> np.ndarray:
        return self.y

    @property
    def base_peak_mass(self) -> float | None:
        if self.is_empty or not np.any(self.y > 0):
            return None
        return float(self.x[int(np.argmax(self.y))])

    @property
    def base_peak_intensity(self) -> float:
        return float(self.y.max()) if self.y.size else 0.0


class DistributionSeries(Mapping[str, _PopulationDistribution]):
    """Read-only named views used by a multi-series CLD/MWD."""

    def __init__(self, values: Mapping[str, _PopulationDistribution]):
        self._values = dict(values)

    def __getitem__(self, name: str) -> _PopulationDistribution:
        return self._values[name]

    def __iter__(self) -> Iterator[str]:
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)


@dataclass(frozen=True)
class MultiDistribution:
    """Several distributions from one snapshot on one common x grid."""

    series: DistributionSeries
    metadata: dict[str, Any]

    @property
    def series_names(self) -> tuple[str, ...]:
        return tuple(self.series)

    @property
    def x(self) -> np.ndarray:
        if not self.series:
            return _readonly([])
        return next(iter(self.series.values())).x

    @property
    def log10_x(self) -> np.ndarray:
        if not self.series:
            return _readonly([])
        return next(iter(self.series.values())).log10_x

    @property
    def y(self):
        raise ValueError("Distribution contains multiple series. Use dist.series['name'].y.")

    @property
    def meta(self) -> dict[str, Any]:
        return self.metadata

    @property
    def snapshot_id(self) -> int | None:
        value = self.metadata.get("snapshot_id")
        return None if value is None else int(value)

    @property
    def time(self) -> float | None:
        return self.metadata.get("t")

    @property
    def kmc_event(self) -> int | None:
        return self.metadata.get("kmc_event")

    @property
    def basis(self) -> str:
        return str(self.metadata["basis"])

    @property
    def method(self) -> str:
        return str(self.metadata["method"])

    @property
    def coordinate(self) -> str:
        return str(self.metadata["coordinate"])

    @property
    def normalization(self) -> str:
        return str(self.metadata["normalization"])

    @property
    def is_empty(self) -> bool:
        return all(value.is_empty for value in self.series.values())

    @property
    def data(self) -> Table:
        x_name = "dp" if self.metadata["kind"] == "CLD" else "mass"
        log_name = "log10_dp" if x_name == "dp" else "log10_mass"
        columns = (x_name, log_name, *self.series_names)
        rows = []
        for i, x in enumerate(self.x):
            rows.append((float(x), math.log10(float(x)) if x > 0 else float("nan"),
                         *(float(self.series[name].y[i]) for name in self.series_names)))
        return Table(columns, rows, name=self.metadata["kind"])

    def as_table(self) -> Table:
        return self.data

    def _metadata_lines(self) -> list[str]:
        lines = ["# pyslimmc_export: 1", "# object: multi_distribution"]
        for key in ("kind", "descriptor", "snapshot_id", "t", "kmc_event", "basis",
                    "method", "coordinate", "output", "normalization", "reference", "sigma",
                    "normalize", "neutral_mass", "mz", "ionization_model"):
            value = self.metadata.get(key)
            if value is not None:
                lines.append(f"# {key}: {value}")
        lines.append("# series: " + ", ".join(self.series_names))
        for name in self.series_names:
            lines.append(f"# integral.{name}: {self.series[name].metadata.get('integral', 0.0):.17g}")
        return lines

    def _write_meta(self, target: Path) -> None:
        target.with_suffix(".meta.json").write_text(
            json.dumps(self.metadata, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8"
        )

    def to_tsv(self, path: str | Path, *, metadata: str = "comments", layout: str = "wide") -> Path:
        if layout not in {"wide", "long"}:
            raise ValueError("layout must be wide or long")
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="") as handle:
            if metadata == "comments":
                handle.write("\n".join(self._metadata_lines()) + "\n")
            if layout == "wide":
                table = self.data
                handle.write("\t".join(table.columns) + "\n")
                for row in table:
                    handle.write("\t".join(str(value) for value in row) + "\n")
            else:
                handle.write("series\tx\tlog10_x\tvalue\n")
                for name in self.series_names:
                    for x, y in zip(self.x, self.series[name].y):
                        handle.write(f"{name}\t{float(x)}\t{math.log10(float(x)) if x > 0 else 'nan'}\t{float(y)}\n")
        self._write_meta(target)
        return target


    def plot(self, ax=None, *, path: str | Path | None = None, dpi: int = 300,
             mode: str = "overlay", xscale: str | None = None, yscale: str = "linear",
             display_normalization: str | None = None, title: str | None = None,
             styles: Mapping[str, Mapping[str, Any]] | None = None,
             style: str = "screen", span: str | None = None):
        from .plotting import apply_axes_style, create_axes, require_owned_geometry, style_kwargs
        if mode not in {"overlay", "stacked", "grouped"}:
            raise ValueError("mode must be overlay, stacked, or grouped")
        if mode == "stacked" and self.normalization == "per_series":
            raise ValueError("stacked mode is not defined for per_series normalization")
        if mode == "stacked" and not self.metadata.get("series_disjoint", False):
            raise ValueError("stacked mode requires disjoint series")
        if display_normalization not in {None, "peak_each", "peak_global", "percent_peak_each"}:
            raise ValueError("unknown display_normalization")
        require_owned_geometry(ax, span)
        if ax is None:
            _, ax = create_axes(style, span=span)
        # See _PopulationDistribution.plot() for the rationale: pre-transform
        # x ourselves when coordinate="log10" and hand matplotlib a plain
        # linear axis, instead of letting it log-scale physical x at draw
        # time.
        if self.coordinate == "log10":
            x_plot = np.log10(self.x)
            if xscale is None:
                xscale = "linear"
            elif xscale != "linear":
                raise ValueError(
                    "xscale must be 'linear' (or omitted) when coordinate='log10' -- "
                    "x is already log10-transformed here, so a matplotlib log-scale "
                    "axis on top of it would double-log. Build with coordinate='linear' "
                    "instead if you want physical x on a log-scaled matplotlib axis."
                )
        else:
            x_plot = self.x
            if xscale is None:
                xscale = "linear"
        arrays = {name: self.series[name].y.copy() for name in self.series_names}
        if display_normalization in {"peak_each", "percent_peak_each"}:
            target = 100.0 if display_normalization == "percent_peak_each" else 1.0
            arrays = {name: values * target / values.max() if values.size and values.max() > 0 else values
                      for name, values in arrays.items()}
        elif display_normalization == "peak_global":
            peak = max((float(values.max()) for values in arrays.values() if values.size), default=0.0)
            if peak > 0:
                arrays = {name: values / peak for name, values in arrays.items()}
        if mode == "stacked":
            from .plotting import get_style
            palette = get_style(style).palette
            ax.stackplot(x_plot, *(arrays[name] for name in self.series_names),
                         labels=self.series_names,
                         colors=[palette[i % len(palette)] for i in range(len(self.series_names))])
        elif mode == "grouped":
            if self.method not in {"sticks", "hist"}:
                raise ValueError("grouped mode is only available for sticks or hist")
            width = (np.min(np.diff(x_plot)) if len(x_plot) > 1 else 1.0) / max(1, len(self.series))
            start = -(len(self.series) - 1) * width / 2.0
            for i, name in enumerate(self.series_names):
                ax.bar(x_plot + start + i * width, arrays[name], width=width,
                       label=name, **{**style_kwargs(style, index=i),
                                      **dict((styles or {}).get(name, {}))})
        else:
            for i, name in enumerate(self.series_names):
                kwargs = {**style_kwargs(style, index=i), **dict((styles or {}).get(name, {}))}
                if self.method == "sticks":
                    ax.vlines(x_plot, 0, arrays[name], label=name, **kwargs)
                else:
                    ax.plot(x_plot, arrays[name], label=name, **kwargs)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        if self.metadata["kind"] == "CLD":
            ax.set_xlabel("log10(DP)" if self.coordinate == "log10" else "DP")
        else:
            ax.set_xlabel("log10(M)" if self.coordinate == "log10" else "M / g mol$^{-1}$")
        ax.set_ylabel(self.metadata["descriptor"])
        if title is not None:
            ax.set_title(title)
        ax.legend()
        apply_axes_style(ax, style)
        if path is not None:
            ax.figure.savefig(path, dpi=dpi)
        return ax


    def info_text(self) -> str:
        return "\n".join([
            f"distribution: {self.metadata['kind']}", f"series: {', '.join(self.series_names)}",
            f"snapshot_id: {self.snapshot_id}", f"basis: {self.basis}", f"method: {self.method}",
            f"coordinate: {self.coordinate}", f"normalization: {self.normalization}",
            "", "Common next steps:", "  spectrum.plot()", "  spectrum.plot(mode=\"overlay\")",
            "  spectrum.to_tsv(\"spectrum.tsv\")",
        ])

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text


class MultiChainLengthDistribution(MultiDistribution):
    pass


class MultiMolarMassDistribution(MultiDistribution):
    pass


class MultiChainMassSpectrum(MultiDistribution):
    def plot(self, *args, mode: str = "overlay", **kwargs):
        if mode == "stacked" and self.metadata.get("normalize") != "count":
            raise ValueError("stacked chain spectra require normalize='count'")
        return super().plot(*args, mode=mode, **kwargs)



def _canonical_y(dist: Distribution, *, coordinate: str, output: str,
                 normalization: str) -> tuple[np.ndarray, float]:
    y = np.asarray(dist.y, dtype=float).copy()
    if dist.method == "kde":
        total = (dist.n if dist.basis == "number" else
                 dist.source_stats.get("total_mass_weight", 0.0))
    else:
        total = float(np.sum(y))
    # 'sticks' is the one representation whose x stays physical/linear
    # regardless of `coordinate` (see _representation(): method="sticks"
    # returns `base` untouched) -- hist/kde/gaussian already carry a
    # genuinely log10-transformed `dist.x` here when coordinate="log10",
    # so only sticks needs this explicit correction. Before this fix,
    # requesting output="density" with coordinate="log10" on a sticks
    # distribution silently returned the exact same numbers as
    # coordinate="linear" (dW/dM, not dW/dlog10(M)) even though the
    # descriptor/metadata claimed otherwise.
    x_axis = np.log10(np.asarray(dist.x, dtype=float)) if (
        dist.method == "sticks" and coordinate == "log10"
    ) else np.asarray(dist.x, dtype=float)
    if output == "density" and dist.method != "kde" and y.size:
        if len(x_axis) > 1:
            widths = np.gradient(x_axis)
            y /= widths
    denominator = total if normalization == "per_series" else 1.0
    if normalization == "per_series" and denominator > 0:
        y /= denominator
    if output == "fraction" and normalization == "absolute" and total > 0:
        y /= total
    integral = (_trapezoid(y, x_axis) if output == "density" and len(y) > 1
                else float(np.sum(y)))
    return _readonly(y), integral


def build_cld(population: ChainPopulation, *, masses=None, mass_model: str | None = None,
              pool: str = "all", basis: str = "number",
              method: str = "sticks", coordinate: str = "linear", output: str = "fraction",
              normalization: str = "per_series", bins: int | None = None,
              bin_width: float | None = None, sigma: float | None = None,
              grid_step: float | None = None) -> ChainLengthDistribution:
    basis, method, coordinate, output, normalization = _validate_common(
        basis=basis, method=method, coordinate=coordinate, output=output, normalization=normalization
    )
    selected = select_pool(population, pool)
    dp_pairs = [(float(dp), float(n), float(dp)) for dp, n in zip(selected.dp, selected.count)]
    meta = dict(_snapshot_meta(selected, pool=pool), kind="CLD", mass_model=mass_model)
    moment_base = Distribution.from_pairs("CLD", "DP", dp_pairs, basis=basis, method="sticks", meta=meta)
    base = moment_base
    if basis == "mass":
        all_masses = np.asarray(masses, dtype=float)
        if len(all_masses) != len(population):
            raise ValueError("mass-basis CLD requires one molar mass per compressed population row")
        if selected is population:
            selected_masses = all_masses
        else:
            raw = selected._raw_arrays()
            if "mass" in raw:
                selected_masses = np.asarray(raw["mass"], dtype=float)
            elif hasattr(selected, "masses"):
                selected_masses = np.asarray(selected.masses(mass_model=mass_model or "repeat_units"), dtype=float)
            else:
                raise ValueError("cannot map molar masses after CLD population selection")
        weight_pairs = [(float(dp), float(n), float(mass))
                        for dp, n, mass in zip(selected.dp, selected.count, selected_masses)]
        base = Distribution.from_pairs("CLD", "DP", weight_pairs, basis=basis, method="sticks", meta=meta)
    represented = _representation(base, method=method, coordinate=coordinate, bins=bins,
                                  bin_width=bin_width, sigma=sigma, grid_step=grid_step)
    if basis == "mass":
        stats = dict(moment_base.source_stats)
        stats["total_mass_weight"] = base.source_stats["total_mass_weight"]
        stats["total_z_weight"] = base.source_stats["total_z_weight"]
        stats["total_z2_weight"] = base.source_stats["total_z2_weight"]
        represented = Distribution(
            represented.name, represented.x_name, represented.x, represented.count,
            represented.mass_weight, represented.z_weight, represented.z2_weight,
            basis=represented.basis, method=represented.method,
            source_stats=stats, meta=represented.meta,
        )
    y, integral = _canonical_y(represented, coordinate=coordinate, output=output,
                               normalization=normalization)
    axis = np.asarray(represented.x, dtype=float)
    physical_x = np.power(10.0, axis) if coordinate == "log10" and method != "sticks" else axis
    descriptor = _build_descriptor(kind="CLD", basis=basis, coordinate=coordinate, output=output)
    metadata = dict(meta)
    metadata.update(represented.meta)
    metadata.update(basis=basis, method=method, coordinate=coordinate, output=output,
                    normalization=normalization, sigma=sigma, descriptor=descriptor, integral=integral)
    return ChainLengthDistribution(represented, _readonly(physical_x), y, metadata)


def build_mwd(population: ChainPopulation, masses, *, pool: str = "all", basis: str = "mass",
              method: str = "gaussian", coordinate: str = "log10", output: str = "density",
              normalization: str = "per_series", bins: int | None = None,
              bin_width: float | None = None, sigma: float | None = None,
              grid_step: float | None = None, mass_model: str | None = None) -> MolarMassDistribution:
    """Build the canonical MWD.

    Default ``method="gaussian"`` (log10-binned histogram + Gaussian
    smoothing), not ``"kde"``. A raw KDE's +-5*sigma tails in log10-space
    extend the reported mass range well beyond the actually-simulated
    chains (on a typical fixture: real range [600.7, 130856.8] Da vs. a
    KDE-reported range of [379.0, 207334.7] Da at the default sigma=0.04
    decades) -- 'gaussian' never leaves the original histogram's bin
    range, so it does not imply chains outside what was actually observed.
    ``method="kde"`` remains fully available for callers who want a smooth
    nonparametric estimate and accept that tradeoff.
    """
    basis, method, coordinate, output, normalization = _validate_common(
        basis=basis, method=method, coordinate=coordinate, output=output, normalization=normalization
    )
    selected = select_pool(population, pool)
    all_masses = np.asarray(masses, dtype=float)
    if len(all_masses) != len(population):
        raise ValueError("masses must contain one value per compressed ChainPopulation row")
    if selected is population:
        selected_masses = all_masses
    else:
        # Resolve after selection from the selected population whenever possible.
        raw = selected._raw_arrays()
        if "mass" in raw:
            selected_masses = np.asarray(raw["mass"], dtype=float)
        elif hasattr(selected, "masses"):
            selected_masses = np.asarray(selected.masses(mass_model=mass_model or "repeat_units"), dtype=float)
        else:
            raise ValueError("cannot map molar masses after population selection")
    pairs = [(float(m), float(n), float(m)) for m, n in zip(selected_masses, selected.count)]
    meta = dict(_snapshot_meta(selected, pool=pool), kind="MWD", mass_model=mass_model)
    base = Distribution.from_pairs("MWD", "M / g mol-1", pairs, basis=basis, method="sticks", meta=meta)
    represented = _representation(base, method=method, coordinate=coordinate, bins=bins,
                                  bin_width=bin_width, sigma=sigma, grid_step=grid_step)
    y, integral = _canonical_y(represented, coordinate=coordinate, output=output,
                               normalization=normalization)
    axis = np.asarray(represented.x, dtype=float)
    physical_x = np.power(10.0, axis) if coordinate == "log10" and method != "sticks" else axis
    descriptor = _build_descriptor(kind="MWD", basis=basis, coordinate=coordinate, output=output)
    metadata = dict(meta)
    metadata.update(represented.meta)
    metadata.update(basis=basis, method=method, coordinate=coordinate, output=output,
                    normalization=normalization, sigma=sigma, descriptor=descriptor, integral=integral)
    return MolarMassDistribution(represented, _readonly(physical_x), y, metadata)


def _resolve_series(population: ChainPopulation, series) -> dict[str, ChainPopulation]:
    if isinstance(series, (tuple, list)):
        series = {str(name): str(name) for name in series}
    if not isinstance(series, Mapping) or not series:
        raise TypeError("series must be a non-empty mapping or tuple of pool selectors")
    resolved: dict[str, ChainPopulation] = {}
    for name, selector in series.items():
        if not isinstance(name, str) or not name:
            raise ValueError("series names must be non-empty strings")
        if isinstance(selector, ChainPopulation):
            selected = selector
        elif isinstance(selector, str):
            selected = select_pool(population, selector)
        else:
            raise TypeError("each series value must be a ChainPopulation or pool selector string")
        if int(selected.snapshot_id) != int(population.snapshot_id):
            raise ValueError("all series must come from the same snapshot_id")
        if selected.t != population.t:
            raise ValueError("all series must come from the same saved time")
        resolved[name] = selected
    return resolved


def _row_signatures(population: ChainPopulation) -> set[tuple[Any, ...]]:
    raw = population._raw_arrays()
    names = tuple(sorted(raw))
    signatures = set()
    for i in range(len(population)):
        values = []
        for name in names:
            value = raw[name][i]
            values.append(value.item() if isinstance(value, np.generic) else value)
        signatures.add(tuple(values))
    return signatures


def _are_disjoint(populations: Mapping[str, ChainPopulation]) -> bool:
    seen: set[tuple[Any, ...]] = set()
    for population in populations.values():
        current = _row_signatures(population)
        if seen.intersection(current):
            return False
        seen.update(current)
    return True


def _common_coordinate_grid(values: Mapping[str, _PopulationDistribution], *, method: str,
                            coordinate: str) -> np.ndarray:
    axes = [np.log10(value.x) if coordinate == "log10" else value.x
            for value in values.values() if value.x.size]
    if not axes:
        return _readonly([])
    if method == "kde":
        steps = [float(np.min(np.diff(axis))) for axis in axes if len(axis) > 1]
        step = min(steps) if steps else 1.0
        lo = min(float(axis.min()) for axis in axes)
        hi = max(float(axis.max()) for axis in axes)
        return _readonly(np.arange(lo, hi + step * 0.5, step))
    return _readonly(sorted(set(float(x) for axis in axes for x in axis)))


def _align_y(value: _PopulationDistribution, grid: np.ndarray, *, method: str,
             coordinate: str) -> np.ndarray:
    if not value.x.size:
        return np.zeros(grid.size, dtype=float)
    axis = np.log10(value.x) if coordinate == "log10" else value.x
    if method == "kde":
        return np.interp(grid, axis, value.y, left=0.0, right=0.0)
    lookup = {round(float(x), 13): float(y) for x, y in zip(axis, value.y)}
    return np.asarray([lookup.get(round(float(x), 13), 0.0) for x in grid], dtype=float)


def _series_total(y: np.ndarray, grid: np.ndarray, *, output: str) -> float:
    if output == "density" and len(grid) > 1:
        return _trapezoid(y, grid)
    return float(np.sum(y))


def _normalize_multi(raw: Mapping[str, _PopulationDistribution], populations: Mapping[str, ChainPopulation],
                     *, kind: str, basis: str, method: str, coordinate: str, output: str,
                     normalization: str, reference: str | None, sigma: float | None):
    if normalization not in {"absolute", "per_series", "combined", "reference"}:
        raise ValueError("normalization must be absolute, per_series, combined, or reference")
    disjoint = _are_disjoint(populations)
    if normalization == "combined" and not disjoint:
        raise ValueError("normalization='combined' requires pairwise-disjoint series")
    if normalization == "reference":
        if reference is None:
            raise ValueError("normalization='reference' requires reference='series_name'")
        if reference not in raw:
            raise KeyError(f"unknown reference series {reference!r}; available: {', '.join(raw)}")
    grid = _common_coordinate_grid(raw, method=method, coordinate=coordinate)
    aligned = {name: _align_y(value, grid, method=method, coordinate=coordinate)
               for name, value in raw.items()}
    standard_partition = False
    if {"all", "live", "dead"} <= set(populations):
        live_rows = _row_signatures(populations["live"])
        dead_rows = _row_signatures(populations["dead"])
        all_rows = _row_signatures(populations["all"])
        standard_partition = not live_rows.intersection(dead_rows) and all_rows == live_rows.union(dead_rows)
        if standard_partition:
            # Preserve the exact physical identity on the common grid; separate
            # KDE grids would otherwise introduce small interpolation residue.
            aligned["all"] = aligned["live"] + aligned["dead"]
    totals = {name: _series_total(y, grid, output=output) for name, y in aligned.items()}
    combined_total = sum(totals.values())
    physical_x = np.power(10.0, grid) if coordinate == "log10" else grid
    built = {}
    for name, value in raw.items():
        if normalization == "absolute":
            denominator = 1.0
        elif normalization == "per_series":
            denominator = totals[name]
        elif normalization == "combined":
            denominator = combined_total
        else:
            denominator = totals[reference]  # type: ignore[index]
        y = aligned[name] / denominator if denominator > 0 else aligned[name]
        integral = _series_total(y, grid, output=output)
        meta = dict(value.metadata, normalization=normalization, reference=reference,
                    integral=integral, series_name=name)
        cls = ChainLengthDistribution if kind == "CLD" else MolarMassDistribution
        built[name] = cls(value._distribution, _readonly(physical_x), _readonly(y), meta)
    descriptor = next(iter(built.values())).metadata["descriptor"] if built else kind
    first = next(iter(populations.values()))
    metadata = dict(
        _snapshot_meta(first, pool="multiple"), kind=kind, descriptor=descriptor,
        basis=basis, method=method, coordinate=coordinate, output=output,
        normalization=normalization, reference=reference, sigma=sigma,
        series_names=tuple(built), series_disjoint=disjoint,
        standard_partition=standard_partition,
        integrals={name: built[name].metadata["integral"] for name in built},
    )
    cls_multi = MultiChainLengthDistribution if kind == "CLD" else MultiMolarMassDistribution
    return cls_multi(DistributionSeries(built), metadata)


def build_cld_series(population: ChainPopulation, *, series, pool: str = "all", basis: str = "number",
                     method: str = "sticks", coordinate: str = "linear",
                     output: str = "fraction", normalization: str = "per_series",
                     reference: str | None = None, bins: int | None = None,
                     bin_width: float | None = None, sigma: float | None = None,
                     grid_step: float | None = None, mass_model: str | None = None):
    if pool != "all":
        raise ValueError("pool= and series= are mutually exclusive; put pool selectors inside series")
    populations = _resolve_series(population, series)
    raw = {}
    raw_output = "density" if method in {"kde"} or output == "density" else "amount"
    for name, selected in populations.items():
        raw[name] = selected.cld(
            pool="all", basis=basis, method=method, coordinate=coordinate,
            output=raw_output, normalization="absolute", bins=bins,
            bin_width=bin_width, sigma=sigma, grid_step=grid_step,
            mass_model=mass_model,
        )
    return _normalize_multi(raw, populations, kind="CLD", basis=basis, method=method,
                            coordinate=coordinate, output=raw_output if output == "density" else output,
                            normalization=normalization, reference=reference, sigma=sigma)


def build_mwd_series(population: ChainPopulation, *, series, pool: str = "all", basis: str = "mass",
                     method: str = "gaussian", coordinate: str = "log10", output: str = "density",
                     normalization: str = "per_series", reference: str | None = None,
                     bins: int | None = None, bin_width: float | None = None,
                     sigma: float | None = None, grid_step: float | None = None,
                     mass_model: str | None = None):
    if pool != "all":
        raise ValueError("pool= and series= are mutually exclusive; put pool selectors inside series")
    populations = _resolve_series(population, series)
    raw = {}
    raw_output = "density" if method in {"kde"} or output == "density" else "amount"
    for name, selected in populations.items():
        raw[name] = selected.mwd(
            pool="all", basis=basis, method=method, coordinate=coordinate,
            output=raw_output, normalization="absolute", bins=bins,
            bin_width=bin_width, sigma=sigma, grid_step=grid_step,
            mass_model=mass_model,
        )
    return _normalize_multi(raw, populations, kind="MWD", basis=basis, method=method,
                            coordinate=coordinate, output=raw_output if output == "density" else output,
                            normalization=normalization, reference=reference, sigma=sigma)


def _spectrum_y(counts: np.ndarray, normalize: str) -> np.ndarray:
    normalize = str(normalize).lower()
    if normalize not in {"count", "fraction", "base_peak"}:
        raise ValueError("normalize must be count, fraction, or base_peak")
    y = np.asarray(counts, dtype=float).copy()
    if normalize == "fraction":
        total = float(y.sum())
        if total > 0:
            y /= total
    elif normalize == "base_peak":
        peak = float(y.max()) if y.size else 0.0
        if peak > 0:
            y /= peak
    return _readonly(y)


def build_chain_mass_spectrum(population: ChainPopulation, masses, *, pool: str = "all",
                         normalize: str = "base_peak", mass_model: str | None = None) -> ChainMassSpectrum:
    selected = select_pool(population, pool)
    all_masses = np.asarray(masses, dtype=float)
    if len(all_masses) != len(population):
        raise ValueError("masses must contain one value per compressed ChainPopulation row")
    if selected is population:
        selected_masses = all_masses
    else:
        raw = selected._raw_arrays()
        if "mass" in raw:
            selected_masses = np.asarray(raw["mass"], dtype=float)
        elif hasattr(selected, "masses"):
            selected_masses = np.asarray(selected.masses(mass_model=mass_model or "repeat_units"), dtype=float)
        else:
            raise ValueError("cannot map molar masses after spectrum population selection")
    pairs = [(float(mass), float(n), float(mass))
             for mass, n in zip(selected_masses, selected.count)]
    meta = dict(
        _snapshot_meta(selected, pool=pool), kind="chain_mass_spectrum",
        descriptor=f"neutral chain mass spectrum ({normalize})", mass_model=mass_model,
        normalize=normalize, basis="number", method="sticks", coordinate="linear",
        output="intensity", normalization=normalize, neutral_mass=True,
        ionization_model=None, mz=False,
    )
    base = Distribution.from_pairs(
        "Chain spectrum", "M / g mol-1", pairs, basis="number",
        method="chain_mass_spectrum", meta=meta,
    )
    y = _spectrum_y(np.asarray(base.count, dtype=float), normalize)
    meta["integral"] = float(y.sum())
    return ChainMassSpectrum(base, _readonly(base.x), y, meta)


def build_chain_mass_spectrum_series(population: ChainPopulation, *, series, pool: str = "all",
                                normalize: str = "base_peak", mass_model: str | None = None):
    if pool != "all":
        raise ValueError("pool= and series= are mutually exclusive; put pool selectors inside series")
    populations = _resolve_series(population, series)
    raw = {}
    for name, selected in populations.items():
        raw_arrays = selected._raw_arrays()
        if "mass" in raw_arrays:
            masses = raw_arrays["mass"]
        elif hasattr(selected, "masses"):
            masses = selected.masses(mass_model=mass_model or "repeat_units")
        else:
            raise ValueError("chain spectrum requires per-chain molar masses")
        raw[name] = build_chain_mass_spectrum(
            selected, masses, pool="all", normalize="count", mass_model=mass_model
        )
    grid = _common_coordinate_grid(raw, method="sticks", coordinate="linear")
    built = {}
    for name, value in raw.items():
        counts = _align_y(value, grid, method="sticks", coordinate="linear")
        y = _spectrum_y(counts, normalize)
        meta = dict(value.metadata, normalize=normalize, normalization=normalize,
                    descriptor=f"neutral chain mass spectrum ({normalize})",
                    integral=float(y.sum()), series_name=name)
        built[name] = ChainMassSpectrum(value._distribution, _readonly(grid), y, meta)
    first = next(iter(populations.values()))
    disjoint = _are_disjoint(populations)
    metadata = dict(
        _snapshot_meta(first, pool="multiple"), kind="chain_mass_spectrum",
        descriptor=f"neutral chain mass spectrum ({normalize})", basis="number", method="sticks",
        coordinate="linear", output="intensity", normalization=normalize,
        normalize=normalize, mass_model=mass_model, neutral_mass=True,
        ionization_model=None, mz=False, series_names=tuple(built),
        series_disjoint=disjoint,
        integrals={name: built[name].metadata["integral"] for name in built},
    )
    return MultiChainMassSpectrum(DistributionSeries(built), metadata)


def _chain_mass_spectrum_from_table(table: Table, *, normalize: str = "base_peak",
                              snapshot_id: int | None = None) -> ChainMassSpectrum:
    required = {"mass", "count"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"oligomer table is missing required columns: {', '.join(sorted(missing))}")
    rows = table.rows()
    if "snapshot_id" in table.columns:
        ids = sorted({int(row["snapshot_id"]) for row in rows})
        if snapshot_id is None:
            if len(ids) > 1:
                raise ValueError(
                    f"oligomer table contains multiple snapshots {ids}; pass snapshot_id= explicitly"
                )
            snapshot_id = ids[0] if ids else None
        rows = [row for row in rows if int(row["snapshot_id"]) == snapshot_id]
    pairs = [(float(row["mass"]), float(row["count"]), float(row["mass"])) for row in rows]
    meta = dict(
        kind="chain_mass_spectrum", descriptor=f"neutral oligomer mass spectrum ({normalize})",
        snapshot_id=snapshot_id, t=None, kmc_event=None, pool="oligomers",
        mass_model=None, normalize=normalize, basis="number", method="sticks",
        coordinate="linear", output="intensity", normalization=normalize,
        neutral_mass=True, ionization_model=None, mz=False,
    )
    base = Distribution.from_pairs(
        "Oligomer spectrum", "M / g mol-1", pairs, basis="number",
        method="chain_mass_spectrum", meta=meta,
    )
    y = _spectrum_y(np.asarray(base.count, dtype=float), normalize)
    meta["integral"] = float(y.sum())
    return ChainMassSpectrum(base, _readonly(base.x), y, meta)
