"""Chain-composition analyses for Slimmc Storage populations."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np

from .plotting import figure_size, get_style, apply_axes_style


def _readonly(values, *, dtype=None):
    out = np.asarray(values, dtype=dtype)
    out.setflags(write=False)
    return out


def _save(fig, path, dpi):
    if path is not None:
        fig.savefig(Path(path), dpi=dpi, bbox_inches="tight")


def _axes(ax, *, style, span, title=None):
    import matplotlib.pyplot as plt
    cfg = get_style(style)
    if ax is None:
        fig, ax = plt.subplots(figsize=figure_size(style, span=span))
    else:
        fig = ax.figure
    if title:
        ax.set_title(title)
    apply_axes_style(ax, style)
    return fig, ax


def _validate_bounds(minimum, maximum, *, name):
    if minimum is not None and maximum is not None and minimum > maximum:
        raise ValueError(f"{name}: min must be <= max")


def _bin_edges(values: np.ndarray, bins, *, integer: bool) -> np.ndarray:
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    if bins is not None:
        if np.isscalar(bins):
            n = int(bins)
            if n < 1:
                raise ValueError("bins must be positive")
            if finite.size == 0:
                return np.linspace(0.0, 1.0, n + 1)
            lo, hi = float(np.min(finite)), float(np.max(finite))
            if hi == lo:
                hi = lo + 1.0
            return np.linspace(lo, hi, n + 1)
        edges = np.asarray(bins, dtype=float)
        if edges.ndim != 1 or len(edges) < 2 or np.any(np.diff(edges) <= 0):
            raise ValueError("bins must be increasing one-dimensional edges")
        return edges
    if finite.size == 0:
        return np.asarray([0.0, 1.0])
    lo, hi = float(np.min(finite)), float(np.max(finite))
    if integer and hi - lo <= 200:
        return np.arange(np.floor(lo), np.ceil(hi) + 2, dtype=float) - 0.5
    if hi == lo:
        return np.asarray([lo - 0.5, hi + 0.5], dtype=float)
    return np.linspace(lo, hi, 51)


def _weighted_quantile(values, weights, q):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(valid):
        return np.nan
    values = values[valid]
    weights = weights[valid]
    order = np.argsort(values, kind="stable")
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    target = float(q) * cumulative[-1]
    return float(values[min(np.searchsorted(cumulative, target, side="left"), len(values) - 1)])


@dataclass(frozen=True)
class CompositionByDP:
    names: tuple[str, ...]
    dp_left: np.ndarray
    dp_right: np.ndarray
    dp_center: np.ndarray
    record_count: np.ndarray
    chain_count: np.ndarray
    mean: Mapping[str, np.ndarray]
    median: Mapping[str, np.ndarray]
    q25: Mapping[str, np.ndarray]
    q75: Mapping[str, np.ndarray]

    def plot(self, path=None, *, statistic="mean", interval=None, style="screen",
             ax=None, span=None, dpi=300, title=None):
        values = getattr(self, statistic, None)
        if values is None or statistic not in {"mean", "median"}:
            raise ValueError("statistic must be 'mean' or 'median'")
        if interval not in {None, "iqr"}:
            raise ValueError("interval must be None or 'iqr'")
        fig, ax = _axes(ax, style=style, span=span, title=title)
        for name in self.names:
            line, = ax.plot(self.dp_center, values[name], label=name)
            if interval == "iqr":
                ax.fill_between(self.dp_center, self.q25[name], self.q75[name],
                                alpha=0.18, color=line.get_color())
        ax.set_xlabel("DP")
        ax.set_ylabel("chain mole fraction")
        ax.set_ylim(0.0, 1.0)
        ax.legend()
        fig.tight_layout()
        _save(fig, path, dpi)
        return ax


@dataclass(frozen=True)
class CompositionMap:
    x_edges: np.ndarray
    y_edges: np.ndarray
    values: np.ndarray
    x_label: str
    y_label: str
    value_label: str = "chain count"

    def plot(self, path=None, *, log=False, style="screen", ax=None, span=None,
             dpi=300, title=None):
        from matplotlib.colors import LogNorm
        fig, ax = _axes(ax, style=style, span=span, title=title)
        data = np.asarray(self.values, dtype=float).T
        norm = None
        if log:
            positive = data[data > 0]
            if positive.size:
                norm = LogNorm(vmin=float(np.min(positive)), vmax=float(np.max(positive)))
        mesh = ax.pcolormesh(self.x_edges, self.y_edges, data, shading="auto", norm=norm)
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        fig.colorbar(mesh, ax=ax, label=self.value_label)
        fig.tight_layout()
        _save(fig, path, dpi)
        return ax


@dataclass(frozen=True)
class ComponentClasses:
    labels: tuple[str, ...]
    record_count: np.ndarray
    chain_count: np.ndarray
    number_fraction: np.ndarray
    mass_fraction: np.ndarray

    def plot(self, path=None, *, value="number_fraction", style="screen", ax=None,
             span=None, dpi=300, title=None):
        if value not in {"record_count", "chain_count", "number_fraction", "mass_fraction"}:
            raise ValueError("unsupported value")
        fig, ax = _axes(ax, style=style, span=span, title=title)
        data = np.asarray(getattr(self, value))
        ax.bar(self.labels, data)
        ax.set_xlabel("components")
        ax.set_ylabel(value.replace("_", " "))
        fig.tight_layout()
        _save(fig, path, dpi)
        return ax


def composition_by_dp(chains, *, bins=None) -> CompositionByDP:
    dp = np.asarray(chains.dp, dtype=float)
    weights = np.asarray(chains.count, dtype=float)
    fractions = np.asarray(chains.composition.fractions.matrix, dtype=float)
    names = tuple(chains.composition.names)
    edges = _bin_edges(dp, bins, integer=True)
    n = len(edges) - 1
    idx = np.searchsorted(edges, dp, side="right") - 1
    idx[dp == edges[-1]] = n - 1
    valid = (idx >= 0) & (idx < n) & np.isfinite(dp) & (weights > 0)
    record_count = np.bincount(idx[valid], minlength=n).astype(np.uint64)
    chain_count = np.bincount(idx[valid], weights=weights[valid], minlength=n)
    mean, median, q25, q75 = {}, {}, {}, {}
    for col, name in enumerate(names):
        vals = fractions[:, col]
        ok = valid & np.isfinite(vals)
        sums = np.bincount(idx[ok], weights=weights[ok] * vals[ok], minlength=n)
        den = np.bincount(idx[ok], weights=weights[ok], minlength=n)
        mean[name] = _readonly(np.divide(sums, den, out=np.full(n, np.nan), where=den > 0))
        med = np.full(n, np.nan); lo = np.full(n, np.nan); hi = np.full(n, np.nan)
        for b in range(n):
            use = ok & (idx == b)
            med[b] = _weighted_quantile(vals[use], weights[use], 0.5)
            lo[b] = _weighted_quantile(vals[use], weights[use], 0.25)
            hi[b] = _weighted_quantile(vals[use], weights[use], 0.75)
        median[name], q25[name], q75[name] = map(_readonly, (med, lo, hi))
    return CompositionByDP(
        names=names,
        dp_left=_readonly(edges[:-1]), dp_right=_readonly(edges[1:]),
        dp_center=_readonly((edges[:-1] + edges[1:]) / 2),
        record_count=_readonly(record_count), chain_count=_readonly(chain_count),
        mean=mean, median=median, q25=q25, q75=q75,
    )


def composition_dp_map(chains, monomer: str, *, dp_bins=None, fraction_bins=None) -> CompositionMap:
    dp = np.asarray(chains.dp, dtype=float)
    frac = np.asarray(chains.composition.fractions[monomer], dtype=float)
    weights = np.asarray(chains.count, dtype=float)
    x_edges = _bin_edges(dp, dp_bins, integer=True)
    y_edges = _bin_edges(frac, np.linspace(0.0, 1.0, 51) if fraction_bins is None else fraction_bins,
                         integer=False)
    valid = np.isfinite(dp) & np.isfinite(frac) & np.isfinite(weights) & (weights > 0)
    values, _, _ = np.histogram2d(dp[valid], frac[valid], bins=(x_edges, y_edges), weights=weights[valid])
    return CompositionMap(_readonly(x_edges), _readonly(y_edges), _readonly(values),
                          "DP", f"chain mole fraction {monomer}")


def composition_mass_map(chains, monomer: str, *, mass_model="with_end_groups",
                         mass_bins=None, fraction_bins=None) -> CompositionMap:
    mass = np.asarray(chains.masses(mass_model=mass_model), dtype=float)
    frac = np.asarray(chains.composition.fractions[monomer], dtype=float)
    weights = np.asarray(chains.count, dtype=float)
    x_edges = _bin_edges(mass, mass_bins, integer=False)
    y_edges = _bin_edges(frac, np.linspace(0.0, 1.0, 51) if fraction_bins is None else fraction_bins,
                         integer=False)
    valid = np.isfinite(mass) & np.isfinite(frac) & np.isfinite(weights) & (weights > 0)
    values, _, _ = np.histogram2d(mass[valid], frac[valid], bins=(x_edges, y_edges), weights=weights[valid])
    return CompositionMap(_readonly(x_edges), _readonly(y_edges), _readonly(values),
                          "molar mass / g mol$^{-1}$", f"chain mole fraction {monomer}")


def composition_map(chains, x: str, y: str, *, bins=None) -> CompositionMap:
    xf = np.asarray(chains.composition.fractions[x], dtype=float)
    yf = np.asarray(chains.composition.fractions[y], dtype=float)
    weights = np.asarray(chains.count, dtype=float)
    if bins is None:
        bins = np.linspace(0.0, 1.0, 51)
    if isinstance(bins, tuple):
        x_edges = _bin_edges(xf, bins[0], integer=False); y_edges = _bin_edges(yf, bins[1], integer=False)
    else:
        x_edges = _bin_edges(xf, bins, integer=False); y_edges = _bin_edges(yf, bins, integer=False)
    valid = np.isfinite(xf) & np.isfinite(yf) & np.isfinite(weights) & (weights > 0)
    values, _, _ = np.histogram2d(xf[valid], yf[valid], bins=(x_edges, y_edges), weights=weights[valid])
    return CompositionMap(_readonly(x_edges), _readonly(y_edges), _readonly(values),
                          f"chain mole fraction {x}", f"chain mole fraction {y}")


def component_classes(chains) -> ComponentClasses:
    matrix = np.asarray(chains.composition.matrix)
    present = matrix > 0
    names = tuple(chains.composition.names)
    weights = np.asarray(chains.count, dtype=float)
    masses = np.asarray(chains.molar_mass, dtype=float)
    labels_by_row = ["".join(name for name, flag in zip(names, row) if flag) or "none" for row in present]
    preferred = []
    for size in range(1, len(names) + 1):
        for mask in range(1, 1 << len(names)):
            if mask.bit_count() == size:
                preferred.append("".join(names[i] for i in range(len(names)) if mask & (1 << i)))
    labels = tuple(label for label in preferred if label in set(labels_by_row))
    rc=[]; cc=[]; mc=[]
    for label in labels:
        use = np.asarray([item == label for item in labels_by_row])
        rc.append(int(np.count_nonzero(use)))
        cc.append(float(np.sum(weights[use])))
        mc.append(float(np.sum(weights[use] * masses[use])))
    rc=np.asarray(rc,dtype=np.uint64); cc=np.asarray(cc); mc=np.asarray(mc)
    nf=np.divide(cc,np.sum(cc),out=np.zeros_like(cc),where=np.sum(cc)>0)
    mf=np.divide(mc,np.sum(mc),out=np.zeros_like(mc),where=np.sum(mc)>0)
    return ComponentClasses(labels,_readonly(rc),_readonly(cc),_readonly(nf),_readonly(mf))

class RunPlotNamespace:
    """Thin plotting shortcuts bound to one run.

    The namespace never replaces analysis objects and their ``.plot()``
    methods.  It only provides discoverable one-call shortcuts for common
    exploratory plots.
    """
    def __init__(self, run):
        self._run = run

    def _chains(self, snapshot):
        return self._run._resolve_chain_snapshot(snapshot).chains

    @staticmethod
    def _split_plot_kwargs(kwargs, extra=()):
        keys = {"path", "style", "ax", "span", "dpi", "title",
                "mode", "display_normalization", *extra}
        return {key: kwargs.pop(key) for key in tuple(kwargs) if key in keys}

    def _x(self, x):
        if x == "time":
            return np.asarray(self._run.t, dtype=float), "time"
        if x == "conversion":
            return np.asarray(self._run.conv.total, dtype=float), "overall conversion"
        raise ValueError("x must be 'time' or 'conversion'")

    def _lines(self, values, *, x="time", names=None, total=None, ylabel,
               path=None, style="screen", ax=None, span=None, dpi=300,
               title=None, ylim=None):
        x_values, xlabel = self._x(x)
        fig, ax = _axes(ax, style=style, span=span, title=title)
        available = tuple(values.keys())
        selected = available if names is None else tuple(names)
        unknown = [name for name in selected if name not in available]
        if unknown:
            raise KeyError(f"unknown series {unknown}; available: {list(available)}")
        cfg = get_style(style)
        for i, name in enumerate(selected):
            ax.plot(x_values, np.asarray(values[name], dtype=float), label=str(name),
                    color=cfg.palette[i % len(cfg.palette)], linewidth=cfg.line_width)
        if total is not None:
            ax.plot(x_values, np.asarray(total, dtype=float), label="total", linestyle="--",
                    color=cfg.palette[len(selected) % len(cfg.palette)], linewidth=cfg.line_width)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        if ylim is not None: ax.set_ylim(*ylim)
        if selected or total is not None: ax.legend()
        apply_axes_style(ax, style); _save(fig, path, dpi)
        return ax

    # Existing distributions -------------------------------------------------
    def mwd(self, *, snapshot="final", **kwargs):
        plot = self._split_plot_kwargs(kwargs)
        return self._run.mwd(snapshot=snapshot, **kwargs).plot(**plot)

    def mass_distribution(self, *, snapshot="final", **kwargs):
        plot = self._split_plot_kwargs(kwargs)
        return self._run.mass_distribution(snapshot=snapshot, **kwargs).plot(**plot)

    def cld(self, *, snapshot="final", **kwargs):
        plot = self._split_plot_kwargs(kwargs)
        return self._run.cld(snapshot=snapshot, **kwargs).plot(**plot)


    def dp_counts(self, *, snapshot="final", pool="all", **plot_kwargs):
        return self._run.dp_counts(snapshot=snapshot, pool=pool).plot(**plot_kwargs)

    def mass_counts(self, *, snapshot="final", pool="all", mass_model=None, **plot_kwargs):
        return self._run.mass_counts(snapshot=snapshot, pool=pool, mass_model=mass_model).plot(**plot_kwargs)

    # General run series ------------------------------------------------------
    def conversion(self, *, x="time", monomers=None, total=True, **plot_kwargs):
        values = {name: self._run.conv[name] for name in self._run.conv.keys()}
        return self._lines(values, x=x, names=monomers,
                           total=self._run.conv.total if total else None,
                           ylabel="conversion", ylim=(0.0, 1.05), **plot_kwargs)

    def concentrations(self, *, x="time", entities=None, **plot_kwargs):
        values = {name: self._run.conc[name] for name in self._run.state.names}
        return self._lines(values, x=x, names=entities, ylabel="concentration (mol/L)", **plot_kwargs)

    def counts(self, *, x="time", entities=None, **plot_kwargs):
        values = {name: self._run.count[name] for name in self._run.state.names}
        return self._lines(values, x=x, names=entities, ylabel="count", **plot_kwargs)

    def moles(self, *, x="time", entities=None, **plot_kwargs):
        values = {name: self._run.moles[name] for name in self._run.state.names}
        return self._lines(values, x=x, names=entities, ylabel="amount (mol)", **plot_kwargs)

    def temperature(self, *, x="time", **plot_kwargs):
        values = {"temperature": np.asarray(self._run.temp, dtype=float)}
        return self._lines(values, x=x, ylabel="temperature", **plot_kwargs)

    def volume(self, *, x="time", **plot_kwargs):
        values = {"volume": np.asarray(self._run.volume, dtype=float)}
        return self._lines(values, x=x, ylabel="volume (L)", **plot_kwargs)

    # Copolymerization series -------------------------------------------------
    def monomer_composition(self, *, x="conversion", **plot_kwargs):
        return self._run.copolymerization.monomer_composition().plot(x=x, **plot_kwargs)

    def incremental_composition(self, *, x="conversion", **plot_kwargs):
        return self._run.copolymerization.incremental_composition().plot(x=x, **plot_kwargs)

    def cumulative_composition(self, *, x="conversion", **plot_kwargs):
        return self._run.copolymerization.cumulative_composition().plot(x=x, **plot_kwargs)

    def composition_drift(self, *, monomer_reference="start", **plot_kwargs):
        return self._run.copolymerization.composition_drift(
            monomer_reference=monomer_reference).plot(**plot_kwargs)

    def mayo_lewis(self, **plot_kwargs):
        return self._run.copolymerization.mayo_lewis().plot(**plot_kwargs)

    def compare_mayo_lewis(self, *, monomer_reference="start",
                           parameter_reference="start", **plot_kwargs):
        return self._run.copolymerization.compare_mayo_lewis(
            monomer_reference=monomer_reference,
            parameter_reference=parameter_reference).plot(**plot_kwargs)

    # Chain-composition analysis ---------------------------------------------
    def composition_by_dp(self, *, snapshot="final", bins=None, **plot_kwargs):
        return self._chains(snapshot).composition_by_dp(bins=bins).plot(**plot_kwargs)

    def composition_dp_map(self, monomer: str, *, snapshot="final", dp_bins=None,
                           fraction_bins=None, **plot_kwargs):
        result = self._chains(snapshot).composition_dp_map(
            monomer, dp_bins=dp_bins, fraction_bins=fraction_bins)
        return result.plot(**plot_kwargs)

    def composition_mass_map(self, monomer: str, *, snapshot="final",
                             mass_model="with_end_groups", mass_bins=None,
                             fraction_bins=None, **plot_kwargs):
        result = self._chains(snapshot).composition_mass_map(
            monomer, mass_model=mass_model, mass_bins=mass_bins,
            fraction_bins=fraction_bins)
        return result.plot(**plot_kwargs)

    def composition_map(self, x: str, y: str, *, snapshot="final", bins=None,
                        **plot_kwargs):
        return self._chains(snapshot).composition_map(x, y, bins=bins).plot(**plot_kwargs)

    def component_classes(self, *, snapshot="final", **plot_kwargs):
        return self._chains(snapshot).component_classes().plot(**plot_kwargs)

    # Full-sequence analysis --------------------------------------------------
    def block_lengths(self, monomer=None, *, snapshot="final", progress=None, **plot_kwargs):
        return self._chains(snapshot).block_lengths(monomer, progress=progress).plot(**plot_kwargs)

    def transition_matrix(self, *, snapshot="final", normalize=None, progress=None, **plot_kwargs):
        return self._chains(snapshot).transition_matrix(normalize=normalize, progress=progress).plot(**plot_kwargs)

    def microstructure_by_dp(self, statistic, *, snapshot="final", monomer=None,
                             bins=None, progress=None, **plot_kwargs):
        result = self._chains(snapshot).microstructure_by_dp(
            statistic, monomer=monomer, bins=bins, progress=progress)
        return result.plot(**plot_kwargs)


    def ngrams(self, n=4, *, snapshot="final", min_count=1, progress=None, **plot_kwargs):
        return self._chains(snapshot).ngrams(n=n, min_count=min_count, progress=progress).plot(**plot_kwargs)

    def position_profile(self, *, snapshot="final", bins=20, progress=None, **plot_kwargs):
        return self._chains(snapshot).position_profile(bins=bins, progress=progress).plot(**plot_kwargs)

    def microstructure_map(self, statistic, *, snapshot="final", monomer=None,
                           dp_bins=None, value_bins=None, progress=None, **plot_kwargs):
        result = self._chains(snapshot).microstructure_map(
            statistic, monomer=monomer, dp_bins=dp_bins, value_bins=value_bins, progress=progress)
        return result.plot(**plot_kwargs)
