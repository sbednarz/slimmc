from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable
import csv
import math

from .table import Table
from .run import UnsupportedFeatureError


def _validate_basis(basis: str) -> str:
    basis = str(basis).lower()
    if basis not in {"number", "mass"}:
        raise ValueError("basis must be 'number' or 'mass'")
    return basis


@dataclass(frozen=True)
class Distribution:
    """Small no-Pandas distribution object used for CLD/MWD (and related
    per-chain distributions: chain-mass spectra, composition histograms).

    Design (see docs/PYSLIMMC.md and the review that prompted this
    rewrite):

    - Always stores three parallel weight arrays -- ``count`` (number),
      ``mass_weight`` (= count * x, i.e. mass-weighted), and ``z_weight``
      (= count * x**2, i.e. second-moment) -- computed once, consistently,
      by :meth:`from_pairs`. ``basis`` only *selects which one* ``.y``
      exposes for plotting; it never changes what's stored. This is what
      makes ``basis="mass"`` actually correct (previously it silently
      returned the mass-weighted result on both engines) and what makes
      ``.mn()``/``.mw()``/``.pdi()`` well-defined regardless of which
      basis a Distribution happened to be constructed with.
    - A fourth array, ``z2_weight`` (= count * x**3), exists purely to
      make ``.mz()`` mathematically correct (critical-analysis P0.1):
      Mz = Σ(N·M³)/Σ(N·M²) genuinely requires the third moment, which the
      original three-array design (count/mass_weight/z_weight) had no
      way to provide -- the previous ``.mz()`` recomputed
      Σ(N·M²)/Σ(N·M), which is exactly the Mw formula, so ``mz() ==
      mw()`` always, silently. ``z2_weight`` is not exposed as a
      plotting ``basis`` (there is no meaningful "z2 distribution" to
      plot) -- it exists solely to feed ``source_stats["mz"]``.
    - ``count``/``mass_weight``/``z_weight``/``z2_weight`` are always
      *true bin/stick totals* (their sum is the real total chain count /
      mass / z-moment / third-moment of the population), never a density
      like dN/bin_width. A log-spaced histogram's bins are still evenly
      spaced in log-x (useful for plotting), but the values in each bin
      are plain totals, not divided by the bin width -- this is what
      makes ``.n``, ``.total_weight``, and ``.as_table()``'s fraction
      columns mean what their names say, on every method
      (sticks/hist/gaussian/kde), not just on "sticks".
    - ``.mn()``/``.mw()``/``.mz()``/``.pdi()`` are computed **once**, at
      :meth:`from_pairs` (i.e. from the original, unbinned, linear-space
      data), and carried forward unchanged by every derived transform
      (:meth:`from_histogram`, :meth:`from_histogram_log`,
      :meth:`gaussian`, :meth:`kde`). This is deliberate: computing "mn"
      from a log10(x)-binned or Gaussian-smoothed array would silently
      give the mean of log10(x), not Mn in the distribution's real units
      -- a real bug found in the pre-rewrite version. Regardless of which
      transform produced the object you're holding, ``.mn()``/``.mw()``/
      ``.mz()``/``.pdi()`` always answer for the *original* population.
    """

    name: str
    x_name: str
    x: tuple[float, ...]
    count: tuple[float, ...]
    mass_weight: tuple[float, ...]
    z_weight: tuple[float, ...]
    z2_weight: tuple[float, ...] = ()
    basis: str = "mass"
    method: str = "sticks"
    source_stats: dict = field(default_factory=dict)
    meta: dict = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.meta is None:
            object.__setattr__(self, "meta", {})

    # -- construction -------------------------------------------------

    @classmethod
    def from_pairs(
        cls,
        name: str,
        x_name: str,
        pairs: Iterable[tuple[float, float, float]],
        *,
        basis: str = "mass",
        method: str = "sticks",
        meta: dict | None = None,
    ) -> "Distribution":
        """Build an exact (unbinned) "sticks" distribution from raw
        ``(x, count, mass)`` triples (summed over repeated ``x`` values).
        ``mass_weight``/``z_weight``/``z2_weight`` are derived
        automatically (``count*x``/``count*x**2``/``count*x**3``) -- this
        is the single place all four weight arrays and the source-level
        moments (including the corrected ``mz``, critical-analysis P0.1)
        are computed from, and every other constructor/transform on this
        class carries these forward rather than recomputing them."""
        basis = _validate_basis(basis)
        data: dict[float, list[float]] = {}
        for x, count, mass in pairs:
            if x not in data:
                data[x] = [0.0, 0.0, 0.0, 0.0]
            data[x][0] += float(count)
            data[x][1] += float(count) * float(mass)
            data[x][2] += float(count) * float(mass) ** 2
            data[x][3] += float(count) * float(mass) ** 3
        xs = tuple(sorted(data))
        counts = tuple(data[v][0] for v in xs)
        mass_weights = tuple(data[v][1] for v in xs)
        z_weights = tuple(data[v][2] for v in xs)
        z2_weights = tuple(data[v][3] for v in xs)
        total_n = sum(counts)
        total_mass = sum(mass_weights)
        total_z = sum(z_weights)
        total_z2 = sum(z2_weights)
        mn = (total_mass / total_n) if total_n > 0 else 0.0
        mw = (total_z / total_mass) if total_mass > 0 else 0.0
        mz = (total_z2 / total_z) if total_z > 0 else 0.0
        stats = {
            "n": total_n, "total_mass_weight": total_mass, "total_z_weight": total_z,
            "total_z2_weight": total_z2,
            "mn": mn, "mw": mw, "mz": mz, "pdi": (mw / mn if mn > 0 else 0.0),
        }
        return cls(name=name, x_name=x_name, x=xs, count=counts, mass_weight=mass_weights,
                   z_weight=z_weights, z2_weight=z2_weights, basis=basis, method=method,
                   source_stats=stats, meta=dict(meta or {}))

    @classmethod
    def from_histogram(
        cls,
        name: str,
        x_name: str,
        pairs: Iterable[tuple[float, float, float]],
        *,
        bins: int | None = 100,
        bin_width: float | None = None,
        basis: str = "mass",
        meta: dict | None = None,
    ) -> "Distribution":
        """Uniform *linear*-width histogram over ``(x, count, mass)``
        triples. Keeps empty bins (a regular grid, consistent with
        :meth:`from_histogram_log` -- the pre-rewrite version silently
        dropped them here, which could make line plots connect
        non-adjacent bins). Bin values are true totals, not densities."""
        basis = _validate_basis(basis)
        base = cls.from_pairs(name, x_name, pairs, basis=basis, method="sticks")
        base_meta = dict(meta or {})
        if not base.x:
            return cls(name=name, x_name=x_name, x=(), count=(), mass_weight=(), z_weight=(), z2_weight=(),
                       basis=basis, method="hist", source_stats=base.source_stats, meta=base_meta)
        if bins is not None and bin_width is not None:
            raise ValueError("bins and bin_width are mutually exclusive")
        xmin, xmax = min(base.x), max(base.x)
        if bin_width is not None:
            width = float(bin_width)
            if width <= 0:
                raise ValueError("bin_width must be > 0")
            lo = math.floor(xmin / width) * width
            hi = math.ceil(xmax / width) * width
            if hi <= lo:
                hi = lo + width
            bins = max(1, int(round((hi - lo) / width)))
        else:
            bins = max(1, int(100 if bins is None else bins))
            if xmin == xmax or bins == 1:
                width = 1.0
                lo = xmin - 0.5
                bins = 1
            else:
                lo = xmin
                width = (xmax - xmin) / bins
        base_meta.update(method="hist", coordinate="linear", bins=bins, bin_width=width)
        counts = [0.0] * bins
        mass_weights = [0.0] * bins
        z_weights = [0.0] * bins
        z2_weights = [0.0] * bins
        for x, c, mw, zw, z2w in zip(base.x, base.count, base.mass_weight, base.z_weight, base.z2_weight):
            idx = min(bins - 1, max(0, int((x - lo) / width)))
            counts[idx] += c
            mass_weights[idx] += mw
            z_weights[idx] += zw
            z2_weights[idx] += z2w
        centers = tuple(lo + (i + 0.5) * width for i in range(bins))  # keeps empty bins
        return cls(name=name, x_name=x_name, x=centers, count=tuple(counts),
                   mass_weight=tuple(mass_weights), z_weight=tuple(z_weights), z2_weight=tuple(z2_weights),
                   basis=basis, method="hist", source_stats=base.source_stats, meta=base_meta)

    @classmethod
    def from_histogram_log(
        cls,
        base: "Distribution",
        *,
        bins: int | None = None,
        bin_width: float | None = None,
    ) -> "Distribution":
        """Histogram on a regular ``log10(x)`` grid.

        ``bins`` and ``bin_width`` are alternatives. ``bin_width`` is in
        decades. If neither is supplied, 0.02 decade is used.
        """
        if bins is not None and bin_width is not None:
            raise ValueError("bins and bin_width are mutually exclusive")
        if any(x <= 0 for x in base.x):
            raise ValueError("log10 histogram requires all x values > 0")
        if not base.x:
            meta = dict(base.meta, method="hist", coordinate="log10")
            return cls(base.name, "log10(" + base.x_name + ")", (), (), (), (), (),
                       basis=base.basis, method="hist", source_stats=base.source_stats, meta=meta)

        logx = [math.log10(v) for v in base.x]
        lo_raw, hi_raw = min(logx), max(logx)
        if bins is not None:
            bins = max(1, int(bins))
            if lo_raw == hi_raw or bins == 1:
                width = 1.0
                lo = lo_raw - 0.5
                n_bins = 1
            else:
                width = (hi_raw - lo_raw) / bins
                lo = lo_raw
                n_bins = bins
        else:
            width = 0.02 if bin_width is None else float(bin_width)
            if width <= 0:
                raise ValueError("bin_width must be > 0")
            lo = math.floor(lo_raw / width) * width
            hi = math.ceil(hi_raw / width) * width
            if hi <= lo:
                hi = lo + width
            n_bins = max(1, int(round((hi - lo) / width)))

        counts = [0.0] * n_bins
        mass_weights = [0.0] * n_bins
        z_weights = [0.0] * n_bins
        z2_weights = [0.0] * n_bins
        for value, c, mw, zw, z2w in zip(logx, base.count, base.mass_weight, base.z_weight, base.z2_weight):
            idx = min(n_bins - 1, max(0, int((value - lo) / width)))
            counts[idx] += c
            mass_weights[idx] += mw
            z_weights[idx] += zw
            z2_weights[idx] += z2w
        centers = tuple(lo + (i + 0.5) * width for i in range(n_bins))
        meta = dict(base.meta, method="hist", coordinate="log10", bin_width=width, bins=n_bins)
        return cls(base.name, "log10(" + base.x_name + ")", centers,
                   tuple(counts), tuple(mass_weights), tuple(z_weights), tuple(z2_weights),
                   basis=base.basis, method="hist", source_stats=base.source_stats, meta=meta)

    def gaussian(self, *, sigma: float) -> "Distribution":
        """Gaussian-smoothed version of a (log-spaced) histogram
        distribution. Ported from classic slimmc's ``mwd_gaussian``,
        generalized to smooth ``count``/``mass_weight``/``z_weight`` in
        parallel with the same kernel matrix, mass-conserving (the sum of
        each array is unchanged by smoothing).

        Requires ``method == "hist"`` on an evenly-spaced axis -- calling
        this on a "sticks" distribution (irregular x-spacing, by
        construction) previously ran anyway and produced a mathematically
        meaningless result; it's a hard error now."""
        import numpy as np

        if self.method != "hist":
            raise ValueError(
                f"gaussian() requires a method='hist' distribution with an evenly-spaced "
                f"axis (call .from_histogram()/.from_histogram_log() first) -- got "
                f"method={self.method!r}, which has no guaranteed regular spacing"
            )
        if sigma <= 0.0:
            raise ValueError("sigma must be > 0")
        if len(self.x) < 2:
            return Distribution(self.name, self.x_name, self.x, self.count, self.mass_weight,
                                self.z_weight, self.z2_weight, basis=self.basis, method="gaussian",
                                source_stats=self.source_stats,
                                meta=dict(self.meta, method="gaussian", sigma=sigma))
        xs = np.asarray(self.x, dtype=float)
        diff = xs[:, None] - xs[None, :]
        kernels = np.exp(-0.5 * (diff / sigma) ** 2)
        kernels /= kernels.sum(axis=0, keepdims=True)
        counts_out = kernels @ np.asarray(self.count, dtype=float)
        mass_out = kernels @ np.asarray(self.mass_weight, dtype=float)
        z_out = kernels @ np.asarray(self.z_weight, dtype=float)
        z2_out = kernels @ np.asarray(self.z2_weight, dtype=float)
        return Distribution(self.name, self.x_name, self.x, tuple(counts_out.tolist()),
                            tuple(mass_out.tolist()), tuple(z_out.tolist()), tuple(z2_out.tolist()),
                            basis=self.basis, method="gaussian", source_stats=self.source_stats,
                            meta=dict(self.meta, method="gaussian", sigma=sigma))

    def kde(self, *, bandwidth: float, transform: str | None = None,
            grid_step: float | None = None) -> "Distribution":
        """Kernel density estimate from the unbinned "sticks" distribution.
        Implements the formal KDE contract (docs/PYSLIMMC.md):

        - ``transform`` selects the axis KDE is evaluated on: ``None``
          (plain x, e.g. DP) or ``"log10"`` (log10(x), the usual choice
          for molecular-weight-like data spanning orders of magnitude).
          The transform is never applied silently -- it's always recorded
          in ``meta["axis_transform"]`` and reflected in ``x_name``.
        - Computes ``count_density``/``weight_density`` in parallel (both
          integrate, over the returned axis, to the source's total count/
          weight -- KDE never changes the source population).
        - ``.mn()``/``.mw()``/``.pdi()`` still come from the untouched
          source-level stats (this was already this class's design from
          the Category 2 rewrite; the KDE contract just makes explicit
          that computing them from a log10(x) density axis would be
          wrong).
        """
        import numpy as np

        if transform not in (None, "log10"):
            raise ValueError(f"transform must be None or 'log10', got {transform!r}")
        if not (math.isfinite(bandwidth) and bandwidth > 0.0):
            raise ValueError("bandwidth must be a finite number > 0")
        if grid_step is not None and not (math.isfinite(grid_step) and grid_step > 0.0):
            raise ValueError("grid_step must be a finite number > 0 when given")
        total_n = sum(self.count)
        if not self.x or total_n <= 0:
            raise ValueError("KDE requires a non-empty distribution with positive total count")
        if transform == "log10" and any(x <= 0 for x in self.x):
            raise ValueError(
                "kde(transform='log10') requires all x values > 0 (log10 of a "
                "non-positive mass/DP is undefined) -- check for zero or "
                "negative values in the source data, or use transform=None"
            )

        step = grid_step if grid_step is not None else bandwidth / 20.0
        u_source = np.log10(np.asarray(self.x, dtype=float)) if transform == "log10" else np.asarray(self.x, dtype=float)
        lo = float(u_source.min()) - 5.0 * bandwidth
        hi = float(u_source.max()) + 5.0 * bandwidth
        grid = np.arange(lo, hi + step * 0.5, step)
        norm = math.sqrt(2.0 * math.pi) * bandwidth
        z = (grid[:, None] - u_source[None, :]) / bandwidth
        kernel = np.exp(-0.5 * z * z) / norm  # each column integrates to 1 over the grid axis
        count_density = kernel @ np.asarray(self.count, dtype=float)
        weight_density = kernel @ np.asarray(self.mass_weight, dtype=float)
        z_density = kernel @ np.asarray(self.z_weight, dtype=float)
        z2_density = kernel @ np.asarray(self.z2_weight, dtype=float)

        x_name = f"log10({self.x_name})" if transform == "log10" else self.x_name
        meta = dict(
            self.meta, method="kde", transform=transform, axis_transform=transform,
            bandwidth=bandwidth, grid_step=step,
            source_method=self.method, source_x_name=self.x_name,
            source_count_total=total_n, source_weight_total=sum(self.mass_weight),
            axis_kind=("log10" if transform == "log10" else "linear"),
            density=True, statistics_source="sticks",
        )
        return Distribution(self.name, x_name, tuple(grid.tolist()), tuple(count_density.tolist()),
                            tuple(weight_density.tolist()), tuple(z_density.tolist()), tuple(z2_density.tolist()),
                            basis=self.basis, method="kde", source_stats=self.source_stats, meta=meta)

    # -- summary statistics --------------------------------------------
    # Always sourced from source_stats (computed once from the original,
    # unbinned, linear-space data) -- never recomputed from self.x/count/
    # etc., which may currently be log10-transformed and/or smoothed.

    @property
    def n(self) -> float:
        return self.source_stats.get("n", sum(self.count))

    @property
    def total_weight(self) -> float:
        """Total for the current public basis: number or mass."""
        if self.basis == "mass":
            return self.source_stats.get("total_mass_weight", sum(self.mass_weight))
        return self.n

    @property
    def y_label(self) -> str:
        return {"number": "count", "mass": "weight"}[self.basis]

    @property
    def y(self) -> tuple[float, ...]:
        if self.basis == "number":
            return self.count
        return self.mass_weight

    @property
    def total_y(self) -> float:
        return sum(self.y)

    def mn(self) -> float:
        """Number-average, always from the original linear-space data
        regardless of the current method/basis (see class docstring)."""
        return self.source_stats.get("mn", 0.0)

    def mw(self) -> float:
        """Weight-average, always from the original linear-space data."""
        return self.source_stats.get("mw", 0.0)

    def mz(self) -> float:
        """z-average, always from the original linear-space data.

        Previously this recomputed ``total_z_weight / total_mass_weight``
        (i.e. Σ(N·M²)/Σ(N·M)) -- which is the formula for Mw, not Mz. The
        class had no stored third moment, so ``mz() == mw()`` always,
        silently. Now sourced from ``source_stats["mz"]`` (computed once
        in ``from_pairs`` from the real third moment, ``z2_weight`` =
        count*x**3), exactly like ``.mn()``/``.mw()``/``.pdi()``. See
        docs/PYSLIMMC.md / critical-analysis P0.1."""
        return self.source_stats.get("mz", 0.0)

    def pdi(self) -> float:
        return self.source_stats.get("pdi", 0.0)

    def _require_cld(self, method_name: str) -> None:
        if self.meta.get("kind") != "CLD":
            raise UnsupportedFeatureError(
                f"{method_name}() is only defined for a CLD distribution "
                f"(meta['kind'] == 'CLD'); this object's meta['kind'] is "
                f"{self.meta.get('kind')!r}. Use mn()/mw()/mz() instead -- "
                f"same numbers, mass-distribution naming."
            )

    def dpn(self) -> float:
        """DPn -- alias for ``mn()``, only defined on a CLD distribution
        (``meta["kind"] == "CLD"``). Same formula (Σ(N·x)/Σ(N)); DPn is
        just the physically correct name when ``x`` is chain length
        rather than mass. See ``info_text()``, which already used this
        naming for display; ``to_dict()`` adds these keys alongside
        ``mn``/``mw``/``mz`` for CLD objects specifically."""
        self._require_cld("dpn")
        return self.mn()

    def dpw(self) -> float:
        self._require_cld("dpw")
        return self.mw()

    def dpz(self) -> float:
        self._require_cld("dpz")
        return self.mz()

    # -- export ---------------------------------------------------------

    def as_table(self) -> Table:
        """Return the current distribution as a small table.

        Critical-analysis P2.4/#13: for a KDE result specifically, the
        density columns get a matching ``area_fraction_*`` column that
        is **area-under-curve=1** (dividing by ``total * grid_step``,
        the same true-integral normalization ``.plot(normalize=True)``
        uses) -- not a point-sum-to-1 convenience view. An earlier
        version of this method exported ``point_fraction_*`` (sum=1
        across grid *points*, ignoring grid spacing) specifically to
        make the difference from ``.plot()`` visible in the column name;
        that turned out to still be a not-quite-physical quantity worth
        removing rather than just renaming -- ``area_fraction_*`` is the
        one number here that means what "fraction of the distribution"
        actually means for a density, consistent with ``.plot()``."""
        total_count = self.n
        is_kde = self.method == "kde"
        y = self.mass_weight
        total_y = self.source_stats.get("total_mass_weight", sum(self.mass_weight))
        weight_name = "weight"
        count_col = "count_density" if is_kde else "count"
        weight_col = f"{weight_name}_density" if is_kde else weight_name
        frac_count_col = "area_fraction_count_density" if is_kde else "fraction_count"
        frac_weight_col = f"area_fraction_{weight_name}_density" if is_kde else f"fraction_{weight_name}"
        sum_count = sum(self.count) or 1.0
        sum_y = sum(y) or 1.0
        if is_kde:
            grid_step = self.meta.get("grid_step")
            area_count = sum_count * grid_step if grid_step else sum_count
            area_y = sum_y * grid_step if grid_step else sum_y
        rows = []
        for x, c, yv in zip(self.x, self.count, y):
            if is_kde:
                rows.append((x, c, yv,
                            c / area_count if area_count > 0 else 0.0,
                            yv / area_y if area_y > 0 else 0.0))
            else:
                rows.append((
                    x, c, yv,
                    c / total_count if total_count > 0 else 0.0,
                    yv / total_y if total_y > 0 else 0.0,
                ))
        return Table([self.x_name, count_col, weight_col, frac_count_col, frac_weight_col], rows, name=self.name)

    def to_dict(self) -> dict:
        """The dict always carries "mn"/"mw"/"mz"/"pdi" (a stable schema
        existing consumers may already parse) -- for a CLD object
        (``meta["kind"] == "CLD"``), these are really DPn/DPw/DPz
        (chain-length averages, not molar mass): ``x`` for a CLD is DP,
        and ``from_pairs`` was given ``mass=DP`` to make DPw/DPz
        computable the same way Mw/Mz are for a real mass distribution.
        For CLD objects specifically, the dict *also* carries
        "dpn"/"dpw"/"dpz" (same numbers, physically correct keys) --
        additive, so nothing that already reads "mn"/"mw"/"mz" breaks."""
        out = {
            "name": self.name, "x_name": self.x_name, "basis": self.basis, "method": self.method,
            "bins": len(self.x), "n": self.n, "total_weight": self.total_weight,
            "mn": self.mn(), "mw": self.mw(), "mz": self.mz(), "pdi": self.pdi(),
            "meta": dict(self.meta),
            "rows": self.as_table().rows(),
        }
        if self.meta.get("kind") == "CLD":
            out["dpn"] = self.dpn()
            out["dpw"] = self.dpw()
            out["dpz"] = self.dpz()
        return out

    def _write_delimited(self, path: str | Path, *, delimiter: str) -> Path:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        table = self.as_table()
        with out.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(list(table.columns))
            for row in table:
                writer.writerow(list(row))
        return out

    def to_tsv(self, path: str | Path) -> Path:
        return self._write_delimited(path, delimiter="\t")



    def info_text(self) -> str:
        """Critical-analysis #25/26: for a CLD object (``meta["kind"] ==
        "CLD"``), ``.mn()``/``.mw()``/``.mz()`` are really DPn/DPw/DPz
        (chain-length averages) -- labeled accordingly here, rather than
        the misleading Mn/Mw/Mz used for an actual mass distribution."""
        is_cld = self.meta.get("kind") == "CLD"
        n_label, w_label, z_label = ("DPn", "DPw", "DPz") if is_cld else ("Mn", "Mw", "Mz")
        return "\n".join([
            f"distribution: {self.name}", f"axis: {self.x_name}", f"basis: {self.basis}",
            f"method: {self.method}", f"bins: {len(self.x)}", f"n: {self.n:g}",
            f"total_weight: {self.total_weight:g}",
            f"{n_label}: {self.mn():.6g}", f"{w_label}: {self.mw():.6g}",
            f"{z_label}: {self.mz():.6g}", f"PDI: {self.pdi():.6g}",
        ])

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text

    def help(self) -> str:
        text = (
            "Distribution helper examples:\n"
            "  d.mn(); d.mw(); d.mz(); d.pdi()   # always from the original linear-space data\n"
            "  d.to_tsv('out.tsv')\n"
            "  d.as_table()\n"
        )
        print(text)
        return text

    def plot(self, ax=None, *, normalize: bool = True, style: str | None = None,
             span: str | None = None, **kwargs):
        """``normalize=True`` (default) plots a shape-normalized view;
        ``normalize=False`` plots the raw count/weight (or density)
        values. For a KDE result, "normalized" means area-under-curve=1
        (accounting for the grid spacing), matching the KDE contract --
        for sticks/hist, it means fractions summing to 1 (the pre-rewrite
        version always did the sum=1 kind, with no way to see absolute
        numbers, and never distinguished the two conventions).

        For ``method in {"sticks", "chain_mass_spectrum"}``, draws real
        impulses (``ax.vlines``) by default -- previously this used
        circle markers (``ax.plot(x, y, "o")``), so a "sticks"
        distribution's defining visual feature (discrete impulses from
        the axis, not a scatter of points) never actually appeared on
        screen. The sticks representation uses true vertical impulses. Pass ``style=...``
        explicitly to override with a plain ``ax.plot()`` style string
        instead (e.g. ``style="o"`` for the old marker look)."""
        from .plotting import (apply_axes_style, available_styles, create_axes,
                               require_owned_geometry, style_kwargs)

        theme = style if style in available_styles() else "screen"
        require_owned_geometry(ax, span)
        if ax is None:
            _, ax = create_axes(theme, span=span)
        if normalize:
            if self.method == "kde":
                grid_step = self.meta.get("grid_step")
                total_area = self.total_y * grid_step if grid_step else self.total_y
                y = [v / total_area if total_area > 0 else 0.0 for v in self.y]
                ylabel = f"{self.y_label}_density (normalized, area=1)"
            else:
                total = self.total_y
                y = [v / total if total > 0 else 0.0 for v in self.y]
                # critical-analysis #15: this used to hardcode
                # "fraction_weight" for any non-"number" basis, so
                # basis="mass" plots were labeled as if they were
                # mass-weighted -- now basis-aware via self.y_label,
                # matching as_table()'s already-fixed column naming (P2.3).
                ylabel = "fraction_count" if self.basis == "number" else f"fraction_{self.y_label}"
        else:
            y = list(self.y)
            ylabel = self.y_label
        line_style = None if style is None or style in available_styles() else style
        defaults = style_kwargs(theme)
        defaults.update(kwargs)
        if line_style is not None:
            ax.plot(list(self.x), y, line_style, **defaults)
        elif self.method in {"sticks", "chain_mass_spectrum"}:
            ax.vlines(list(self.x), 0, y, **defaults)
        else:
            ax.plot(list(self.x), y, "-", **defaults)
        ax.set_xlabel(self.x_name)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{self.name} ({self.basis}, {self.method})")
        apply_axes_style(ax, theme)
        return ax

    def __repr__(self) -> str:
        return (f"Distribution(name={self.name!r}, x_name={self.x_name!r}, "
                f"bins={len(self.x)}, basis={self.basis!r}, method={self.method!r})")
