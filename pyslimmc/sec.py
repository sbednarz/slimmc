from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


def _readonly(values, *, dtype=float):
    a = np.asarray(values, dtype=dtype)
    a.flags.writeable = False
    return a


def _validate_positive_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and strictly positive")
    return value


@dataclass(frozen=True)
class SECDistribution:
    _x: np.ndarray
    _y: np.ndarray
    sigma_log10M: float
    mass_model: str
    total_chains: int
    snapshot_id: int | None
    t: float | None
    mn: float
    mw: float
    mz: float

    def __post_init__(self):
        self._x.flags.writeable = False
        self._y.flags.writeable = False

    @property
    def x(self):
        """log10 molar-mass coordinate."""
        return self._x

    @property
    def log10_mass(self):
        return self._x

    @property
    def mass(self):
        result = np.power(10.0, self._x)
        result.flags.writeable = False
        return result

    @property
    def y(self):
        """Apparent polymer-mass density dW_app/dlog10(M)."""
        return self._y

    @property
    def dispersity(self):
        return self.mw / self.mn if self.mn else float("nan")

    @property
    def is_empty(self):
        return self._x.size == 0

    @property
    def metadata(self):
        return {
            "kind": "SEC",
            "representation": "continuous",
            "coordinate": "log10_mass",
            "ordinate": "dW_app/dlog10M",
            "sigma_log10M": self.sigma_log10M,
            # Buback et al. (1996): for log10(M) = a - b v, the same
            # broadening width is written as b*sigma_v.
            "b_sigma_v_equivalent": self.sigma_log10M,
            "mass_model": self.mass_model,
        }

    @property
    def meta(self):
        return self.metadata

    def to_tsv(self, path: str | Path):
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="") as handle:
            handle.write("# representation: continuous SEC response\n")
            handle.write(f"# sigma_log10M: {self.sigma_log10M:.17g}\n")
            handle.write(f"# mass_model: {self.mass_model}\n")
            handle.write("log10_mass\tdW_app_dlog10M\n")
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
        ax.set_ylabel("Apparent polymer mass density")
        return ax

    def info(self):
        return (
            f"{type(self).__name__}\n"
            f"  representation: continuous SEC response\n"
            f"  coordinate: log10(M)\n"
            f"  ordinate: dW_app/dlog10(M)\n"
            f"  sigma_log10M: {self.sigma_log10M:.6g}\n"
            f"  mass_model: {self.mass_model}\n"
            f"  total_chains: {self.total_chains}"
        )

    def help(self):
        return self.info()


def build_sec(
    population,
    *,
    sigma_log10M: float,
    mass_model: str | None = None,
    step_log10M: float | None = None,
):
    """Apply the Buback-style Gaussian SEC response in log10(M).

    Exact chain masses are converted to exact polymer-mass fractions ``w_i``
    and treated as a discrete true log-MWD measure.  The apparent SEC trace is

        S(u) = sum_i w_i G_sigma(u - log10(M_i))

    with ``u = log10(M)`` and a normalized Gaussian ``G_sigma``.

    ``sigma_log10M`` is the Gaussian standard deviation in log10 molar-mass
    units.  For the linear SEC calibration used by Buback et al. (1996),
    ``log10(M) = a - b*v``, it is exactly the literature product ``b*sigma_v``.
    """
    sigma = _validate_positive_finite("sigma_log10M", sigma_log10M)

    counts = population.mass_counts(mass_model=mass_model)
    mass = np.asarray(counts.mass, dtype=float)
    count = np.asarray(counts.count, dtype=float)
    if mass.size == 0 or not np.sum(count):
        return SECDistribution(
            _readonly([]), _readonly([]), sigma, counts.mass_model,
            int(np.sum(count)), getattr(population, "snapshot_id", None),
            getattr(population, "t", None), float("nan"), float("nan"), float("nan"),
        )

    weighted_mass = mass * count
    total_mass = float(np.sum(weighted_mass))
    if not np.isfinite(total_mass) or total_mass <= 0.0:
        raise ValueError("SEC requires a finite, strictly positive total polymer mass")
    weights = weighted_mass / total_mass
    support = np.log10(mass)

    if step_log10M is None:
        step = sigma / 20.0
    else:
        step = _validate_positive_finite("step_log10M", step_log10M)

    lo = float(np.min(support) - 6.0 * sigma)
    hi = float(np.max(support) + 6.0 * sigma)
    n_points = max(2, int(np.ceil((hi - lo) / step)) + 1)
    if n_points > 1_000_000:
        raise ValueError(
            "SEC output grid would exceed 1,000,000 points; provide a larger step_log10M"
        )
    x = np.linspace(lo, hi, n_points, dtype=float)

    delta = (x[:, None] - support[None, :]) / sigma
    kernels = np.exp(-0.5 * delta * delta) / (sigma * np.sqrt(2.0 * np.pi))
    y = kernels @ weights

    # Do not renormalize the sampled curve.  The Gaussian mixture is already
    # normalized analytically; forcing numerical area to one would hide an
    # inadequate output grid.  Six-sigma padding makes truncation negligible
    # for the default grid.
    area = float(np.trapezoid(y, x))
    if not np.isfinite(area) or area <= 0.0:
        raise ValueError("SEC numerical response has a non-positive area")
    if abs(area - 1.0) > 5.0e-4:
        raise ValueError(
            "SEC output grid is too coarse to represent the Gaussian response "
            "accurately; use a smaller step_log10M"
        )

    # Moments remain properties of the exact source population, not of the
    # broadened apparent distribution.
    source_moments = population.moments(mass_model=counts.mass_model)

    return SECDistribution(
        _readonly(x), _readonly(y), sigma, counts.mass_model,
        source_moments.total_chains, int(population.snapshot_id), population.t,
        source_moments.mn, source_moments.mw, source_moments.mz,
    )
