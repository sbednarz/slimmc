from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .mass_model import record_masses


def _safe_ratio(num: float, den: float) -> float:
    return float(num / den) if den else float("nan")


@dataclass(frozen=True)
class PopulationMoments:
    """Exact moments for one selected chain population.

    DP moments and molar-mass moments are kept distinct.  ``dpz`` is NaN
    when the backing storage contains only aggregate moments without the
    third DP moment required to reconstruct it.
    """

    total_chains: int
    dpn: float
    dpw: float
    dpz: float
    mn: float
    mw: float
    mz: float
    dp_dispersity: float
    mass_dispersity: float
    mass_model: str
    snapshot_id: int | None = None
    t: float | None = None
    source: str = "chains"

    @property
    def has_dpz(self) -> bool:
        return math.isfinite(self.dpz)

    def info(self) -> str:
        dpz = f"{self.dpz:.8g}" if self.has_dpz else "unavailable"
        text = (
            "PopulationMoments\n"
            f"  source: {self.source}\n"
            f"  total_chains: {self.total_chains}\n"
            f"  mass_model: {self.mass_model}\n"
            f"  DPn / DPw / DPz: {self.dpn:.8g} / {self.dpw:.8g} / {dpz}\n"
            f"  Mn / Mw / Mz: {self.mn:.8g} / {self.mw:.8g} / {self.mz:.8g}\n"
            f"  DP dispersity: {self.dp_dispersity:.8g}\n"
            f"  mass dispersity: {self.mass_dispersity:.8g}"
        )
        print(text)
        return text

    def help(self) -> str:
        return self.info()


def calculate_population_moments(population, *, mass_model: str | None = None) -> PopulationMoments:
    """Calculate exact DP and molar-mass moments from selected chain records."""
    count = np.asarray(population.count, dtype=float)
    dp = np.asarray(population.dp, dtype=float)
    if count.shape != dp.shape:
        raise ValueError("chain dp and count arrays must have the same shape")
    if np.any(count < 0) or np.any(~np.isfinite(count)):
        raise ValueError("chain counts must be finite and non-negative")
    if np.any(dp <= 0) or np.any(~np.isfinite(dp)):
        raise ValueError("degree of polymerization must be finite and strictly positive")

    mass, model = record_masses(population, mass_model)

    n0 = float(np.sum(count))
    d1 = float(np.dot(count, dp))
    d2 = float(np.dot(count, dp * dp))
    d3 = float(np.dot(count, dp * dp * dp))
    m1 = float(np.dot(count, mass))
    m2 = float(np.dot(count, mass * mass))
    m3 = float(np.dot(count, mass * mass * mass))

    dpn = _safe_ratio(d1, n0)
    dpw = _safe_ratio(d2, d1)
    dpz = _safe_ratio(d3, d2)
    mn = _safe_ratio(m1, n0)
    mw = _safe_ratio(m2, m1)
    mz = _safe_ratio(m3, m2)

    return PopulationMoments(
        total_chains=int(round(n0)),
        dpn=dpn,
        dpw=dpw,
        dpz=dpz,
        mn=mn,
        mw=mw,
        mz=mz,
        dp_dispersity=_safe_ratio(dpw, dpn),
        mass_dispersity=_safe_ratio(mw, mn),
        mass_model=model,
        snapshot_id=int(population.snapshot_id) if getattr(population, "snapshot_id", None) is not None else None,
        t=float(population.t) if getattr(population, "t", None) is not None else None,
        source="chains",
    )
