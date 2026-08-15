from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .chains import ChainPopulation


def _readonly(values, *, dtype=None) -> np.ndarray:
    result = np.asarray(values, dtype=dtype)
    result.flags.writeable = False
    return result


@dataclass(frozen=True)
class DPCounts:
    """Exact, unnormalized chain count grouped by degree of polymerization."""

    dp: np.ndarray
    count: np.ndarray
    snapshot_id: int
    t: float | None
    pool: str = "selected"

    def __post_init__(self) -> None:
        self.dp.flags.writeable = False
        self.count.flags.writeable = False

    @classmethod
    def from_population(cls, population: ChainPopulation, *, pool: str = "selected") -> "DPCounts":
        dp = np.asarray(population.dp, dtype=np.int64)
        counts = np.asarray(population.count, dtype=np.int64)
        if dp.size:
            unique, inverse = np.unique(dp, return_inverse=True)
            totals = np.zeros(unique.size, dtype=np.int64)
            np.add.at(totals, inverse, counts)
        else:
            unique = np.empty(0, dtype=np.int64)
            totals = np.empty(0, dtype=np.int64)
        return cls(
            _readonly(unique, dtype=np.int64),
            _readonly(totals, dtype=np.int64),
            int(population.snapshot_id),
            population.t,
            pool,
        )

    @property
    def total_chains(self) -> int:
        return int(np.sum(self.count, dtype=np.int64))

    @property
    def total_repeat_units(self) -> int:
        return int(np.dot(self.dp.astype(np.int64), self.count.astype(np.int64)))

    @property
    def min_dp(self) -> int | None:
        return int(self.dp.min()) if self.dp.size else None

    @property
    def max_dp(self) -> int | None:
        return int(self.dp.max()) if self.dp.size else None

    @property
    def is_empty(self) -> bool:
        return self.dp.size == 0

    def to_tsv(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="") as handle:
            handle.write(f"# snapshot_id: {self.snapshot_id}\n# pool: {self.pool}\n")
            handle.write("dp\tcount\n")
            for dp, count in zip(self.dp, self.count):
                handle.write(f"{int(dp)}\t{int(count)}\n")
        return target

    def plot(self, *, ax=None, path: str | Path | None = None, dpi: int = 300,
             style: str = "screen", span: str | None = None, **plot_kwargs):
        try:
            import matplotlib.pyplot as plt  # noqa: F401
        except ImportError as exc:
            raise ImportError("DPCounts.plot() requires optional dependency matplotlib") from exc
        from .plotting import apply_axes_style, create_axes, require_owned_geometry, style_kwargs
        require_owned_geometry(ax, span)
        if ax is None:
            _, ax = create_axes(style, span=span)
        kwargs = style_kwargs(style)
        kwargs.update(plot_kwargs)
        ax.vlines(self.dp, 0, self.count, **kwargs)
        ax.set_xlabel("DP")
        ax.set_ylabel("chain count")
        apply_axes_style(ax, style)
        if path is not None:
            ax.figure.savefig(path, dpi=dpi)
        return ax


@dataclass(frozen=True)
class MassCounts:
    """Exact, unnormalized chain count grouped by neutral chain molar mass."""

    mass: np.ndarray
    count: np.ndarray
    snapshot_id: int
    t: float | None
    mass_model: str
    pool: str = "selected"

    def __post_init__(self) -> None:
        self.mass.flags.writeable = False
        self.count.flags.writeable = False

    @classmethod
    def from_population(
        cls,
        population: ChainPopulation,
        masses,
        *,
        mass_model: str,
        pool: str = "selected",
    ) -> "MassCounts":
        masses = np.asarray(masses, dtype=float)
        counts = np.asarray(population.count, dtype=np.int64)
        if masses.shape != counts.shape:
            raise ValueError("per-chain masses and counts must have the same shape")
        if masses.size:
            unique, inverse = np.unique(masses, return_inverse=True)
            totals = np.zeros(unique.size, dtype=np.int64)
            np.add.at(totals, inverse, counts)
        else:
            unique = np.empty(0, dtype=float)
            totals = np.empty(0, dtype=np.int64)
        return cls(
            _readonly(unique, dtype=float),
            _readonly(totals, dtype=np.int64),
            int(population.snapshot_id),
            population.t,
            mass_model,
            pool,
        )

    @property
    def total_chains(self) -> int:
        return int(np.sum(self.count, dtype=np.int64))

    @property
    def min_mass(self) -> float | None:
        return float(self.mass.min()) if self.mass.size else None

    @property
    def max_mass(self) -> float | None:
        return float(self.mass.max()) if self.mass.size else None

    @property
    def is_empty(self) -> bool:
        return self.mass.size == 0

    def to_tsv(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", encoding="utf-8", newline="") as handle:
            handle.write(
                f"# snapshot_id: {self.snapshot_id}\n"
                f"# pool: {self.pool}\n"
                f"# mass_model: {self.mass_model}\n"
            )
            handle.write("mass\tcount\n")
            for mass, count in zip(self.mass, self.count):
                handle.write(f"{float(mass):.17g}\t{int(count)}\n")
        return target

    def plot(self, *, ax=None, path: str | Path | None = None, dpi: int = 300,
             style: str = "screen", span: str | None = None, **plot_kwargs):
        try:
            import matplotlib.pyplot as plt  # noqa: F401
        except ImportError as exc:
            raise ImportError("MassCounts.plot() requires optional dependency matplotlib") from exc
        from .plotting import apply_axes_style, create_axes, require_owned_geometry, style_kwargs
        require_owned_geometry(ax, span)
        if ax is None:
            _, ax = create_axes(style, span=span)
        kwargs = style_kwargs(style)
        kwargs.update(plot_kwargs)
        ax.vlines(self.mass, 0, self.count, **kwargs)
        ax.set_xlabel("Molar mass, g mol$^{-1}$")
        ax.set_ylabel("chain count")
        apply_axes_style(ax, style)
        if path is not None:
            ax.figure.savefig(path, dpi=dpi)
        return ax
