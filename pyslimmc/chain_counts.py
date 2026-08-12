from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .chains import ChainPopulation
from .table import Table


def _readonly(values, *, dtype=None) -> np.ndarray:
    result = np.asarray(values, dtype=dtype)
    result.flags.writeable = False
    return result


@dataclass(frozen=True)
class ChainCounts:
    """Exact, unnormalized chain count grouped by degree of polymerization."""

    dp: np.ndarray
    count: np.ndarray
    snapshot_id: int
    t: float | None
    pool: str

    def __post_init__(self) -> None:
        self.dp.flags.writeable = False
        self.count.flags.writeable = False

    @classmethod
    def from_population(cls, population: ChainPopulation, *, pool: str) -> "ChainCounts":
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
            _readonly(unique, dtype=np.int64), _readonly(totals, dtype=np.int64),
            int(population.snapshot_id), population.t, pool,
        )

    @property
    def x(self) -> np.ndarray:
        return self.dp

    @property
    def y(self) -> np.ndarray:
        return self.count

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

    def as_table(self) -> Table:
        return Table(("dp", "count"), zip(self.dp.tolist(), self.count.tolist()), name="chain_counts")

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
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError("ChainCounts.plot() requires optional dependency matplotlib") from exc
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


class ChainCountsGroup(Mapping[str, ChainCounts]):
    def __init__(self, spectra: Mapping[str, ChainCounts]):
        self._spectra = dict(spectra)

    def __getitem__(self, key: str) -> ChainCounts:
        return self._spectra[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._spectra)

    def __len__(self) -> int:
        return len(self._spectra)

    def plot(self, *, ax=None, path: str | Path | None = None, dpi: int = 300,
             style: str = "screen", span: str | None = None, **plot_kwargs):
        try:
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError("ChainCountsGroup.plot() requires optional dependency matplotlib") from exc
        from .plotting import apply_axes_style, create_axes, require_owned_geometry, style_kwargs
        require_owned_geometry(ax, span)
        if ax is None:
            _, ax = create_axes(style, span=span)
        for i, (name, spectrum) in enumerate(self._spectra.items()):
            kwargs: dict[str, Any] = {**style_kwargs(style, index=i), "label": name}
            kwargs.update(plot_kwargs)
            ax.vlines(spectrum.dp, 0, spectrum.count, **kwargs)
        ax.set_xlabel("DP")
        ax.set_ylabel("chain count")
        ax.legend()
        apply_axes_style(ax, style)
        if path is not None:
            ax.figure.savefig(path, dpi=dpi)
        return ax
