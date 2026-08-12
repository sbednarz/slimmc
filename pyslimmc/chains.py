from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any

import numpy as np

from .core import InvalidOutputError
from .operations import analysis_operation, MWD_HELP, CLD_HELP, SPECTRUM_HELP


def _readonly(values, *, dtype=None) -> np.ndarray:
    result = np.asarray(values, dtype=dtype)
    result.flags.writeable = False
    return result


class ChainCountsView(Mapping[str, Any]):
    """Per-chain repeat-unit counts, keyed by real monomer name."""

    def __init__(self, values: Mapping[str, Any]):
        self._values = MappingProxyType(dict(values))

    def __getitem__(self, name: str) -> Any:
        return self._values[name]

    def __getattr__(self, name: str) -> Any:
        try:
            return self._values[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __iter__(self) -> Iterator[str]:
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)

    @property
    def array(self) -> np.ndarray:
        if not self._values:
            return _readonly(np.empty((0,), dtype=np.int64))
        arrays = [np.asarray(value) for value in self._values.values()]
        result = np.stack(arrays, axis=-1)
        result.flags.writeable = False
        return result

    def __repr__(self) -> str:
        return f"ChainCountsView({', '.join(self._values)})"


@dataclass(frozen=True)
class ChainEndgroups:
    left: Any
    right: Any
    _population: "ChainPopulation | None" = None

    def __call__(self, left: str, right: str, *, ordered: bool = False) -> "ChainPopulation":
        if self._population is None:
            raise TypeError("end groups of a single ChainRow are scalar values, not a population selector")
        return self._population._select_endgroups(left, right, ordered=ordered)

    def summary(self):
        if self._population is None:
            raise TypeError("end groups of a single ChainRow do not have a population summary")
        from .table import Table
        groups = {}
        for left, right, count in zip(self.left, self.right, self._population.count):
            key = (str(left), str(right))
            rec = groups.setdefault(key, [0, 0])
            rec[0] += 1
            rec[1] += int(count)
        total = sum(v[1] for v in groups.values())
        rows = [(l, r, nr, nc, nc/total if total else float("nan")) for (l, r), (nr, nc) in sorted(groups.items())]
        return Table(["left_end", "right_end", "n_records", "n_chains", "fraction"], rows, name="endgroups_summary")


@dataclass(frozen=True)
class ChainRow(Mapping[str, Any]):
    """One compressed structural chain type, not one physical molecule.

    ``count`` tells how many physical chains this row represents.
    ``counts[name]`` is the number of repeat units of monomer ``name`` in
    each represented chain.
    """

    _population: "ChainPopulation"
    _index: int

    @property
    def dp(self) -> int:
        return int(self._population.dp[self._index])

    @property
    def count(self) -> int:
        return int(self._population.count[self._index])

    @property
    def counts(self) -> ChainCountsView:
        return ChainCountsView({name: int(values[self._index]) for name, values in self._population.counts.items()})

    def counts_total(self) -> int:
        return int(sum(self.counts.values()))

    @property
    def endgroups(self) -> ChainEndgroups:
        ends = self._population.endgroups
        return ChainEndgroups(str(ends.left[self._index]), str(ends.right[self._index]))

    def as_dict(self) -> dict[str, Any]:
        result = self._population._row_dict(self._index)
        result.setdefault("dp", self.dp)
        if "count" in self._population._raw_arrays():
            result.setdefault("count", self.count)
        return result

    def __getitem__(self, key: str) -> Any:
        return self.as_dict()[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self.as_dict())

    def __len__(self) -> int:
        return len(self.as_dict())

    def __getattr__(self, name: str) -> Any:
        data = self.as_dict()
        if name in data:
            return data[name]
        raise AttributeError(name)


class ChainPopulation:
    """Shared read-only contract for one exact saved chain population.

    Concrete homo/copo snapshot classes retain their engine-specific fields
    and filters, while this base supplies the chemically common surface.
    """

    def _raw_arrays(self) -> Mapping[str, np.ndarray]:
        raise NotImplementedError

    def _count_arrays(self) -> Mapping[str, np.ndarray]:
        raise NotImplementedError

    @property
    def count(self) -> np.ndarray:
        try:
            result = self._raw_arrays()["count"]
        except KeyError as exc:
            raise InvalidOutputError("chain population has no compressed 'count' column") from exc
        result.flags.writeable = False
        return result

    @property
    def counts(self) -> ChainCountsView:
        return ChainCountsView(self._count_arrays())

    def counts_total(self) -> np.ndarray:
        arrays = tuple(self._count_arrays().values())
        if not arrays:
            result = np.zeros(len(self), dtype=np.int64)
        else:
            result = np.sum(np.stack(arrays, axis=0), axis=0, dtype=np.int64)
        result = _readonly(result, dtype=np.int64)
        if not np.array_equal(result, np.asarray(self.dp, dtype=np.int64)):
            raise InvalidOutputError("chain composition counts do not sum to dp")
        return result

    @property
    def endgroups(self) -> ChainEndgroups:
        raw = self._raw_arrays()
        left = raw.get("left_end", raw.get("eg1"))
        right = raw.get("right_end", raw.get("eg2"))
        return ChainEndgroups(left, right, self)

    def _select_endgroups(self, left: str, right: str, *, ordered: bool = False) -> "ChainPopulation":
        ends = self.endgroups
        if ordered or left == right:
            mask = (ends.left == left) & (ends.right == right)
        else:
            mask = ((ends.left == left) & (ends.right == right)) | (
                (ends.left == right) & (ends.right == left)
            )
        return self._with_mask(mask)

    @property
    def total_chains(self) -> int:
        return int(np.sum(self.count, dtype=np.int64))

    @property
    def total_repeat_units(self) -> int:
        return int(np.dot(np.asarray(self.dp, dtype=np.int64), self.count.astype(np.int64)))

    @property
    def compressed_rows(self) -> int:
        return len(self)

    def all(self) -> "ChainPopulation":
        return self

    @analysis_operation(CLD_HELP)
    def cld(self, *, mass_model: str | None = None, series=None, **kwargs):
        """Build a CLD directly from this exact population snapshot."""
        if series is not None:
            from .spectra import build_cld_series
            return build_cld_series(self, series=series, mass_model=mass_model, **kwargs)
        kwargs.pop("reference", None)
        basis = str(kwargs.get("basis", "number")).lower()
        masses = None
        if basis in {"mass", "weight"}:
            raw = self._raw_arrays()
            if "mass" in raw:
                masses = raw["mass"]
            else:
                calculator = getattr(self, "masses", None)
                if calculator is None:
                    raise InvalidOutputError("mass-basis CLD requires per-chain molar masses")
                masses = calculator(mass_model=mass_model or "repeat_units")
        from .spectra import build_cld
        return build_cld(self, masses=masses, mass_model=mass_model, **kwargs)

    @analysis_operation(MWD_HELP)
    def mwd(self, *, mass_model: str | None = None, series=None, **kwargs):
        """Build an MWD directly from this exact population snapshot.

        Copolymer snapshots carry engine-computed per-chain masses.  Classic
        snapshots compute them from the selected mass model.
        """
        if series is not None:
            from .spectra import build_mwd_series
            return build_mwd_series(self, series=series, mass_model=mass_model, **kwargs)
        kwargs.pop("reference", None)
        raw = self._raw_arrays()
        if "mass" in raw:
            masses = raw["mass"]
        else:
            calculator = getattr(self, "masses", None)
            if calculator is None:
                raise InvalidOutputError("chain population has no molar-mass data")
            masses = calculator(mass_model=mass_model or "repeat_units")
        from .spectra import build_mwd
        return build_mwd(self, masses, mass_model=mass_model, **kwargs)

    @analysis_operation(SPECTRUM_HELP)
    def chain_mass_spectrum(self, *, mass_model: str | None = None, series=None, normalize: str = "count", **kwargs):
        """Exact neutral-chain mass spectrum (not m/z or ionization)."""
        if series is not None:
            from .spectra import build_chain_mass_spectrum_series
            return build_chain_mass_spectrum_series(self, series=series, mass_model=mass_model, normalize=normalize, **kwargs)
        raw = self._raw_arrays()
        if "mass" in raw:
            masses = raw["mass"]
        else:
            calculator = getattr(self, "masses", None)
            if calculator is None:
                raise InvalidOutputError("chain population has no molar-mass data")
            masses = calculator(mass_model=mass_model or "repeat_units")
        from .spectra import build_chain_mass_spectrum
        return build_chain_mass_spectrum(self, masses, mass_model=mass_model, normalize=normalize, **kwargs)

    def row(self, index: int) -> ChainRow:
        if not isinstance(index, int) or isinstance(index, bool):
            raise TypeError("chain row index must be an integer")
        if index < 0:
            index += len(self)
        if index < 0 or index >= len(self):
            raise IndexError(index)
        return ChainRow(self, index)

    def rows(self) -> tuple[ChainRow, ...]:
        return tuple(self.row(i) for i in range(len(self)))

    def _row_dict(self, index: int) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for name, values in self._raw_arrays().items():
            value = values[index]
            result[name] = value.item() if isinstance(value, np.generic) else value
        return result

    def __getattr__(self, name: str):
        if name.startswith("counts_"):
            monomer = name[len("counts_"):]
            try:
                return self.counts[monomer]
            except KeyError as exc:
                raise AttributeError(name) from exc
        raise AttributeError(name)

    def __dir__(self) -> list[str]:
        standard = set(super().__dir__())
        dynamic = {
            f"counts_{name}" for name in self._count_arrays()
            if name.isidentifier()
        }
        return sorted(standard | dynamic)
