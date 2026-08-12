from __future__ import annotations

from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

from .column_family import _is_public_identifier


def _readonly(values, *, dtype=None) -> np.ndarray:
    array = np.asarray(values, dtype=dtype)
    array.flags.writeable = False
    return array


class NamedValues:
    """Small immutable mapping used by analysis result objects."""

    def __init__(self, values: Mapping[str, Any]):
        self._values = MappingProxyType(dict(values))

    def __getitem__(self, name: str):
        return self._values[name]

    def __getattr__(self, name: str):
        try:
            return self._values[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __iter__(self):
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)

    def keys(self) -> tuple[str, ...]:
        return tuple(self._values)

    def items(self):
        return tuple(self._values.items())

    @property
    def array(self) -> np.ndarray:
        if not self._values:
            return _readonly(np.empty((0,), dtype=float))
        return _readonly(np.stack([np.asarray(value) for value in self._values.values()], axis=-1))

    def __repr__(self) -> str:
        return f"NamedValues({', '.join(self._values)})"

    def __dir__(self) -> list[str]:
        standard = set(super().__dir__())
        dynamic = {name for name in self._values if _is_public_identifier(name) and name not in standard}
        return sorted(standard | dynamic)
