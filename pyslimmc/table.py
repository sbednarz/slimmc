from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Iterator, Sequence

TABLE_HELP = """# A pyslimmc.Table is a small no-Pandas table reader.

table.shape
list(table.columns)
table.columns.tolist()
table['t'].iloc[0]
table.t.iloc[0]
table.row(0)
table.rows()
table.head(5)
table.tail(5)
table.to_numpy()
"""


class ColumnNames(tuple):
    """Tuple-like column-name container with a convenient tolist() helper."""

    def tolist(self) -> list[str]:
        return list(self)


class _ILoc:
    def __init__(self, values: Sequence[Any]):
        self._values = values

    def __getitem__(self, idx: int) -> Any:
        return self._values[idx]


@dataclass(frozen=True)
class Column:
    """A lightweight column view returned by ``Table['name']``."""

    name: str
    values: tuple[Any, ...]

    @property
    def iloc(self) -> _ILoc:
        return _ILoc(self.values)

    def tolist(self) -> list[Any]:
        return list(self.values)

    def to_numpy(self):
        import numpy as np

        result = np.array(self.values)
        result.flags.writeable = False
        return result

    def __len__(self) -> int:
        return len(self.values)

    def __iter__(self) -> Iterator[Any]:
        return iter(self.values)

    def __getitem__(self, idx: int) -> Any:
        return self.values[idx]

    def __repr__(self) -> str:
        return f"Column(name={self.name!r}, len={len(self.values)})"


class Table:
    """Small text table table used by pyslimmc.

    The core reader intentionally uses only the Python standard library.  It
    provides the small API needed for interactive inspection and tests: ``shape``,
    ``columns.tolist()``, column access, attribute access for identifier-like
    column names, ``equals()``, ``head()``, ``tail()`` and optional ``to_numpy()``.
    """

    def __init__(
        self,
        columns: Iterable[str],
        rows: Iterable[Iterable[Any]],
        name: str | None = None,
        source: str | None = None,
    ):
        self.columns = ColumnNames(str(c) for c in columns)
        self._rows = tuple(tuple(row) for row in rows)
        self.name = name or "table"
        self.source = source
        width = len(self.columns)
        for i, row in enumerate(self._rows):
            if len(row) != width:
                raise ValueError(f"row {i} has {len(row)} values, expected {width}")
        seen: set[str] = set()
        duplicates = {c for c in self.columns if c in seen or seen.add(c)}
        if duplicates:
            raise ValueError(
                f"duplicate column name(s) {sorted(duplicates)} in table {name or 'table'!r} "
                f"-- the earlier column(s) would silently become unreachable by name"
            )
        self._index = {name: i for i, name in enumerate(self.columns)}

    @property
    def shape(self) -> tuple[int, int]:
        return (len(self._rows), len(self.columns))

    def __len__(self) -> int:
        return len(self._rows)

    def __iter__(self) -> Iterator[tuple[Any, ...]]:
        return iter(self._rows)

    def __getitem__(self, key: str | int) -> Column | tuple[Any, ...]:
        if isinstance(key, int):
            return self._rows[key]
        if key not in self._index:
            raise KeyError(key)
        j = self._index[key]
        return Column(key, tuple(row[j] for row in self._rows))

    def __getattr__(self, name: str) -> Column:
        # Allow run.state.t and run.state.conv_A without committing to Pandas.
        if name in self._index and name.isidentifier():
            return self[name]  # type: ignore[return-value]
        raise AttributeError(name)

    def __dir__(self) -> list[str]:
        standard = set(super().__dir__())
        dynamic = {c for c in self.columns if c.isidentifier() and c not in standard}
        return sorted(standard | dynamic)

    def row(self, idx: int) -> dict[str, Any]:
        return dict(zip(self.columns, self._rows[idx]))

    def rows(self) -> list[dict[str, Any]]:
        return [self.row(i) for i in range(len(self._rows))]

    def head(self, n: int = 5) -> "Table":
        return Table(self.columns, self._rows[:n], name=self.name, source=self.source)

    def tail(self, n: int = 5) -> "Table":
        return Table(self.columns, self._rows[-n:], name=self.name, source=self.source)

    def equals(self, other: object) -> bool:
        return isinstance(other, Table) and self.columns == other.columns and self._rows == other._rows

    def to_numpy(self):
        import numpy as np

        result = np.array(self._rows, dtype=object)
        result.flags.writeable = False
        return result

    def info_text(self) -> str:
        lines = [
            f"table: {self.name}",
            f"shape: {self.shape[0]} rows x {self.shape[1]} columns",
        ]
        if self.source:
            lines.append(f"source: {self.source}")
        lines.append("columns: " + (", ".join(self.columns) if self.columns else "none"))
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text

    def help(self) -> str:
        print(TABLE_HELP)
        return TABLE_HELP

    def __repr__(self) -> str:
        nrow, ncol = self.shape
        preview = ", ".join(self.columns[:6])
        if ncol > 6:
            preview += ", ..."
        return f"Table(name={self.name!r}, shape=({nrow}, {ncol}), columns=[{preview}])"
