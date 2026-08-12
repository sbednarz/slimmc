from __future__ import annotations

import keyword
from typing import Any


def _is_public_identifier(name: str) -> bool:
    return bool(name) and name.isidentifier() and not keyword.iskeyword(name) and not name.startswith("_")


class ColumnFamily:
    """Named-column-group accessor shared by homo's ``History`` and
    copo's ``StateTable``: ``family.LABEL`` looks up
    ``table[pattern.format(label=LABEL)]`` -- e.g. ``history.c.M`` ->
    ``table["c_M"]`` (homo, prefix-only pattern ``"c_{label}"``),
    ``state.conv.A`` -> ``table["conv_A"]`` (copo, prefix pattern), or
    ``state.balance.A`` -> ``table["balance_A_count"]`` (copo, a pattern
    with both a prefix and a suffix around ``{label}``).

    Previously two separate, near-identical implementations (homo's
    ``_PrefixView``, copo's ``ColumnFamily``) doing the same job. Unified
    into one: ``pattern`` takes the general prefix+``{label}``+suffix
    shape (a bare prefix like ``"c_{label}"`` is just the suffix-empty
    case), and labels are always discovered dynamically by scanning the
    underlying table's columns for ones matching the pattern shape -- no
    static label list needs to be precomputed or kept in sync (this was
    already homo's more general approach; copo's previous version needed
    a precomputed ``labels`` tuple passed in at construction, which is no
    longer necessary).

    Works with any ``table`` that supports ``__getitem__(column_name)``
    raising ``KeyError`` for a missing column, and a ``.columns``
    sequence of column names -- both ``History`` (returns raw numpy
    arrays per column) and ``Table``/``StateTable`` (returns ``Column``
    objects) satisfy this without any special-casing here."""

    def __init__(self, table: Any, pattern: str, *, restrict_to: tuple[str, ...] | None = None):
        if pattern.count("{label}") != 1:
            raise ValueError(f"ColumnFamily pattern must contain exactly one '{{label}}', got {pattern!r}")
        self._table = table
        self._pattern = pattern
        self._prefix, self._suffix = pattern.split("{label}")
        # Needed when this family's prefix/suffix combination isn't
        # unambiguous on its own -- e.g. copo's "{label}_count" (empty
        # prefix, suffix "_count") would otherwise also match
        # "poly_A_count"/"balance_A_count" (which end in "_count" too),
        # silently producing spurious labels "poly_A"/"balance_A"
        # alongside the real "A"/"B". Pass the known-good label set
        # (e.g. StateTable.monomers) to filter the dynamic scan against,
        # instead of trusting prefix/suffix matching alone in these
        # cases.
        self._restrict_to = set(restrict_to) if restrict_to is not None else None

    def column_name(self, label: str) -> str:
        return self._pattern.format(label=label)

    def __getitem__(self, label: str) -> Any:
        return self._table[self.column_name(label)]

    def __getattr__(self, label: str) -> Any:
        if label.startswith("_"):
            raise AttributeError(label)
        try:
            return self[label]
        except KeyError as exc:
            raise AttributeError(f"no column {self.column_name(label)!r}") from exc

    def keys(self) -> tuple[str, ...]:
        tokens: list[str] = []
        for name in self._table.columns:
            if not name.startswith(self._prefix):
                continue
            if self._suffix and not name.endswith(self._suffix):
                continue
            end = len(name) - len(self._suffix) if self._suffix else len(name)
            token = name[len(self._prefix):end]
            if token:
                tokens.append(token)
        if self._restrict_to is not None:
            tokens = [t for t in tokens if t in self._restrict_to]
        return tuple(sorted(set(tokens)))

    def items(self) -> list[tuple[str, Any]]:
        return [(label, self[label]) for label in self.keys()]

    def final(self) -> dict[str, Any]:
        """Last-row value per label -- works uniformly for both a raw
        numpy array (homo's History columns) and a Column object (copo's
        Table columns), since both support plain ``[-1]`` indexing."""
        if len(self._table) == 0:
            return {label: None for label in self.keys()}
        return {label: self[label][-1] for label in self.keys()}

    def __dir__(self) -> list[str]:
        standard = set(super().__dir__())
        dynamic = {name for name in self.keys() if _is_public_identifier(name) and name not in standard}
        return sorted(standard | dynamic)

    def __repr__(self) -> str:
        return f"ColumnFamily(pattern={self._pattern!r}, keys={list(self.keys())!r})"
