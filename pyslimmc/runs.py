from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator
from fnmatch import fnmatchcase
import json
import re

from .table import Table
from .run import Run, DataConsistencyError
from .core import InvalidOutputError
from .column_family import _is_public_identifier
from .operations import analysis_operation

METADATA_FILES = ("run_metadata.json",)

_NUM_RE = re.compile(r"(?<![A-Za-z_])[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")
_QUOTED_RE = re.compile(r'"(?:[^"\\]|\\.)*"|\'(?:[^\'\\]|\\.)*\'')


class SelectionError(KeyError):
    """Raised when a requested run selection is ambiguous or invalid."""


def _extract_star(text: str, pattern: str, *, context: str) -> str:
    """Extract the substring matched by the single ``*`` in ``pattern``."""
    stars = pattern.count("*")
    if stars != 1:
        raise ValueError(f"{context} pattern must contain exactly one '*'; got {pattern!r}")
    prefix, suffix = pattern.split("*", 1)
    if not text.startswith(prefix) or not text.endswith(suffix):
        raise SelectionError(
            f"{context} pattern {pattern!r} does not match run_id {text!r}"
        )
    start = len(prefix)
    end = len(text) - len(suffix) if suffix else len(text)
    if end < start:
        raise SelectionError(
            f"{context} pattern {pattern!r} does not match run_id {text!r}"
        )
    return text[start:end]


# --------------------------------------------------------------------------
# discovery
# --------------------------------------------------------------------------

def _unique_run_dirs(root: Path, recursive: bool = True) -> list[Path]:
    """Find directories containing a complete Slimmc Storage identity."""
    if root.is_file():
        root = root.parent
    dirs: set[Path] = set()
    candidates = root.rglob("run_metadata.json") if recursive else root.glob("*/run_metadata.json")
    if (root / "run_metadata.json").is_file() and (root / "schema.jsonl").is_file():
        dirs.add(root)
    for metadata in candidates:
        directory = metadata.parent
        if (directory / "schema.jsonl").is_file():
            dirs.add(directory)
    return sorted(dirs, key=str)


def scan(path: str | Path = ".", *, recursive: bool = True, skip_bad: bool = False,
         engine_filter: str | None = None) -> "Runs":
    """Scan a directory tree for Slimmc Storage runs only."""
    root = Path(path).expanduser().resolve()
    found: dict[str, Run] = {}
    from ._storage import open_storage
    for run_dir in _unique_run_dirs(root, recursive=recursive):
        try:
            # Collections represent the lifecycle of a study, including
            # diagnostic failed/interrupted runs. Direct ``pyslimmc.open``
            # remains strict unless the caller passes allow_incomplete=True.
            run = open_storage(run_dir, allow_incomplete=True)
            if engine_filter is not None and run.engine != engine_filter:
                continue
        except (FileNotFoundError, OSError, ValueError, KeyError, DataConsistencyError) as exc:
            if skip_bad:
                print(f"[pyslimmc.scan] skipping {run_dir}: {exc}")
                continue
            raise InvalidOutputError(f"cannot open Slimmc run {run_dir}: {exc}") from exc
        try:
            relative_dir = run_dir.relative_to(root).as_posix() or "."
        except ValueError:
            relative_dir = "."
        run.relative_dir = relative_dir
        found[str(run.path)] = run
    return Runs(root, found)


# --------------------------------------------------------------------------
# model_diff -- copo's normalized-line algorithm, used for both engines
# (objectively better than classic slimmc's old raw-line-presence version:
# it recognizes "same directive, different value" instead of treating
# every numeric variant as an unrelated line -- see discrepancy table).
# --------------------------------------------------------------------------

def _model_path(run: Run) -> Path:
    return run.path / "input.model"


def _strip_comment(line: str) -> str:
    quote = None
    escaped = False
    out = []
    for ch in line:
        if escaped:
            out.append(ch)
            escaped = False
            continue
        if ch == "\\" and quote:
            out.append(ch)
            escaped = True
            continue
        if ch in {'"', "'"}:
            if quote is None:
                quote = ch
            elif quote == ch:
                quote = None
            out.append(ch)
            continue
        if ch == "#" and quote is None:
            break
        out.append(ch)
    return "".join(out).strip()


def _model_lines(run: Run) -> list[str]:
    model_path = _model_path(run)
    if not model_path.exists():
        return []
    lines = []
    for raw in model_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = _strip_comment(raw)
        if line:
            lines.append(" ".join(line.split()))
    return lines


def _line_key(line: str) -> str:
    collapsed = _QUOTED_RE.sub('"<str>"', line)
    collapsed = _NUM_RE.sub("<num>", collapsed)
    return collapsed


def _make_model_diff(runs: "Runs", *, include_same: bool = False) -> Table:
    run_list = list(runs)
    run_ids = [run.run_id for run in run_list]
    per_run: list[dict[str, str]] = []
    all_keys: list[str] = []
    seen_global: set[str] = set()

    for run in run_list:
        counts: dict[str, int] = {}
        mapping: dict[str, str] = {}
        for line in _model_lines(run):
            base = _line_key(line)
            counts[base] = counts.get(base, 0) + 1
            key = base if counts[base] == 1 else f"{base} [{counts[base]}]"
            mapping[key] = line
            if key not in seen_global:
                all_keys.append(key)
                seen_global.add(key)
        per_run.append(mapping)

    rows = []
    for key in all_keys:
        values = [m.get(key, "") for m in per_run]
        if include_same or len(set(values)) > 1:
            rows.append([key, *values])

    return Table(["item", *run_ids], rows, name="model_diff", source=str(runs.root))


def _is_float_like(value: object) -> bool:
    try:
        float(value)  # type: ignore[arg-type]
        return True
    except (TypeError, ValueError):
        return False


def _value_matches(actual: object, wanted: object) -> bool:
    if isinstance(actual, (int, float)) and isinstance(wanted, (int, float)):
        import math
        return math.isclose(float(actual), float(wanted))
    return actual == wanted


def _variable(run: Run, name: str):
    try:
        return run.var[name]
    except KeyError:
        return None


# --------------------------------------------------------------------------
# dynamic, tab-completable, "never error, always a subcollection" indices
# --------------------------------------------------------------------------

class _DynamicSubcollectionIndex:
    """Base for .var / .prefix / .run_id -- ``index[key]`` always returns a
    ``Runs`` subcollection (0, 1, or many matches), never raises for an
    ambiguous match (see the discrepancy table's Runs/Sweep resolution:
    ambiguity is completely normal here, e.g. the same file prefix reused
    across two separate experiments)."""

    _field: str

    def __init__(self, collection: "Runs"):
        self._collection = collection

    def _match(self, run: Run, key: str) -> bool:
        raise NotImplementedError

    def keys(self) -> tuple[str, ...]:
        raise NotImplementedError

    def __getitem__(self, key: str) -> "Runs":
        selected = {str(r.path): r for r in self._collection if self._match(r, key)}
        return Runs(self._collection.root, selected, _sweep_variables=self._collection._sweep_variables, _order_paths=tuple(selected))

    def __getattr__(self, name: str):
        if name.startswith("_"):
            raise AttributeError(name)
        return self[name]

    def __dir__(self) -> list[str]:
        standard = set(super().__dir__())
        dynamic = {k for k in self.keys() if _is_public_identifier(k) and k not in standard}
        return sorted(standard | dynamic)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(keys={list(self.keys())!r})"


class PrefixIndex(_DynamicSubcollectionIndex):
    """``runs.prefix["case_A"]`` -- matches on the raw model-file prefix."""

    def _match(self, run: Run, key: str) -> bool:
        return run.prefix == key

    def keys(self) -> tuple[str, ...]:
        return tuple(sorted({r.prefix for r in self._collection}))


class RunIdIndex(_DynamicSubcollectionIndex):
    """``runs.run_id["case_A"]`` -- matches on the engine-reported run_id."""

    def _match(self, run: Run, key: str) -> bool:
        return run.run_id == key

    def keys(self) -> tuple[str, ...]:
        return tuple(sorted({r.run_id for r in self._collection}))


class VarIndex(_DynamicSubcollectionIndex):
    """``runs.var["kp_ab"]`` -- matches runs that declared this variable
    name via ``var <kind> <name> <unit>``. ``[value]`` on the result (a
    Runs subcollection) is not directly supported here -- for a single,
    ordered, value->run mapping, use ``runs.sweep("kp_ab")`` instead."""

    def _match(self, run: Run, key: str) -> bool:
        return key in run.var

    def keys(self) -> tuple[str, ...]:
        return tuple(sorted({name for r in self._collection for name in r.var.keys()}))


@dataclass
class Runs:
    """Collection of Run objects returned by ``scan()``, spanning both
    engines transparently.

    Indexing:
    - ``runs[path_string]`` -- exact match on the collection's underlying
      key (``str(run.path)``), the only value guaranteed unique here.
    - ``runs[0]`` / ``runs[1:3]`` -- position in a stable, deterministic
      order: sorted ascending alphabetically by ``run_id`` (see the
      discrepancy table's Runs/Sweep resolution -- filesystem discovery
      order is not deterministic across platforms, so this is sorted once
      here rather than left to rglob()'s arbitrary order).
    - ``runs.prefix[...]`` / ``runs.run_id[...]`` / ``runs.var[...]`` --
      always return a subcollection, never raise on ambiguity (a repeated
      prefix/run_id/var name across separate experiments is normal, not
      an error condition).
    """

    root: Path
    _runs: dict[str, Run]  # keyed by str(run.path) -- naturally unique
    _run_id_index_cache: dict[str, tuple[Run, ...]] | None = field(
        default=None, init=False, repr=False, compare=False
    )
    _sweep_variables: tuple[str, ...] = field(default_factory=tuple, repr=False, compare=False)
    _order_paths: tuple[str, ...] | None = field(default=None, repr=False, compare=False)
    _selection_note: str | None = field(default=None, repr=False, compare=False)

    def _run_id_index(self) -> dict[str, tuple[Run, ...]]:
        """Build the lightweight run_id index only on first interactive use."""
        if self._run_id_index_cache is None:
            grouped: dict[str, list[Run]] = {}
            for run in self._runs.values():
                grouped.setdefault(run.run_id, []).append(run)
            self._run_id_index_cache = {
                run_id: tuple(sorted(items, key=lambda r: str(r.path)))
                for run_id, items in grouped.items()
            }
        return self._run_id_index_cache

    def __getattr__(self, name: str) -> Run:
        """Lazy interactive access: ``runs.<run_id>`` returns one Run.

        Real attributes and methods are resolved by Python before this fallback,
        so the public Runs API always wins on a name collision.
        """
        if name.startswith("_"):
            raise AttributeError(name)
        matches = self._run_id_index().get(name, ())
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            paths = "\n  ".join(str(run.path) for run in matches)
            raise SelectionError(
                f"run_id {name!r} is not unique across {len(matches)} runs:\n  {paths}\n"
                f"Use runs.run_id[{name!r}] for the subcollection or runs.one(...) with more filters."
            )
        raise AttributeError(f"{type(self).__name__!s} has no attribute or unique run_id {name!r}")

    def __dir__(self) -> list[str]:
        standard = set(super().__dir__())
        dynamic = {
            run_id
            for run_id, matches in self._run_id_index().items()
            if len(matches) == 1
            and _is_public_identifier(run_id)
            and run_id not in standard
        }
        return sorted(standard | dynamic)

    def _ordered(self) -> list[Run]:
        if self._order_paths is not None:
            ordered = [self._runs[path] for path in self._order_paths if path in self._runs]
            known = set(self._order_paths)
            ordered.extend(
                sorted(
                    (run for path, run in self._runs.items() if path not in known),
                    key=lambda run: (run.run_id, str(run.path)),
                )
            )
            return ordered
        # Default deterministic order for ordinary scan collections.
        return sorted(self._runs.values(), key=lambda r: (r.run_id, str(r.path)))

    def __len__(self) -> int:
        return len(self._runs)

    def __iter__(self) -> Iterator[Run]:
        return iter(self._ordered())

    def __getitem__(self, key):
        if isinstance(key, slice):
            selected = self._ordered()[key]
            return Runs(self.root, {str(r.path): r for r in selected}, _sweep_variables=self._sweep_variables, _order_paths=tuple(str(r.path) for r in selected), _selection_note=self._selection_note)
        if isinstance(key, int):
            return self._ordered()[key]
        run_id = str(key)
        matches = self._run_id_index().get(run_id, ())
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            paths = "\n  ".join(str(run.path) for run in matches)
            raise SelectionError(f"run_id {run_id!r} is not unique across {len(matches)} runs:\n  {paths}")
        raise KeyError(f"unknown run_id {run_id!r} in this collection")

    def by_path(self, path: str | Path) -> Run:
        """Return the run at an exact filesystem path."""
        resolved = str(Path(path).expanduser().resolve())
        try:
            return self._runs[resolved]
        except KeyError as exc:
            raise KeyError(f"unknown run path {path!r} in this collection") from exc

    @property
    def paths(self) -> list[Path]:
        return [r.path for r in self]

    @property
    def engines(self) -> list[str]:
        return sorted({r.engine for r in self._runs.values()})

    @property
    def schemas(self) -> list[str]:
        return sorted({r.schema for r in self._runs.values()})

    @property
    def var(self) -> VarIndex:
        return VarIndex(self)

    @property
    def prefix(self) -> PrefixIndex:
        return PrefixIndex(self)

    @property
    def run_id(self) -> RunIdIndex:
        return RunIdIndex(self)

    @analysis_operation("""Order runs as a one- or multidimensional parameter sweep.

Usage:
    sweep = runs.sweep("IA")
    sweep = runs.sweep("IA", "temperature")
    sweep = runs.match("f5a_*").sweep("IA", "temperature")

The variables are read only from run.var / run_metadata.json["variables"].
Every run must declare every requested variable. Runs are sorted
lexicographically by the real numeric variable values. Duplicate parameter
points are allowed (for example, replicate seeds) and remain in stable order.
The result is another Runs collection, so slicing, match(), filter(), pack(),
model_diff(), info(), and normal iteration continue to work.
""")
    def sweep(self, *variables: str) -> "Runs":
        """Return a Runs collection ordered by one or more declared variables."""
        if not variables:
            raise TypeError("sweep() requires at least one variable name")
        if any(not isinstance(name, str) or not name for name in variables):
            raise TypeError("sweep variable names must be non-empty strings")
        if len(set(variables)) != len(variables):
            raise ValueError("sweep variable names must be unique")

        missing: dict[str, list[str]] = {name: [] for name in variables}
        for run in self:
            for name in variables:
                if name not in run.var:
                    missing[name].append(run.run_id)
        missing = {name: ids for name, ids in missing.items() if ids}
        if missing:
            details = "; ".join(
                f"{name!r} missing in {len(ids)} run(s): {', '.join(ids[:5])}"
                + (" ..." if len(ids) > 5 else "")
                for name, ids in missing.items()
            )
            raise SelectionError(f"cannot build sweep: {details}")

        ordered = sorted(
            self,
            key=lambda run: (
                tuple(float(run.var[name].value) for name in variables),
                run.run_id,
                str(run.path),
            ),
        )
        return Runs(
            self.root,
            {str(run.path): run for run in ordered},
            _sweep_variables=tuple(variables),
            _order_paths=tuple(str(run.path) for run in ordered),
            _selection_note=self._selection_note,
        )

    @property
    def sweep_variables(self) -> tuple[str, ...]:
        """Variables defining this ordered sweep, or an empty tuple."""
        return self._sweep_variables

    def _sweep_summary(self) -> dict[str, Any]:
        variables = self._sweep_variables
        if not variables:
            return {}
        axes = {
            name: tuple(sorted({float(run.var[name].value) for run in self}))
            for name in variables
        }
        points = [tuple(float(run.var[name].value) for name in variables) for run in self]
        unique_points = set(points)
        expected = 1
        for values in axes.values():
            expected *= len(values)
        return {
            "variables": variables,
            "axes": axes,
            "point_count": len(points),
            "unique_point_count": len(unique_points),
            "duplicate_runs": len(points) - len(unique_points),
            "missing_points": max(0, expected - len(unique_points)),
            "complete": len(unique_points) == expected,
        }

    @analysis_operation("""Select runs whose full run_id matches a shell-style glob.

Usage:
    subset = runs.match("f5a_DBI*_I3")
    subset = runs.match("run_00?")
    subset = runs.match("run_[0-3]*")
    subset = runs.match("run_[!0]*")

Glob syntax:
    *       any sequence of characters
    ?       exactly one character
    [abc]   one character from a set
    [a-z]   one character from a range
    [!abc]  one character not in a set

Matching is case-sensitive, applies only to the complete run_id, preserves
the collection order, and returns a new Runs collection. No match returns an
empty Runs collection. Heavy run data are not loaded.
""")
    def match(self, pattern: str) -> "Runs":
        """Return runs whose complete ``run_id`` matches ``pattern``."""
        if not isinstance(pattern, str):
            raise TypeError("pattern must be a string")
        selected = {str(run.path): run for run in self if fnmatchcase(run.run_id, pattern)}
        return Runs(
            self.root, selected,
            _sweep_variables=self._sweep_variables,
            _order_paths=tuple(selected),
            _selection_note=f"match({pattern!r})",
        )


    @analysis_operation("""Pack runs into an ordered dictionary with convenient keys and optional user fields.

Usage:
    packed = runs.pack()
    packed = runs.pack(key="f5a_*_I3")
    packed = runs.pack(
        key="f5a_*_I3",
        label="f5a_*",
        color=["tab:blue", "tab:orange", "tab:green"],
        offset=[0.00, 0.05, 0.10],
        linewidth=2,
    )

Rules:
    - key=None uses the complete run_id.
    - key="prefix*suffix" extracts the dictionary key from the one ``*``.
    - A string field containing one ``*`` extracts its value from run_id.
    - A list or tuple assigns values positionally in Runs order.
    - Any other value is copied to every record.
    - Every record always contains ``{"run": Run}``; ``run`` is reserved.
    - Patterns must match every run, keys must be unique, and positional
      fields must have exactly the same length as the collection.

The returned object is a normal insertion-ordered Python dict.
""")
    def pack(self, *, key: str | None = None, **fields: Any) -> dict[str, dict[str, Any]]:
        """Return an ordered ``dict`` of run records and optional user data."""
        if key is not None and not isinstance(key, str):
            raise TypeError("key must be a string pattern or None")
        if "run" in fields:
            raise ValueError("field name 'run' is reserved by pack()")

        ordered = self._ordered()
        if key is None:
            keys = [run.run_id for run in ordered]
        else:
            keys = [_extract_star(run.run_id, key, context="key") for run in ordered]

        seen: dict[str, str] = {}
        for packed_key, run in zip(keys, ordered):
            if packed_key in seen:
                raise SelectionError(
                    f"pack key {packed_key!r} is not unique for run_ids "
                    f"{seen[packed_key]!r} and {run.run_id!r}"
                )
            seen[packed_key] = run.run_id

        resolved: dict[str, list[Any]] = {}
        for name, value in fields.items():
            if isinstance(value, str) and "*" in value:
                resolved[name] = [
                    _extract_star(run.run_id, value, context=f"field {name!r}")
                    for run in ordered
                ]
            elif isinstance(value, (list, tuple)):
                if len(value) != len(ordered):
                    raise ValueError(
                        f"field {name!r} has {len(value)} values, but Runs contains "
                        f"{len(ordered)} runs"
                    )
                resolved[name] = list(value)
            else:
                resolved[name] = [value] * len(ordered)

        packed: dict[str, dict[str, Any]] = {}
        for index, (packed_key, run) in enumerate(zip(keys, ordered)):
            record: dict[str, Any] = {"run": run}
            for name, values in resolved.items():
                record[name] = values[index]
            packed[packed_key] = record
        return packed

    def filter(
        self, *, engine: str | None = None, schema: str | None = None, version: str | None = None,
        model_class: str | None = None, has_output: str | None = None, path: str | None = None,
        prefix: str | None = None, var_name: str | None = None, var_value: Any = None,
        run_id: str | None = None, status: str | set[str] | tuple[str, ...] | list[str] | None = None,
    ) -> "Runs":
        """Return a filtered Runs collection without modifying this one."""
        path_token = str(path).replace("\\", "/").casefold() if path is not None else None
        prefix_token = str(prefix).casefold() if prefix is not None else None
        selected: dict[str, Run] = {}
        statuses = None if status is None else ({status} if isinstance(status, str) else set(status))
        for key, run in self._runs.items():
            if engine is not None and run.engine != engine:
                continue
            if schema is not None and run.schema != schema:
                continue
            if version is not None and run.version != version:
                continue
            if model_class is not None and getattr(run, "model_class", "") != model_class:
                continue
            if has_output is not None and has_output not in run.available_outputs():
                continue
            if path_token is not None and path_token not in getattr(run, "relative_dir", ".").casefold():
                continue
            if prefix_token is not None and prefix_token not in run.prefix.casefold():
                continue
            if var_name is not None and var_name not in run.var:
                continue
            if var_value is not None:
                if var_name is None:
                    raise ValueError("var_value requires var_name")
                if not _value_matches(run.var[var_name].value, var_value):
                    continue
            if run_id is not None and run.run_id != run_id:
                continue
            if statuses is not None and run.status not in statuses:
                continue
            selected[key] = run
        return Runs(self.root, selected, _sweep_variables=self._sweep_variables, _order_paths=tuple(selected), _selection_note=self._selection_note)

    @property
    def completed(self) -> "Runs":
        return self.filter(status="completed")

    @property
    def failed(self) -> "Runs":
        return self.filter(status="failed")

    @property
    def interrupted(self) -> "Runs":
        return self.filter(status="interrupted")

    def one(self, **filters: Any) -> Run:
        matches = self.filter(**filters) if filters else self
        if len(matches) == 1:
            return next(iter(matches))
        candidates = "\n  ".join(str(r.path) for r in matches) or "none"
        raise SelectionError(f"expected exactly one run for {filters!r}; found {len(matches)}:\n  {candidates}")

    def first(self) -> Run:
        if not self._runs:
            raise IndexError("empty Runs collection")
        return self._ordered()[0]

    def as_table(self) -> Table:
        rows = []
        for run in self:
            rows.append([
                str(run.path), run.run_id, run.engine, run.version, run.schema,
                getattr(run, "model_class", ""),
                ",".join(run.var.keys()),
                ",".join(sorted(run.available_outputs())),
            ])
        return Table(
            ["path", "run_id", "engine", "version", "schema", "model_class",
             "variables", "outputs"],
            rows, name="runs", source=str(self.root),
        )

    def model_diff(self, *, include_same: bool = False) -> Table:
        """Compare .model files across this collection. Numeric literals
        and quoted strings are normalized so e.g. ``rate kp_ab const 120``
        and ``rate kp_ab const 240`` land on the same comparison row."""
        return _make_model_diff(self, include_same=include_same)

    def info_text(self, max_rows: int = 10) -> str:
        if self._run_id_index_cache is not None:
            grouped = self._run_id_index_cache
        else:
            temporary: dict[str, list[Run]] = {}
            for run in self:
                temporary.setdefault(run.run_id, []).append(run)
            grouped = {run_id: tuple(items) for run_id, items in temporary.items()}
        api_names = set(super().__dir__())
        interactive = sorted(
            run_id for run_id, matches in grouped.items()
            if len(matches) == 1 and _is_public_identifier(run_id) and run_id not in api_names
        )
        collisions = sorted(
            run_id for run_id, matches in grouped.items()
            if len(matches) == 1 and run_id in api_names
        )
        duplicates = sorted(run_id for run_id, matches in grouped.items() if len(matches) > 1)

        lines = [
            "Runs",
            f"  root: {self.root.name}/" if self.root.name else "  root: ./",
            f"  count: {len(self)}",
            f"  unique run_id: {len(grouped)}",
            f"  interactive attributes: {len(interactive)}",
            f"  run-id index: {'cached' if self._run_id_index_cache is not None else 'lazy'}",
            "  engines: " + (", ".join(self.engines) if self.engines else "none"),
            "  declared variables: " + (", ".join(self.var.keys()) if self.var.keys() else "none"),
        ]
        if self._selection_note:
            lines.append(f"  selection: {self._selection_note}")
        if collisions:
            lines.append("  API-name collisions: " + ", ".join(collisions))
        if duplicates:
            lines.append("  duplicate run_id: " + ", ".join(duplicates))

        summary = self._sweep_summary()
        if summary:
            lines += [
                "  sweep dimensions: " + ", ".join(summary["variables"]),
                "  sweep grid: " + ("complete" if summary["complete"] else "incomplete"),
                f"  duplicate parameter runs: {summary['duplicate_runs']}",
                f"  missing parameter points: {summary['missing_points']}",
            ]
            for name, values in summary["axes"].items():
                lines.append(f"  {name} values: " + ", ".join(f"{value:g}" for value in values))
        lines.append("")
        if self:
            lines.append("  path                                      status       X_final       Mw          D")
        for run in self._ordered()[:max_rows]:
            rel = run.relative_dir.rstrip("/") + "/"
            status = str(run.status or "unknown")
            try:
                x_text = f"{float(run.conv.total[-1]):.4g}"
            except Exception:
                x_text = "-"
            try:
                mw_text = f"{float(run.mw[-1]):.6g}"
            except Exception:
                mw_text = "-"
            try:
                d_text = f"{float(run.dispersity[-1]):.4g}"
            except Exception:
                d_text = "-"
            lines.append(f"  {rel:<40} {status:<12} {x_text:>9} {mw_text:>11} {d_text:>8}")
        if len(self) > max_rows:
            lines.append(f"  ... {len(self) - max_rows} more")
        lines.extend([
            "",
            "Common next steps:",
            "  run = runs.<run_id>",
            '  subset = runs.match("f5a_DBI*_I3")',
            '  packed = subset.pack(key="f5a_*_I3")',
            '  sweep = runs.sweep("IA", "temperature")',
            "  runs.match.help(); runs.pack.help(); runs.sweep.help()",
            "  runs.model_diff()",
        ])
        return "\n".join(lines)

    def info(self, max_rows: int = 10) -> str:
        text = self.info_text(max_rows=max_rows)
        print(text)
        return text

    def help(self) -> str:
        text = (
            "Runs\n"
            "\nInteractive access:\n"
            "  run = runs.<run_id>\n"
            '  run = runs.one(run_id="f5a_DBI30_I3")\n'
            '  runs.run_id["f5a_DBI30_I3"]   # explicit subcollection\n'
            "\nSelection:\n"
            '  subset = runs.match("f5a_DBI*_I3")\n'
            '  subset = runs.filter(var_name="IA", var_value=0.1)\n'
            "  runs.match.help()\n"
            "\nPacking:\n"
            '  packed = subset.pack(key="f5a_*_I3", color=[...])\n'
            "  runs.pack.help()\n"
            "\nParameter studies:\n"
            '  sweep = runs.sweep("IA", "temperature")\n'
            "  runs.sweep.help()\n"
            "\nInspection:\n"
            "  runs.info(); runs.as_table(); runs.model_diff()\n"
            "\nNotes:\n"
            "  Direct access and completion are lazy. Public Runs API names win\n"
            "  on collisions; use runs.one(run_id=...) for such run IDs.\n"
        )
        print(text)
        return text

    def __repr__(self) -> str:
        return f"Runs(n={len(self)}, root={str(self.root)!r})"
