from __future__ import annotations

import copy
import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Mapping

from .core import FeatureUnavailableError, InvalidOutputError


class UnsupportedFeatureError(AttributeError, FeatureUnavailableError):
    """Raised when a caller asks a Run for an attribute that exists on the
    other engine's subclass but not on this one.

    This is deliberately an ``AttributeError`` subclass (not a plain
    ``Exception``) -- see docs/PYSLIMMC.md's note on
    interactive convenience: ``hasattr()`` and IPython/Jupyter's tab
    completion rely on catching plain ``AttributeError`` to test whether an
    attribute exists. Raising anything else here would make completion
    machinery misbehave instead of just quietly not offering the
    unavailable name.
    """


class DataConsistencyError(InvalidOutputError):
    """Raised when a run's own output files disagree with each other in a
    way that suggests corrupted/inconsistent data -- e.g. copo's
    the Storage channel-events table "total_fires" column not matching the sum of its own
    "fires_*" columns. Deliberately loud rather than silently trusting
    one source over the other."""


class UnknownColumnError(KeyError):
    """Raised when a selector/argument names a column that does not exist
    in the underlying table -- e.g. ``Oligomers.top(by="cout")`` (typo for
    "count"). Carries ``.available`` so callers can show/inspect the valid
    options instead of getting a silent empty result or a bare KeyError.
    See docs/PYSLIMMC.md."""

    def __init__(self, name: str, context: str, available: list[str]):
        available_text = ", ".join(available) if available else "none"
        msg = f"Unknown column {name!r} for {context}. Available: {available_text}"
        super().__init__(msg)
        self.name = name
        self.context = context
        self.available = list(available)


class UnsupportedChainSchema(DataConsistencyError):
    """Raised when a chain-population table (the Storage chains table) is missing a
    column that the requested operation structurally depends on -- e.g.
    no ``dp`` column at all. Distinct from a plain missing/typo'd
    *filter* value: this means the file itself doesn't have the shape the
    operation needs. Shared across engines; see docs/PYSLIMMC.md,
    points 7 and 16."""


class UnknownMonomerError(KeyError):
    """Raised when a composition selector (e.g. ``by_composition(AAa=1)``)
    names a monomer that is not part of this run's declared monomer set.
    See docs/PYSLIMMC.md."""

    def __init__(self, name: str, available: list[str]):
        available_text = ", ".join(available) if available else "none"
        msg = f"Unknown monomer {name!r}. Available composition keys: {available_text}"
        super().__init__(msg)
        self.name = name
        self.available = list(available)


class _DictView:
    """Read-only dotted-attribute view over a nested dict section of
    ``run_metadata.json`` (e.g. ``run.metadata.var``), with ``.raw()`` to get
    the plain dict back when needed."""

    __slots__ = ("_value",)

    def __init__(self, value: dict[str, Any]):
        object.__setattr__(self, "_value", value)

    def __getitem__(self, key: str) -> Any:
        value = self._value[key]
        return _DictView(value) if isinstance(value, dict) else value

    def __getattr__(self, key: str):
        try:
            return self[key]
        except KeyError as exc:
            raise AttributeError(key) from exc

    def __dir__(self) -> list[str]:
        standard = set(super().__dir__())
        dynamic = {k for k in self._value if isinstance(k, str) and k.isidentifier() and k not in standard}
        return sorted(standard | dynamic)

    def __repr__(self) -> str:
        return f"_DictView({self._value!r})"

    def raw(self) -> dict[str, Any]:
        # Returning the internal dict would make this nominally read-only
        # view mutable through an alias.  A deep copy also protects nested
        # model/execution sections.
        return copy.deepcopy(self._value)




@dataclass(frozen=True)
class ReproducibilityReport:
    input_status: str
    storage_manifest_status: str
    storage_hash_status: str
    binary_status: str
    overall: str
    messages: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return self.overall == "VERIFIED"

    def info_text(self) -> str:
        lines = [
            "Reproducibility verification",
            f"  input file:       {self.input_status}",
            f"  storage manifest: {self.storage_manifest_status}",
            f"  storage hash:     {self.storage_hash_status}",
            f"  binary:           {self.binary_status}",
            f"  overall:          {self.overall}",
        ]
        lines.extend(f"  note: {m}" for m in self.messages)
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text


class Reproducibility:
    __slots__ = ("_run",)

    def __init__(self, run):
        self._run = run

    @staticmethod
    def _sha256(path: Path) -> str:
        h = hashlib.sha256()
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    @property
    def input_hash(self): return getattr(self._run.input, "hash", None)
    @property
    def model_hash(self): return getattr(self._run.model, "hash", None)
    @property
    def binary_hash(self): return getattr(self._run.execution, "binary_hash", None)
    @property
    def storage_hash(self): return getattr(self._run.storage, "hash", None)
    @property
    def git_commit(self): return getattr(self._run.execution, "git_commit", None)
    @property
    def git_dirty(self): return getattr(self._run.execution, "git_dirty", None)

    def verify(self, binary: str | Path | None = None) -> ReproducibilityReport:
        messages = []
        input_path = self._run.path / str(getattr(self._run.input, "file", "input.model"))
        expected_input = self.input_hash
        if not expected_input:
            input_status = "NOT RECORDED"
        elif not input_path.is_file():
            input_status = "MISSING"
        else:
            input_status = "PASS" if self._sha256(input_path) == expected_input else "FAIL"

        manifest_name = getattr(self._run.storage, "manifest_file", "checksums.sha256")
        manifest_path = self._run.path / str(manifest_name)
        manifest_status = "NOT RECORDED"
        storage_hash_status = "NOT RECORDED"
        if manifest_path.is_file():
            manifest_status = "PASS"
            for line_no, line in enumerate(manifest_path.read_text(encoding="utf-8").splitlines(), 1):
                if not line.strip():
                    continue
                try:
                    expected, rel = line.split("  ", 1)
                except ValueError:
                    manifest_status = "FAIL"
                    if len(messages) < 100:
                        messages.append(f"invalid manifest line {line_no}")
                    continue
                target = self._run.path / Path(rel)
                if not target.is_file() or self._sha256(target) != expected:
                    manifest_status = "FAIL"
                    if len(messages) < 100:
                        messages.append(f"mismatch: {rel}")
            expected_storage = self.storage_hash
            if expected_storage:
                payload = b"slimmc-storage-hash-v1\n" + manifest_path.read_bytes()
                storage_hash_status = "PASS" if hashlib.sha256(payload).hexdigest() == expected_storage else "FAIL"
        elif self.storage_hash:
            manifest_status = "MISSING"

        expected_binary = self.binary_hash
        if not expected_binary:
            binary_status = "NOT RECORDED"
        elif binary is None:
            binary_status = "NOT CHECKED"
        else:
            bp = Path(binary)
            binary_status = "PASS" if bp.is_file() and self._sha256(bp) == expected_binary else "FAIL"

        required = (input_status, manifest_status, storage_hash_status)
        if any(x in {"FAIL", "MISSING"} for x in required) or binary_status == "FAIL": overall = "FAILED"
        elif all(x == "PASS" for x in required) and binary_status in {"PASS", "NOT CHECKED", "NOT RECORDED"}: overall = "VERIFIED"
        else: overall = "PARTIAL"
        return ReproducibilityReport(input_status, manifest_status, storage_hash_status, binary_status, overall, tuple(messages))

    def compare(self, other) -> dict[str, str]:
        rhs = other.reproducibility if hasattr(other, "reproducibility") else other
        def state(a, b):
            if a is None or b is None: return "NOT AVAILABLE"
            return "IDENTICAL" if a == b else "DIFFERENT"
        return {
            "input": state(self.input_hash, rhs.input_hash),
            "model": state(self.model_hash, rhs.model_hash),
            "binary": state(self.binary_hash, rhs.binary_hash),
            "storage": state(self.storage_hash, rhs.storage_hash),
        }

    def info(self) -> str:
        text = "\n".join([
            "Reproducibility",
            f"  input hash:   {self.input_hash}",
            f"  model hash:   {self.model_hash}",
            f"  binary hash:  {self.binary_hash}",
            f"  storage hash: {self.storage_hash}",
            f"  git commit:   {self.git_commit}",
            f"  git dirty:    {self.git_dirty}",
        ])
        print(text); return text


@dataclass(frozen=True)
class Variable:
    """One resolved ``var KIND NAME UNIT`` declaration from run metadata."""

    kind: str
    name: str
    value: float
    unit: str


class Variables(Mapping[str, Variable]):
    """Read-only collection of variables declared by one model run."""

    __slots__ = ("_items",)

    def __init__(self, records: list[dict[str, Any]] | tuple[dict[str, Any], ...] | None = None):
        items: dict[str, Variable] = {}
        for index, record in enumerate(records or ()):
            if not isinstance(record, dict):
                raise DataConsistencyError(f"variables[{index}] must be an object")
            try:
                kind = str(record["kind"])
                name = str(record["name"])
                value = float(record["value"])
                unit = str(record["unit"])
            except KeyError as exc:
                raise DataConsistencyError(f"variables[{index}] is missing {exc.args[0]!r}") from exc
            except (TypeError, ValueError) as exc:
                raise DataConsistencyError(f"variables[{index}].value must be numeric") from exc
            if name in items:
                raise DataConsistencyError(f"duplicate variable name {name!r} in run metadata")
            items[name] = Variable(kind=kind, name=name, value=value, unit=unit)
        self._items = items

    def __getitem__(self, name: str) -> Variable:
        return self._items[name]

    def __iter__(self) -> Iterator[str]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __getattr__(self, name: str) -> Variable:
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return self._items[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __dir__(self) -> list[str]:
        standard = set(super().__dir__())
        dynamic = {name for name in self._items if name.isidentifier() and name not in standard}
        return sorted(standard | dynamic)

    def info_text(self) -> str:
        lines = ["Variables", f"  count: {len(self)}"]
        if self._items:
            lines += ["", "  name                 kind       value          unit"]
            for variable in self._items.values():
                lines.append(f"  {variable.name:<20} {variable.kind:<10} {variable.value:<14.8g} {variable.unit}")
        else:
            lines.append("  none declared")
        lines += ["", "Common next steps:", '  run.var["name"].value', "  run.var.keys()"]
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text

    def help(self) -> str:
        text = """Variables

Access:
  variable = run.var["name"]
  variable.kind / variable.name / variable.value / variable.unit
  run.var.keys(); run.var.info()

Values are read from run_metadata.json["variables"]."""
        print(text)
        return text

    def __repr__(self) -> str:
        return f"Variables({list(self._items)!r})"


@dataclass(frozen=True)
class MassEntry:
    """One declared-or-not mass value (a monomer or an end group).

    Distinguishes "no data" (``mw=None``) from "genuinely declared as
    zero" (``mw=0.0, declared=True``) -- collapsing these with something
    like ``float(x or 0.0)`` silently destroys the difference between an
    intentional zero and a missing value. See docs/PYSLIMMC.md,
    point 11."""

    name: str
    mw: float | None
    declared: bool
    source: str = "missing"  # one of: "model", "default", "inferred", "missing"


@dataclass
class MassAuditResult:
    """Shared return type for both engines' ``mass_audit()``.

    Replaces homo's old standalone ``MassAudit`` dataclass and copo's raw
    ``Table`` -- the same method name on both engines previously returned
    genuinely different types, which defeats the point of sharing the
    name. See docs/PYSLIMMC.md.

    ``details`` carries whatever richer, engine-specific data was
    available (e.g. copo's per-chain expected-vs-actual mass table) for
    callers that want the full picture; the fields above are the common,
    guaranteed-present summary.
    """

    ok: bool
    mass_model: str
    checked_records: int = 0
    checked_chains: int = 0
    missing_monomers: tuple[str, ...] = ()
    missing_endgroups: tuple[str, ...] = ()
    implicit_zero_monomers: tuple[str, ...] = ()
    implicit_zero_endgroups: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    details: Any = None
    entries: tuple[MassEntry, ...] = ()

    def __bool__(self) -> bool:
        return self.ok

    def info(self) -> str:
        lines = [
            "MassAuditResult:",
            f"  ok: {self.ok}",
            f"  mass_model: {self.mass_model}",
            f"  checked records: {self.checked_records}",
            f"  checked chains: {self.checked_chains}",
        ]
        if self.missing_monomers:
            lines.append(f"  missing monomer masses: {', '.join(self.missing_monomers)}")
        if self.missing_endgroups:
            lines.append(f"  missing endgroup masses: {', '.join(self.missing_endgroups)}")
        if self.implicit_zero_monomers:
            lines.append(f"  monomers treated as zero (permissive): {', '.join(self.implicit_zero_monomers)}")
        if self.implicit_zero_endgroups:
            lines.append(f"  endgroups treated as zero (permissive): {', '.join(self.implicit_zero_endgroups)}")
        for w in self.warnings:
            lines.append(f"  warning: {w}")
        for e in self.entries:
            status = f"{e.mw:g}" if e.declared else "undeclared"
            lines.append(f"  {e.name}: mw={status} (source={e.source})")
        return "\n".join(lines)

    def as_table(self):
        from .table import Table
        if self.details is not None and hasattr(self.details, "columns"):
            return self.details
        return Table(
            ["field", "value"],
            [
                ["ok", self.ok],
                ["mass_model", self.mass_model],
                ["checked_records", self.checked_records],
                ["checked_chains", self.checked_chains],
                ["missing_monomers", ", ".join(self.missing_monomers)],
                ["missing_endgroups", ", ".join(self.missing_endgroups)],
                ["implicit_zero_monomers", ", ".join(self.implicit_zero_monomers)],
                ["implicit_zero_endgroups", ", ".join(self.implicit_zero_endgroups)],
            ],
            name="mass_audit_summary",
        )

    def raise_if_failed(self) -> None:
        """Raise ``DataConsistencyError`` if this audit did not pass.
        Called automatically by ``mwd(strict=True)`` (the default on both
        engines, see docs/PYSLIMMC.md) -- also usable
        directly for a one-line "check and raise" on a standalone
        ``mass_audit()`` result."""
        if self.ok:
            return
        reason = self.warnings[0] if self.warnings else (
            f"missing monomer mass(es): {list(self.missing_monomers)}; "
            f"missing endgroup mass(es): {list(self.missing_endgroups)}"
        )
        raise DataConsistencyError(f"mass_audit failed: {reason}")

    def entry(self, name: str) -> MassEntry | None:
        """Look up a single monomer/endgroup's MassEntry by name, or
        None if this audit didn't track one under that name (e.g. an
        engine/mass_model combination that doesn't populate .entries
        yet). Prefer this over scanning .entries by hand."""
        for e in self.entries:
            if e.name == name:
                return e
        return None

    def __repr__(self) -> str:
        return f"MassAuditResult(ok={self.ok!r}, mass_model={self.mass_model!r})"


@dataclass
class Run:
    """Native public base for one Slimmc Storage run."""

    path: Path
    _metadata: dict[str, Any] = field(default_factory=dict, repr=False)
    run_id: str = ""
    relative_dir: str = "."
    _prefix: str = field(default="", repr=False)

    def __post_init__(self) -> None:
        if not self.run_id:
            self.run_id = str(self._metadata.get("run_id") or self.path.name)

    @property
    def prefix(self) -> str:
        return self._prefix or self.run_id

    @property
    def engine(self) -> str:
        return str(self._metadata.get("engine", ""))

    @property
    def engine_family(self) -> str:
        return str(self._metadata.get("engine_family", ""))

    @property
    def version(self) -> str:
        return str(self._metadata.get("engine_version", self._metadata.get("version", "")))

    @property
    def schema(self) -> str:
        return str(self._metadata.get("storage", "slimmc-storage"))

    @property
    def status(self) -> str:
        return str(self._metadata.get("run_status", self._metadata.get("status", "")))

    @property
    def termination_reason(self) -> str:
        return str(self._metadata.get("termination_reason", ""))

    @property
    def metadata(self) -> _DictView:
        return _DictView(self._metadata)

    @property
    def model(self) -> _DictView:
        return _DictView(self._metadata.get("model", {}))

    @property
    def input(self) -> _DictView:
        section = self._metadata.get("input")
        if not isinstance(section, dict):
            section = {
                "file": self._metadata.get("input_model_file", "input.model"),
                "source_name": self._metadata.get("source_model_name", ""),
                "hash": self._metadata.get("input_model_sha256"),
                "hash_algorithm": "sha256",
            }
        return _DictView(section)

    @property
    def execution(self) -> _DictView:
        section = self._metadata.get("execution")
        if not isinstance(section, dict):
            keys = (
                "engine", "engine_version", "cli_version", "git_commit", "git_dirty", "git_tag",
                "binary_name", "binary_hash", "binary_hash_algorithm", "build_timestamp_utc",
                "started_at_utc", "finished_at_utc", "wall_time_s", "run_status",
                "exit_code", "platform", "threads", "compiler", "build_mode",
            )
            section = {key: self._metadata[key] for key in keys if key in self._metadata}
            if "run_status" in section:
                section["status"] = section["run_status"]
        return _DictView(section)

    @property
    def storage(self) -> _DictView:
        section = self._metadata.get("storage_info")
        if not isinstance(section, dict):
            section = {
                "name": self._metadata.get("storage", "slimmc-storage"),
                "format_version": self._metadata.get("storage_format_version"),
                "hash": self._metadata.get("storage_hash"),
                "hash_algorithm": self._metadata.get("storage_hash_algorithm", "sha256"),
                "hash_schema": self._metadata.get("storage_hash_schema"),
                "complete": (self.path / "RESULTS_COMPLETE").is_file(),
            }
        return _DictView(section)

    @property
    def reproducibility(self) -> Reproducibility:
        return Reproducibility(self)

    @property
    def var(self) -> Variables:
        raise NotImplementedError

    @property
    def monomers(self) -> dict[str, dict[str, Any]]:
        raise NotImplementedError

    @property
    def endgroups(self) -> dict[str, dict[str, Any]]:
        raise NotImplementedError

    def monomer_mw(self, name: str) -> float:
        try:
            value = self.monomers[name].get("molar_mass_increment", self.monomers[name].get("molar_mass"))
        except KeyError as exc:
            raise KeyError(f"unknown monomer {name!r}; declared: {sorted(self.monomers)}") from exc
        if value is None:
            raise KeyError(f"monomer {name!r} has no declared molar mass")
        return float(value)

    def endgroup_mw(self, name: str) -> float:
        try:
            value = self.endgroups[name].get("molar_mass_contribution", self.endgroups[name].get("molar_mass"))
        except KeyError as exc:
            raise KeyError(f"unknown endgroup {name!r}; declared: {sorted(self.endgroups)}") from exc
        if value is None:
            raise KeyError(f"endgroup {name!r} has no declared mass")
        return float(value)

    def help(self) -> str:
        text = """Run

Obtain a Run:
  run = sl.open(path)
  run = runs.<run_id>
  run = runs.one(run_id="...")
  run = runs.match("...")[0]
  run = runs.pack(...)[key]["run"]

Snapshots and axes:
  run.first / run.last / run.final / run.at_snapshot(id) / run.at_time(t) / run.at_event(n)
  run.t / run.event / run.sid

State and moments:
  run.count / run.moles / run.conc / run.conv
  run.mn / run.mw / run.mz / run.dispersity
  run.temp / run.k / run.channels / run.firings / run.actions

Chains and analyses:
  run.chains / run.chain_counts() / run.mass_audit()
  run.mwd.help(); run.cld.help(); run.chain_mass_spectrum.help()
  run.copolymerization / run.microstructure

Variables and diagnostics:
  run.var / run.var["name"].value / run.var.info() / run.var.help()
  run.raw / run.diagnostics / run.validate() / run.refresh()
  run.summary(); run.info()"""
        print(text)
        return text

    def __repr__(self) -> str:
        return (f"<{self.__class__.__name__} engine={self.engine!r} "
                f"version={self.version!r} run_id={self.run_id!r} path={str(self.path)!r}>")
