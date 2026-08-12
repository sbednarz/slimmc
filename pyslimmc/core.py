from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


class PyslimmcError(Exception):
    """Base class for domain-level pyslimmc failures."""


class FeatureUnavailableError(PyslimmcError):
    pass


class FeatureVariantNotImplementedError(FeatureUnavailableError):
    pass


class ChemicalDomainError(PyslimmcError):
    pass


class UndefinedChemicalQuantityError(ChemicalDomainError):
    pass


class ChemicalAnalysisNotApplicableError(ChemicalDomainError):
    pass


class AnalysisNotApplicableError(ChemicalAnalysisNotApplicableError):
    """Requested physical analysis is invalid or undefined for this run."""


class ChemicalModelIncompatibleError(ChemicalDomainError):
    pass


class DataUnavailableError(PyslimmcError):
    pass


class IncompleteSequenceDataError(DataUnavailableError):
    pass


class InvalidOutputError(PyslimmcError):
    pass


class ValidationFailedError(PyslimmcError):
    pass


class NumericalAnalysisError(PyslimmcError):
    pass


@dataclass(frozen=True)
class OutputStatus:
    available: tuple[str, ...]
    missing: tuple[str, ...]
    invalid: tuple[str, ...]
    run_completed: bool
    has_final_snapshot: bool

    @property
    def complete(self) -> bool:
        return (
            self.run_completed
            and self.has_final_snapshot
            and not self.missing
            and not self.invalid
        )

    @property
    def partial(self) -> bool:
        return not self.complete

    def info_text(self) -> str:
        state = "complete" if self.complete else "partial"
        lines = [f"output status: {state}", f"available: {', '.join(self.available) or 'none'}"]
        if self.missing:
            lines.append("missing: " + ", ".join(self.missing))
        if self.invalid:
            lines.append("invalid: " + ", ".join(self.invalid))
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text


@dataclass(frozen=True)
class ValidationReport:
    is_valid: bool
    is_complete: bool
    missing_outputs: tuple[str, ...] = ()
    invalid_outputs: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    failed_checks: tuple[str, ...] = ()
    first_failure: dict | None = None
    failure_counts: dict[str, int] | None = None
    details: tuple[dict, ...] = ()

    def raise_if_failed(self) -> None:
        if not self.is_valid:
            details = list(self.invalid_outputs) + list(self.warnings)
            raise ValidationFailedError("; ".join(details) or "run validation failed")

    def info_text(self) -> str:
        lines = [
            f"validation: {'valid' if self.is_valid else 'invalid'}",
            f"complete: {self.is_complete}",
        ]
        if self.missing_outputs:
            lines.append("missing outputs: " + ", ".join(self.missing_outputs))
        if self.invalid_outputs:
            lines.append("invalid outputs: " + ", ".join(self.invalid_outputs))
        if self.failed_checks:
            lines.append("failed checks: " + ", ".join(self.failed_checks))
        if self.first_failure:
            f = self.first_failure
            lines.extend([
                "first failure:",
                f"  subcheck: {f.get('subcheck', '')}",
                f"  physical line: {f.get('physical_line', '')}",
                f"  snapshot_id: {f.get('snapshot_id', '')}",
                f"  chain_record_id: {f.get('chain_record_id', '')}",
                f"  field: {f.get('field', '')}",
                f"  actual: {f.get('actual', '')}",
                f"  expected: {f.get('expected', '')}",
                f"  difference: {f.get('difference', '')}",
                f"  tolerance: {f.get('tolerance', '')}",
            ])
        if self.failure_counts:
            lines.append("failure counts: " + ", ".join(
                f"{key}={value}" for key, value in sorted(self.failure_counts.items())
            ))
        if self.warnings:
            lines.append("warnings: " + "; ".join(self.warnings))
        return "\n".join(lines)

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text


@dataclass(frozen=True)
class Diagnostics:
    output_status: OutputStatus
    validation: ValidationReport

    @property
    def warnings(self) -> tuple[str, ...]:
        return self.validation.warnings

    def info_text(self) -> str:
        return self.output_status.info_text() + "\n" + self.validation.info_text()

    def info(self) -> str:
        text = self.info_text()
        print(text)
        return text


def require_file(path: Path, *, output: str) -> Path:
    if not path.is_file():
        raise DataUnavailableError(f"required output {output!r} is unavailable: {path}")
    return path


class SnapshotUnavailableError(DataUnavailableError):
    """Requested snapshot does not exist."""


class FinalSnapshotUnavailableError(SnapshotUnavailableError):
    """A completed final snapshot is unavailable."""


class MassModelUnavailableError(DataUnavailableError):
    """The requested chain-mass model cannot be evaluated."""


class InvalidDistributionConfigurationError(InvalidOutputError):
    """A CLD/MWD configuration contains incompatible options."""
