"""pyslimmc -- read-only NumPy/Matplotlib analysis for Slimmc Storage.

Use ``pyslimmc.open(path)`` for one run and ``pyslimmc.scan(path)`` for a
collection. The public interface is documented in ``docs/reference/PYSLIMMC_API.md``.
"""
from __future__ import annotations

from pathlib import Path

from .core import (
    PyslimmcError, FeatureUnavailableError, ChemicalAnalysisNotApplicableError, AnalysisNotApplicableError,
    ChemicalModelIncompatibleError, DataUnavailableError, IncompleteSequenceDataError,
    InvalidOutputError, ValidationFailedError, NumericalAnalysisError,
    SnapshotUnavailableError, FinalSnapshotUnavailableError,
    MassModelUnavailableError, InvalidDistributionConfigurationError,
)
from .run import (
    Run, Variable, Variables, DataConsistencyError, MassAuditResult, UnknownColumnError,
    UnknownMonomerError, UnsupportedChainSchema,
)
from .runs import Runs, SelectionError, scan as _scan
from .plotting import PlotStyle, available_styles, get_style, figure_size
from .report import Report, report
from ._version import __version__


def open(path: str | Path, *, allow_incomplete: bool = False) -> Run:
    """Open one Slimmc Storage run directory."""
    path_obj = Path(path)
    if not (path_obj / "run_metadata.json").is_file() or not (path_obj / "schema.jsonl").is_file():
        raise InvalidOutputError(
            "Not a Slimmc Storage run: run_metadata.json or schema.jsonl is missing"
        )
    from ._storage import open_storage
    return open_storage(path_obj, allow_incomplete=allow_incomplete)


def scan(path: str | Path = ".", *, recursive: bool = True, skip_bad: bool = False) -> Runs:
    """Scan a directory tree for Slimmc Storage runs."""
    return _scan(path, recursive=recursive, skip_bad=skip_bad)


def help() -> str:
    """Print and return the shortest current pyslimmc workflow."""
    text = (
        "pyslimmc " + __version__ + ":\n"
        "  run = pyslimmc.open(\"results/run_000001\")\n"
        "  run.info(); snap = run.final\n"
        "  run.t; run.conc[\"A\"]; run.conv[\"A\"]; run.conv.total\n"
        "  run.mn; run.mw; run.dispersity; run.temp; run.k\n"
        "  snap.chains; snap.mwd(); snap.cld(); snap.chain_mass_spectrum()\n"
        "  run.channels; run.firings; run.actions; run.kinetics\n"
        "  run.copolymerization; run.microstructure\n"
        "  run.validate(strict=True); run.mass_audit()\n"
        "  run.raw; run.diagnostics; run.summary()\n"
        "  runs = pyslimmc.scan(\"results\")"
    )
    print(text)
    return text


__all__ = [
    "__version__", "open", "scan", "help",
    "Run", "Variable", "Variables", "Runs", "SelectionError",
    "Report", "report",
    "PlotStyle", "available_styles", "get_style", "figure_size",
    "MassAuditResult",
    "PyslimmcError", "FeatureUnavailableError",
    "ChemicalAnalysisNotApplicableError", "AnalysisNotApplicableError", "ChemicalModelIncompatibleError",
    "DataUnavailableError", "IncompleteSequenceDataError", "InvalidOutputError",
    "ValidationFailedError", "NumericalAnalysisError", "DataConsistencyError",
    "UnknownColumnError", "UnknownMonomerError", "UnsupportedChainSchema",
    "SnapshotUnavailableError", "FinalSnapshotUnavailableError",
    "MassModelUnavailableError", "InvalidDistributionConfigurationError",
]

# Install the uniform public help/info contract after the public modules are loaded.
from .helpinfo import install as _install_helpinfo
_install_helpinfo()
del _install_helpinfo

from .options import options
