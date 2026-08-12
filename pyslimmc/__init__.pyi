from pathlib import Path
from .run import (
    Run, Variable, Variables, MassAuditResult, DataConsistencyError,
    UnknownColumnError, UnknownMonomerError, UnsupportedChainSchema,
)
from .runs import Runs, SelectionError
from .report import Report
from .plotting import PlotStyle, available_styles, get_style, figure_size
from .options import Options, options
from .core import (
    PyslimmcError, FeatureUnavailableError, ChemicalAnalysisNotApplicableError, AnalysisNotApplicableError,
    ChemicalModelIncompatibleError, DataUnavailableError, IncompleteSequenceDataError,
    InvalidOutputError, ValidationFailedError, NumericalAnalysisError,
    SnapshotUnavailableError, FinalSnapshotUnavailableError,
    MassModelUnavailableError, InvalidDistributionConfigurationError,
)

__version__: str

def open(path: str | Path, *, allow_incomplete: bool = ...) -> Run: ...
def scan(path: str | Path = ..., *, recursive: bool = ..., skip_bad: bool = ...) -> Runs: ...
def help() -> str: ...
def report(title: str = ...) -> Report: ...
