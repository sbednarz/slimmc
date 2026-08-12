from __future__ import annotations

from dataclasses import dataclass
from functools import update_wrapper
from inspect import signature
from typing import Any, Callable, Generic, TypeVar

T = TypeVar("T")
R = TypeVar("R")


@dataclass(frozen=True)
class _BoundOperation(Generic[R]):
    _owner: Any
    _func: Callable[..., R]
    _help_text: str

    def __call__(self, *args: Any, **kwargs: Any) -> R:
        return self._func(self._owner, *args, **kwargs)

    def help(self) -> str:
        print(self._help_text)
        return self._help_text

    @property
    def __doc__(self) -> str:
        return self._func.__doc__ or self._help_text

    @property
    def __name__(self) -> str:
        return self._func.__name__

    @property
    def __signature__(self):
        """Expose the wrapped public signature without the bound ``self``."""
        original = signature(self._func)
        parameters = tuple(original.parameters.values())
        if parameters and parameters[0].name == "self":
            parameters = parameters[1:]
        return original.replace(parameters=parameters)


class analysis_operation(Generic[R]):
    """Descriptor for an analysis that is both callable and has ``.help()``."""

    def __init__(self, help_text: str):
        self.help_text = help_text
        self.func: Callable[..., R] | None = None

    def __call__(self, func: Callable[..., R]) -> "analysis_operation[R]":
        self.func = func
        update_wrapper(self, func)
        return self

    def __get__(self, instance: Any, owner: type | None = None):
        if instance is None:
            return self
        if self.func is None:
            raise RuntimeError("analysis operation is not bound to a function")
        return _BoundOperation(instance, self.func, self.help_text)

    def help(self) -> str:
        print(self.help_text)
        return self.help_text


MWD_HELP = """Molar-mass distribution: mwd(...)

Start here:
    mwd = run.mwd()
    mwd.info()
    mwd.plot()

Defaults:
    snapshot='final', pool='all', mass_model='repeat_units'
    method='gaussian', basis='mass', coordinate='log10'
    output='density', normalization='per_series'

Methods:
    sticks    exact discrete masses
    hist      histogram without smoothing
    gaussian  Gaussian-smoothed histogram
    kde       kernel-density estimate

Axes:
    mwd.x        physical molar mass in g/mol
    mwd.log10_x  log10(mwd.x)
    mwd.y        distribution in the selected coordinate

SEC/GPC-like example:
    run.mwd(method='gaussian', basis='mass', coordinate='log10',
            output='density', normalization='per_series',
            bin_width=0.01, sigma=0.04)

Exact chain counts by mass:
    run.mwd(method='sticks', basis='number', output='amount',
            normalization='absolute')
"""

CLD_HELP = """Chain-length distribution: cld(...)

Start here:
    cld = run.cld()
    cld.info()
    cld.plot()

Defaults:
    snapshot='final', pool='all'
    method='sticks', basis='number', coordinate='linear'
    output='fraction', normalization='per_series'

Interpretation:
    cld.x        degree of polymerization (DP)
    cld.log10_x  log10(cld.x)
    basis='number' weights chains by count
    basis='mass' weights DP classes by chain mass
"""

SPECTRUM_HELP = """Neutral chain-mass spectrum: chain_mass_spectrum(...)

Start here:
    spectrum = run.chain_mass_spectrum()
    spectrum.info()
    spectrum.plot()

Arguments:
    snapshot='final', pool='all', mass_model='repeat_units'
    normalize='count' | 'fraction' | 'base_peak'

Fields:
    spectrum.mass
    spectrum.intensity
    spectrum.base_peak_mass
    spectrum.base_peak_intensity

This is a neutral-chain mass spectrum, not an m/z spectrum. Charges,
isotopes, adducts, fragmentation and detector response are not modelled.
"""
