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


MWD_HELP = """Molar-mass density: mwd(...)

Start here:
    mwd = run.mwd()
    mwd.info()
    mwd.plot()

Defaults:
    snapshot='final', pool='all', mass_model=None

Meaning:
    mass-weighted density dW/dlog10(M), normalized to unit area.
    Reconstruction follows mcPolymer-style linear interpolation in log10(M).
    Homopolymer affine M(DP) lattices are zero-filled at missing integer DP;
    general/copolymer populations interpolate occupied exact-mass support.

Exact source representations:
    run.mass_counts()
    run.mass_distribution()

MWD is a derived density representation. It is not the exact discrete mass
population and it is not an SEC instrument-response curve.
"""

CLD_HELP = """Chain-length distribution: cld(...)

Start here:
    cld = run.cld()
    cld.info()
    cld.plot()

Defaults:
    snapshot='final', pool='all', weighting='number', mass_model=None

Weighting:
    number  discrete chain-number fraction by integer DP
    mass    discrete polymer-mass fraction grouped by integer DP
    z       discrete z-weighted fraction by integer DP

Exact chain counts by DP:
    run.dp_counts()

CLD is exact and discrete. A logarithmic plotting axis does not change it into
an independent logarithmic density.
"""

SEC_HELP = """SEC broadening: sec(...)

Start here:
    sec = run.sec(sigma_log10M=0.05)
    sec.info()
    sec.plot()

Required:
    sigma_log10M   standard deviation of the Gaussian instrument response
                   in log10(molar-mass) units

Optional:
    snapshot='final', pool='all', mass_model=None, step_log10M=None

SEC operates directly on exact chain masses and polymer-mass fractions.
It returns a continuous apparent density dW_app/dlog10(M).  It is an
instrument-response model, not generic smoothing of an MWD curve.
"""
