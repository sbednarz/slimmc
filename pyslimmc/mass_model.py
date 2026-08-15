from __future__ import annotations

import numpy as np

from .core import InvalidOutputError

_ALLOWED_MASS_MODELS = {"repeat_units", "with_end_groups"}


def _validate_mass_model(model: str) -> str:
    if model not in _ALLOWED_MASS_MODELS:
        raise ValueError("mass_model must be 'repeat_units' or 'with_end_groups'")
    return model


def resolve_run_mass_model(run, mass_model: str | None = None) -> str:
    """Resolve an explicit or canonical mass model for one run."""
    if mass_model is not None:
        return _validate_mass_model(str(mass_model))
    metadata = getattr(run, "_metadata", {})
    model_meta = metadata.get("model", {}) if isinstance(metadata, dict) else {}
    params = model_meta.get("parameters", {}) if isinstance(model_meta, dict) else {}
    # Legacy storage did not record the model. Its stored neutral chain mass
    # corresponds to the end-group-aware representation.
    return _validate_mass_model(str(params.get("mass_model", "with_end_groups")))


def resolve_mass_model(population, mass_model: str | None = None) -> str:
    """Resolve an explicit or canonical mass model for a selected population."""
    run = getattr(population, "run", None)
    if run is None and mass_model is None:
        return "with_end_groups"
    return resolve_run_mass_model(run, mass_model)


def record_masses(population, mass_model: str | None = None) -> tuple[np.ndarray, str]:
    """Return exact per-record neutral molar masses and the resolved model."""
    model = resolve_mass_model(population, mass_model)
    calculator = getattr(population, "masses", None)
    if calculator is not None:
        mass = np.asarray(calculator(mass_model=model), dtype=float)
    else:
        raw = population._raw_arrays()
        if "mass" not in raw:
            raise InvalidOutputError("chain population has no molar-mass data")
        mass = np.asarray(raw["mass"], dtype=float)
    count = np.asarray(population.count)
    if mass.shape != count.shape:
        raise ValueError("per-chain masses and counts must have the same shape")
    if np.any(~np.isfinite(mass)) or np.any(mass <= 0):
        raise ValueError("molar masses must be finite and strictly positive")
    mass.flags.writeable = False
    return mass, model
