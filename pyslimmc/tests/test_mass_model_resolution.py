from types import SimpleNamespace

import pytest

from pyslimmc.mass_model import resolve_mass_model


def _population(metadata):
    return SimpleNamespace(run=SimpleNamespace(_metadata=metadata))


def test_explicit_mass_model_wins():
    pop = _population({"model": {"parameters": {"mass_model": "with_end_groups"}}})
    assert resolve_mass_model(pop, "repeat_units") == "repeat_units"


def test_mass_model_defaults_to_recorded_model_parameter():
    pop = _population({"model": {"parameters": {"mass_model": "repeat_units"}}})
    assert resolve_mass_model(pop, None) == "repeat_units"


def test_legacy_mass_model_defaults_to_stored_neutral_mass_basis():
    pop = _population({})
    assert resolve_mass_model(pop, None) == "with_end_groups"


def test_invalid_mass_model_is_rejected():
    pop = _population({})
    with pytest.raises(ValueError):
        resolve_mass_model(pop, "mystery")
