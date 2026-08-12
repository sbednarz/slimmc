from __future__ import annotations

import pytest

from pyslimmc.run import DataConsistencyError, Variable, Variables


def test_variables_collection_and_interactive_access():
    variables = Variables([
        {"kind": "monomer", "name": "IA", "value": 0.1, "unit": "mol_L"},
        {"kind": "param", "name": "temperature", "value": 343.15, "unit": "K"},
    ])
    assert list(variables) == ["IA", "temperature"]
    assert variables["IA"] == Variable("monomer", "IA", 0.1, "mol_L")
    assert variables.temperature.value == 343.15
    assert "temperature" in dir(variables)
    assert variables.info_text().startswith("Variables")


def test_variables_reject_duplicate_name():
    with pytest.raises(DataConsistencyError, match="duplicate variable name"):
        Variables([
            {"kind": "rate", "name": "kp", "value": 1.0, "unit": "L_mol_s"},
            {"kind": "rate", "name": "kp", "value": 2.0, "unit": "L_mol_s"},
        ])


def test_variables_reject_missing_or_non_numeric_value():
    with pytest.raises(DataConsistencyError, match="missing 'unit'"):
        Variables([{"kind": "param", "name": "temperature", "value": 343.15}])
    with pytest.raises(DataConsistencyError, match="must be numeric"):
        Variables([{"kind": "param", "name": "temperature", "value": "hot", "unit": "K"}])
