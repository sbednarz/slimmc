from __future__ import annotations

import numpy as np
import pytest


def test_penultimate_topology_has_twelve_propagation_channels(copo_penultimate_run):
    propagation = [name for name in copo_penultimate_run.event_counts.names if name.startswith("prop_")]
    assert len(propagation) == 12


def test_penultimate_parameter_analysis_is_available(copo_penultimate_run):
    result = copo_penultimate_run.copolymerization.penultimate_parameters()
    assert len(result) == len(copo_penultimate_run.t)


def test_penultimate_rejects_terminal_mayo_lewis(copo_penultimate_run):
    import pyslimmc

    with pytest.raises(pyslimmc.ChemicalModelIncompatibleError):
        copo_penultimate_run.copolymerization.mayo_lewis()


def test_terpolymer_round_trips_three_monomers(copo_terpolymer_run):
    assert copo_terpolymer_run.monomer_names == ("A", "B", "C")
    fractions = np.vstack([copo_terpolymer_run.F.cum[name] for name in copo_terpolymer_run.monomer_names])
    mask = np.all(np.isfinite(fractions), axis=0)
    assert mask.any() and np.allclose(fractions[:, mask].sum(axis=0), 1.0)


def test_terpolymer_binary_analysis_is_not_applicable(copo_terpolymer_run):
    import pyslimmc

    with pytest.raises(pyslimmc.ChemicalAnalysisNotApplicableError):
        copo_terpolymer_run.copolymerization.mayo_lewis()


def test_composition_mode_retains_aggregate_microstructure(copo_composition_run):
    chains = copo_composition_run.last_with_chains.chains
    assert not chains.has_sequences
    assert len(copo_composition_run.microstructure.dyads()) == 4


def test_composition_mode_rejects_sequence_only_analysis(copo_composition_run):
    import pyslimmc

    with pytest.raises(pyslimmc.DataUnavailableError, match="sequences"):
        copo_composition_run.last_with_chains.chains.sequence_stats()

