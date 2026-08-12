from __future__ import annotations

import numpy as np
import pytest


HOMO_CHANNEL_KINDS = (
    "rxn_", "macro_init_", "macro_prop_", "macro_term_c_", "macro_term_d_",
    "macro_term_x_", "macro_transfer_", "macro_transfer_m_", "macro_deprop_",
)

COPO_CHANNELS = {
    "rxn_0", "init_R_A", "init_R_B", "prop_PA_A", "prop_PA_B",
    "prop_PB_A", "prop_PB_B", "term_c_PA_PB", "term_d_PA_PB",
    "term_x_PA_CAP", "transfer_PA_CTA", "init_Rcta_A",
    "transfer_m_PA_A", "deprop_PA_A",
}

ACTION_TYPES = {
    "print", "print_info", "print_memory", "save", "save_chains", "set_k",
    "add_k", "set_temp", "add_temp", "set_c", "add_c", "feed",
}


def test_homo_all_reaction_kinds_reach_storage(homo_mechanisms_run):
    names = homo_mechanisms_run.event_counts.names
    assert all(any(name.startswith(prefix) for name in names) for prefix in HOMO_CHANNEL_KINDS)


def test_homo_all_declared_channels_fire(homo_mechanisms_run):
    final_counts = [int(homo_mechanisms_run.event_counts[name][-1]) for name in homo_mechanisms_run.event_counts.names]
    assert all(count > 0 for count in final_counts)


def test_homo_action_vocabulary_is_recorded(homo_mechanisms_run):
    assert ACTION_TYPES.issubset({action.type for action in homo_mechanisms_run.actions})


def test_homo_mutable_rates_and_temperature(homo_mechanisms_run):
    assert homo_mechanisms_run.k["kp"][-1] == pytest.approx(220000.0)
    assert homo_mechanisms_run.temp[-1] == pytest.approx(353.15)
    assert np.ptp(homo_mechanisms_run.k["kd"]) > 0


def test_homo_var_targets_round_trip(homo_mechanisms_run):
    expected = {"kp", "temperature", "CTA", "M", "CAP"}
    assert set(homo_mechanisms_run.var) == expected
    assert homo_mechanisms_run.var.temperature.unit == "K"


def test_homo_set_c_invalidates_only_requested_balance(homo_mechanisms_run):
    import pyslimmc

    with pytest.raises(pyslimmc.AnalysisNotApplicableError, match="set_c"):
        _ = homo_mechanisms_run.balance.total["CAP"]
    assert np.all(np.isfinite(homo_mechanisms_run.balance.total["M"]))


def test_homo_with_end_groups_mass_audit(homo_mechanisms_run):
    audit = homo_mechanisms_run.mass_audit(mass_model="with_end_groups")
    assert audit.ok and audit.checked_chains > 0


def test_homo_conditional_stop_is_clean(homo_stop_run):
    assert homo_stop_run.status == "completed"
    assert homo_stop_run.termination_reason == "stop_condition"
    assert [action.type for action in homo_stop_run.actions] == ["save", "save_chains", "stop"]


def test_copo_all_advanced_channels_reach_storage(copo_mechanisms_run):
    assert COPO_CHANNELS == set(copo_mechanisms_run.event_counts.names)


def test_copo_all_advanced_channels_fire(copo_mechanisms_run):
    assert all(
        int(copo_mechanisms_run.event_counts[name][-1]) > 0
        for name in copo_mechanisms_run.event_counts.names
    )


def test_copo_action_vocabulary_is_recorded(copo_mechanisms_run):
    assert ACTION_TYPES.issubset({action.type for action in copo_mechanisms_run.actions})


def test_copo_arrhenius_and_mutable_parameters(copo_mechanisms_run):
    assert copo_mechanisms_run.k["kp_aa"][-1] == pytest.approx(220000.0)
    assert copo_mechanisms_run.temp[-1] == pytest.approx(353.15)
    assert np.ptp(copo_mechanisms_run.k["kd"]) > 0


def test_copo_feed_and_full_sequences(copo_mechanisms_run):
    assert copo_mechanisms_run.feed_events.n_events == 1
    assert copo_mechanisms_run.last_with_chains.chains.has_sequences


def test_copo_with_end_groups_mass_audit(copo_mechanisms_run):
    audit = copo_mechanisms_run.mass_audit(mass_model="with_end_groups")
    assert audit.ok and audit.checked_chains > 0

