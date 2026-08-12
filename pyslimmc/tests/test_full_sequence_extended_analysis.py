from __future__ import annotations
import tempfile
from pathlib import Path
import numpy as np
import pyslimmc
from pyslimmc.tests.test_l2_4_chains import build


def _run():
    td = tempfile.TemporaryDirectory()
    root = Path(td.name) / "run"; root.mkdir(); build(root, full=True)
    return td, pyslimmc.open(root)


def test_motif_counts_and_ngrams_are_multiplicity_weighted():
    td, run = _run()
    try:
        c = run.final.chains
        motif = c.motif_counts("ABA")
        assert motif.occurrence_count.tolist() == [0, 1, 0]
        assert motif.normalized_frequency.tolist() == [0.0, 1.0, 0.0]
        grams = c.ngrams(2)
        data = {"".join(k): v for k, v in zip(grams.motifs, grams.count)}
        assert data == {"AA": 2.0, "AB": 3.0, "BA": 1.0, "BB": 4.0}
        assert np.isclose(np.sum(grams.fraction), 1.0)
    finally:
        td.cleanup()


def test_position_profile_and_microstructure_map():
    td, run = _run()
    try:
        c = run.final.chains
        profile = c.position_profile(bins=3)
        assert profile.names == ("A", "B")
        assert np.allclose(profile.fraction["A"] + profile.fraction["B"], 1.0)
        assert profile.plot().get_xlabel() == "relative chain position"
        result = c.microstructure_map("transition_fraction", dp_bins=[1.5,2.5,3.5], value_bins=[0,0.5,1.01])
        assert result.values.shape == (2,2)
        assert np.sum(result.values) == 7
        assert result.plot().get_xlabel() == "DP"
    finally:
        td.cleanup()


def test_run_plot_extended_full_shortcuts():
    td, run = _run()
    try:
        assert run.plot.ngrams(2).get_xlabel() == "2-gram"
        assert run.plot.position_profile(bins=3).get_xlabel() == "relative chain position"
        assert run.plot.microstructure_map("transition_fraction", dp_bins=2, value_bins=2).get_xlabel() == "DP"
    finally:
        td.cleanup()


def test_sequence_stats_cache_is_shared_with_filtered_views():
    td, run = _run()
    try:
        chains = run.final.chains
        root_stats = chains.sequence_stats(progress=False)
        assert chains.sequence_stats(progress=False) is root_stats

        subset = chains.where(dp_min=3)
        subset_stats = subset.sequence_stats(progress=False)
        expected = np.flatnonzero(np.asarray(chains.dp) >= 3)
        assert subset._analysis_root is chains
        assert np.array_equal(subset_stats.transition_count, root_stats.transition_count[expected])
        assert np.array_equal(
            subset_stats.max_block_length["A"],
            root_stats.max_block_length["A"][expected],
        )
        assert not subset_stats.transition_count.flags.writeable
    finally:
        td.cleanup()
