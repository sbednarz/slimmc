"""Mathematical oracle for the pyslimmc distribution API redesign.

This file is intentionally contract-first (TDD): it targets the proposed vNext
API and serves as an independent mathematical oracle for the redesigned
distribution API.

The fixture is deliberately chosen so that two distinct chain records have the
same DP but different molar masses.  Therefore a correct implementation must
project CLD by DP and MWD independently by actual chain mass.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import pyslimmc


# Exact source population used by every test:
#
# record   DP   repeat-unit mass   count   composition [A, B]
# A         1        100             2          [1, 0]
# B         2        200             3          [2, 0]
# C         2        250             1          [1, 1]
# D         3        400             4          [1, 2]
#
# Monomer mass increments: A=100, B=150 g/mol.
# Total chains = 10; total polymer mass = 2650 (common proportional unit).

DP = np.array([1.0, 2.0, 2.0, 3.0])
MASS = np.array([100.0, 200.0, 250.0, 400.0])
COUNT = np.array([2, 3, 1, 4], dtype=np.uint64)

EXPECTED_DP = np.array([1.0, 2.0, 3.0])
EXPECTED_DP_COUNT = np.array([2, 4, 4], dtype=np.uint64)
EXPECTED_MASS = np.array([100.0, 200.0, 250.0, 400.0])
EXPECTED_MASS_COUNT = np.array([2, 3, 1, 4], dtype=np.uint64)

EXPECTED_CLD_NUMBER = np.array([0.2, 0.4, 0.4])
EXPECTED_CLD_MASS = np.array([200.0, 850.0, 1600.0]) / 2650.0
EXPECTED_CLD_Z = np.array([2.0, 16.0, 36.0]) / 54.0

EXPECTED_MWD_NUMBER = np.array([0.2, 0.3, 0.1, 0.4])
EXPECTED_MWD_MASS = np.array([200.0, 600.0, 250.0, 1600.0]) / 2650.0
EXPECTED_MWD_Z = np.array([
    2.0 * 100.0**2,
    3.0 * 200.0**2,
    1.0 * 250.0**2,
    4.0 * 400.0**2,
])
EXPECTED_MWD_Z /= EXPECTED_MWD_Z.sum()

EXPECTED_DPN = 2.2
EXPECTED_DPW = 2.4545454545454546
EXPECTED_DPZ = 2.6296296296296298

EXPECTED_MN = 265.0
EXPECTED_MW = 317.92452830188677
EXPECTED_MZ = 353.2640949554896
EXPECTED_MASS_DISPERSITY = EXPECTED_MW / EXPECTED_MN
EXPECTED_DP_DISPERSITY = EXPECTED_DPW / EXPECTED_DPN


def _write_table(root: Path, name: str, columns: dict[str, np.ndarray]) -> None:
    table = root / name
    table.mkdir(parents=True)
    for column, values in columns.items():
        np.save(table / f"{column}.npy", values, allow_pickle=False)


def _build_oracle_storage(root: Path) -> None:
    """Create the smallest storage fixture that preserves the oracle population."""
    metadata = {
        "run_id": "distribution-math-oracle",
        "storage": "slimmc-storage",
        "storage_format_version": "1.2.0",
        "run_status": "completed",
        "validation_error_count": 0,
        "engine": "slimmc-copo",
        "kinetic_model": "copo",
        "sequence_mode": "composition",
    }
    (root / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    schema: list[dict[str, object]] = [
        {
            "record_type": "schema_header",
            "schema_name": "slimmc-storage",
            "schema_version": "1.2.0",
        }
    ]
    for table in (
        "snapshots",
        "state",
        "chains",
        "chain_composition",
        "sequences",
        "moments",
    ):
        schema.append({"record_type": "table", "name": table, "required": True})

    for i, (name, mass) in enumerate((('A', 100.0), ('B', 150.0))):
        schema.extend(
            [
                {
                    "record_type": "dictionary_entry",
                    "dictionary": "monomers",
                    "id": i,
                    "name": name,
                    "molar_mass_increment": mass,
                },
                {
                    "record_type": "dictionary_entry",
                    "dictionary": "state_entities",
                    "id": i,
                    "name": name,
                    "kind": "monomer",
                },
            ]
        )

    dictionaries = {
        "chain_populations": ("live", "dead"),
        "chain_pools": ("not_applicable", "terminal_A"),
        "chain_origins": ("unknown", "init"),
        "chain_end_types": ("not_applicable", "unknown"),
        "population_scope": ("all", "live", "dead"),
        "mass_bases": ("repeat_units", "with_end_groups"),
    }
    for dictionary, names in dictionaries.items():
        for i, name in enumerate(names):
            schema.append(
                {
                    "record_type": "dictionary_entry",
                    "dictionary": dictionary,
                    "id": i,
                    "name": name,
                }
            )

    (root / "schema.jsonl").write_text(
        "".join(json.dumps(row, separators=(",", ":")) + "\n" for row in schema)
    )

    _write_table(
        root,
        "snapshots",
        {
            "snapshot_id": np.array([0, 1], dtype="<u8"),
            "time": np.array([0.0, 1.0]),
            "kmc_event": np.array([0, 10], dtype="<u8"),
            "snapshot_reason_id": np.array([0, 4], dtype="<u4"),
            "is_final": np.array([False, True]),
            "has_chains": np.array([False, True]),
            "has_sequences": np.array([False, False]),
            "kinetic_parameter_set_id": np.zeros(2, dtype="<u8"),
        },
    )

    _write_table(
        root,
        "state",
        {
            "snapshot_id": np.repeat(np.arange(2, dtype="<u8"), 2),
            "entity_id": np.tile(np.arange(2, dtype="<u4"), 2),
            "count": np.array([10, 10, 4, 6], dtype="<u8"),
            "moles": np.array([10.0, 10.0, 4.0, 6.0]),
            "concentration": np.array([10.0, 10.0, 4.0, 6.0]),
        },
    )

    # All four records are dead.  Stored molar_mass equals repeat-unit mass in
    # this oracle so both supported mass models intentionally give the same
    # expected projection.
    _write_table(
        root,
        "chains",
        {
            "chain_record_id": np.arange(4, dtype="<u8"),
            "snapshot_id": np.ones(4, dtype="<u8"),
            "population_id": np.ones(4, dtype="<u4"),
            "pool_id": np.zeros(4, dtype="<u4"),
            "origin_id": np.ones(4, dtype="<u4"),
            "dp": DP.astype("<u8"),
            "molar_mass": MASS,
            "count": COUNT,
            "moles": COUNT.astype(float),
            "concentration": COUNT.astype(float),
            "left_end_id": np.ones(4, dtype="<u4"),
            "right_end_id": np.ones(4, dtype="<u4"),
            "has_first_monomer": np.ones(4, dtype=bool),
            "first_monomer_id": np.array([0, 0, 0, 0], dtype="<u4"),
            "has_penultimate_monomer": np.ones(4, dtype=bool),
            "penultimate_monomer_id": np.array([0, 0, 0, 1], dtype="<u4"),
            "has_last_monomer": np.ones(4, dtype=bool),
            "last_monomer_id": np.array([0, 0, 1, 1], dtype="<u4"),
            "has_sequence": np.zeros(4, dtype=bool),
            "sequence_offset": np.zeros(4, dtype="<u8"),
            "sequence_length": np.zeros(4, dtype="<u8"),
        },
    )

    # [A, B] composition gives masses 100, 200, 250, 400 g/mol.
    composition = np.array(
        [
            [1, 0],
            [2, 0],
            [1, 1],
            [1, 2],
        ],
        dtype="<u8",
    )
    _write_table(
        root,
        "chain_composition",
        {
            "chain_record_id": np.repeat(np.arange(4, dtype="<u8"), 2),
            "monomer_id": np.tile(np.arange(2, dtype="<u4"), 4),
            "unit_count": composition.reshape(-1),
        },
    )

    _write_table(root, "sequences", {"symbols": np.array([], dtype="<u4")})

    # Storage moments are internally consistent with the oracle.  The new API
    # tests below nevertheless require distribution result moments to be
    # invariant across form, rather than derived from displayed y values.
    rows = 4  # all/dead x repeat_units/with_end_groups
    _write_table(
        root,
        "moments",
        {
            "snapshot_id": np.ones(rows, dtype="<u8"),
            "population_scope_id": np.array([0, 0, 2, 2], dtype="<u4"),
            "mass_basis_id": np.array([0, 1, 0, 1], dtype="<u4"),
            "chain_count": np.full(rows, 10, dtype="<u8"),
            "sum_dp": np.full(rows, 22.0),
            "sum_dp2": np.full(rows, 54.0),
            "dp_n": np.full(rows, EXPECTED_DPN),
            "dp_w": np.full(rows, EXPECTED_DPW),
            "sum_molar_mass": np.full(rows, 2650.0),
            "sum_molar_mass2": np.full(rows, 842500.0),
            "sum_molar_mass3": np.full(rows, 297625000.0),
            "mn": np.full(rows, EXPECTED_MN),
            "mw": np.full(rows, EXPECTED_MW),
            "mz": np.full(rows, EXPECTED_MZ),
            "dispersity": np.full(rows, EXPECTED_MASS_DISPERSITY),
        },
    )

    (root / "RESULTS_COMPLETE").write_text("slimmc-storage-v1\n")


@pytest.fixture
def oracle_run(tmp_path: Path):
    root = tmp_path / "run"
    root.mkdir()
    _build_oracle_storage(root)
    return pyslimmc.open(root)


def test_dp_counts_is_exact_projection(oracle_run):
    c = oracle_run.dp_counts(pool="dead")
    np.testing.assert_array_equal(c.dp, EXPECTED_DP)
    np.testing.assert_array_equal(c.count, EXPECTED_DP_COUNT)
    assert c.total_chains == 10
    assert c.total_repeat_units == 22


def test_mass_counts_is_exact_projection(oracle_run):
    c = oracle_run.mass_counts(pool="dead", mass_model="repeat_units")
    np.testing.assert_allclose(c.mass, EXPECTED_MASS)
    np.testing.assert_array_equal(c.count, EXPECTED_MASS_COUNT)
    assert c.total_chains == 10


def test_same_dp_can_map_to_multiple_masses(oracle_run):
    """Core architecture invariant: MWD must not be reconstructed from CLD."""
    dp = oracle_run.dp_counts(pool="dead")
    mass = oracle_run.mass_counts(pool="dead", mass_model="repeat_units")
    assert len(dp.dp) == 3
    assert len(mass.mass) == 4


@pytest.mark.parametrize(
    ("form", "expected"),
    [
        ("number", EXPECTED_CLD_NUMBER),
        ("mass", EXPECTED_CLD_MASS),
        ("z", EXPECTED_CLD_Z),
    ],
)
def test_cld_discrete_forms(oracle_run, form, expected):
    d = oracle_run.cld(pool="dead", form=form, mass_model="repeat_units")
    np.testing.assert_allclose(d.x, EXPECTED_DP)
    np.testing.assert_allclose(d.y, expected)
    assert d.form == form
    assert d.y.sum() == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("form", "expected"),
    [
        ("number", EXPECTED_MWD_NUMBER),
        ("mass", EXPECTED_MWD_MASS),
        ("z", EXPECTED_MWD_Z),
    ],
)
def test_mwd_discrete_forms(oracle_run, form, expected):
    d = oracle_run.mwd(pool="dead", form=form, mass_model="repeat_units")
    np.testing.assert_allclose(d.x, EXPECTED_MASS)
    np.testing.assert_allclose(d.y, expected)
    assert d.form == form
    assert d.y.sum() == pytest.approx(1.0)


def test_cld_mass_uses_actual_chain_mass_within_dp_class(oracle_run):
    """DP=2 contains both M=200 and M=250 and must aggregate both masses."""
    cld = oracle_run.cld(pool="dead", form="mass", mass_model="repeat_units")
    mwd = oracle_run.mwd(pool="dead", form="mass", mass_model="repeat_units")

    cld_dp2 = cld.y[np.flatnonzero(cld.x == 2.0)[0]]
    mwd_dp2_mass = mwd.y[np.isin(mwd.x, [200.0, 250.0])].sum()

    assert cld_dp2 == pytest.approx(850.0 / 2650.0)
    assert cld_dp2 == pytest.approx(mwd_dp2_mass)


def test_log_mwd_changes_support_not_atomic_mass_fractions(oracle_run):
    linear = oracle_run.mwd(pool="dead", form="mass", mass_model="repeat_units")
    log = oracle_run.mwd(pool="dead", form="log", mass_model="repeat_units")

    np.testing.assert_allclose(log.x, np.log10(linear.x))
    np.testing.assert_allclose(log.y, linear.y)
    assert log.form == "log"
    assert log.y.sum() == pytest.approx(1.0)


def test_log_cld_changes_support_not_atomic_mass_fractions(oracle_run):
    linear = oracle_run.cld(pool="dead", form="mass", mass_model="repeat_units")
    log = oracle_run.cld(pool="dead", form="log", mass_model="repeat_units")

    np.testing.assert_allclose(log.x, np.log10(linear.x))
    np.testing.assert_allclose(log.y, linear.y)
    assert log.form == "log"
    assert log.y.sum() == pytest.approx(1.0)


@pytest.mark.parametrize("form", ["number", "mass", "z", "log"])
def test_mwd_exact_moments_are_invariant_to_form(oracle_run, form):
    d = oracle_run.mwd(pool="dead", form=form, mass_model="repeat_units")
    assert d.mn == pytest.approx(EXPECTED_MN)
    assert d.mw == pytest.approx(EXPECTED_MW)
    assert d.mz == pytest.approx(EXPECTED_MZ)
    assert d.dispersity == pytest.approx(EXPECTED_MASS_DISPERSITY)


@pytest.mark.parametrize("form", ["number", "mass", "z", "log"])
def test_cld_exact_moments_are_invariant_to_form(oracle_run, form):
    d = oracle_run.cld(pool="dead", form=form, mass_model="repeat_units")
    assert d.dpn == pytest.approx(EXPECTED_DPN)
    assert d.dpw == pytest.approx(EXPECTED_DPW)
    assert d.dpz == pytest.approx(EXPECTED_DPZ)
    assert d.dispersity == pytest.approx(EXPECTED_DP_DISPERSITY)


def test_cld_does_not_expose_molar_mass_moment_aliases(oracle_run):
    d = oracle_run.cld(pool="dead")
    assert not hasattr(d, "mn")
    assert not hasattr(d, "mw")
    assert not hasattr(d, "mz")


def test_exact_count_objects_do_not_expose_generic_xy(oracle_run):
    dp = oracle_run.dp_counts(pool="dead")
    mass = oracle_run.mass_counts(pool="dead", mass_model="repeat_units")
    assert not hasattr(dp, "x")
    assert not hasattr(dp, "y")
    assert not hasattr(mass, "x")
    assert not hasattr(mass, "y")


def test_population_moments_separates_dp_and_mass_dispersity(oracle_run):
    m = oracle_run.final.chains.dead.moments(mass_model="repeat_units")
    assert m.total_chains == 10
    assert m.dpn == pytest.approx(EXPECTED_DPN)
    assert m.dpw == pytest.approx(EXPECTED_DPW)
    assert m.dpz == pytest.approx(EXPECTED_DPZ)
    assert m.mn == pytest.approx(EXPECTED_MN)
    assert m.mw == pytest.approx(EXPECTED_MW)
    assert m.mz == pytest.approx(EXPECTED_MZ)
    assert m.dp_dispersity == pytest.approx(EXPECTED_DP_DISPERSITY)
    assert m.mass_dispersity == pytest.approx(EXPECTED_MASS_DISPERSITY)
    assert m.dp_dispersity != pytest.approx(m.mass_dispersity)
    assert m.mass_model == "repeat_units"
    assert m.source == "chains"
    assert m.has_dpz
    assert not hasattr(m, "dispersity")


def test_run_moments_callable_delegates_to_selected_population(oracle_run):
    m = oracle_run.moments(
        snapshot="final", population="dead", mass_model="repeat_units"
    )
    direct = oracle_run.final.chains.dead.moments(mass_model="repeat_units")
    assert m == direct


def test_snapshot_moments_callable(oracle_run):
    m = oracle_run.final.moments(population="dead", mass_model="repeat_units")
    assert m.dpn == pytest.approx(EXPECTED_DPN)
    assert m.dpz == pytest.approx(EXPECTED_DPZ)
    assert m.mn == pytest.approx(EXPECTED_MN)


def test_legacy_moments_tree_is_removed(oracle_run):
    assert not hasattr(oracle_run.moments, "all")
    assert not hasattr(oracle_run.moments, "live")
    assert not hasattr(oracle_run.moments, "dead")
    assert not hasattr(oracle_run.moments, "default")
    assert not hasattr(oracle_run.moments, "select")


def test_run_moments_falls_back_to_stored_aggregates_without_chains(tmp_path):
    root = tmp_path / "aggregate_only"
    root.mkdir()
    _build_oracle_storage(root)
    # Force the final snapshot to advertise no stored chains while retaining
    # the exact aggregate moments table.
    np.save(root / "snapshots" / "has_chains.npy", np.array([False, False], dtype=bool), allow_pickle=False)
    run = pyslimmc.open(root)
    m = run.moments(snapshot="final", population="dead", mass_model="repeat_units")
    assert m.source == "stored_aggregate"
    assert m.total_chains == 10
    assert m.dpn == pytest.approx(EXPECTED_DPN)
    assert m.dpw == pytest.approx(EXPECTED_DPW)
    assert not m.has_dpz
    assert np.isnan(m.dpz)
    assert m.mn == pytest.approx(EXPECTED_MN)
    assert m.mw == pytest.approx(EXPECTED_MW)
    assert m.mz == pytest.approx(EXPECTED_MZ)



def _oracle_partition(oracle_run):
    chains = oracle_run.final.chains.dead
    left = chains._with_mask(np.array([True, True, False, False]))
    right = chains._with_mask(np.array([False, False, True, True]))
    return chains, left, right


def test_mwd_series_per_series_keeps_each_series_normalized(oracle_run):
    _, left, right = _oracle_partition(oracle_run)
    group = oracle_run.mwd_series(
        series={"left": left, "right": right},
        form="mass",
        normalization="per_series",
        mass_model="repeat_units",
    )
    assert group.series_names == ("left", "right")
    assert group.normalization == "per_series"
    assert group.series_disjoint is True
    assert group["left"].y.sum() == pytest.approx(1.0)
    assert group["right"].y.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(group["left"].mass, [100.0, 200.0])
    np.testing.assert_allclose(group["right"].mass, [250.0, 400.0])


def test_mwd_series_combined_preserves_exact_mass_contributions(oracle_run):
    _, left, right = _oracle_partition(oracle_run)
    group = oracle_run.mwd_series(
        series={"left": left, "right": right},
        form="mass",
        normalization="combined",
        mass_model="repeat_units",
    )
    assert group["left"].y.sum() == pytest.approx(800.0 / 2650.0)
    assert group["right"].y.sum() == pytest.approx(1850.0 / 2650.0)
    assert sum(d.y.sum() for d in group.series.values()) == pytest.approx(1.0)


def test_mwd_series_combined_uses_form_specific_exact_source_total(oracle_run):
    _, left, right = _oracle_partition(oracle_run)
    number = oracle_run.mwd_series(
        series={"left": left, "right": right}, form="number", normalization="combined"
    )
    assert number["left"].y.sum() == pytest.approx(5.0 / 10.0)
    assert number["right"].y.sum() == pytest.approx(5.0 / 10.0)

    z = oracle_run.mwd_series(
        series={"left": left, "right": right},
        form="z", normalization="combined", mass_model="repeat_units"
    )
    assert z["left"].y.sum() == pytest.approx(140000.0 / 842500.0)
    assert z["right"].y.sum() == pytest.approx(702500.0 / 842500.0)


def test_cld_series_combined_uses_dp_z_total_for_z_form(oracle_run):
    _, left, right = _oracle_partition(oracle_run)
    group = oracle_run.cld_series(
        series={"left": left, "right": right}, form="z", normalization="combined"
    )
    assert group["left"].y.sum() == pytest.approx(14.0 / 54.0)
    assert group["right"].y.sum() == pytest.approx(40.0 / 54.0)


def test_combined_rejects_overlapping_series(oracle_run):
    all_chains, left, _ = _oracle_partition(oracle_run)
    with pytest.raises(ValueError, match="pairwise-disjoint"):
        oracle_run.mwd_series(
            series={"all": all_chains, "left": left},
            form="mass", normalization="combined", mass_model="repeat_units"
        )


def test_per_series_allows_overlapping_series(oracle_run):
    all_chains, left, _ = _oracle_partition(oracle_run)
    group = oracle_run.mwd_series(
        series={"all": all_chains, "left": left},
        form="mass", normalization="per_series", mass_model="repeat_units"
    )
    assert group.series_disjoint is False
    assert group["all"].y.sum() == pytest.approx(1.0)
    assert group["left"].y.sum() == pytest.approx(1.0)


def test_series_normalization_rejects_removed_modes(oracle_run):
    _, left, right = _oracle_partition(oracle_run)
    for mode in ("absolute", "reference"):
        with pytest.raises(ValueError, match="per_series.*combined"):
            oracle_run.mwd_series(
                series={"left": left, "right": right}, normalization=mode
            )


def test_combined_scaling_does_not_change_source_moments(oracle_run):
    _, left, right = _oracle_partition(oracle_run)
    per = oracle_run.mwd_series(
        series={"left": left, "right": right}, form="mass", normalization="per_series",
        mass_model="repeat_units",
    )
    combined = oracle_run.mwd_series(
        series={"left": left, "right": right}, form="mass", normalization="combined",
        mass_model="repeat_units",
    )
    for name in ("left", "right"):
        assert combined[name].mn == pytest.approx(per[name].mn)
        assert combined[name].mw == pytest.approx(per[name].mw)
        assert combined[name].mz == pytest.approx(per[name].mz)


def test_chain_population_selection_is_explicit(oracle_run):
    chains = oracle_run.final.chains
    assert not hasattr(chains, "select")
    with pytest.raises(TypeError):
        chains.mwd(pool="dead")
    with pytest.raises(TypeError):
        chains.cld(pool="dead")
    with pytest.raises(TypeError):
        chains.sec(pool="dead", sigma_log10M=0.05)

    # Convenience pool selection remains on Snapshot/Run only.
    assert oracle_run.mwd(pool="dead").total_chains == chains.dead.total_chains
    assert oracle_run.final.cld(pool="dead").total_chains == chains.dead.total_chains
