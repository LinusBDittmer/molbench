import numpy as np
import pytest

from molbench.comparison import Comparison
from molbench.molecule import Datapoint, Molecule, MoleculeList
from molbench.statistics import (
    Statistics,
    mae,
    median_se,
    mse,
    register_as_error_measure,
    rmsd,
    sde,
)
from molbench.statistics import (
    max as stat_max,
)
from molbench.statistics import (
    min as stat_min,
)

INTEREST = {"method": "HF"}
REFERENCE = {"method": "TBE"}
PROPTYPE = "energy"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _get_errors(comparison, interest=INTEREST, reference=REFERENCE):
    return Statistics(comparison).compare(interest, reference)


def _assign(proptype=PROPTYPE):
    return Statistics.assign_by_proptype(proptype)


# ---------------------------------------------------------------------------
# identify
# ---------------------------------------------------------------------------


def test_identify_reference(known_comparison):
    stats = Statistics(known_comparison)
    ident = stats.identify(INTEREST, REFERENCE)
    # A key-path matching reference method=TBE
    seps = ("water", "cc-pvdz", "TBE", "energy", "ref")
    assert ident(seps) == "reference"


def test_identify_interest(known_comparison):
    stats = Statistics(known_comparison)
    ident = stats.identify(INTEREST, REFERENCE)
    seps = ("water", "cc-pvdz", "HF", "energy", "computed")
    assert ident(seps) == "interest"


def test_identify_neither(known_comparison):
    stats = Statistics(known_comparison)
    ident = stats.identify(INTEREST, REFERENCE)
    seps = ("water", "cc-pvdz", "MP2", "energy", "other")
    assert ident(seps) is None


# ---------------------------------------------------------------------------
# compare — absolute errors
# ---------------------------------------------------------------------------


def test_compare_absolute_error(known_comparison):
    errors = _get_errors(known_comparison)
    assert len(errors) == 1
    for interest_dict in errors.values():
        for se in interest_dict.values():
            assert se.value == pytest.approx(0.1)


def test_compare_signed_positive(known_comparison):
    # interest (-75.9) - reference (-76.0) = +0.1
    errors = _get_errors(known_comparison)
    for id_ in errors.values():
        for se in id_.values():
            assert se.value > 0


def test_compare_empty_interest():
    ref = Molecule(
        "water",
        "ref",
        {},
        {
            "gs": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": Datapoint(-76.0, "au")},
            }
        },
    )
    c = Comparison()
    c.add(MoleculeList([ref]))
    errors = Statistics(c).compare({"method": "HF"}, {"method": "TBE"})
    assert errors == {}


# ---------------------------------------------------------------------------
# compare — relative errors
# ---------------------------------------------------------------------------


def test_compare_relative_error(known_comparison):
    errors = Statistics(known_comparison).compare(INTEREST, REFERENCE, relative=True)
    for id_ in errors.values():
        for se in id_.values():
            # (−75.9 − (−76.0)) / |−76.0| = 0.1/76.0
            assert se.value == pytest.approx(0.1 / 76.0, rel=1e-5)


def test_compare_relative_with_damping(known_comparison):
    errors = Statistics(known_comparison).compare(
        INTEREST, REFERENCE, relative=True, relative_damping=1.0
    )
    for id_ in errors.values():
        for se in id_.values():
            # (−75.9 − (−76.0)) / (|−76.0| + 1.0) = 0.1/77.0
            assert se.value == pytest.approx(0.1 / 77.0, rel=1e-5)


# ---------------------------------------------------------------------------
# evaluate
# ---------------------------------------------------------------------------


def test_evaluate_mse(known_comparison):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).evaluate(errors, "mse", proptype=PROPTYPE)
    val, count = result["mse"]
    assert val == pytest.approx(0.1)
    assert count == 1


def test_evaluate_mae(known_comparison):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).evaluate(errors, "mae", proptype=PROPTYPE)
    val, _ = result["mae"]
    assert val == pytest.approx(0.1)


def test_evaluate_rmsd(known_comparison):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).evaluate(errors, "rmsd", proptype=PROPTYPE)
    val, _ = result["rmsd"]
    assert val == pytest.approx(0.1)


def test_evaluate_sde_single_point(known_comparison):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).evaluate(errors, "sde", proptype=PROPTYPE)
    val, _ = result["sde"]
    # std of a single value is 0
    assert val == pytest.approx(0.0)


def test_evaluate_min_max(known_comparison):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).evaluate(
        errors, "min", "max", proptype=PROPTYPE
    )
    assert result["min"][0] == pytest.approx(0.1)
    assert result["max"][0] == pytest.approx(0.1)


def test_evaluate_median_se(known_comparison):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).evaluate(
        errors, "median_se", proptype=PROPTYPE
    )
    assert result["median_se"][0] == pytest.approx(0.1)


def test_evaluate_all_keyword(known_comparison):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).evaluate(errors, "all", proptype=PROPTYPE)
    for key in ("mse", "mae", "rmsd", "sde", "min", "max", "median_se"):
        assert key in result


def test_evaluate_unknown_measure_logs_error(known_comparison, capfd):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).evaluate(
        errors, "unknown_measure", proptype=PROPTYPE
    )
    assert "unknown_measure" not in result


def test_evaluate_no_assign_no_proptype_returns_none(known_comparison):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).evaluate(errors, "mse")
    assert result is None


# ---------------------------------------------------------------------------
# two-molecule: known MSE, MAE, RMSD
# ---------------------------------------------------------------------------


def test_mse_two_molecules(two_molecule_comparison):
    errors = _get_errors(two_molecule_comparison)
    result = Statistics(two_molecule_comparison).evaluate(
        errors, "mse", proptype=PROPTYPE
    )
    val, count = result["mse"]
    assert count == 2
    assert val == pytest.approx(0.2)  # (0.1 + 0.3) / 2


def test_mae_two_molecules(two_molecule_comparison):
    errors = _get_errors(two_molecule_comparison)
    result = Statistics(two_molecule_comparison).evaluate(
        errors, "mae", proptype=PROPTYPE
    )
    val, _ = result["mae"]
    assert val == pytest.approx(0.2)


def test_rmsd_two_molecules(two_molecule_comparison):
    errors = _get_errors(two_molecule_comparison)
    result = Statistics(two_molecule_comparison).evaluate(
        errors, "rmsd", proptype=PROPTYPE
    )
    val, _ = result["rmsd"]
    expected = np.sqrt((0.1**2 + 0.3**2) / 2)
    assert val == pytest.approx(expected)


def test_min_max_two_molecules(two_molecule_comparison):
    errors = _get_errors(two_molecule_comparison)
    result = Statistics(two_molecule_comparison).evaluate(
        errors, "min", "max", proptype=PROPTYPE
    )
    assert result["min"][0] == pytest.approx(0.1)
    assert result["max"][0] == pytest.approx(0.3)


# ---------------------------------------------------------------------------
# mae on empty errors
# ---------------------------------------------------------------------------


def test_mae_empty_errors(known_comparison):
    assign = _assign()
    val, count = mae({}, assign)
    assert val == pytest.approx(0.0)
    assert count == 0


# ---------------------------------------------------------------------------
# register_as_error_measure decorator
# ---------------------------------------------------------------------------


def test_register_as_error_measure():
    @register_as_error_measure
    def my_custom_measure(signed_errors, assign):
        return (42.0, 0)

    assert "my_custom_measure" in Statistics.available_error_measures
    # Clean up
    del Statistics.available_error_measures["my_custom_measure"]


# ---------------------------------------------------------------------------
# assign_by_proptype
# ---------------------------------------------------------------------------


def test_assign_by_proptype_matches():
    assign = Statistics.assign_by_proptype("energy")
    ref_keys = ("water", "cc-pvdz", "TBE", "energy", "ref")
    int_keys = ("water", "cc-pvdz", "HF", "energy", "computed")
    assert assign(ref_keys, int_keys) is True


def test_assign_by_proptype_no_match():
    assign = Statistics.assign_by_proptype("energy")
    ref_keys = ("water", "cc-pvdz", "TBE", "dipole", "ref")
    int_keys = ("water", "cc-pvdz", "HF", "energy", "computed")
    assert assign(ref_keys, int_keys) is False


def test_assign_by_proptype_different_ref_int():
    assign = Statistics.assign_by_proptype("dipole", refproptype="energy")
    ref_keys = ("water", "cc-pvdz", "TBE", "energy", "ref")
    int_keys = ("water", "cc-pvdz", "HF", "dipole", "computed")
    assert assign(ref_keys, int_keys) is True


# ---------------------------------------------------------------------------
# __init__ type guard
# ---------------------------------------------------------------------------


def test_init_wrong_type_exits_cleanly():
    # Must hit the intended log.critical() path immediately, not log a
    # non-fatal error and then crash later on an unrelated AttributeError.
    with pytest.raises(SystemExit):
        Statistics("not a comparison")


# ---------------------------------------------------------------------------
# relative error with a zero reference value
# ---------------------------------------------------------------------------


def test_compare_relative_zero_reference_skips_pair(caplog):
    ref = Molecule(
        "water",
        "ref",
        {},
        {
            "gs": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": Datapoint(0.0, "au")},
            }
        },
    )
    interest = Molecule(
        "water",
        "computed",
        {},
        {
            "gs": {
                "basis": "cc-pvdz",
                "method": "HF",
                "data": {"energy": Datapoint(1.0, "au")},
            }
        },
    )
    c = Comparison()
    c.add(MoleculeList([ref, interest]))
    with caplog.at_level("WARNING", logger="molbench"):
        errors = Statistics(c).compare(INTEREST, REFERENCE, relative=True)
    assert errors == {}
    assert any("is zero" in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# extreme_error_keys
# ---------------------------------------------------------------------------


def test_extreme_error_keys_two_molecules(two_molecule_comparison):
    errors = _get_errors(two_molecule_comparison)
    result = Statistics(two_molecule_comparison).extreme_error_keys(
        errors, proptype=PROPTYPE
    )
    assert result["min"]["value"] == pytest.approx(0.1)
    assert result["min"]["reference"][0] == "water"
    assert result["max"]["value"] == pytest.approx(0.3)
    assert result["max"]["reference"][0] == "benzene"


def test_extreme_error_keys_empty_returns_empty_dict(known_comparison):
    result = Statistics(known_comparison).extreme_error_keys({}, proptype=PROPTYPE)
    assert result == {}


def test_extreme_error_keys_no_assign_no_proptype_logs_error(known_comparison):
    errors = _get_errors(known_comparison)
    result = Statistics(known_comparison).extreme_error_keys(errors)
    assert result == {}


# ---------------------------------------------------------------------------
# empty-input behavior of the built-in error measures
# ---------------------------------------------------------------------------


def test_empty_errors_all_measures_no_crash():
    assign = _assign()
    for measure in (mse, sde, stat_min, stat_max, median_se, rmsd):
        val, count = measure({}, assign)
        assert count == 0
        assert np.isnan(val)


def test_evaluate_all_on_empty_errors_no_crash(known_comparison):
    # A proptype filter matching zero pairs used to crash "all" on min/max.
    result = Statistics(known_comparison).evaluate(
        {}, "all", proptype="nonexistent_proptype"
    )
    for key in ("mse", "mae", "rmsd", "sde", "min", "max", "median_se"):
        assert key in result
        assert result[key][1] == 0
