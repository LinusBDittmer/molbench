"""Integration: Statistics.compare → evaluate with known analytical results."""

import numpy as np
import pytest

from molbench.comparison import Comparison
from molbench.molecule import Datapoint, Molecule, MoleculeList
from molbench.statistics import Statistics

INTEREST = {"method": "HF"}
REFERENCE = {"method": "TBE"}
PROPTYPE = "energy"


def _build_comparison(pairs):
    """Build a Comparison from (name, ref_energy, int_energy) triples (au)."""
    mols = MoleculeList()
    for name, ref_e, int_e in pairs:
        mols.append(
            Molecule(
                name,
                "ref",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "TBE",
                        "data": {"energy": Datapoint(ref_e, "au")},
                    }
                },
            )
        )
        mols.append(
            Molecule(
                name,
                "computed",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "HF",
                        "data": {"energy": Datapoint(int_e, "au")},
                    }
                },
            )
        )
    c = Comparison()
    c.add(mols)
    return c


def _eval(c, *measures):
    stats = Statistics(c)
    errors = stats.compare(INTEREST, REFERENCE)
    return stats.evaluate(errors, *measures, proptype=PROPTYPE)


# ---------------------------------------------------------------------------
# Single molecule: delta = +0.1
# ---------------------------------------------------------------------------


def test_mse_single_delta():
    c = _build_comparison([("water", -76.0, -75.9)])
    result = _eval(c, "mse")
    assert result["mse"][0] == pytest.approx(0.1)
    assert result["mse"][1] == 1


def test_mae_single_delta():
    c = _build_comparison([("water", -76.0, -75.9)])
    result = _eval(c, "mae")
    assert result["mae"][0] == pytest.approx(0.1)


def test_rmsd_single_delta():
    c = _build_comparison([("water", -76.0, -75.9)])
    result = _eval(c, "rmsd")
    assert result["rmsd"][0] == pytest.approx(0.1)


def test_sde_single_point():
    c = _build_comparison([("water", -76.0, -75.9)])
    result = _eval(c, "sde")
    assert result["sde"][0] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Two molecules: errors +0.1 and +0.3
# ---------------------------------------------------------------------------


def test_mse_two_molecules():
    c = _build_comparison([("water", -76.0, -75.9), ("benzene", -10.0, -9.7)])
    result = _eval(c, "mse")
    assert result["mse"][0] == pytest.approx(0.2)
    assert result["mse"][1] == 2


def test_mae_two_molecules():
    c = _build_comparison([("water", -76.0, -75.9), ("benzene", -10.0, -9.7)])
    result = _eval(c, "mae")
    assert result["mae"][0] == pytest.approx(0.2)


def test_rmsd_two_molecules():
    c = _build_comparison([("water", -76.0, -75.9), ("benzene", -10.0, -9.7)])
    result = _eval(c, "rmsd")
    expected = np.sqrt((0.1**2 + 0.3**2) / 2)
    assert result["rmsd"][0] == pytest.approx(expected)


def test_min_max_two_molecules():
    c = _build_comparison([("water", -76.0, -75.9), ("benzene", -10.0, -9.7)])
    result = _eval(c, "min", "max")
    assert result["min"][0] == pytest.approx(0.1)
    assert result["max"][0] == pytest.approx(0.3)


def test_median_two_molecules():
    c = _build_comparison([("water", -76.0, -75.9), ("benzene", -10.0, -9.7)])
    result = _eval(c, "median_se")
    assert result["median_se"][0] == pytest.approx(0.2)


def test_evaluate_all_keyword():
    c = _build_comparison([("water", -76.0, -75.9)])
    result = _eval(c, "all")
    for key in ("mse", "mae", "rmsd", "sde", "min", "max", "median_se"):
        assert key in result


# ---------------------------------------------------------------------------
# Relative errors
# ---------------------------------------------------------------------------


def test_relative_error():
    c = _build_comparison([("water", -76.0, -75.9)])
    stats = Statistics(c)
    errors = stats.compare(INTEREST, REFERENCE, relative=True)
    result = stats.evaluate(errors, "mse", proptype=PROPTYPE)
    expected = 0.1 / 76.0
    assert result["mse"][0] == pytest.approx(expected, rel=1e-5)


def test_relative_error_with_damping():
    c = _build_comparison([("water", -76.0, -75.9)])
    stats = Statistics(c)
    errors = stats.compare(INTEREST, REFERENCE, relative=True, relative_damping=1.0)
    result = stats.evaluate(errors, "mse", proptype=PROPTYPE)
    expected = 0.1 / 77.0
    assert result["mse"][0] == pytest.approx(expected, rel=1e-5)


# ---------------------------------------------------------------------------
# Empty comparison
# ---------------------------------------------------------------------------


def test_empty_interest_empty_errors():
    c = _build_comparison([("water", -76.0, -75.9)])
    stats = Statistics(c)
    # Ask for method "MP2" as interest — doesn't exist
    errors = stats.compare({"method": "MP2"}, {"method": "TBE"})
    assert errors == {}


def test_evaluate_empty_errors():
    c = _build_comparison([("water", -76.0, -75.9)])
    stats = Statistics(c)
    errors = {}
    result = stats.evaluate(errors, "mae", proptype=PROPTYPE)
    val, count = result["mae"]
    assert count == 0
    assert val == pytest.approx(0.0)
