import pytest
import numpy
from molbench.molecule import Molecule, MoleculeList, Datapoint
from molbench.comparison import Comparison


def _mol(name, data_id, basis, method, energy, unit="au", transition_id=None):
    state = {
        "basis": basis,
        "method": method,
        "data": {"energy": Datapoint(energy, unit)},
    }
    if transition_id:
        state["transition_id"] = transition_id
    return Molecule(name, data_id, {}, {"gs": state})


# ---------------------------------------------------------------------------
# Initialization
# ---------------------------------------------------------------------------

def test_default_separators():
    c = Comparison()
    assert c.data_separators == ("basis", "method")


def test_custom_separators():
    c = Comparison("basis", "method", "transition_id")
    assert "transition_id" in c.data_separators


def test_forbidden_separator_names_filtered():
    c = Comparison("name", "proptype", "data_id", "basis")
    assert "name" not in c.data_separators
    assert "proptype" not in c.data_separators
    assert "data_id" not in c.data_separators
    assert "basis" in c.data_separators


def test_structure_property():
    c = Comparison("basis", "method")
    assert c.structure == ("name", "basis", "method", "proptype", "data_id")


# ---------------------------------------------------------------------------
# add_molecule / add
# ---------------------------------------------------------------------------

def test_add_molecule_inserts_at_correct_path():
    c = Comparison()
    mol = _mol("water", "bench", "cc-pvdz", "HF", -76.0)
    c.add_molecule(mol)
    assert "water" in c
    assert "cc-pvdz" in c["water"]
    assert "HF" in c["water"]["cc-pvdz"]
    assert "energy" in c["water"]["cc-pvdz"]["HF"]


def test_add_molecule_value_is_datapoint():
    c = Comparison()
    mol = _mol("water", "bench", "cc-pvdz", "HF", -76.0)
    c.add_molecule(mol)
    val = c["water"]["cc-pvdz"]["HF"]["energy"]["bench"]
    assert isinstance(val, Datapoint)
    assert val.value == pytest.approx(-76.0)


def test_add_molecule_list():
    c = Comparison()
    ml = MoleculeList([
        _mol("water", "bench", "cc-pvdz", "HF", -76.0),
        _mol("benzene", "bench", "cc-pvdz", "HF", -10.0),
    ])
    c.add(ml)
    assert "water" in c and "benzene" in c


def test_add_single_molecule_via_add():
    c = Comparison()
    mol = _mol("water", "bench", "cc-pvdz", "HF", -76.0)
    c.add(mol)
    assert "water" in c


def test_add_molecule_duplicate_data_id_warns(recwarn):
    c = Comparison()
    mol = _mol("water", "bench", "cc-pvdz", "HF", -76.0)
    c.add_molecule(mol)
    c.add_molecule(mol)  # same data_id → warning
    # No exception, just a warning
    assert "water" in c


def test_add_molecule_assigned_transition_id_overrides():
    mol = Molecule(
        "mol", "bench", {},
        {"gs": {
            "basis": "cc-pvdz",
            "method": "HF",
            "data": {"energy": Datapoint(-1.0, "au")},
            "transition_id": "old_tid",
            "assigned_transition_id": "new_tid",
        }}
    )
    c = Comparison("basis", "method", "transition_id")
    c.add_molecule(mol)
    # The "transition_id" separator picks up "new_tid" not "old_tid"
    assert "new_tid" in c["mol"]["cc-pvdz"]["HF"]


def test_add_molecule_none_separator_skipped():
    # State is missing a separator key → silently skipped
    mol = Molecule(
        "mol", "bench", {},
        {"gs": {"method": "HF",
                "data": {"energy": Datapoint(-1.0, "au")}}}
        # no "basis" key → separator value is None
    )
    c = Comparison()  # separators: basis, method
    c.add_molecule(mol)
    assert "mol" not in c


# ---------------------------------------------------------------------------
# walk_by_key / walk_values
# ---------------------------------------------------------------------------

def test_walk_by_key_finds_energy(simple_comparison):
    results = list(simple_comparison.walk_by_key("energy"))
    assert len(results) > 0
    for keys, val in results:
        assert "energy" in keys


def test_walk_values_yields_datapoints(simple_comparison):
    for _, val in simple_comparison.walk_values():
        assert isinstance(val, Datapoint)


def test_walk_values_count(simple_comparison):
    # simple_comparison has 2 molecules, each with 1 energy value
    values = list(simple_comparison.walk_values())
    assert len(values) == 2


# ---------------------------------------------------------------------------
# _import_value
# ---------------------------------------------------------------------------

def test_import_value_scalar():
    c = Comparison()
    assert c._import_value(1) == 1
    assert c._import_value(1.5) == 1.5


def test_import_value_datapoint():
    c = Comparison()
    dp = Datapoint(1.0, "au")
    result = c._import_value(dp)
    assert result is dp


def test_import_value_dict_with_value_unit():
    c = Comparison()
    result = c._import_value({"value": 2.0, "unit": "eV"})
    assert isinstance(result, Datapoint)
    assert result.value == 2.0
    assert result.unit == "eV"


def test_import_value_list_to_ndarray():
    c = Comparison()
    result = c._import_value([1.0, 2.0, 3.0])
    assert isinstance(result, numpy.ndarray)


def test_import_value_malformed_dict_logs_error(caplog):
    # A dict that isn't exactly {"value", "unit"} must be reported, not
    # silently swallowed into None.
    c = Comparison()
    with caplog.at_level("ERROR", logger="molbench"):
        result = c._import_value({"value": 2.0, "unit": "eV", "note": "extra"})
    assert result is None
    assert any("Could not interpret" in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# add_molecule — invalid input type
# ---------------------------------------------------------------------------

def test_add_molecule_wrong_type_logs_error_and_returns(caplog):
    c = Comparison()
    with caplog.at_level("ERROR", logger="molbench"):
        c.add_molecule({"not": "a molecule"})  # must not raise
    assert c == {}
    assert any("Can't add data of type" in rec.message for rec in caplog.records)
