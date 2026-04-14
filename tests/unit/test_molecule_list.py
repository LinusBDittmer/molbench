import pytest
from molbench.molecule import Molecule, MoleculeList, Datapoint


def _make_mol(name, data_id="bench", charge=0, basis="cc-pvdz", method="HF",
              energy=-1.0, unit="au"):
    return Molecule(
        name=name, data_id=data_id,
        system_data={"xyz": "H 0 0 0", "charge": charge, "multiplicity": 1},
        state_data={
            "gs": {
                "basis": basis,
                "method": method,
                "data": {"energy": Datapoint(energy, unit)},
            }
        },
    )


@pytest.fixture
def mixed_list():
    ml = MoleculeList()
    ml.extend([
        _make_mol("water", charge=0, basis="cc-pvdz", method="HF", energy=-76.0),
        _make_mol("benzene", charge=0, basis="cc-pvtz", method="HF", energy=-10.0),
        _make_mol("methane", charge=1, basis="cc-pvdz", method="MP2", energy=-5.0),
    ])
    return ml


# ---------------------------------------------------------------------------
# filter / remove
# ---------------------------------------------------------------------------

def test_filter_by_name(mixed_list):
    result = mixed_list.filter("name", "water")
    assert len(result) == 1
    assert result[0].name == "water"


def test_filter_by_data_id():
    ml = MoleculeList([
        _make_mol("a", data_id="bench1"),
        _make_mol("b", data_id="bench2"),
    ])
    result = ml.filter("data_id", "bench1")
    assert len(result) == 1
    assert result[0].name == "a"


def test_filter_by_system_key(mixed_list):
    result = mixed_list.filter("charge", 0)
    assert all(m.system_data["charge"] == 0 for m in result)
    assert len(result) == 2


def test_filter_by_state_key_basis(mixed_list):
    result = mixed_list.filter("basis", "cc-pvdz")
    names = {m.name for m in result}
    assert "benzene" not in names  # benzene has cc-pvtz
    assert "water" in names


def test_filter_removes_non_matching_states():
    mol = Molecule(
        "mol", "bench", {"charge": 0},
        {
            "s1": {"basis": "cc-pvdz", "method": "HF",
                   "data": {"energy": Datapoint(-1.0, "au")}},
            "s2": {"basis": "cc-pvtz", "method": "HF",
                   "data": {"energy": Datapoint(-2.0, "au")}},
        }
    )
    ml = MoleculeList([mol])
    result = ml.filter("basis", "cc-pvdz")
    assert len(result) == 1
    assert "s2" not in result[0].state_data
    assert "s1" in result[0].state_data


def test_filter_drops_molecule_when_all_states_removed(mixed_list):
    result = mixed_list.filter("method", "MP2")
    assert len(result) == 1
    assert result[0].name == "methane"


def test_remove_is_inverse_of_filter(mixed_list):
    kept = mixed_list.filter("name", "water")
    removed = mixed_list.remove("name", "water")
    assert len(kept) + len(removed) == len(mixed_list)
    assert all(m.name != "water" for m in removed)


def test_filter_returns_molecule_list(mixed_list):
    result = mixed_list.filter("name", "water")
    assert isinstance(result, MoleculeList)


# ---------------------------------------------------------------------------
# filter_by_range
# ---------------------------------------------------------------------------

def test_filter_by_range_both(mixed_list):
    result = mixed_list.filter_by_range("charge", min=0, max=0)
    assert all(m.system_data["charge"] == 0 for m in result)


def test_filter_by_range_min_only(mixed_list):
    result = mixed_list.filter_by_range("charge", min=1)
    assert all(m.system_data["charge"] >= 1 for m in result)


def test_filter_by_range_max_only(mixed_list):
    result = mixed_list.filter_by_range("charge", max=0)
    assert all(m.system_data["charge"] <= 0 for m in result)


def test_filter_by_range_none_none(mixed_list):
    result = mixed_list.filter_by_range("charge")
    assert len(result) == len(mixed_list)


# ---------------------------------------------------------------------------
# filter_by_vec_norm
# ---------------------------------------------------------------------------

def _mol_with_vec_norm(name, vec):
    return Molecule(
        name, "bench",
        {"vec_norm": vec},
        {"gs": {"basis": "b", "method": "m",
                "data": {"e": Datapoint(0.0, "au")}}}
    )


def test_filter_by_vec_norm_passes_in_range():
    ml = MoleculeList([_mol_with_vec_norm("a", [0.5, 0.5])])
    result = ml.filter_by_vec_norm("vec_norm", min=[0.0, 0.0], max=[1.0, 1.0])
    assert len(result) == 1


def test_filter_by_vec_norm_excludes_out_of_range():
    ml = MoleculeList([_mol_with_vec_norm("a", [1.5, 0.5])])
    result = ml.filter_by_vec_norm("vec_norm", min=[0.0, 0.0], max=[1.0, 1.0])
    assert len(result) == 0


def test_filter_by_vec_norm_scalar_promoted():
    ml = MoleculeList([_mol_with_vec_norm("a", [0.5])])
    result = ml.filter_by_vec_norm("vec_norm", min=0.0, max=1.0)
    assert len(result) == 1


def test_filter_by_vec_norm_dict_value():
    mol = Molecule(
        "m", "b", {"vec_norm": {"x": 0.5, "y": 0.3}},
        {"gs": {"basis": "b", "method": "m",
                "data": {"e": Datapoint(0.0, "au")}}}
    )
    result = MoleculeList([mol]).filter_by_vec_norm(
        "vec_norm", min=[0.0, 0.0], max=[1.0, 1.0]
    )
    assert len(result) == 1


def test_filter_by_vec_norm_none_none():
    ml = MoleculeList([_mol_with_vec_norm("a", [5.0])])
    result = ml.filter_by_vec_norm("vec_norm")
    assert len(result) == 1


# ---------------------------------------------------------------------------
# apply_stochiometry
# ---------------------------------------------------------------------------

def test_apply_stochiometry_basic():
    mol_a = _make_mol("molA", energy=-10.0)
    mol_b = _make_mol("molB", energy=-5.0)
    ml = MoleculeList([mol_a, mol_b])
    stochiometry = {
        "combined": {
            "molecules": ["molA", "molB"],
            "factors": [1.0, -1.0],
        }
    }
    result = ml.apply_stochiometry(stochiometry)
    assert len(result) == 1
    assert result[0].name == "combined"


def test_apply_stochiometry_missing_molecule_logged(capfd):
    mol_a = _make_mol("molA")
    ml = MoleculeList([mol_a])
    stochiometry = {
        "combined": {
            "molecules": ["molA", "nonexistent"],
            "factors": [1.0, -1.0],
        }
    }
    result = ml.apply_stochiometry(stochiometry)
    # Missing molecule → error logged, entry skipped
    assert len(result) == 0


def test_apply_stochiometry_system_data_merged():
    mol_a = _make_mol("molA", energy=-10.0)
    mol_b = _make_mol("molB", energy=-5.0)
    ml = MoleculeList([mol_a, mol_b])
    stochiometry = {
        "combined": {
            "molecules": ["molA", "molB"],
            "factors": [1.0, 1.0],
        }
    }
    result = ml.apply_stochiometry(stochiometry)
    # system_data keys become _list variants
    assert "xyz_list" in result[0].system_data
    assert len(result[0].system_data["xyz_list"]) == 2
