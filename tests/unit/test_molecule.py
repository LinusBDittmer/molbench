import pytest
from molbench.molecule import Molecule, Datapoint


# ---------------------------------------------------------------------------
# from_benchmark
# ---------------------------------------------------------------------------

def _entry_single(name=None):
    entry = {
        "charge": 0,
        "multiplicity": 1,
        "n_atoms": 1,
        "xyz": ["H 0 0 0"],
        "properties": {
            "gs": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": {"value": -0.5, "unit": "au"}},
            }
        },
    }
    if name:
        entry["name"] = name
    return entry


def _entry_multi():
    return {
        "charge_list": [0, 1],
        "multiplicity_list": [1, 2],
        "n_atoms_list": [1, 1],
        "xyz_list": [["H 0 0 0"], ["H 0 1 0"]],
        "properties": {},
    }


def test_from_benchmark_single_geo_xyz_joined():
    mol = Molecule.from_benchmark(_entry_single(), "bench", "molA")
    assert mol.system_data["xyz"] == "H 0 0 0"


def test_from_benchmark_xyz_already_string():
    entry = _entry_single()
    entry["xyz"] = "H 0 0 0"  # already a string
    mol = Molecule.from_benchmark(entry, "bench", "molA")
    assert mol.system_data["xyz"] == "H 0 0 0"


def test_from_benchmark_multi_geo_xyz_joined():
    mol = Molecule.from_benchmark(_entry_multi(), "bench", "rel")
    xyz_list = mol.system_data["xyz_list"]
    assert isinstance(xyz_list, list)
    assert all(isinstance(x, str) for x in xyz_list)


def test_from_benchmark_name_from_entry():
    entry = _entry_single(name="from_entry")
    mol = Molecule.from_benchmark(entry, "bench", "fallback")
    assert mol.name == "from_entry"
    assert "name" not in mol.system_data  # removed from system_data


def test_from_benchmark_name_fallback_to_molname():
    mol = Molecule.from_benchmark(_entry_single(), "bench", "fallback")
    assert mol.name == "fallback"


def test_from_benchmark_no_name_exits():
    with pytest.raises(SystemExit):
        Molecule.from_benchmark(_entry_single(), "bench", None)


def test_from_benchmark_data_id():
    mol = Molecule.from_benchmark(_entry_single(), "my_bench", "mol")
    assert mol.data_id == "my_bench"


def test_from_benchmark_no_properties():
    entry = {
        "charge": 0, "multiplicity": 1,
        "xyz": "H 0 0 0",
    }
    mol = Molecule.from_benchmark(entry, "bench", "mol")
    assert mol.state_data == {}


def test_from_benchmark_system_data_excludes_properties():
    mol = Molecule.from_benchmark(_entry_single(), "bench", "mol")
    assert "properties" not in mol.system_data


# ---------------------------------------------------------------------------
# from_external
# ---------------------------------------------------------------------------

def _valid_state_data():
    return {
        "gs": {
            "basis": "cc-pvdz",
            "method": "HF",
            "data": {"energy": {"value": -0.5, "unit": "au"}},
        }
    }


def test_from_external_valid():
    mol = Molecule.from_external(
        {"xyz": "H 0 0 0"}, _valid_state_data(), "/path/file.out", "molA"
    )
    assert mol.name == "molA"
    assert isinstance(mol.state_data["gs"]["data"]["energy"], Datapoint)


def test_from_external_both_empty_exits():
    with pytest.raises(SystemExit):
        Molecule.from_external({}, {}, "f.out", "mol")


def test_from_external_missing_method_exits():
    sd = {"gs": {"basis": "cc-pvdz", "data": {"energy": {"value": 0, "unit": "au"}}}}
    with pytest.raises(SystemExit):
        Molecule.from_external({}, sd, "f.out", "mol")


def test_from_external_missing_basis_exits():
    sd = {"gs": {"method": "HF", "data": {"energy": {"value": 0, "unit": "au"}}}}
    with pytest.raises(SystemExit):
        Molecule.from_external({}, sd, "f.out", "mol")


def test_from_external_missing_data_exits():
    sd = {"gs": {"method": "HF", "basis": "cc-pvdz"}}
    with pytest.raises(SystemExit):
        Molecule.from_external({}, sd, "f.out", "mol")


def test_from_external_malformed_datapoint_exits():
    # data dict entry must have exactly "value" and "unit"
    sd = {"gs": {"method": "HF", "basis": "cc-pvdz",
                 "data": {"energy": {"value": 0}}}}  # missing "unit"
    with pytest.raises(SystemExit):
        Molecule.from_external({}, sd, "f.out", "mol")


def test_from_external_non_dict_state_exits():
    sd = {"gs": "not a dict"}
    with pytest.raises(SystemExit):
        Molecule.from_external({}, sd, "f.out", "mol")


def test_from_external_non_string_state_key_exits():
    sd = {123: {"method": "HF", "basis": "cc-pvdz",
                "data": {"energy": {"value": 0, "unit": "au"}}}}
    with pytest.raises(SystemExit):
        Molecule.from_external({}, sd, "f.out", "mol")


def test_from_external_system_data_only():
    mol = Molecule.from_external({"xyz": "H 0 0 0"}, {}, "f.out", "mol")
    assert mol.system_data == {"xyz": "H 0 0 0"}
    assert mol.state_data == {}


# ---------------------------------------------------------------------------
# add_assignments
# ---------------------------------------------------------------------------

def _mol_with_transitions():
    return Molecule(
        "mol", "bench", {},
        {
            "s1": {"method": "HF", "basis": "cc-pvdz", "data": {},
                   "transition_id": "s0->s1"},
            "s2": {"method": "HF", "basis": "cc-pvdz", "data": {},
                   "transition_id": "s0->s2"},
        }
    )


def test_add_assignments_assigns_correctly():
    mol = _mol_with_transitions()
    mol.add_assignments({"s0->s1": "state_001", "s0->s2": "state_002"})
    assert mol.state_data["s1"]["assigned_transition_id"] == "state_001"
    assert mol.state_data["s2"]["assigned_transition_id"] == "state_002"


def test_add_assignments_removes_unassigned():
    mol = _mol_with_transitions()
    mol.add_assignments({"s0->s1": "state_001"})  # s2 not in assignments
    assert "s2" not in mol.state_data
    assert "s1" in mol.state_data


def test_add_assignments_no_transition_id_skipped():
    mol = Molecule(
        "mol", "bench", {},
        {"gs": {"method": "HF", "basis": "cc-pvdz", "data": {}}}  # no transition_id
    )
    mol.add_assignments({"some_key": "val"})
    # state without transition_id is unchanged and not removed
    assert "gs" in mol.state_data


def test_add_assignments_custom_keys():
    mol = Molecule(
        "mol", "bench", {},
        {"s1": {"method": "HF", "basis": "cc-pvdz", "data": {},
                "my_tid": "A"}}
    )
    mol.add_assignments({"A": "B"},
                        old_transition_id_key="my_tid",
                        new_transition_id_key="result_tid")
    assert mol.state_data["s1"]["result_tid"] == "B"
