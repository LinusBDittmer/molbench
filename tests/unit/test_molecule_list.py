import pytest

from molbench.molecule import Datapoint, Molecule, MoleculeList


def _make_mol(
    name,
    data_id="bench",
    charge=0,
    basis="cc-pvdz",
    method="HF",
    energy=-1.0,
    unit="au",
):
    return Molecule(
        name=name,
        data_id=data_id,
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
    ml.extend(
        [
            _make_mol("water", charge=0, basis="cc-pvdz", method="HF", energy=-76.0),
            _make_mol("benzene", charge=0, basis="cc-pvtz", method="HF", energy=-10.0),
            _make_mol("methane", charge=1, basis="cc-pvdz", method="MP2", energy=-5.0),
        ]
    )
    return ml


# ---------------------------------------------------------------------------
# filter / remove
# ---------------------------------------------------------------------------


def test_filter_by_name(mixed_list):
    result = mixed_list.filter("name", "water")
    assert len(result) == 1
    assert result[0].name == "water"


def test_filter_by_data_id():
    ml = MoleculeList(
        [
            _make_mol("a", data_id="bench1"),
            _make_mol("b", data_id="bench2"),
        ]
    )
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
        "mol",
        "bench",
        {"charge": 0},
        {
            "s1": {
                "basis": "cc-pvdz",
                "method": "HF",
                "data": {"energy": Datapoint(-1.0, "au")},
            },
            "s2": {
                "basis": "cc-pvtz",
                "method": "HF",
                "data": {"energy": Datapoint(-2.0, "au")},
            },
        },
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


def test_filter_by_name_and_data_id_delegate(mixed_list):
    assert [m.name for m in mixed_list.filter("name", "water")] == [
        m.name for m in mixed_list.filter_names("water")
    ]
    assert [m.name for m in mixed_list.remove("name", "water")] == [
        m.name for m in mixed_list.remove_names("water")
    ]
    assert [m.name for m in mixed_list.filter("data_id", "bench")] == [
        m.name for m in mixed_list.filter_data_ids("bench")
    ]
    assert [m.name for m in mixed_list.remove("data_id", "bench")] == [
        m.name for m in mixed_list.remove_data_ids("bench")
    ]


# ---------------------------------------------------------------------------
# filter_properties / remove_properties
# ---------------------------------------------------------------------------


@pytest.fixture
def property_list():
    return MoleculeList(
        [
            Molecule(
                name="multi",
                data_id="bench",
                system_data={"charge": 0},
                state_data={
                    "s1": {
                        "basis": "cc-pvdz",
                        "method": "adc2",
                        "data": {
                            "excitation_energy": Datapoint(1.0, "eV"),
                            "oscillator_strength": Datapoint(0.1, "au"),
                        },
                    },
                    "s2": {
                        "basis": "cc-pvdz",
                        "method": "adc2",
                        "data": {"oscillator_strength": Datapoint(0.2, "au")},
                    },
                },
            ),
            Molecule(
                name="osc_only",
                data_id="bench",
                system_data={"charge": 0},
                state_data={
                    "s1": {
                        "basis": "cc-pvdz",
                        "method": "adc2",
                        "data": {"oscillator_strength": Datapoint(0.3, "au")},
                    },
                },
            ),
        ]
    )


def test_filter_properties(property_list):
    # ensure that osc_only is dropped and multi only keeps the s1 energy
    result = property_list.filter_properties("excitation_energy")
    assert len(result) == 1
    assert result[0].name == "multi"
    # the state without excitation energy is dropped
    assert list(result[0].state_data) == ["s1"]
    assert list(result[0].state_data["s1"]["data"]) == ["excitation_energy"]
    # the remaining state keeps its other entries
    assert result[0].state_data["s1"]["method"] == "adc2"
    # Ensure that nothing is dropped
    result = property_list.filter_properties("excitation_energy", "oscillator_strength")
    assert len(result) == 2
    assert set(result[0].state_data) == {"s1", "s2"}
    assert set(result[0].state_data["s1"]["data"]) == {
        "excitation_energy",
        "oscillator_strength",
    }
    assert list(result[1].state_data["s1"]["data"]) == ["oscillator_strength"]
    # ensure that multi looses its energy
    removed = property_list.remove_properties("excitation_energy")
    assert len(removed) == 2
    assert removed[0].name == "multi"
    assert set(removed[0].state_data) == {"s1", "s2"}
    assert removed[1].name == "osc_only"
    assert list(removed[1].state_data) == ["s1"]
    assert all(
        ptype == "oscillator_strength"
        for m in removed
        for _, data in m.state_data.items()
        for ptype in data["data"]
    )
    # ensure that everything is dropped
    result = property_list.remove_properties("excitation_energy", "oscillator_strength")
    assert len(result) == 0
    # ensure that a state without data is dropped
    no_data = MoleculeList()
    no_data.append(
        Molecule(
            name="no_data",
            data_id="test",
            system_data={},
            state_data={"s1": {"method": "adc2", "basis": "cc-pvdz"}},
        )
    )
    res = no_data.filter_properties("bla")
    assert len(res) == 0
    res = no_data.remove_properties("bla")
    assert len(res) == 0


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
        name,
        "bench",
        {"vec_norm": vec},
        {"gs": {"basis": "b", "method": "m", "data": {"e": Datapoint(0.0, "au")}}},
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
        "m",
        "b",
        {"vec_norm": {"x": 0.5, "y": 0.3}},
        {"gs": {"basis": "b", "method": "m", "data": {"e": Datapoint(0.0, "au")}}},
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


def test_apply_stochiometry_mismatched_lengths_exits():
    mol_a = _make_mol("molA")
    mol_b = _make_mol("molB")
    ml = MoleculeList([mol_a, mol_b])
    stochiometry = {
        "combined": {
            "molecules": ["molA", "molB"],
            "factors": [1.0],  # only one factor for two molecules
        }
    }
    with pytest.raises(SystemExit):
        ml.apply_stochiometry(stochiometry)


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
