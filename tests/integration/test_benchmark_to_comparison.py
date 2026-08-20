"""Integration: load benchmark → filter → add to Comparison → walk."""
import pytest
from molbench.benchmark_parser import JSONBenchmarkParser
from molbench.molecule import Datapoint
from molbench.comparison import Comparison


@pytest.fixture(scope="module")
def ascdb():
    return JSONBenchmarkParser().load("ascdb", benchmark_id="TBE")


@pytest.fixture(scope="module")
def questdb():
    return JSONBenchmarkParser().load("questdb", benchmark_id="TBE")


# ---------------------------------------------------------------------------
# ascdb
# ---------------------------------------------------------------------------

def test_ascdb_loads_and_populates_comparison(ascdb):
    c = Comparison()
    c.add(ascdb)
    values = list(c.walk_values())
    assert len(values) > 0
    for _, v in values:
        assert isinstance(v, Datapoint)


def test_ascdb_all_molecules_in_comparison(ascdb):
    c = Comparison()
    c.add(ascdb)
    for mol in ascdb:
        assert mol.name in c


def test_ascdb_filter_then_add(ascdb):
    filtered = ascdb.filter("method", "TBE")
    c = Comparison()
    c.add(filtered)
    assert len(c) > 0
    # Only TBE molecules
    for name, basis_dict in c.items():
        for basis, method_dict in basis_dict.items():
            assert "TBE" in method_dict or list(method_dict.keys())


# ---------------------------------------------------------------------------
# questdb
# ---------------------------------------------------------------------------

def test_questdb_loads(questdb):
    assert len(questdb) > 0


def test_questdb_transition_id_separator(questdb):
    c = Comparison("basis", "method", "transition_id")
    c.add(questdb)
    values = list(c.walk_values())
    assert len(values) > 0


def test_questdb_populates_comparison(questdb):
    c = Comparison()
    c.add(questdb)
    assert len(c) > 0


def test_questdb_energy_values_are_datapoints(questdb):
    c = Comparison()
    c.add(questdb)
    energy_entries = list(c.walk_by_key("excitation_energy"))
    assert len(energy_entries) > 0
    for _, val_dict in energy_entries:
        for _, v in val_dict.items():
            assert isinstance(v, Datapoint)


# ---------------------------------------------------------------------------
# stochiometry pipeline (ascdb has stochiometry entries)
# ---------------------------------------------------------------------------

def test_ascdb_stochiometry_entries_exist(ascdb):
    """Some ascdb entries use stochiometry (multi-geometry)."""
    has_xyz_list = any(
        "xyz_list" in mol.system_data for mol in ascdb
    )
    assert has_xyz_list
