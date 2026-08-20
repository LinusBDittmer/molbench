import pytest

from molbench.benchmark_parser import BenchmarkParser, JSONBenchmarkParser
from molbench.molecule import Molecule, MoleculeList


@pytest.fixture
def parser():
    return JSONBenchmarkParser()


# ---------------------------------------------------------------------------
# Built-in benchmarks
# ---------------------------------------------------------------------------


def test_load_ascdb(parser):
    ml = parser.load("ascdb")
    assert isinstance(ml, MoleculeList)
    assert len(ml) > 0


def test_load_questdb(parser):
    ml = parser.load("questdb")
    assert len(ml) > 0


def test_load_jacquemin18(parser):
    ml = parser.load("jacquemin18")
    assert len(ml) > 0


def test_ascdb_all_molecules_instances(parser):
    ml = parser.load("ascdb")
    assert all(isinstance(m, Molecule) for m in ml)


def test_questdb_has_transition_ids(parser):
    ml = parser.load("questdb")
    found = False
    for mol in ml:
        for state in mol.state_data.values():
            if "transition_id" in state:
                found = True
                break
    assert found


# ---------------------------------------------------------------------------
# Custom file
# ---------------------------------------------------------------------------


def test_load_from_file_path(parser, minimal_benchmark_file):
    ml = parser.load(minimal_benchmark_file, benchmark_id="test_bench")
    assert len(ml) == 2


def test_benchmark_id_propagated(parser, minimal_benchmark_file):
    ml = parser.load(minimal_benchmark_file, benchmark_id="my_id")
    assert all(m.data_id == "my_id" for m in ml)


def test_benchmark_id_defaults_to_benchmark_arg(parser, minimal_benchmark_file):
    ml = parser.load(minimal_benchmark_file)
    assert all(m.data_id == minimal_benchmark_file for m in ml)


def test_molecule_names_are_keys(parser, minimal_benchmark_file):
    ml = parser.load(minimal_benchmark_file)
    names = {m.name for m in ml}
    assert "molA" in names
    assert "molB" in names


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------


def test_load_nonexistent_exits(parser):
    with pytest.raises(SystemExit):
        parser.load("/nonexistent/path/bench.json")


def test_load_malformed_json_exits(parser, tmp_path):
    f = tmp_path / "bad.json"
    f.write_text("{not valid json")
    with pytest.raises(SystemExit):
        parser.load(str(f))


def test_load_empty_content_exits(parser, tmp_path):
    f = tmp_path / "empty.json"
    f.write_text("{}")
    with pytest.raises(SystemExit):
        parser.load(str(f))


def test_use_local_benchmark_skips_premade(parser, minimal_benchmark_file):
    ml = parser.load(minimal_benchmark_file, use_local_benchmark=True)
    assert len(ml) > 0


# ---------------------------------------------------------------------------
# premade_benchmarks discovery
# ---------------------------------------------------------------------------


def test_premade_benchmarks_discovered(parser):
    BenchmarkParser._collect_premade_benchmarks()
    assert "ascdb" in parser.premade_benchmarks
    assert "questdb" in parser.premade_benchmarks


def test_premade_benchmarks_are_json_paths(parser):
    for path in parser.premade_benchmarks.values():
        assert path.endswith(".json")


# ---------------------------------------------------------------------------
# multi-geometry benchmark
# ---------------------------------------------------------------------------


def test_load_xyz_list_benchmark(parser, minimal_benchmark_list_file):
    ml = parser.load(minimal_benchmark_list_file)
    assert len(ml) == 1
    assert "xyz_list" in ml[0].system_data
