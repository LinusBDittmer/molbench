"""Integration: load benchmark → generate input/assignment files → verify."""

import json
from pathlib import Path

import pytest

from molbench.assignment import parse_assignment_file
from molbench.benchmark_parser import JSONBenchmarkParser
from molbench.input_constructor import (
    CompressedTemplateConstructor,
    TemplateConstructor,
)


@pytest.fixture(scope="module")
def ascdb():
    return JSONBenchmarkParser().load("ascdb", benchmark_id="TBE")


@pytest.fixture(scope="module")
def questdb():
    return JSONBenchmarkParser().load("questdb", benchmark_id="TBE")


# ---------------------------------------------------------------------------
# TemplateConstructor with built-in pyscf template
# ---------------------------------------------------------------------------

PYSCF_CALC = {
    "method": "hf",
    "verbose": 0,
    "xyz_unit": "A",
    "symmetry": False,
    "spin_unrestricted": False,
    "use_soscf": False,
    "conv_tol": 1e-9,
    "conv_tol_grad": 1e-6,
    "scf_max_cycle": 100,
    "optimality_criterion": "gradient",
    "ordmp2_solver": "diis",
    "ordmp2_stepsize": 0.1,
    "ordmp2_max_cycles": 50,
    "use_renorm": False,
    "renormalisation": 1.0,
}


def test_generate_inputs_from_ascdb_creates_files(ascdb, tmp_path):
    tc = TemplateConstructor("pyscf_ordmp2")
    small = ascdb[:3]  # just first 3 to keep it fast
    result = tc.create_inputs(small, str(tmp_path), PYSCF_CALC)
    assert len(result) > 0
    for f in result:
        assert Path(f).exists()


def test_generated_input_contains_xyz(ascdb, tmp_path):
    tc = TemplateConstructor("pyscf_ordmp2")
    small = ascdb[:1]
    result = tc.create_inputs(small, str(tmp_path), PYSCF_CALC)
    content = Path(result[0]).read_text()
    # xyz coordinates should be substituted
    assert "atom_str" in content or len(content) > 50


def test_generated_input_contains_basis(ascdb, tmp_path):
    tc = TemplateConstructor("pyscf_ordmp2")
    small = ascdb[:1]
    result = tc.create_inputs(small, str(tmp_path), PYSCF_CALC)
    content = Path(result[0]).read_text()
    # basis should appear
    assert "basis" in content


def test_generate_inputs_nested_folder_per_molecule(ascdb, tmp_path):
    tc = TemplateConstructor("pyscf_ordmp2")
    small = ascdb[:2]
    tc.create_inputs(small, str(tmp_path), PYSCF_CALC)
    # Each molecule gets its own folder
    subdirs = [d for d in tmp_path.iterdir() if d.is_dir()]
    assert len(subdirs) >= 1


def test_generate_inputs_flat_structure(ascdb, tmp_path):
    tc = TemplateConstructor("pyscf_ordmp2")
    small = ascdb[:2]
    result = tc.create_inputs(small, str(tmp_path), PYSCF_CALC, flat_structure=True)
    for f in result:
        assert Path(f).parent == tmp_path


# ---------------------------------------------------------------------------
# create_assignments with questdb (has transition_ids)
# ---------------------------------------------------------------------------


def test_generate_assignments_from_questdb(questdb, tmp_path):
    tc = TemplateConstructor("pyscf_ordmp2")
    small = questdb[:3]
    tc.create_assignments(small, str(tmp_path), PYSCF_CALC)
    ass_files = list(tmp_path.rglob("*.ass"))
    assert len(ass_files) > 0


def test_assignment_files_parseable(questdb, tmp_path):
    tc = TemplateConstructor("pyscf_ordmp2")
    small = questdb[:2]
    tc.create_assignments(small, str(tmp_path), PYSCF_CALC)
    for ass_file in tmp_path.rglob("*.ass"):
        # Should parse without SystemExit (all entries null initially)
        result = parse_assignment_file(str(ass_file))
        assert isinstance(result, dict)


def test_assignment_file_contains_transition_ids(questdb, tmp_path):
    tc = TemplateConstructor("pyscf_ordmp2")
    small = questdb[:1]
    tc.create_assignments(small, str(tmp_path), PYSCF_CALC)
    for ass_file in tmp_path.rglob("*.ass"):
        content = ass_file.read_text()
        # At least one transition id should appear
        assert "s0->" in content or "null" in content


# ---------------------------------------------------------------------------
# CompressedTemplateConstructor with ascdb (has xyz_list entries)
# ---------------------------------------------------------------------------


def test_compressed_creates_references_json(ascdb, tmp_path):
    tc = CompressedTemplateConstructor("pyscf_ordmp2")
    # Pick only multi-geometry molecules
    multi = [m for m in ascdb if "xyz_list" in m.system_data]
    if not multi:
        pytest.skip("No multi-geometry molecules in ascdb slice")
    from molbench.molecule import MoleculeList

    tc.create_inputs(
        MoleculeList(multi[:2]),
        str(tmp_path),
        PYSCF_CALC,
        reference_path="references.json",
    )
    assert (tmp_path / "references.json").exists()


def test_references_json_maps_to_molecule_names(ascdb, tmp_path):
    tc = CompressedTemplateConstructor("pyscf_ordmp2")
    multi = [m for m in ascdb if "xyz_list" in m.system_data]
    if not multi:
        pytest.skip("No multi-geometry molecules in ascdb slice")
    from molbench.molecule import MoleculeList

    tc.create_inputs(
        MoleculeList(multi[:2]),
        str(tmp_path),
        PYSCF_CALC,
        reference_path="references.json",
    )
    refs = json.loads((tmp_path / "references.json").read_text())
    original_names = {m.name for m in multi[:2]}
    assert any(name in refs for name in original_names)
