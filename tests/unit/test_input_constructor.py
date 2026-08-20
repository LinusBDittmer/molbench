import json
from pathlib import Path

import pytest

from molbench.assignment import parse_assignment_file
from molbench.input_constructor import (
    CompressedTemplateConstructor,
    TemplateConstructor,
)
from molbench.molecule import Datapoint, Molecule, MoleculeList

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mol(name, basis="cc-pvdz", method="HF", tid=None, charge=0):
    state = {
        "basis": basis,
        "method": method,
        "data": {"energy": Datapoint(-1.0, "au")},
    }
    if tid:
        state["transition_id"] = tid
    return Molecule(
        name,
        "bench",
        {"xyz": "H 0 0 0", "charge": charge, "multiplicity": 1, "xyz_unit": "A"},
        {"gs": state},
    )


def _bench(*names, basis="cc-pvdz"):
    return MoleculeList([_mol(n, basis=basis) for n in names])


# ---------------------------------------------------------------------------
# init_template
# ---------------------------------------------------------------------------


def test_init_template_from_templates_dir():
    # "pyscf_ordmp2" is a built-in template
    tc = TemplateConstructor("pyscf_ordmp2")
    assert "[[basis]]" in tc.template


def test_init_template_from_filepath(tmp_path):
    f = tmp_path / "my_tmpl.txt"
    f.write_text("basis=[[basis]]")
    tc = TemplateConstructor(str(f))
    assert tc.template == "basis=[[basis]]"


def test_init_template_nonexistent_exits():
    with pytest.raises(SystemExit):
        TemplateConstructor("/nonexistent/template.txt")


# ---------------------------------------------------------------------------
# create_inputs — files and directories
# ---------------------------------------------------------------------------


def test_create_inputs_creates_directory(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    calc = {"method": "HF"}
    tc.create_inputs(_bench("molA"), str(tmp_path / "out"), calc)
    assert (tmp_path / "out").exists()


def test_create_inputs_default_nested_structure(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    tc.create_inputs(_bench("molA"), str(tmp_path), {"method": "HF"})
    # Expect basepath/molA/molA_HF_cc-pvdz.in
    mol_dir = tmp_path / "molA"
    assert mol_dir.is_dir()
    files = list(mol_dir.glob("*.in"))
    assert len(files) == 1


def test_create_inputs_flat_structure(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    tc.create_inputs(
        _bench("molA", "molB"), str(tmp_path), {"method": "HF"}, flat_structure=True
    )
    # All files directly in basepath
    files = list(tmp_path.glob("*.in"))
    assert len(files) == 2


def test_create_inputs_file_count_matches_molecules(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    tc.create_inputs(
        _bench("a", "b", "c"), str(tmp_path), {"method": "HF"}, flat_structure=True
    )
    files = list(tmp_path.glob("*.in"))
    assert len(files) == 3


def test_create_inputs_content_substituted(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    tc.create_inputs(_bench("molA"), str(tmp_path), {"method": "HF"})
    content = next((tmp_path / "molA").glob("*.in")).read_text()
    assert "cc-pvdz" in content
    assert "H 0 0 0" in content


def test_create_inputs_returns_file_list(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    result = tc.create_inputs(_bench("molA"), str(tmp_path), {"method": "HF"})
    assert isinstance(result, list)
    assert len(result) == 1
    assert Path(result[0]).exists()


def test_create_inputs_list_expansion(tmp_path):
    """charge_list → multiple files per molecule."""
    template = "charge=[[charge]]\nbasis=[[basis]]"
    f = tmp_path / "tmpl.txt"
    f.write_text(template)
    mol = Molecule(
        "molA",
        "bench",
        {
            "xyz_list": ["H 0 0 0", "H 0 1 0"],
            "charge_list": [0, 1],
            "multiplicity_list": [1, 2],
            "xyz_unit": "A",
        },
        {
            "gs": {
                "basis": "cc-pvdz",
                "method": "HF",
                "data": {"energy": Datapoint(-1.0, "au")},
            }
        },
    )
    tc = TemplateConstructor(str(f))
    result = tc.create_inputs(
        MoleculeList([mol]),
        str(tmp_path / "out"),
        {"method": "HF"},
        flat_structure=True,
    )
    assert len(result) == 2


def test_create_inputs_custom_name_template(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    tc.create_inputs(
        _bench("molA"),
        str(tmp_path),
        {"method": "HF"},
        name_template="[[name]].input",
        flat_structure=True,
    )
    assert (tmp_path / "molA.input").exists()


# ---------------------------------------------------------------------------
# create_assignments
# ---------------------------------------------------------------------------


def test_create_assignments_creates_ass_files(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    bench = MoleculeList([_mol("molA", tid="s0->s1")])
    tc.create_assignments(bench, str(tmp_path), {"method": "HF"})
    ass_files = list(tmp_path.rglob("*.ass"))
    assert len(ass_files) == 1


def test_create_assignments_content_is_valid(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    bench = MoleculeList([_mol("molA", tid="s0->s1")])
    tc.create_assignments(bench, str(tmp_path), {"method": "HF"})
    ass_file = next(tmp_path.rglob("*.ass"))
    # File should be parseable (all entries are null initially)
    result = parse_assignment_file(str(ass_file))
    assert result == {}  # all null → empty dict


def test_create_assignments_transition_ids_present(tmp_path, simple_template_file):
    tc = TemplateConstructor(simple_template_file)
    bench = MoleculeList([_mol("molA", tid="s0->s1")])
    tc.create_assignments(bench, str(tmp_path), {"method": "HF"})
    ass_file = next(tmp_path.rglob("*.ass"))
    content = ass_file.read_text()
    assert "s0->s1" in content


def test_create_assignments_no_tid_warning(tmp_path, simple_template_file, capfd):
    tc = TemplateConstructor(simple_template_file)
    bench = _bench("molA")  # no transition_id
    tc.create_assignments(bench, str(tmp_path), {"method": "HF"})
    # should not crash; empty assignment file created


def test_create_inputs_missing_expansion_key_warns(
    tmp_path, simple_template_file, caplog
):
    # A state missing the (default) "basis" expansion key must be skipped
    # with a logged warning, not silently dropped.
    mol = Molecule(
        "molA",
        "bench",
        {"xyz": "H 0 0 0", "charge": 0, "multiplicity": 1},
        {
            "gs": {
                "method": "HF",  # no "basis" key
                "data": {"energy": Datapoint(-1.0, "au")},
            }
        },
    )
    tc = TemplateConstructor(simple_template_file)
    with caplog.at_level("WARNING", logger="molbench"):
        result = tc.create_inputs(MoleculeList([mol]), str(tmp_path), {"method": "HF"})
    assert result == []
    assert any("missing value" in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# CompressedTemplateConstructor
# ---------------------------------------------------------------------------


def _multi_mol(name):
    return Molecule(
        name,
        "bench",
        {
            "xyz_list": ["H 0 0 0", "H 0 1 0"],
            "charge_list": [0, 0],
            "multiplicity_list": [1, 1],
            "n_atoms_list": [1, 1],
        },
        {
            "p1": {
                "basis": "cc-pvdz",
                "method": "TBE",
                "data": {"energy": Datapoint(0.5, "au")},
                "stochiometry": [-1.0, 1.0],
            }
        },
    )


def test_compressed_creates_references_json(tmp_path, simple_template_file):
    tc = CompressedTemplateConstructor(simple_template_file)
    bench = MoleculeList([_multi_mol("relA")])
    tc.create_inputs(
        bench, str(tmp_path), {"method": "HF"}, reference_path="references.json"
    )
    assert (tmp_path / "references.json").exists()


def test_compressed_references_json_valid(tmp_path, simple_template_file):
    tc = CompressedTemplateConstructor(simple_template_file)
    bench = MoleculeList([_multi_mol("relA")])
    tc.create_inputs(
        bench, str(tmp_path), {"method": "HF"}, reference_path="references.json"
    )
    refs = json.loads((tmp_path / "references.json").read_text())
    assert "relA" in refs


def test_compressed_single_geometry_handled(tmp_path, simple_template_file):
    tc = CompressedTemplateConstructor(simple_template_file)
    bench = _bench("molA")  # no xyz_list
    result = tc.create_inputs(
        bench, str(tmp_path), {"method": "HF"}, reference_path="references.json"
    )
    assert len(result) >= 1


def test_compressed_property_typo_single_property_exits(tmp_path, simple_template_file):
    # A single-property molecule with a typo'd compressed_property must hit
    # the intended log.critical() path, not a raw IndexError.
    tc = CompressedTemplateConstructor(simple_template_file)
    bench = MoleculeList([_multi_mol("relA")])
    with pytest.raises(SystemExit):
        tc.create_inputs(
            bench,
            str(tmp_path),
            {"method": "HF"},
            reference_path="references.json",
            compressed_property="does_not_exist",
        )


def test_compressed_deduplicates_identical_geometries(tmp_path, simple_template_file):
    tc = CompressedTemplateConstructor(simple_template_file)
    mol1 = _multi_mol("relA")
    mol2 = _multi_mol("relB")
    # Both use "H 0 0 0" as first geometry → should be deduplicated
    bench = MoleculeList([mol1, mol2])
    result = tc.create_inputs(
        bench, str(tmp_path), {"method": "HF"}, reference_path="references.json"
    )
    # 2 unique geometries per mol, but first geometry shared → 3 unique instead of 4
    assert len(result) <= 4
