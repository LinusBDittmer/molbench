import json
import pytest
from pathlib import Path
from molbench.external_parser import ExternalParser
from molbench.molecule import MoleculeList


# ---------------------------------------------------------------------------
# Minimal mock parsers
# ---------------------------------------------------------------------------

def mock_parser_1param(filepath):
    return (
        Path(filepath).stem,
        {"xyz": "H 0 0 0", "charge": 0, "multiplicity": 1},
        {"gs": {"basis": "cc-pvdz", "method": "HF",
                "data": {"energy": {"value": -0.5, "unit": "au"}}}},
    )


def mock_parser_2param(filepath, name):
    return (
        name or Path(filepath).stem,
        {"xyz": "H 0 0 0", "charge": 0, "multiplicity": 1},
        {"gs": {"basis": "cc-pvdz", "method": "HF",
                "data": {"energy": {"value": -0.5, "unit": "au"}}}},
    )


def bad_parser_3param(filepath, name, extra):
    return (Path(filepath).stem, {}, {})


def bad_parser_wrong_return(filepath):
    return ("only_name", {})  # 2-tuple, not 3


# ---------------------------------------------------------------------------
# _fetch_all_outfiles
# ---------------------------------------------------------------------------

def test_fetch_all_outfiles_flat(tmp_path):
    (tmp_path / "a.out").write_text("output")
    (tmp_path / "b.out").write_text("output")
    ep = ExternalParser()
    found = ep._fetch_all_outfiles(str(tmp_path), ".out")
    assert len(found) == 2


def test_fetch_all_outfiles_recursive(tmp_path):
    sub = tmp_path / "sub"
    sub.mkdir()
    (tmp_path / "a.out").write_text("")
    (sub / "b.out").write_text("")
    ep = ExternalParser()
    found = ep._fetch_all_outfiles(str(tmp_path), ".out")
    assert len(found) == 2


def test_fetch_all_outfiles_ignores_wrong_suffix(tmp_path):
    (tmp_path / "a.txt").write_text("")
    (tmp_path / "b.out").write_text("")
    ep = ExternalParser()
    found = ep._fetch_all_outfiles(str(tmp_path), ".out")
    assert len(found) == 1
    assert found[0].endswith(".out")


def test_fetch_all_outfiles_empty_dir(tmp_path):
    ep = ExternalParser()
    found = ep._fetch_all_outfiles(str(tmp_path), ".out")
    assert found == []


# ---------------------------------------------------------------------------
# _assignmentfile_from_outfile
# ---------------------------------------------------------------------------

def test_assignmentfile_exists(tmp_path):
    out = tmp_path / "mol.out"
    ass = tmp_path / "mol.ass"
    out.write_text("")
    ass.write_text("")
    ep = ExternalParser()
    result = ep._assignmentfile_from_outfile(str(out), ".ass")
    assert result == str(ass)


def test_assignmentfile_absent(tmp_path):
    out = tmp_path / "mol.out"
    out.write_text("")
    ep = ExternalParser()
    result = ep._assignmentfile_from_outfile(str(out), ".ass")
    assert result is None


# ---------------------------------------------------------------------------
# load — basic
# ---------------------------------------------------------------------------

def test_load_single_file(tmp_path):
    (tmp_path / "molA.out").write_text("fake")
    ml = ExternalParser().load(str(tmp_path), parser=mock_parser_1param)
    assert isinstance(ml, MoleculeList)
    assert len(ml) == 1


def test_load_multiple_files(tmp_path):
    for name in ("a", "b", "c"):
        (tmp_path / f"{name}.out").write_text("fake")
    ml = ExternalParser().load(str(tmp_path), parser=mock_parser_1param)
    assert len(ml) == 3


def test_load_returns_molecules_with_correct_names(tmp_path):
    (tmp_path / "molA.out").write_text("")
    ml = ExternalParser().load(str(tmp_path), parser=mock_parser_1param)
    assert ml[0].name == "molA"


def test_load_1param_parser(tmp_path):
    (tmp_path / "mol.out").write_text("")
    ml = ExternalParser().load(str(tmp_path), parser=mock_parser_1param)
    assert len(ml) == 1


def test_load_2param_parser(tmp_path):
    (tmp_path / "mol.out").write_text("")
    ml = ExternalParser().load(str(tmp_path), parser=mock_parser_2param)
    assert len(ml) == 1


def test_load_wrong_param_count_exits(tmp_path):
    (tmp_path / "mol.out").write_text("")
    with pytest.raises(SystemExit):
        ExternalParser().load(str(tmp_path), parser=bad_parser_3param)


def test_load_wrong_return_length_exits(tmp_path):
    (tmp_path / "mol.out").write_text("")
    with pytest.raises(SystemExit):
        ExternalParser().load(str(tmp_path), parser=bad_parser_wrong_return)


def test_load_parser_none_exits_cleanly(tmp_path):
    # Must hit the intended log.critical() path, not a raw TypeError from
    # inspect.signature(None).
    (tmp_path / "mol.out").write_text("")
    with pytest.raises(SystemExit):
        ExternalParser().load(str(tmp_path), parser=None)


def test_load_parser_returns_none_exits_cleanly(tmp_path):
    # A parser that forgets its return statement returns None, which has no
    # len(); must hit the intended log.critical() path, not a raw TypeError.
    (tmp_path / "mol.out").write_text("")

    def forgetful_parser(filepath):
        pass

    with pytest.raises(SystemExit):
        ExternalParser().load(str(tmp_path), parser=forgetful_parser)


# ---------------------------------------------------------------------------
# load — with assignment files
# ---------------------------------------------------------------------------

def test_load_with_assignment_file(tmp_path):
    # Molecule states have external computation ids as transition_id.
    # Assignment file maps: ref_id ==> external_id  (returns {external: ref}).
    # add_assignments looks up transition_id (external) as key in that dict.
    (tmp_path / "mol.out").write_text("")

    def parser_with_tid(filepath):
        return (
            "mol",
            {},
            {"s1": {"basis": "cc-pvdz", "method": "HF",
                    "data": {"energy": {"value": -1.0, "unit": "au"}},
                    "transition_id": "state_001"},   # external computation id
             "s2": {"basis": "cc-pvdz", "method": "HF",
                    "data": {"energy": {"value": -2.0, "unit": "au"}},
                    "transition_id": "state_002"}},  # external computation id
        )

    # File format: ref_id ==> external_id; state_002 left as null → removed
    ass_content = "s0->s1 ==> state_001\ns0->s2 ==> null\n"
    (tmp_path / "mol.ass").write_text(ass_content)

    ml = ExternalParser().load(str(tmp_path), parser=parser_with_tid)
    assert len(ml) == 1
    # s2 had null assignment → state removed by add_assignments
    assert "s2" not in ml[0].state_data
    assert "s1" in ml[0].state_data


def test_load_custom_out_suffix(tmp_path):
    (tmp_path / "mol.log").write_text("")
    ml = ExternalParser().load(str(tmp_path), parser=mock_parser_1param, out_suffix=".log")
    assert len(ml) == 1


def test_load_custom_assignment_suffix(tmp_path):
    (tmp_path / "mol.out").write_text("")

    def parser_with_tid(filepath):
        return ("mol", {},
                {"s1": {"basis": "b", "method": "m",
                        "data": {"energy": {"value": -1.0, "unit": "au"}},
                        "transition_id": "ref_001"}})  # external computation id

    # File format: ref_id ==> external_id; parse returns {external: ref}
    # add_assignments looks up transition_id (ref_001) → assigned = "s0->s1"
    (tmp_path / "mol.asgn").write_text("s0->s1 ==> ref_001\n")
    ml = ExternalParser().load(str(tmp_path), parser=parser_with_tid,
                               assignment_suffix=".asgn")
    assert ml[0].state_data["s1"].get("assigned_transition_id") == "s0->s1"
