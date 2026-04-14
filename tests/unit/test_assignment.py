import pytest
from molbench.assignment import new_assignment_file, parse_assignment_file


# ---------------------------------------------------------------------------
# new_assignment_file
# ---------------------------------------------------------------------------

def test_new_assignment_file_header():
    content = new_assignment_file([])
    assert content.startswith("# ref_state_id ==>")


def test_new_assignment_file_entries():
    content = new_assignment_file(["s1", "s2"])
    assert "s1 ==> null" in content
    assert "s2 ==> null" in content


def test_new_assignment_file_empty():
    content = new_assignment_file([])
    lines = [l for l in content.splitlines() if not l.startswith("#")]
    assert lines == []


def test_new_assignment_file_format():
    content = new_assignment_file(["t1"])
    lines = content.splitlines()
    entry_lines = [l for l in lines if not l.startswith("#")]
    assert entry_lines == ["t1 ==> null"]


# ---------------------------------------------------------------------------
# parse_assignment_file
# ---------------------------------------------------------------------------

def _write_ass(tmp_path, content, filename="test.ass"):
    f = tmp_path / filename
    f.write_text(content)
    return str(f)


def test_parse_basic(tmp_path):
    # File format: ref_id ==> ext_id; returns {ext_id: ref_id}
    f = _write_ass(tmp_path, "ref1 ==> ext1\nref2 ==> ext2\n")
    result = parse_assignment_file(f)
    assert result == {"ext1": "ref1", "ext2": "ref2"}


def test_parse_skips_comment_lines(tmp_path):
    f = _write_ass(tmp_path, "# this is a comment\nref ==> ext\n")
    result = parse_assignment_file(f)
    assert "# this is a comment" not in result
    assert result == {"ext": "ref"}


def test_parse_skips_empty_lines(tmp_path):
    f = _write_ass(tmp_path, "\nref ==> ext\n\n")
    result = parse_assignment_file(f)
    assert result == {"ext": "ref"}


def test_parse_skips_null_token(tmp_path):
    f = _write_ass(tmp_path, "ref1 ==> null\n")
    result = parse_assignment_file(f)
    assert result == {}


def test_parse_null_ref_skipped(tmp_path):
    f = _write_ass(tmp_path, "null ==> ext1\n")
    result = parse_assignment_file(f)
    assert result == {}


def test_parse_inline_comment_stripped(tmp_path):
    # File format: ref_id ==> ext_id; inline comment after ext_id is stripped
    f = _write_ass(tmp_path, "ref ==> ext # inline comment\n")
    result = parse_assignment_file(f)
    assert result == {"ext": "ref"}


def test_parse_duplicate_external_warns(tmp_path, recwarn):
    # Two different refs pointing to the same external id — last wins
    f = _write_ass(tmp_path, "ref1 ==> ext\nref2 ==> ext\n")
    result = parse_assignment_file(f)
    # Last value wins
    assert result["ext"] == "ref2"


def test_parse_bad_separator_exits(tmp_path):
    f = _write_ass(tmp_path, "ext ref\n")
    with pytest.raises(SystemExit):
        parse_assignment_file(f)


def test_parse_import_external_callback(tmp_path):
    # ref ==> ext; import_external transforms the RIGHT side (ext)
    f = _write_ass(tmp_path, "ref ==> ext\n")
    result = parse_assignment_file(f, import_external=lambda x: x.upper())
    assert "EXT" in result


def test_parse_import_ref_callback(tmp_path):
    # ref ==> ext; import_ref transforms the LEFT side (ref)
    f = _write_ass(tmp_path, "ref ==> ext\n")
    result = parse_assignment_file(f, import_ref=lambda x: x.upper())
    assert result.get("ext") == "REF"


def test_roundtrip(tmp_path):
    """new_assignment_file content that is partially filled can be read back."""
    ids = ["s0->s1", "s0->t1"]
    raw = new_assignment_file(ids)
    # Replace null with actual external (computation) state ids
    raw = raw.replace("s0->s1 ==> null", "s0->s1 ==> state_001")
    raw = raw.replace("s0->t1 ==> null", "s0->t1 ==> state_002")
    f = tmp_path / "rt.ass"
    f.write_text(raw)
    result = parse_assignment_file(str(f))
    # parse returns {external: ref} = {state_001: s0->s1, ...}
    assert result == {"state_001": "s0->s1", "state_002": "s0->t1"}
