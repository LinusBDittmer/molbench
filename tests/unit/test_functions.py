import pytest
from molbench.functions import (
    substitute_template,
    _substitute_single_template,
    default_name_template,
    walk_dict_by_key,
    walk_dict_values,
)


class TestSubstituteTemplate:
    def test_simple_substitution(self):
        result = substitute_template("charge=[[charge]]", {"charge": 0})
        assert result == ("charge=0",)

    def test_multiple_keys(self):
        result = substitute_template("[[name]] [[basis]]", {"name": "mol", "basis": "cc-pvdz"})
        assert result == ("mol cc-pvdz",)

    def test_no_placeholders(self):
        result = substitute_template("no placeholders", {"x": 1})
        assert result == ("no placeholders",)

    def test_list_expansion(self):
        subvals = {"charge_list": [0, 1], "multiplicity_list": [1, 2]}
        tmpl = "charge=[[charge]] mult=[[multiplicity]]"
        result = substitute_template(tmpl, subvals)
        assert len(result) == 2
        assert "charge=0 mult=1" in result
        assert "charge=1 mult=2" in result

    def test_list_expansion_with_common_key(self):
        subvals = {"charge_list": [0, 1], "basis": "cc-pvdz"}
        tmpl = "[[charge]] [[basis]]"
        result = substitute_template(tmpl, subvals)
        assert len(result) == 2
        assert all("cc-pvdz" in r for r in result)

    def test_index_syntax(self):
        result = substitute_template("[[xyz->1]]", {"xyz": ["a", "b", "c"]})
        assert result == ("b",)

    def test_index_syntax_zero(self):
        result = substitute_template("[[xyz->0]]", {"xyz": ["first", "second"]})
        assert result == ("first",)

    def test_mismatched_list_lengths_exits(self):
        subvals = {"a_list": [1, 2], "b_list": [1, 2, 3]}
        with pytest.raises(SystemExit):
            substitute_template("[[a]] [[b]]", subvals)

    def test_missing_key_exits(self):
        with pytest.raises(SystemExit):
            substitute_template("[[missing]]", {})

    def test_non_numeric_index_exits(self):
        with pytest.raises(SystemExit):
            substitute_template("[[xyz->abc]]", {"xyz": [1, 2]})

    def test_index_on_scalar_exits(self):
        with pytest.raises(SystemExit):
            substitute_template("[[val->0]]", {"val": 42})

    def test_returns_tuple(self):
        result = substitute_template("[[x]]", {"x": 1})
        assert isinstance(result, tuple)


class TestDefaultNameTemplate:
    def test_basic(self):
        result = default_name_template(("basis",), ".in")
        assert "[[name]]" in result
        assert "[[method]]" in result
        assert "[[basis]]" in result
        assert result.endswith(".in")

    def test_no_extra_keys(self):
        # name and method are always included; listing them again shouldn't duplicate
        result = default_name_template(("name", "method"), ".in")
        assert result.count("[[name]]") == 1
        assert result.count("[[method]]") == 1

    def test_multiple_expansion_keys(self):
        result = default_name_template(("basis", "charge"), ".ass")
        assert "[[basis]]" in result
        assert "[[charge]]" in result
        assert result.endswith(".ass")


class TestWalkDictByKey:
    def test_shallow(self):
        d = {"a": 1, "b": 2}
        results = list(walk_dict_by_key(d, "a"))
        assert len(results) == 1
        assert results[0][1] == 1

    def test_deep(self):
        d = {"x": {"y": {"target": 42}}}
        results = list(walk_dict_by_key(d, "target"))
        assert len(results) == 1
        keys, val = results[0]
        assert val == 42
        assert "target" in keys

    def test_full_path_returned(self):
        d = {"x": {"target": 99}}
        results = list(walk_dict_by_key(d, "target"))
        keys, val = results[0]
        assert keys == ("x", "target")

    def test_absent_key(self):
        d = {"a": {"b": 1}}
        results = list(walk_dict_by_key(d, "missing"))
        assert results == []

    def test_multiple_occurrences(self):
        d = {"a": {"key": 1}, "b": {"key": 2}}
        results = list(walk_dict_by_key(d, "key"))
        assert len(results) == 2
        values = [v for _, v in results]
        assert 1 in values and 2 in values


class TestWalkDictValues:
    def test_shallow(self):
        d = {"a": 1, "b": 2}
        results = list(walk_dict_values(d))
        values = [v for _, v in results]
        assert 1 in values and 2 in values

    def test_deep(self):
        d = {"x": {"y": 42, "z": 99}}
        results = list(walk_dict_values(d))
        values = [v for _, v in results]
        assert 42 in values and 99 in values

    def test_nested_dict_not_yielded_as_value(self):
        d = {"outer": {"inner": 5}}
        results = list(walk_dict_values(d))
        values = [v for _, v in results]
        assert isinstance(results[0][1], int)
        assert len(results) == 1

    def test_key_path_returned(self):
        d = {"a": {"b": 7}}
        results = list(walk_dict_values(d))
        keys, val = results[0]
        assert keys == ("a", "b")
        assert val == 7
