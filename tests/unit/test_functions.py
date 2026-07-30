import pytest
from molbench.functions import (
    substitute_template,
    _substitute_single_template,
    default_name_template,
    determine_basis_cardinality,
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

    def test_list_key_with_scalar_value_not_expanded(self):
        # A "_list"-suffixed key whose value isn't actually a list/tuple must
        # not be treated as something to expand (previously crashed with a
        # raw IndexError since to_expand ended up empty).
        result = substitute_template("charge=[[charge_list]]",
                                     {"charge_list": 5})
        assert result == ("charge=5",)

    def test_out_of_range_index_exits(self):
        with pytest.raises(SystemExit):
            substitute_template("[[coords->5]]", {"coords": [1, 2, 3]})

    def test_empty_list_expansion_logs_warning_and_returns_empty(self, caplog):
        with caplog.at_level("WARNING", logger="molbench"):
            result = substitute_template("[[name]]", {"name_list": []})
        assert result == ()
        assert any("empty list" in rec.message for rec in caplog.records)

    def test_unclosed_placeholder_left_as_literal_and_warns(self, caplog):
        with caplog.at_level("ERROR", logger="molbench"):
            result = substitute_template("charge=[[chargex] end", {"charge": 0})
        assert result == ("charge=[[chargex] end",)
        assert any("unclosed placeholder" in rec.message for rec in caplog.records)

    def test_stray_closing_bracket_before_real_placeholder(self):
        # A stray "]]" earlier in the template must not be mistaken for the
        # closing bracket of the real "[[charge]]" placeholder that follows.
        result = substitute_template("stray ]] bracket [[charge]]", {"charge": 0})
        assert result == ("stray ]] bracket 0",)

    def test_bash_bracket_collision_exits_with_helpful_message(self, caplog):
        # bash's own "[[ ... ]]" test syntax collides with placeholder
        # delimiters and cannot be resolved automatically; the error message
        # should hint at the real cause instead of just naming a garbage key.
        tpl = "if [[ -f file.txt ]]; then echo ok; fi"
        with pytest.raises(SystemExit):
            substitute_template(tpl, {"charge": 0})


class TestDetermineBasisCardinality:
    def test_dunning_valid(self):
        assert determine_basis_cardinality("cc-pvtz") == 3

    def test_karlsruhe_valid(self):
        assert determine_basis_cardinality("def2-tzvp") == 3

    def test_malformed_dunning_returns_zero(self, caplog):
        with caplog.at_level("ERROR", logger="molbench"):
            result = determine_basis_cardinality("cc-p")
        assert result == 0
        assert any("could not be identified" in rec.message
                   for rec in caplog.records)

    def test_malformed_def2_returns_zero(self, caplog):
        with caplog.at_level("ERROR", logger="molbench"):
            result = determine_basis_cardinality("def2")
        assert result == 0
        assert any("could not be identified" in rec.message
                   for rec in caplog.records)


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
