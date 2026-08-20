import pytest
from molbench.formatting import StdFormatter, LatexFormatter
from molbench.export import TableExporter


class TestStdFormatter:
    def setup_method(self):
        self.fmt = StdFormatter(n_decimals=3)

    def test_format_int(self):
        assert self.fmt.format_datapoint(1) == "1"

    def test_format_float_rounding(self):
        assert self.fmt.format_datapoint(3.14159) == "3.142"

    def test_format_float_exact(self):
        assert self.fmt.format_datapoint(1.0) == "1.0"

    def test_format_list(self):
        result = self.fmt.format_datapoint([1.0, 2.0])
        assert result == "1.0, 2.0"

    def test_format_tuple(self):
        result = self.fmt.format_datapoint((1.5, 2.5))
        assert "1.5" in result and "2.5" in result

    def test_format_none_default_empty(self):
        assert self.fmt.format_datapoint(None) == ""

    def test_format_none_custom_empty_field(self):
        fmt = StdFormatter(empty_field="N/A")
        assert fmt.format_datapoint(None) == "N/A"

    def test_format_non_iterable_fallback(self):
        # Non-numeric, non-iterable objects fall through to str()
        class Custom:
            def __str__(self): return "custom_val"
        assert self.fmt.format_datapoint(Custom()) == "custom_val"

    def test_custom_delimiter(self):
        fmt = StdFormatter(value_delimiter=" | ")
        result = fmt.format_datapoint([1.0, 2.0])
        assert result == "1.0 | 2.0"

    def test_format_bare_string_no_recursion(self):
        # Strings are iterable, so without a dedicated check they get
        # treated as a sequence of themselves and recurse forever.
        assert self.fmt.format_datapoint("closed-shell") == "closed-shell"

    def test_format_single_char_string_no_recursion(self):
        assert self.fmt.format_datapoint("x") == "x"


def test_table_exporter_rejects_incompatible_formatter():
    # StdFormatter doesn't implement the table-structure methods
    # TableExporter needs; must fail fast at construction with a clear
    # message instead of crashing mid-export with an AttributeError.
    with pytest.raises(SystemExit):
        TableExporter(StdFormatter())


class TestLatexFormatter:
    def setup_method(self):
        self.fmt = LatexFormatter(n_decimals=3)

    def test_init_table_contains_begin(self):
        result = self.fmt.init_table(n_additional_cols=1, n_columns=2)
        assert r"\begin{table}" in result
        assert r"\begin{tabular}" in result

    def test_init_table_alignment_string(self):
        result = self.fmt.init_table(n_additional_cols=2, n_columns=3)
        # 2 additional (l) + | + 3 columns (c)
        assert "ll|ccc" in result

    def test_finalize_table(self):
        result = self.fmt.finalize_table()
        assert r"\end{tabular}" in result
        assert r"\end{table}" in result

    def test_table_header_ends_with_hline(self):
        result = self.fmt.table_header([["A", "B"]])
        assert r"\hline" in result

    def test_multicolumn(self):
        result = self.fmt.multicolumn(3, "label")
        assert r"\multicolumn{3}{c}{label}" == result

    def test_multirow(self):
        result = self.fmt.multirow(2, "label")
        assert r"\multirow{2}{*}{label}" == result

    def test_join_labels_default_delimiter(self):
        result = self.fmt.join_labels(iter(["a", "b", "c"]))
        assert result == "a/b/c"

    def test_join_columns(self):
        result = self.fmt.join_columns(["a", "b", "c"])
        assert result == "a & b & c"

    def test_join_rows(self):
        result = self.fmt.join_rows(["row1", "row2"])
        assert "row1" in result and "row2" in result
        assert r"\\" in result

    def test_table_content_not_empty(self):
        result = self.fmt.table_content([["val1", "val2"]])
        assert "val1" in result and "val2" in result
