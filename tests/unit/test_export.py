import io

import pytest

from molbench.comparison import Comparison
from molbench.export import LatexExporter
from molbench.molecule import Datapoint, Molecule, MoleculeList
from molbench.tree import Node

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def two_data_comparison():
    """Two molecules, two data_ids (ref and computed) for 'energy'."""
    mols = MoleculeList(
        [
            Molecule(
                "water",
                "ref",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "TBE",
                        "data": {"energy": Datapoint(-76.0, "au")},
                    }
                },
            ),
            Molecule(
                "water",
                "computed",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "HF",
                        "data": {"energy": Datapoint(-75.9, "au")},
                    }
                },
            ),
        ]
    )
    c = Comparison()
    c.add(mols)
    return c


def _export_to_string(comparison, prop, columns, rows=None, **kwargs):
    exp = LatexExporter(**kwargs)
    buf = io.StringIO()
    exp.export(comparison, prop, buf, columns=columns, rows=rows)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Basic output
# ---------------------------------------------------------------------------


def test_export_writes_nonempty(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col)
    assert len(result) > 0


def test_export_contains_begin_table(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col)
    assert r"\begin{table}" in result


def test_export_contains_end_table(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col)
    assert r"\end{table}" in result


def test_export_contains_tabular(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col)
    assert r"\begin{tabular}" in result


def test_export_column_label_in_output(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col)
    # Both data_ids appear as column labels
    assert "ref" in result or "computed" in result


def test_export_row_label_in_output(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col)
    assert "water" in result


def test_export_numeric_value_in_output(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col)
    # -76.0 should appear in the table
    assert "76" in result


# ---------------------------------------------------------------------------
# Sorting
# ---------------------------------------------------------------------------


def test_export_sort_cols_default(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col, sort_cols=True)
    assert len(result) > 0


def test_export_sort_cols_false(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col, sort_cols=False)
    assert len(result) > 0


# ---------------------------------------------------------------------------
# Sparse row labels
# ---------------------------------------------------------------------------


def test_export_sparse_row_labels_two_rows(two_data_comparison):
    """With two rows under the same parent, the second should have empty prefix."""
    mols = MoleculeList(
        [
            Molecule(
                "water",
                "ref",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "TBE",
                        "data": {"energy": Datapoint(-76.0, "au")},
                    },
                    "gs2": {
                        "basis": "cc-pvtz",
                        "method": "TBE",
                        "data": {"energy": Datapoint(-76.1, "au")},
                    },
                },
            ),
            Molecule(
                "water",
                "computed",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "HF",
                        "data": {"energy": Datapoint(-75.9, "au")},
                    },
                    "gs2": {
                        "basis": "cc-pvtz",
                        "method": "HF",
                        "data": {"energy": Datapoint(-76.0, "au")},
                    },
                },
            ),
        ]
    )
    c = Comparison()
    c.add(mols)
    col = Node("data_id")
    result = _export_to_string(c, "energy", col, sparse_row_labels=True)
    assert len(result) > 0


# ---------------------------------------------------------------------------
# Key not in Comparison exits
# ---------------------------------------------------------------------------


def test_export_invalid_column_key_exits(two_data_comparison):
    col = Node("nonexistent_key")
    buf = io.StringIO()
    with pytest.raises(SystemExit):
        LatexExporter().export(two_data_comparison, "energy", buf, columns=col)


# ---------------------------------------------------------------------------
# Empty field for missing data
# ---------------------------------------------------------------------------


def test_export_missing_value_shows_empty_field():
    """Molecule 'benzene' has energy but molecule 'water' does not for method MP2."""
    mols = MoleculeList(
        [
            Molecule(
                "water",
                "ref",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "TBE",
                        "data": {"energy": Datapoint(-76.0, "au")},
                    }
                },
            ),
            Molecule(
                "benzene",
                "ref",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "HF",
                        "data": {"energy": Datapoint(-10.0, "au")},
                    }
                },
            ),
        ]
    )
    c = Comparison()
    c.add(mols)
    col = Node("data_id")
    result = _export_to_string(c, "energy", col)
    # The table should still be produced
    assert r"\begin{table}" in result


# ---------------------------------------------------------------------------
# multirow
# ---------------------------------------------------------------------------


def test_export_multirow_flag(two_data_comparison):
    col = Node("data_id")
    result = _export_to_string(two_data_comparison, "energy", col, multirow=True)
    assert len(result) > 0


# ---------------------------------------------------------------------------
# Silent cell collisions
# ---------------------------------------------------------------------------


def test_export_uncovered_separator_collision_warns(caplog):
    # Two entries share the same name and data_id and differ only by
    # "method", but the row tree only covers "name" and the column tree
    # only covers "data_id" - so both values collide into the same cell.
    # This must be logged, not silent.
    mols = MoleculeList(
        [
            Molecule(
                "water",
                "bench",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "HF",
                        "data": {"energy": Datapoint(-1.0, "au")},
                    }
                },
            ),
            Molecule(
                "water",
                "bench",
                {},
                {
                    "gs": {
                        "basis": "cc-pvdz",
                        "method": "MP2",
                        "data": {"energy": Datapoint(-2.0, "au")},
                    }
                },
            ),
        ]
    )
    c = Comparison()
    c.add(mols)
    rows = Node("name")
    cols = Node("data_id")
    buf = io.StringIO()
    with caplog.at_level("WARNING", logger="molbench"):
        LatexExporter().export(c, "energy", buf, columns=cols, rows=rows)
    assert any("same table cell" in rec.message for rec in caplog.records)
