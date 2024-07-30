from .comparison import Comparison
from .formatting import LatexFormatter, Formatter
from .tree import Node, DummyNode
from . import logger as log
from itertools import chain
import typing


class Exporter:
    """
    Parent class for a comparison between external data and a benchmark set.

    This class is used to perform comparisons between external data and a
    benchmark set. Any class that is supposed to perform such comparisons must
    inherit from this class.

    Methods
    -------
    compare(benchmark: dict, external_data: dict, properties: tuple) -> str
        Compare benchmark data with external data and return file contents of
        the comparison.

    Attributes
    ----------
    None
    """

    def export(self, *args, **kwargs) -> str:
        raise NotImplementedError("Export function has to be implemented on "
                                  "the child classes.")


class TableExporter(Exporter):
    def __init__(self, formatter: Formatter, sort_cols: bool = True,
                 sort_rows: bool = True, sparse_row_labels: bool = True,
                 multirow: bool = False):
        self.formatter = formatter
        self.sort_cols = sort_cols
        self.sort_rows = sort_rows
        self.sparse_row_labels = sparse_row_labels
        self.multirow = multirow

    def export(self, data: Comparison, property, outfile: typing.IO,
               columns, rows=None):
        # no row labels provided -> use all available keys
        if rows is None:
            rows = DummyNode()
            for key in data.structure:
                Node(key, rows)
        # TODO: add some tree constructors
        if not isinstance(columns, Node):
            raise NotImplementedError
        if not isinstance(rows, Node):
            raise NotImplementedError

        # prepare the data for the export by constructing row and column labels
        # for all data points
        prepared_data, col_label_tree, row_label_tree = (
            self._prepare_data(data, columns, rows, property)
        )
        # sort the column and row trees according to the labels
        if self.sort_cols:
            col_label_tree.sort(key=lambda n: n.value)
        if self.sort_rows:
            row_label_tree.sort(key=lambda n: n.value)
        # collect all unique column labels
        column_labels = []
        for leave in col_label_tree.walk_leaves():
            column_labels.append(
                tuple(node.value for node in reversed(leave.path_to_root()))
            )
        # get the labels for the additional columns from the row structure
        additional_col_labels = tuple(
            self.formatter.join_labels(n.value for n in generation)
            for generation in rows.traverse_generations()
        )
        # init the table
        preamble = self.formatter.init_table(len(additional_col_labels),
                                             len(column_labels))
        # Build the header of the table (column labels)
        header = self._prepare_table_header(col_label_tree,
                                            additional_col_labels)
        header = self.formatter.table_header(header)
        # Fill the body of the table
        content = self._prepare_content(prepared_data, row_label_tree,
                                        column_labels)
        content = self.formatter.table_content(content)
        # finalize the table (close environments etc.)
        finish = self.formatter.finalize_table()
        table = "\n".join((preamble, header, content, finish))
        outfile.write(table)

    def _prepare_data(self, comparison: Comparison, column_tree: Node,
                      row_tree: Node, prop):
        # Ensure that all nodes are present in comparison!
        data_structure = comparison.structure
        column_nodes = tuple(column_tree.traverse_generations())
        row_nodes = tuple(row_tree.traverse_generations())
        if any(n.value not in data_structure for n in
               chain.from_iterable(chain(column_nodes, row_nodes))):
            log.critical("A key is not available in the provided Comparison.", "Export")

        col_node_cache = {"root": DummyNode()}
        row_node_cache = {"root": DummyNode()}
        data: dict[tuple, dict[tuple, list]] = {}
        for keys, propdata in comparison.walk_by_key(prop):
            for data_id, value in propdata.items():
                column_l, row_l = self._build_row_column_label(
                    column_nodes, row_nodes, data_structure, keys + (data_id,)
                )
                self._add_to_label_tree(column_l, col_node_cache)
                self._add_to_label_tree(row_l, row_node_cache)
                if row_l not in data:
                    data[row_l] = {}
                if column_l not in data[row_l]:
                    data[row_l][column_l] = []
                data[row_l][column_l].append(value)
        col_label_tree = col_node_cache["root"]
        row_label_tree = row_node_cache["root"]
        return data, col_label_tree, row_label_tree

    def _build_row_column_label(self, columns: tuple, rows: tuple,
                                sep_names: tuple,
                                sep_vals: tuple) -> tuple[tuple, tuple]:
        # construct the row and column label for a given entry
        row_l, column_l = {}, {}
        for name, val in zip(sep_names, sep_vals):
            assigned = False
            for gen, cols in enumerate(columns):
                # also include the index of the node to maintain the label
                # order as given in the tree
                for i, node in enumerate(cols):
                    if name == node.value:
                        if gen not in column_l:
                            column_l[gen] = []
                        column_l[gen].append((i, node.to_string(val)))
                        assigned = True
                        break
                if assigned:
                    break
            if assigned:
                continue
            for gen, row in enumerate(rows):
                for i, node in enumerate(row):
                    if name == node.value:
                        if gen not in row_l:
                            row_l[gen] = []
                        row_l[gen].append((i, node.to_string(val)))
                        assigned = True
                        break
                if assigned:
                    break
        # sort the labels for each generation to maintain the order
        # as given in the input tree
        for gen, labels in column_l.items():
            column_l[gen] = self.formatter.join_labels(
                v for _, v in sorted(labels)
            )
        for gen, labels in row_l.items():
            row_l[gen] = self.formatter.join_labels(
                v for _, v in sorted(labels)
            )
        # finally sort the labels such that the upper entries of the tree
        # appear first
        return (
            tuple(v for _, v in sorted(column_l.items())),
            tuple(v for _, v in sorted(row_l.items()))
        )

    def _add_to_label_tree(self, label: tuple[str], node_cache: dict) -> None:
        # Helper function for building the label trees
        # constructs all nodes for a given entry and inserts them in the tree
        parent = node_cache["root"]
        for value in label:
            key = (value, parent)
            node = node_cache.get(key, None)
            if node is None:  # construct the node and store it
                node = Node(value, parent)
                node_cache[key] = node
            parent = node

    def _prepare_table_header(self, col_label_tree: Node,
                              additional_cols: tuple[str]) -> list:
        # Build the table header as nested list of strings
        rows = []
        prefix = tuple("" for _ in range(len(additional_cols)))
        nodes_by_generation = tuple(col_label_tree.traverse_generations())
        for i, generation in enumerate(nodes_by_generation):
            if i == len(nodes_by_generation) - 1:  # last row
                prefix = additional_cols
            row = [*prefix]
            for node in generation:
                if (width := node.width() - 1) > 1:
                    label = self.formatter.multicolumn(width, node.value)
                else:
                    label = node.value
                row.append(label)
            rows.append(row)
        return rows

    def _prepare_content(self, data: dict, row_label_tree: Node,
                         column_labels: tuple[str]) -> list:
        # Build the content of the table as nested list
        content: list[list[str]] = []
        prev_row_label = None
        for leave in row_label_tree.walk_leaves():
            row = []
            # write the row labels for the values in the current row
            tree_path = tuple(n for n in reversed(leave.path_to_root()))
            row_label = tuple(n.value for n in tree_path)
            if prev_row_label is None:
                write_label = [True for _ in range(len(row_label))]
            else:
                write_label = [label != prev for label, prev in
                               zip(row_label, prev_row_label)]
            for label, write, node in zip(row_label, write_label, tree_path):
                # possibly skip repeating row labels
                if self.sparse_row_labels and not write:
                    row.append("")
                    continue
                # we have to write the row label
                # -> either as multirow or simply as entry
                if self.multirow and (width := node.width() - 1) > 1:
                    row.append(self.formatter.multirow(width, label))
                else:
                    row.append(label)
            prev_row_label = row_label
            # add the values to the row
            for label in column_labels:
                value = data[row_label].get(label, None)
                row.append(self.formatter.format_datapoint(value))
            content.append(row)
        return content


class LatexExporter(TableExporter):
    def __init__(self, formatter: Formatter = None, sort_cols: bool = True,
                 sort_rows: bool = True, sparse_row_labels: bool = True,
                 multirow: bool = False):
        if formatter is None:
            formatter = LatexFormatter()
        super().__init__(formatter, sort_cols=sort_cols, sort_rows=sort_rows,
                         sparse_row_labels=sparse_row_labels,
                         multirow=multirow)
