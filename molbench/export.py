from .comparison import Comparison
from .formatting import StdFormatter, Formatter
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
    def _build_row_column_label(self, columns: tuple, rows: tuple,
                                sep_names: tuple,
                                sep_vals: tuple) -> tuple[tuple, tuple]:
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
        # TODO: custom separator
        column_l = {key: "/".join(v for _, v in sorted(val))
                    for key, val in column_l.items()}
        row_l = {key: "/".join(v for _, v in sorted(val))
                 for key, val in row_l.items()}
        return (
            tuple(v for _, v in sorted(column_l.items())),
            tuple(v for _, v in sorted(row_l.items()))
        )

    def _add_to_label_tree(self, label, node_cache) -> None:
        parent = node_cache["root"]
        for value in label:
            key = (value, parent)
            node = node_cache.get(key, None)
            if node is None:  # construct the node and store it
                node = Node(value, parent)
                node_cache[key] = node
            parent = node

    def prepare_export(self, comparison: Comparison, column_tree: Node,
                       row_tree: Node, prop):
        # Ensure that all nodes are present in comparison!
        data_structure = comparison.structure
        column_nodes = tuple(column_tree.traverse_generations())
        row_nodes = tuple(row_tree.traverse_generations())
        if any(n.value not in data_structure for n in
               chain.from_iterable(chain(column_nodes, row_nodes))):
            log.critical("A key is not available in the provided Comparison.")

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
        # sort the column label tree according to the values
        col_label_tree = col_node_cache["root"]
        col_label_tree.sort(key=lambda n: n.value)
        # build all unique column labels in the right order.
        column_labels = []
        for leave in col_label_tree.walk_leaves():
            path = leave.path_to_root()
            column_labels.append(tuple(node.value for node in reversed(path)))

        # also sort the row labels
        # TODO: if we don't use multirow but instead just don't repeat
        # row labels that are the same as in the previous row,
        # sorting the row labels is optional!
        row_label_tree = row_node_cache["root"]
        row_label_tree.sort(key=lambda n: n.value)

        # get the names for the additional columns that are required due to
        # the structure of the rows.
        # TODO: custom separator
        additional_col_labels = tuple("/".join(n.value for n in gen)
                                      for gen in row_nodes)
        return (
            data, column_labels, col_label_tree, additional_col_labels,
            row_label_tree
        )

    def _prepare_table_header(self, col_label_tree: Node,
                              additional_cols: tuple[str]) -> str:
        rows = []
        prefix = tuple("" for _ in range(len(additional_cols)))
        nodes_by_generation = tuple(col_label_tree.traverse_generations())
        for i, generation in enumerate(nodes_by_generation):
            if i == len(nodes_by_generation) - 1:  # last row
                prefix = additional_cols
            row = [*prefix]
            for node in generation:
                if (width := node.width() - 1) > 1:
                    label = self._multicolumn(width, node.value)
                else:
                    label = node.value
                row.append(label)
            rows.append(row)
        return rows

    def _build_row(self, row_data: dict, row_label: tuple,
                   column_labels: tuple,
                   prev_row_label: tuple = None) -> list:
        if prev_row_label is None:  # first row -> write all row labels
            row = [*row_label]
        else:
            # TODO: do we want to use multirow here?
            # only write row labels that differ from the previous label
            row = ["" if label == prev else label
                   for label, prev in zip(row_label, prev_row_label)]
        for label in column_labels:
            row.append(
                self.formatter.format_datapoint(row_data.get(label, None))
            )
        return row


class LatexExporter(TableExporter):
    def __init__(self, formatter: Formatter = None):
        # TODO: allow to modify the style of the latex table
        # also move formatter to init and use it also for things like
        # linebreaks, column speparators, etc.
        if formatter is None:
            formatter = StdFormatter()
        self.formatter = formatter

    def export(self, data: Comparison, prop, outfile: typing.IO,
               columns, rows) -> None:
        # if no row labels are provided, use all labels
        # -> will put everything in row_label that is not defined as column
        if rows is None:
            rows = DummyNode()
            for key in data.structure:
                Node(key, rows)
        # TODO: add some functions to create trees from certain data structures
        if not isinstance(columns, Node):
            raise NotImplementedError
        if not isinstance(rows, Node):
            raise NotImplementedError

        export_data = self.prepare_export(data, columns, rows, prop)
        prepared_data = export_data[0]
        column_labels = export_data[1]
        col_label_tree = export_data[2]
        additional_cols = export_data[3]
        row_label_tree = export_data[4]
        print(prepared_data)
        print(column_labels)
        print(col_label_tree)
        print(additional_cols)
        print(row_label_tree)
        print()

        # generate the header of the table
        header = self._table_header(col_label_tree, additional_cols)
        print(header)
        print()
        content = self._table_content(prepared_data, row_label_tree,
                                      column_labels)
        print(content)
        exit()
        outfile.write(header)
        outfile.write("\n")
        outfile.write(content)

    def _table_header(self, col_label_tree: Node,
                      additional_cols: tuple[str]) -> str:
        rows = self._prepare_table_header(col_label_tree, additional_cols)
        return (
            "\\\\\n".join(" & ".join(r for r in row) for row in rows) +
            r"\\ \hline"
        )

    def _multicolumn(self, width: int, value: str) -> str:
        return r"\multicolumn{" + str(width) + r"}{c}{" + value + "}"

    def _table_content(self, data: dict, row_label_tree: Node,
                       column_labels: tuple):
        content: list[str] = []
        prev_row_label = None
        for leave_node in row_label_tree.walk_leaves():
            path = leave_node.path_to_root()
            row_label = tuple(n.value for n in reversed(path))
            row = self._build_row(data[row_label], row_label, column_labels,
                                  prev_row_label)
            # TODO: custom separator
            content.append(" & ".join(row))
            prev_row_label = row_label
        return "\\\\\n".join(content)
