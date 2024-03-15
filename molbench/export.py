from .comparison import Comparison
from .formatting import StdFormatter, Formatter
from collections import defaultdict
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

    def export(self, comparison: Comparison, *args, **kwargs) -> str:
        raise NotImplementedError("Export function has to be implemented on "
                                  "the child classes.")


class TableExporter(Exporter):
    def build_row_column_label(self, columns: tuple, rows: tuple,
                               sep_names: tuple,
                               sep_vals: tuple) -> tuple[tuple, tuple]:
        row_l, column_l = [], []
        for name, val in zip(sep_names, sep_vals):
            if name in columns:
                column_l.append(val)
            elif name in rows:
                row_l.append(val)
        return "/".join(row_l), "/".join(column_l)

    def prepare_export(self, comparison: Comparison, columns: tuple,
                       rows: tuple, prop):
        column_labels = set()
        data = defaultdict(lambda: defaultdict(list))
        for keys, propdata in comparison.walk_by_key(prop):
            for data_id, value in propdata.items():
                row_l, column_l = self.build_row_column_label(
                    columns, rows, comparison.structure, keys + (data_id,)
                )
                column_labels.add(column_l)
                data[row_l][column_l].append(value)
        return data, column_labels


class LatexExporter(TableExporter):
    def __init__(self):
        # TODO: allow to modify the style of the latex table
        pass

    def export(self, data: Comparison, prop, outfile: typing.IO,
               columns: tuple, rows: tuple = None,
               formatter: Formatter = None) -> None:
        # if no row labels are provided, use all labels
        # -> will put everything in row_label that is not defined as column
        if rows is None:
            rows = data.structure
        if formatter is None:
            formatter = StdFormatter()

        prepared_data, column_labels = self.prepare_export(data, columns, rows,
                                                           prop)
        column_labels = sorted(column_labels)
        # TODO: set up table
        out = []
        out.append(" & ".join(("", *column_labels)))
        for row_l, row_data in prepared_data.items():
            row = [row_l,
                   *(formatter.format_datapoint(row_data.get(column_l, None))
                     for column_l in column_labels)]
            out.append(" & ".join(row))
        out = r"\\\n".join(out)
        # TODO: finalize table
        outfile.write(out)
