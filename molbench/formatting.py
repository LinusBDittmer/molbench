

class Formatter:
    def format_datapoint(self, value):
        raise NotImplementedError("Method not implemented on the base class.")


class StdFormatter(Formatter):

    def __init__(self, n_decimals: int = 5, empty_field: str = "",
                 value_delimiter: str = ", ") -> None:
        self.n_decimals = n_decimals
        self.empty_field = empty_field
        self.value_delimiter = value_delimiter

    def format_datapoint(self, value) -> str:
        if isinstance(value, (int, float, complex)):
            return str(round(value, self.n_decimals))
        elif hasattr(value, '__iter__'):  # dict, set, list, tuple, ...
            return self.value_delimiter.join(
                self.format_datapoint(v) for v in value
            )
        elif value is None:
            return self.empty_field
        else:
            return str(value)


class LatexFormatter(StdFormatter):
    def __init__(self, n_decimals: int = 5, empty_field: str = "",
                 value_delimiter: str = ", ",
                 label_delimiter: str = "/",
                 column_delimiter: str = " & ",
                 row_delimiter: str = "\\\\ \n",
                 column_allignment: str = "c",
                 additional_column_allignment: str = "l",
                 multicol_allignment: str = "c",
                 multirow_width: str = "*") -> None:
        super().__init__(n_decimals, empty_field, value_delimiter)
        self.label_delimiter = label_delimiter
        self.column_delimiter = column_delimiter
        self.row_delimiter = row_delimiter
        self.column_allignment = column_allignment
        self.additional_column_allignment = additional_column_allignment
        self.multicol_alignment = multicol_allignment
        self.multirow_width = multirow_width

    def init_table(self, n_additional_cols: int, n_columns: int,) -> str:
        allignment = (self.additional_column_allignment * n_additional_cols +
                      "|" + self.column_allignment * n_columns)
        return r"\begin{table}" + "\n" + r"\begin{tabular}{" + allignment + "}"

    def finalize_table(self):
        return r"\end{tabular}" + "\n" + r"\end{table}"

    def table_header(self, labels: list[list[str]]) -> str:
        return (
            self.join_rows(self.join_columns(row) for row in labels) +
            r"\\ \hline"
        )

    def table_content(self, content: list[list[str]]) -> str:
        return (
            self.join_rows(self.join_columns(row) for row in content)
        )

    def join_labels(self, labels: tuple[str]) -> str:
        return self.label_delimiter.join(labels)

    def join_columns(self, columns: tuple[str]) -> str:
        return self.column_delimiter.join(columns)

    def join_rows(self, rows: tuple[str]) -> str:
        return self.row_delimiter.join(rows)

    def multicolumn(self, width: int, value: str) -> str:
        return (r"\multicolumn{" + str(width) + "}{" +
                self.multicol_alignment + "}{" + value + "}")

    def multirow(self, heigth: int, value: str) -> str:
        return (r"\multirow{" + str(heigth) + "}{" +
                self.multirow_width + "}{" + value + "}")
