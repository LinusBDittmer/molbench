

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
                 value_delimiter: str = ", ") -> None:
        # TODO: construct the boilerplate table setup and customize it
        super().__init__(n_decimals, empty_field, value_delimiter)
        self.label_delimiter = "/"
        self.column_delimiter = " & "
        self.row_delimiter = r"\\" + "\n"

    def table_header(self, labels: list[list[str]]) -> str:
        # TODO: customize the layout -> more/fewer hlines
        return (
            self.join_rows(self.join_columns(row) for row in labels) +
            r"\\ \hline"
        )

    def table_content(self, content: list[list[str]]) -> str:
        # TODO: customize layout -> more/fewer hlines
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
        # TODO: custom allignment
        return r"\multicolumn{" + str(width) + r"}{c}{" + value + "}"

    def multirow(self, heigth: int, value: str) -> str:
        # TODO: custom width of the cell
        return r"\multirow{" + str(heigth) + r"}{*}{" + value + "}"
