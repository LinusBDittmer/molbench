from typing import Any


class Formatter:
    """A base class for formatting datapoints.
    """
    def format_datapoint(self, value: Any) -> str:
        """Returns a formatted datapoint value

        Parameters
        ----------
        value : Any
            The datapoint to be formatted. Can be either a number or a
            list/tuple

        Returns
        -------
        str
            The formatted datapoint as a string

        Raises
        ------
        NotImplementedError
            The base class does not implement a formatting
        """
        raise NotImplementedError("Method not implemented on the base class.")


class StdFormatter(Formatter):

    def __init__(self, n_decimals: int = 5, empty_field: str = "",
                 value_delimiter: str = ", ") -> None:
        """The standard formatter for text.

        Parameters
        ----------
        n_decimals : int, optional
            Number of decimal places to print, by default 5
        empty_field : str, optional
            Dummy value for an empty datapoint, by default ""
        value_delimiter : str, optional
            Delimiter between two datapoint values, by default ", "
        """
        self.n_decimals: int = n_decimals
        self.empty_field: str = empty_field
        self.value_delimiter: str = value_delimiter

    def format_datapoint(self, value: Any) -> str:
        """Returns a formatted datapoint value.
        The standard formatter returns a number (int, float) as a number
        rounded to the given number of decimal places. An iterable type is
        formatted as a series of individual datapoints concatenated by the
        class attribute "value_delimiter".

        Parameters
        ----------
        value : Any
            The datapoint to be formatted

        Returns
        -------
        str
            The formatted datapoint
        """
        if isinstance(value, (int, float)):
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
                 column_alignment: str = "c",
                 additional_column_alignment: str = "l",
                 multicol_alignment: str = "c",
                 multirow_width: str = "*") -> None:
        """Used to create a LaTeX table of all datapoints

        Parameters
        ----------
        n_decimals : int, optional
            Number of decimal digits for rounding, by default 5
        empty_field : str, optional
            How to represent an empty datapoint, by default ""
        value_delimiter : str, optional
            Delimiter for iterable datapoints, by default ", "
        label_delimiter : str, optional
            Delimiter for labels in the header, by default "/"
        column_delimiter : str, optional
            The delimiter for different columns, by default " & "
        row_delimiter : str, optional
            The delimiter for different rows, by default "\\ \n"
        column_alignment : str, optional
            Alignment for the columns, by default "c"
        additional_column_alignment : str, optional
            Alignment for the additional columns, by default "l"
        multicol_alignment : str, optional
            Aligment for multiple joined columns, by default "c"
        multirow_width : str, optional
            Width of joined rows, by default "*"
        """
        super().__init__(n_decimals, empty_field, value_delimiter)
        self.label_delimiter = label_delimiter
        self.column_delimiter = column_delimiter
        self.row_delimiter = row_delimiter
        self.column_alignment = column_alignment
        self.additional_column_alignment = additional_column_alignment
        self.multicol_alignment = multicol_alignment
        self.multirow_width = multirow_width

    def init_table(self, n_additional_cols: int, n_columns: int,) -> str:
        alignment = (self.additional_column_alignment * n_additional_cols +
                     "|" + self.column_alignment * n_columns)
        return r"\begin{table}" + "\n" + r"\begin{tabular}{" + alignment + "}"

    def finalize_table(self):
        return r"\end{tabular}" + "\n" + r"\end{table}"

    def table_header(self, labels: tuple[tuple[str, ...], ...]) -> str:
        return (
            self.join_rows(tuple(self.join_columns(row) for row in labels)) +
            r"\\ \hline"
        )

    def table_content(self, content: tuple[tuple[str, ...], ...]) -> str:
        return (
            self.join_rows(tuple(self.join_columns(row) for row in content))
        )

    def join_labels(self, labels: tuple[str, ...]) -> str:
        return self.label_delimiter.join(labels)

    def join_columns(self, columns: tuple[str, ...]) -> str:
        return self.column_delimiter.join(columns)

    def join_rows(self, rows: tuple[str, ...]) -> str:
        return self.row_delimiter.join(rows)

    def multicolumn(self, width: int, value: str) -> str:
        return (r"\multicolumn{" + str(width) + "}{" +
                self.multicol_alignment + "}{" + value + "}")

    def multirow(self, heigth: int, value: str) -> str:
        return (r"\multirow{" + str(heigth) + "}{" +
                self.multirow_width + "}{" + value + "}")
