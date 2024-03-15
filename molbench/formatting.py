""" Formatting of datapoints """


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
