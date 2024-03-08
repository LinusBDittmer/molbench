"""

Comparison

Structure as follows

name -> basis -> method -> property -> id/path

"""

from . import logger as log
import numpy


class Comparison(dict):
    """
    structure of the nested comparison dict:
    - name / molkey
    - arbitrary number of user defined sort keys used in the provided order
    - type of the property (energy, ...)
    - data_id to avoid overwriting data (benchmark_id or the key used by the
      external parser to uniquely identify the individual read out files)
    """
    _data_separators = ("basis", "method")

    def __init__(self, *data_separators: str) -> None:
        if data_separators:
            self._data_separators = data_separators
        super().__init__()

    @property
    def data_separators(self):
        return self._data_separators

    @property
    def structure(self):
        return ("name", *self.data_separators, "proptype", "data_id")

    def _import_value(self, value):
        if isinstance(value, (int, float, complex, str)):
            return value
        else:
            return numpy.array(value)

    def walk_property(self, property):
        return self._walk_key(self, desired_key=property)

    def walk_values(self):
        """walk all values that are no dicts"""
        return self._walk_values(self)

    @staticmethod
    def _walk_values(indict: dict, prev_keys: list = None):
        """Walk the arbitrarily nested dictionary."""
        if prev_keys is None:
            prev_keys = []
        if isinstance(indict, dict):
            for key, val in indict.items():
                if isinstance(val, dict):
                    for d in Comparison._walk_values(val, prev_keys + [key]):
                        yield d
                else:
                    yield prev_keys + [key], val

    @staticmethod
    def _walk_key(indict: dict, desired_key, prev_keys: list = None):
        """Walk the arbitrarily nested dictionary until the desired key is
           found. Return the sequence of keys to reach the value that
           corresponds to the desired key."""
        if prev_keys is None:
            prev_keys = []
        if isinstance(indict, dict):
            for key, val in indict.items():
                if key == desired_key:
                    yield prev_keys + [desired_key], val
                elif isinstance(val, dict):
                    for d in Comparison._walk_key(val, desired_key,
                                                  prev_keys + [key]):
                        yield d

    def add_benchmark(self, benchmark: dict, benchmark_id: str) -> None:
        """
        Read in a benchmark of the following structure:
        {name: {..., 'properties': {_: {
            'basis': val,..., 'type': proptype1, value: 42
        }}}}
        """
        # XXX: 'type' and 'value' are benchmark specific
        # -> could add both as argument to the function
        if isinstance(benchmark, Comparison):
            log.error("Cannot parse a Comparison as a benchmark.")
            return
        for name, moldict in benchmark.items():
            properties = moldict.get("properties", None)
            if not properties:  # key does not exist or prop dict is empty
                continue
            for prop in properties.values():
                separators = [prop.get(key, None)
                              for key in self.data_separators]
                proptype = prop.get("type", None)
                value = prop.get("value", None)
                if proptype is None or any(v is None for v in separators) or \
                        value is None:
                    continue
                if name not in self:
                    self[name] = {}
                d = self[name]
                for separator in separators:
                    if separator not in d:
                        d[separator] = {}
                    d = d[separator]
                if proptype not in d:
                    d[proptype] = {}
                if benchmark_id in (d := d[proptype]):
                    log.warning("Benchmark ID is not unique. Found conflicting"
                                f" entry for {name}, {separators} and "
                                f"{proptype}. Overwriting the exisiting value",
                                Comparison)
                d[benchmark_id] = self._import_value(value)

    def add_external(self, external: dict) -> None:
        """
        Read in external data of the following structure:
        {outfile: {_: {
            data_separator1: val, data_separator2: val, ...,
            'data': {proptype1: val, proptype2: val, ...}
        }}}
        """
        # XXX: 'name' and 'data' are external specific
        # -> could add them as arguments to the function
        if isinstance(external, Comparison):
            log.error("Cannot parse a Comparison as external data.")
            return
        for outfile, dataset in external.items():
            for metadata in dataset.values():
                name = metadata.get("name", None)
                data = metadata.get("data", None)
                separators = [metadata.get(key, None)
                              for key in self.data_separators]
                if name is None or data is None or \
                        any(s is None for s in separators):
                    continue
                if name not in self:
                    self[name] = {}
                d = self[name]
                for separator in separators:
                    if separator not in d:
                        d[separator] = {}
                    d = d[separator]
                for proptype, value in data.items():
                    if proptype not in d:
                        d[proptype] = {}
                    if outfile in d[proptype]:
                        log.warning("Overwriting already existing value for "
                                    f"{name}, {separators} and {proptype}.",
                                    Comparison)
                    d[proptype][outfile] = self._import_value(value)
