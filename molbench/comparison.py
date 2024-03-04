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
        """
        Initialize the Comparison class.

        Parameters
        ----------
        *data_separators : str
            Data separators used to split data into nested dictionaries.
        
        """
        if data_separators:
            self._data_separators = data_separators
        super().__init__()

    @property
    def data_separators(self):
        """
        Get the data separators.

        Returns
        -------
        tuple of str
            Data separators used to split data into nested dictionaries.
        
        """
        return self._data_separators

    @property
    def structure(self):
        """
        Get the structure of the nested comparison dictionary.

        Returns
        -------
        tuple of str
            Structure of the nested comparison dictionary.
        
        """
        return ("name", *self.data_separators, "proptype", "data_id")

    def _import_value(self, value):
        """
        Import value into numpy array if it is not a scalar type.

        Parameters
        ----------
        value : scalar or array-like
            Value to be imported.

        Returns
        -------
        scalar or numpy.ndarray
            Imported value.
        
        """
        if isinstance(value, (int, float, complex, str)):
            return value
        else:
            return numpy.array(value)

    def walk_property(self, property):
        """
        Walk through properties in the comparison dictionary.

        Parameters
        ----------
        property : str
            Property to walk through.

        Returns
        -------
        generator
            Generator yielding paths and corresponding values for the specified property.
        
        """
        return self._walk_key(self, desired_key=property)

    def walk_values(self):
        """
        Walk through all values in the comparison dictionary that are not dictionaries.

        Returns
        -------
        generator
            Generator yielding paths and corresponding values for all values that are not
            dictionaries.
        
        """
        return self._walk_values(self)

    @staticmethod
    def _walk_values(indict: dict, prev_keys: list = None):
        """
        Walk through all values in the nested dictionary.

        Parameters
        ----------
        indict : dict
            Nested dictionary to walk through.
        prev_keys : list of str, optional
            List of previous keys leading to the current dictionary. Default is None.

        Returns
        -------
        generator
            Generator yielding paths and corresponding values. 
        """
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
        """
        Walk through nested dictionary until the desired key is found.

        Parameters
        ----------
        indict : dict
            Nested dictionary to walk through.
        desired_key : str
            Key to search for.
        prev_keys : list of str, optional
            List of previous keys leading to the current dictionary. Default is None.

        Returns
        -------
        generator
            Generator yielding paths and corresponding values for the desired key.

        """
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
        Add benchmark data to the comparison dictionary. Read in a benchmark of the following 
        structure:
        {name: {..., 'properties': {_: {
            'basis': val,..., 'type': proptype1, value: 42
        }}}}


        Parameters
        ----------
        benchmark : dict
            Benchmark data to be added.
        benchmark_id : str
            Identifier for the benchmark data.

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
                                f" entry for {separators} and {proptype}."
                                "Overwriting the exisiting value", Comparison)
                d[benchmark_id] = self._import_value(value)

    def add_external(self, external: dict) -> None:
        """
        Add external data to the comparison dictionary. Reads in external data of the following 
        structure:
        {outfile: {_: {
            data_separator1: val, data_separator2: val, ...,
            'data': {proptype1: val, proptype2: val, ...}
        }}}


        Parameters
        ----------
        external : dict
            External data to be added.

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
