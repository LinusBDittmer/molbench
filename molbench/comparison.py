"""

Comparison

Structure as follows

name -> basis -> method -> property -> id/path

"""

from . import logger as log
from .molecule import Molecule
from .functions import walk_dict_by_key, walk_dict_values
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
            # remove the special separators (used anyway)
            data_separators = tuple(
                separator for separator in data_separators
                if separator not in ["name", "proptype", "data_id"]
            )
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

    def add(self, dataset: tuple[Molecule]) -> None:
        for data in dataset:
            self.add_molecule(data)

    def add_molecule(self, data: Molecule) -> None:
        if not isinstance(data, Molecule):
            log.error(f"Can't add data of type {type(data)}.",
                      "Comparison.add_molecule")
        for prop in data.state_data.values():
            separators = [prop.get(key, None) for key in self.data_separators]
            proptype = prop.get("type", None)
            value = prop.get("value", None)
            if proptype is None or value is None or any(v is None
                                                        for v in separators):
                continue
            # move into and establish the nested dict structure
            if data.name not in self:  # special separator: name
                self[data.name] = {}
            d = self[data.name]
            for sep in separators:
                if sep not in d:
                    d[sep] = {}
                d = d[sep]
            if proptype not in d:  # special separator: proptype
                d[proptype] = {}
            d = d[proptype]
            if data.data_id in d:
                log.warning(f"data_id {data.data_id}is not unique. Found "
                            f"conflicting entry for {data.name}, {separators} "
                            f"and {proptype}. Overwriting the exisiting value",
                            "Comparison.add_molecule")
            d[data.data_id] = self._import_value(value)

    def walk_property(self, property):
        return walk_dict_by_key(self, desired_key=property)

    def walk_values(self):
        """walk all values that are no dicts"""
        return walk_dict_values(self)
