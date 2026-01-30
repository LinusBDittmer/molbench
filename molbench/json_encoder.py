"""

JSON Encoder for comparisons and Molecules

"""

from . import logger as log
from .molecule import Molecule, MoleculeList, Datapoint
from dataclasses import asdict
import json
import numpy

class MolbenchJSONEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, Datapoint):
            return {"value": obj.value, "unit": obj.unit}

        if isinstance(obj, Molecule):
            moldict = {obj.name: obj.system_data}
            moldict[obj.name]["properties"] = obj.state_data
            return moldict

        if isinstance(obj, numpy.ndarray):
            return list(obj)

        return obj
