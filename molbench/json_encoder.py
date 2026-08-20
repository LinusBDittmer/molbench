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
            return {obj.name: dict(obj.system_data, properties=obj.state_data)}

        if isinstance(obj, numpy.ndarray):
            return list(obj)

        if isinstance(obj, numpy.integer):
            return int(obj)

        if isinstance(obj, numpy.floating):
            return float(obj)

        return super().default(obj)
