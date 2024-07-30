from . import logger as log
from .functions import walk_dict_by_key


class Molecule:
    def __init__(self, name, data_id, system_data: dict = None,
                 state_data: dict = None, geometry_num: int = 1) -> None:
        self.name = name
        self.data_id = data_id
        # Whether this molecule encodes relative energies through
        # multiple molecules
        self.geometry_num = geometry_num
        # dict that contains all the information regarding the system:
        # xyz_coords, charge, ...
        # -> information that is mostly required to build input files
        self.system_data: dict = {} if system_data is None else system_data
        # dict that contains all the information regarding the data points
        # for each state. Is of the form:
        # {state: {basis: _, method: _, ... type: _, value: _}
        self.state_data: dict = {} if state_data is None else state_data
        # Compute relative properties from given data
        self._compute_relative_properties()

    @classmethod
    def from_benchmark(cls, benchmark_entry: dict,
                       benchmark_id, molname=None) -> 'Molecule':
        # molname is only used as backup if name is not defined in the
        # benchmark
        system_data = {k: v for k, v in benchmark_entry.items()
                       if k != "properties"}
        # Check if this molecule is a multi-Molecule entry,
        # e. g. through relative energies
        # If so, it contains the usual entries suffixed by "_list"
        # i. e. "xyz_list", "multiplicity_list", "n_atoms_list" etc.
        geometry_num: int = len(system_data["xyz_list"]) if "xyz_list" in system_data else 1 
        if "xyz_list" in system_data:
            if not isinstance(system_data["xyz_list"][0], str):
                system_data["xyz_list"] = ["\n".join(s) for s in system_data["xyz_list"]]
        # ensure that xyz coordinates are a string
        if "xyz" in system_data and not isinstance(system_data["xyz"], str):
            system_data["xyz"] = "\n".join(system_data["xyz"])
        # get a name for the molecule
        if "name" in system_data:
            name = system_data["name"]
            del system_data["name"]
        elif molname is None:
            log.critical("Name not specified in benchmark entry and not "
                         "provided as argument to the method.",
                         "Molecule: from_benchmark")
        else:
            name = molname

        properties = benchmark_entry.get("properties", None)
        return cls(name, benchmark_id, system_data, properties, geometry_num)

    @classmethod
    def from_external(cls, external: dict, data_id,
                      molname=None) -> 'Molecule':
        # possibly the name of the molecule is included in external
        # -> extract it
        # if it is not given -> use molname as backup
        # multiple properties may be listed in external['data']
        # -> flatten the structure by adding a counter to the state
        name = None
        state_data = {}
        for state_id, metadata in external.items():
            n = metadata.get("name", None)
            if n is not None:
                if name is None:
                    name = n
                elif name != n:
                    log.critical("Ambiguous name definition in external data: "
                                 f"{name} and {n}", "Molecule: from_external")
            data = metadata.get("data", None)
            metadata = {k: v for k, v in metadata.items()
                        if k not in ["name", "data"]}
            if data is None:
                state_data[state_id] = metadata
                continue

            if "type" in metadata or "value" in metadata:
                log.warning("The keys 'type' and 'value' in external data "
                            "will be overwritten when the structure is "
                            "flattened.", "Molecule: from_external")
            i = 0
            for proptype, value in data.items():
                prop = metadata.copy()
                prop["type"] = proptype
                prop["value"] = value
                # ensure that we don't overwrite any data!
                key = f"{state_id}_{i}"
                while key in external or key in state_data:
                    i += 1
                    key = f"{state_id}_{i}"
                state_data[key] = prop
                i += 1

        if name is None:
            if molname is None:
                log.critical("Name not specified in external data and not "
                             "provided as argument to the method.",
                             "Molecule: from_external")
            name = molname
        return cls(name, data_id, None, state_data)

    def add_assignments(self, assignments: dict,
                        state_id_key: str = "state_id") -> None:
        """
        Add the state assignment to the state data, i.e., try to add an
        assignment for each property in the state data.

        Parameters
        ----------
        assignments : dict
            The assignments to add.
        state_id_key : str, optional
            The key under which the assignment should be placed
            (default: 'state_id').
        """
        # during import the external_ids might be suffixed with '_number'
        # first check for an exact match in the assignment
        # if not found remove the extension and check if we find an assignment
        # for the shorter key.
        for external_id, state_data in self.state_data.items():
            if external_id not in assignments:
                ext_id = external_id.split("_")
                if not ext_id[-1].isdigit():  # unexpected extension
                    continue
                ext_id = "_".join(ext_id[:-1])
                if ext_id not in assignments:
                    continue
                external_id = ext_id

            # check if the key is already in use (None as value is ok)
            current_id = state_data.get(state_id_key, None)
            if current_id is not None:
                log.warning(f"The state id key {state_id_key} is already "
                            f"used in the state data of property {external_id}"
                            f" of molecule {self.name}. Overwriting the "
                            "existing value.", "Molecule: add_assignments")
            state_data[state_id_key] = assignments[external_id]

    def _compute_relative_properties(self):
        # List of relevant keys in state_data
        relative_keys: list[str] = [k for k in self.state_data if \
                                    "stochiometry" in self.state_data[k]]
        # Leave away all properties with precomputed stochiometry
        relative_keys = [k for k in relative_keys if 
                         isinstance(self.state_data[k]["stochiometry"], dict)]

        if len(relative_keys) == 0:
            return
        for relkey in relative_keys:
            # First mentioned component energy dict
            p0: str = tuple(self.state_data[relkey]["stochiometry"].keys())[0]
            # Method and basis are not necessarily set so we define them here
            for subkey in ("method", "basis"):
                if subkey not in self.state_data[relkey]:
                    self.state_data[relkey][subkey] = self.state_data[p0][subkey]

            # Value of p0 to get the type correct
            v0 = self.state_data[p0]["value"]
            relative_value = None
            if isinstance(v0, (float, int, complex)):
                relative_value = 0.0
            elif isinstance(v0, (tuple, list)):
                relative_value = list()

            for stoch_key, stoch_val in self.state_data[relkey]["stochiometry"].items():
                if isinstance(relative_value, list):
                    if len(relative_value) == 0:
                        for i in range(len(self.state_data[stoch_key]["value"])):
                            relative_value.append(0.0)
                    for idx in range(len(relative_value)):
                        relative_value[idx] += self.state_data[stoch_key]["value"][idx] * stoch_val
                else:
                    relative_value += self.state_data[stoch_key]["value"] * stoch_val
                
                if not self.state_data[stoch_key]["type"].startswith("component"):
                    self.state_data[stoch_key]["type"] = "component " + self.state_data[stoch_key]["type"]
            
            self.state_data[relkey]["value"] = relative_value


class MoleculeList(list):
    def filter(self, key, *values) -> 'MoleculeList':
        return self._filter(key, lambda v: v in values)

    def remove(self, key, *values) -> 'MoleculeList':
        return self._filter(key, lambda v: v not in values)

    def filter_by_range(self, key, min=None, max=None) -> 'MoleculeList':
        # here we don't know how to add elements to min and max
        # and we also don't know how anything about the value we want to
        # compare
        # -> user has to provide good structured input that can be directly
        # compared to the value in the data
        if min is None and max is None:
            return self

        def _filter_range(value) -> bool:
            if min is None:
                return max >= value
            elif max is None:
                return min <= value
            else:
                return min <= value and max >= value

        return self._filter(key, _filter_range)

    def filter_by_vec_norm(self, key, min=None, max=None) -> 'MoleculeList':
        # filter according to the vector norm and only keep states with a norm
        # min <= norm <= max
        # min and max should be provided as list of numbers of arbitrary length
        # but a single number is also possible
        # if more or only min values are provided max is filled with 1, i.e.,
        # a normalized vector norm is assumed (similarly min is filled with 0)
        if isinstance(min, (int, float)):
            min = (min,)
        if isinstance(max, (int, float)):
            max = (max,)

        if min is None and max is None:
            return self
        elif min is None:
            min = tuple(0 for _ in range(len(max)))
        elif max is None:
            max = tuple(1 for _ in range(len(min)))
        elif len(max) != len(min):
            n_missing = len(min) - len(max)
            if n_missing < 0:  # min has not fewer elements -> fill with 0
                min = (*min, *(0 for _ in range(abs(n_missing))))
            else:  # max has fewer elements than min -> fill with 1
                max = (*max, *(1 for _ in range(n_missing)))
        norm_range = tuple((lower, upper) for lower, upper in zip(min, max))

        def _filter_vec_norm(vec_norm) -> bool:
            # assumptions:
            # - if vec norm is a list/tuple it has to be sorted already
            # - if it is a dict it should be possible to obtain the correct
            #   order by sorting the keys
            if isinstance(vec_norm, dict):
                vec_norm = [norm for _, norm in sorted(vec_norm.items())]
            # ensure that vec_norm is at least as long as norm_range
            # -> fill up with 0
            if len(vec_norm) < len(norm_range):
                vec_norm = [
                    *vec_norm,
                    *(0 for _ in range(len(norm_range) - len(vec_norm)))
                ]
            # norm_range might be shorter than vec_norm
            for i, norm in enumerate(vec_norm):
                if len(norm_range) > i:  # have a user input for the norm
                    lower, upper = norm_range[i]
                    if lower > norm or upper < norm:
                        return False
                else:  # user input is exhausted -> allow all norms
                    break
            return True
        return self._filter(key, _filter_vec_norm)

    def _filter(self: list[Molecule], key, callback: callable):
        # XXX: currently it is not possible to filter according to state names!
        if key == "name":  # filter according to molecule names
            filtered = MoleculeList(m for m in self if callback(m.name))
        elif key == "data_id":
            filtered = MoleculeList(m for m in self if callback(m.data_id))
        else:  # check system_data and state_data for the key
            filtered = MoleculeList()
            for molecule in self:
                # start by checking system_data -> possibly drop the molecule
                if not all(callback(val) for _, val in
                           walk_dict_by_key(molecule.system_data, key)):
                    continue
                # now check the state_data -> possibly drop multiple states
                # if no states are left we drop the whole molecule
                remaining_states = [state for state, data in
                                    molecule.state_data.items()
                                    if all(callback(val) for _, val
                                           in walk_dict_by_key(data, key))]
                if not remaining_states:  # no state left -> drop molecule
                    continue
                # no state was removed
                if len(remaining_states) == len(molecule.state_data.keys()):
                    filtered.append(molecule)
                else:  # at least 1 state was removed
                    state_data = {state: molecule.state_data[state]
                                  for state in remaining_states}
                    filtered.append(Molecule(
                        molecule.name, molecule.data_id, molecule.system_data,
                        state_data
                    ))
        return filtered
