from . import logger as log
from .functions import walk_dict_by_key
from typing import Any, Callable


class Molecule:

    __slots__ = ("name", "data_id", "system_data", "state_data")
    name: str
    data_id: str
    system_data: dict[str, Any]
    state_data: dict[str, Any]

    def __init__(self, name, data_id, system_data: dict | None = None,
                 state_data: dict | None = None) -> None:
        self.name = name
        self.data_id = data_id
        # dict that contains all the information regarding the system:
        # xyz_coords, charge, ...
        # -> information that is mostly required to build input files
        self.system_data = {} if system_data is None else system_data
        # dict that contains all the information regarding the data points
        # for each state. Is of the form:
        # {state: {basis: _, method: _, ... type: _, value: _}
        self.state_data = {} if state_data is None else state_data

    def __repr__(self):
        return f"{self.name}: {self.data_id}"

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
        if "xyz_list" in system_data and \
                isinstance(system_data["xyz_list"], (list, tuple)):
            system_data["xyz_list"] = [
                "\n".join(xyz) if not isinstance(xyz, str) else xyz
                for xyz in system_data["xyz_list"]
            ]
        # ensure that xyz coordinates are a string
        if "xyz" in system_data and not isinstance(system_data["xyz"], str):
            system_data["xyz"] = "\n".join(system_data["xyz"])
        # get a name for the molecule
        name = molname
        if "name" in system_data:
            name = system_data["name"]
            del system_data["name"]
        if name is None:
            log.critical("Name not specified in benchmark entry and not "
                         "provided as argument to the method.",
                         "Molecule: from_benchmark")

        properties = benchmark_entry.get("properties", None)
        return cls(name, benchmark_id, system_data, properties)

    @classmethod
    def from_external(cls, external_system_data: dict[str, Any],
                      external_state_data: dict[str, Any],
                      data_id: str, molname: str) -> 'Molecule':
        # Either system data or state data must exist
        if (not external_state_data) and (not external_system_data):
            log.critical("Both state and system data dicts are empty. "
                         f"Molecule {molname}", "Molecule: from_external")
        # We now define variables that are filled with the parsed data
        system_data: dict[str, Any] | None = None
        state_data: dict[str, Any] | None = None
        if external_system_data:
            system_data = external_system_data
        if external_state_data:
            state_data = external_state_data
        # Next we check if the state_data dictionary is correctly set up
        # We expect that the top level contains each state in the format
        # [str]: [dict]
        if state_data is not None:
            for state_key, state_values in external_state_data.items():
                # Make sure that all keys are strings
                if not isinstance(state_key, str):
                    log.critical("All state data keys must be strings. "
                                 f"{state_key} is of type {type(state_key)} "
                                 f"in Molecule {molname}",
                                 "Molecule: from_external")
                if not isinstance(state_values, dict):
                    log.critical("Incorrect state data structre. Found "
                                 f"{type(state_values)} where there "
                                 f"should be a dict in molecule {molname}",
                                 "Molecule: from_external")
                # State_values must have the keywords method, basis and data
                if "method" not in state_values:
                    log.critical("\"method\" keyword is required in state "
                                 f"data. State Key: {state_key}",
                                 "Molecule: from_external")
                if "basis" not in state_values:
                    log.critical("\"basis\" keyword is required in state "
                                 f"data. State Key: {state_key}",
                                 "Molecule: from_external")
                if "data" not in state_values:
                    log.critical("\"data\" keyword is required in state "
                                 f"data. State Key: {state_key}",
                                 "Molecule: from_external")
                # Lastly, we assert that the "data" dict is correctly set up
                # It should contain a series of entries, each containing a pair
                # of entries labeled "value" and "unit"
                if not isinstance(state_values["data"], dict):
                    log.critical("\"data\" keyword in state dictionary must "
                                 f"contain a dictionary (state: {state_key}, "
                                 f"molecule: {molname})",
                                 "Molecule: from_external")
                for dkey, dvals in state_values["data"].items():
                    # only value, unit allowed
                    if len(dvals) != 2:
                        log.critical("Incorrect datapoint specification: "
                                     f"{dvals}", "Molecule: from_external")
                    if "unit" not in dvals:
                        log.critical("Incorrect datapoint specification: "
                                     f"{dvals}", "Molecule: from_external")
                    if "value" not in dvals:
                        log.critical("Incorrect datapoint specification: "
                                     f"{dvals}", "Molecule: from_external")
                    dpoint: Datapoint = Datapoint(dvals["value"],
                                                  str(dvals["unit"]))
                    state_data[state_key]["data"][dkey] = dpoint

        return cls(molname, data_id, system_data, state_data)

    def add_assignments(self, assignments: dict,
                        old_transition_id_key: str = "transition_id",
                        new_transition_id_key: str = "assigned_transition_id") -> None:
        """
        Add the state assignment to the state data, i.e., try to add an
        assignment for each property in the state data.

        Parameters
        ----------
        assignments : dict
            The assignments to add, sorted as external: ref
        old_transition_id_key : str, optional
            The original identifier for the transition
            (default: 'transition_id').
        new_transition_id_key : str, optional
            The identifier under which the assigned transition id should pe placed
            (default: 'assigned_transition_id').
        """
        # We create a list to keep track of previously assigned states
        prev_assigned = list()
        # We also create a list of unassigned states to pop
        unassigned = list()

        # Next, we iterate over all states and look for the correct transition id
        for state_key, state_data in self.state_data.items():
            # If the property does not have a transition id, it cannot be assigned
            if old_transition_id_key not in state_data:
                continue
            # If the property is not found, it cannot be assigned
            tid = state_data[old_transition_id_key]
            if tid not in assignments:
                unassigned.append(state_key)
                continue
            # The assigned transition_id
            ass_tid = assignments[tid]
            if ass_tid in prev_assigned:
                log.warning(f"The transition id key {tid} is assigned to multiple"
                            "transitions. Overwriting.", "Molecule: add_assignments")
            state_data[new_transition_id_key] = ass_tid

        for key in unassigned:
            del self.state_data[key]

class MoleculeList(list[Molecule]):

    def __repr__(self):
        return super().__repr__()

    def __str__(self):
        return super().__str__()

    def filter(self, key, *values) -> 'MoleculeList':
        return self._filter(key, lambda v: v in values)

    def remove(self, key, *values) -> 'MoleculeList':
        return self._filter(key, lambda v: v not in values)

    def apply_stochiometry(self, stochiometry: dict) -> 'MoleculeList':
        # TODO: some explanation: e.g., what is the expected form of the
        # stochiometry dict?
        combined_list = MoleculeList()

        def find_mol(name):
            for mol in self:
                if mol.name == name:
                    return mol
            return None

        for c_name, c_data in stochiometry.items():
            relevant_mol_names = c_data["molecules"]
            data_id = " / ".join(relevant_mol_names)
            factors = c_data["factors"]
            relevant_mols = [find_mol(name) for name in relevant_mol_names]

            if any([x is None for x in relevant_mols]):
                log.error(f"Could not find all molecules for {c_name}",
                          "MoleculeList")
                continue

            combined_mol = Molecule(c_name, data_id, dict(), dict())

            for molidx, rmol in enumerate(relevant_mols):
                rmol: Molecule
                self._join_system_data(combined_mol.system_data,
                                       rmol.system_data, molidx)
                self._join_state_data(combined_mol.state_data, rmol.state_data,
                                      factors[molidx])

            combined_list.append(combined_mol)

        return combined_list

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

    def _filter(self: list[Molecule], key, callback: Callable):
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

    def _join_system_data(self, dst_sys, src_sys, idx):
        """Joins the system data dictionary of a desination system and a
        source system, where the destination system data dict is overwritten
        and the source system data dict is inserted at position idx.

        Args:
            dst_sys (dict): Destiation system data dict
            src_sys (dict): Source system data dict
            idx (int): Index at which the source dict is inserted
        """
        for key, val in src_sys.items():
            # In the inidividual molecules, each key is given "by itself"
            # In the joined Molecule, they are appended by "_list"
            lkey = str(key)
            if not lkey.endswith("_list"):
                lkey += "_list"
            if lkey not in dst_sys and idx > 0:
                # Add the key if it is not yet in dst_sys and pad with
                # the appropriate number of missing entries
                dst_sys[lkey] = [None] * idx
            elif lkey not in dst_sys:
                # If this is the first Molecule added into the combined
                # Molecule, then we do not need to pad.
                dst_sys[lkey] = []
            # Add the value to the combined Molecule
            dst_sys[lkey].append(val)

    def _join_state_data(self, dst_state, src_state, factor):
        """Joins the state data dictionary of a destination system and a
        source system, where the destination state data dict is overwritten
        and the source system data dict is merged under the application
        of a formally stochiometric factor.

        The decision, which properties to merge is performed by the
        individual properties type argument, i. e. the keys given at
        state_data[state_id]["data"].

        Args:
            dst_state (dict): Destination state data dict
            src_state (dict): Source state data dict
            factor (float): Formal stochiometric factor
        """
        # Create and cache a dictionary, which types of properties
        # are handled by which property keys in the destination state
        # property dictionary for easy lookup in the loop
        keydict: dict[str, str] = dict()
        for state_id, state_props in dst_state.items():
            all_props = list(state_props["data"].keys())
            keydict.update(
                {k: state_id for k in all_props}
            )
        for key, val in src_state.items():
            # Extract the property value of the source state data
            # dictionary for out-of-place modification since the
            # same source state data dictionary can influence
            # multiple combined molecules
            for datakey, datapoint in val["data"].items():
                datakey: str
                datapoint: Datapoint
                property_value = datapoint.value

                if isinstance(property_value, (list, tuple)):
                    property_value = [p * factor for p in property_value]
                else:
                    property_value *= factor

                if datakey not in keydict:
                    dst_state[datakey] = dict(val)
                    dst_state[datakey]["data"] = {
                        datakey: Datapoint(datapoint.value, datapoint.unit)
                    }
                else:
                    dst_dp = dst_state[keydict[datakey]]["data"][datakey]
                    if isinstance(property_value, (list, tuple)):
                        assert isinstance(dst_dp.value, (list, tuple))
                        dst_dp.value = [p0 + p1 for p0, p1 in
                                        zip(property_value, dst_dp.value)]
                    else:
                        dst_dp.value += property_value


class Datapoint:

    __slots__ = ("value", "unit")
    value: Any
    unit: str

    def __init__(self, value: Any, unit: str) -> None:
        self.value = value
        self.unit = unit

    def __repr__(self) -> str:
        return f"{self.value} {self.unit}"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Datapoint):
            return False
        return (self.value == other.value) and (self.unit.lower()
                                                == other.unit.lower())

    def __add__(self, other: Datapoint):
        if other.unit != self.unit:
            raise ValueError("Only datapoints with identical units can be compared")
        return Datapoint(self.value + other.value, self.unit)

    def __sub__(self, other: Datapoint):
        if other.unit != self.unit:
            raise ValueError("Only datapoints with identical units can be compared")
        return Datapoint(self.value - other.value, self.unit)

    def __abs__(self):
        return Datapoint(abs(self.value), self.unit)

    def __mul__(self, other):
        if isinstance(other, (float, int, complex)):
            return Datapoint(self.value * other, self.unit)
        raise ValueError("Datapoints can only multiplied by dimensionless quantities")

    def __truediv__(self, other):
        if isinstance(other, (float, int, complex)):
            return Datapoint(self.value / other, self.unit)
        raise ValueError("Datapoints can only divided by dimensionless quantities")

    def __floordiv__(self, other):
        if isinstance(other, (float, int, complex)):
            return Datapoint(self.value // other, self.unit)
        raise ValueError("Datapoints can only divided by dimensionless quantities")
