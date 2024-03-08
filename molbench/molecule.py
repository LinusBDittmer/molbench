from . import logger as log


class Molecule:
    def __init__(self, name, data_id, system_data: dict = None,
                 state_data: dict = None) -> None:
        self.name = name
        self.data_id = data_id
        # dict that contains all the information regarding the system:
        # xyz_coords, charge, ...
        # -> information that is mostly required to build input files
        self.system_data: dict = {} if system_data is None else system_data
        # dict that contains all the information regarding the data points
        # for each state. Is of the form:
        # {state: {basis: _, method: _, ... type: _, value: _}
        self.state_data: dict = {} if state_data is None else state_data

    @classmethod
    def from_benchmark(cls, benchmark_entry: dict,
                       benchmark_id, molname=None) -> 'Molecule':
        # molname is only used as backup if name is not defined in the
        # benchmark
        system_data = {k: v for k, v in benchmark_entry.items()
                       if k != "properties"}
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
                         Molecule.from_benchmark)
        else:
            name = molname

        properties = benchmark_entry.get("properties", None)
        return cls(name, benchmark_id, system_data, properties)

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
                                 f"{name} and {n}", Molecule.from_external)
            data = metadata.get("data", None)
            metadata = {k: v for k, v in metadata.items()
                        if k not in ["name", "data"]}
            if data is None:
                state_data[state_id] = metadata
                continue

            if "type" in metadata or "value" in metadata:
                log.warning("The keys 'type' and 'value' in external data "
                            "will be overwritten when the structure is "
                            "flattened.", Molecule.from_external)
            for i, (proptype, value) in enumerate(data.items()):
                prop = metadata.copy()
                prop["type"] = proptype
                prop["value"] = value
                state_data[f"{state_id}_{i}"] = prop

        if name is None:
            if molname is None:
                log.critical("Name not specified in external data and not "
                             "provided as argument to the method.",
                             Molecule.from_external)
            name = molname
        return cls(name, data_id, None, state_data)
