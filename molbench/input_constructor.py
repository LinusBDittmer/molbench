"""Python file for input constructors.

"""

import json
from pathlib import Path
from collections import Counter, defaultdict
from .assignment import new_assignment_file
from . import logger as log
from .configuration import config
from .functions import substitute_template, default_name_template
from .molecule import MoleculeList, Molecule
from .tree import Node, DummyNode


class InputConstructor:
    """
    Parent class for an Inputfile-API.

    This class is used to write input files for a benchmark.
    Any class that is supposed to write input files for a benchmark must
    inherit from this class.

    Methods
    -------
    create(benchmark: dict, filepath: str, flat_structure: bool = False,
           name_template: str = '[[name]]_[[method]]_[[basis]].in') -> None
        Create input files for the given benchmark at the specified filepath.

    Attributes
    ----------
    None
    """

    def create_inputs(self, *args, **kwargs) -> list:
        raise NotImplementedError("The 'create_inputs' method is only "
                                  "implemented on child classes.")

    def create_assignments(self, *args, **kwargs) -> list:
        raise NotImplementedError("The 'create_assignments' method is only "
                                  "implemented on child classes.")

    def _create_files(self, data_iterable, basepath: str,
                      file_name_generator: callable,
                      file_content_generator: callable,
                      folder_structure_generator: callable):

        basepath: Path = Path(basepath).resolve()
        if not basepath.exists():
            basepath.mkdir(parents=True, exist_ok=True)
        log.debug("Path created successfully.")

        # - data_iterable should return 1 data set for each data point, e.g.,
        #   each molecule, that is to be created
        # - file_name_generator is a callable that takes the data set and
        #   returns the names of all files to create for the given data point
        # - file_content_generator is a callable that takes the data set and
        #   returns the content of all files of the given data point
        # - folder_structure_generator is a callable that produces the
        #   relative path of folders where the files are placed:
        #   basepath / folder_path / file_name
        created_files = []
        for data in data_iterable:
            name_list: tuple[str, ...] = file_name_generator(data)
            content_list: tuple[str, ...] = file_content_generator(data)
            folders: Path = folder_structure_generator(data)
            assert len(name_list) == len(content_list)

            # create the folder
            path = basepath / folders
            if not path.exists():
                path.mkdir(parents=True, exist_ok=True)

            # create the actual files
            for content, name in zip(content_list, name_list):
                file = path / name
                if file.is_file():
                    log.warning(f"Overwriting existing file {file}.",
                                "Input Constructor")
                # write the file
                with open(file, "w") as f:
                    f.write(content)
                created_files.append(file)
        return created_files

    def _folders_from_tree(self, root: Node):
        def _folder_structure(variant_data: dict):
            path = Path("")
            for generation in root.traverse_generations():
                folder_name = []
                for node in generation:
                    val = variant_data.get(node.value, None)
                    if val is None:
                        log.critical("Failed to resolve folder path for "
                                     + f"{variant_data}.", "Input Constructor")
                    folder_name.append(node.to_string(val))
                path /= "_".join(folder_name)
            return path
        return _folder_structure


class TemplateConstructor(InputConstructor):
    """
    Constructor that creates input files by substituting from a template.
    """

    def __init__(self, template: str):
        self.init_template(template)

    def init_template(self, template: str):
        # - try to load a template file from the template folder
        template_dir = Path(__file__).parent.resolve() / "templates"
        template_file = None
        for template_f in template_dir.rglob("*.txt"):
            if not template_f.is_file():
                continue
            elif template_f.stem == template:
                template_file = template_f
                break
        # - try to load a local template if we did not find a template yet
        if template_file is None:
            template_file = Path(template).resolve()
        # actually try to read the file
        try:
            with open(template_file, 'r') as f:
                self.template = f.read()
        except Exception:
            log.critical(f"Template {template} could not be loaded.",
                         "Template Constructor")

    def create_inputs(self, benchmark: MoleculeList[Molecule], basepath: str,
                      calc_details: dict,
                      file_expansion_keys: tuple = ("basis",),
                      flat_structure: bool = False,
                      name_template: str | None = None) -> list:
        """
        Create inputs files for the provided set of Molecules by filling
        in the placeholders in the input template with data from the
        benchmark set itself, additional user input provided through
        'calc_details' and the global config.

        Parameters
        ----------
        benchmark: MoleculeList
            The set of Molecules to generate input files for.
        basepath: str
            The path where all the generated files are placed.
        calc_details: dict
            Additional information used to resolve the placeholders in the
            template
        file_expansion_keys: tuple, optional
            Possibly multiple input files need to be generated for one
            molecule that correspond to different data point represented by
            e.g. different basis sets and/or multiplicities.
            This parameter defines the keys for which the Molecules
            are checked. Each uninque combination of the corresponding
            values defines a data point for which input files are generated.
            (default: ('basis',))
        flat_structure: bool, optional
            True: place all generated files in the same folder
            False: generate a folder for each Molecule
            (default: False)
        name_template: str, optional
            Template to define the names of the generated files. By default,
            a template based on the file_expansion_keys is used that ensures
            that each data point has a unique name.
        """
        if name_template is None:
            name_template = default_name_template(file_expansion_keys, ".in")
        variant_data_iterator = self._molecule_variants_data_iter(
            benchmark, calc_details, file_expansion_keys
        )
        file_name_generator = self._gen_file_names(name_template)
        file_content_generator = self._substitute_template(self.template)

        # create a tree representing the folder structure
        # TODO: allow more complex input for more complex folder
        #       structure?
        if flat_structure:
            tree = DummyNode()
        else:
            tree = Node("name")

        folder_structure_generator = self._gen_folder_structure(tree)

        return self._create_files(variant_data_iterator, basepath,
                                  file_name_generator, file_content_generator,
                                  folder_structure_generator)

    def create_assignments(self, benchmark: MoleculeList[Molecule],
                           basepath: str, calc_details: dict,
                           file_expansion_keys: tuple = ("basis",),
                           flat_structure: bool = False,
                           name_template: str | None = None,
                           state_id_key="state_id") -> list:
        """
        Create assignment files for the provided set of Molecules.
        Note: This does currently not work for relative properties!

        Parameters
        ----------
        benchmark: MoleculeList
            The set of Molecules to generate assignment files for.
        basepath: str
            The path where all the generated files are placed.
        calc_details: dict
            Additional information to generate the filenames.
        file_expansion_keys: tuple, optional
            Possibly multiple assignment files need to be generated for one
            molecule, i.e, there are multiple data points representing e.g.,
            different basis sets and/or multiplicities.
            This parameter defines the keys for which to check the Molecules.
            Each unique combination defines a data point and for which
            assignment files are generated.
            (default: ('basis',))
        flat_structure: bool, optional
            True: place all generated files in the same folder
            False: generate a folder for each Molecule
            (default: False)
        name_template: str, optional
            Template to define the names of the generated files. By default,
            a template based on the file_expansion_keys is used that ensures
            that all data points have a unique name.
        state_id_key
            The key under which the state_ids an be found in the
            properties (the state_data) of the molecules.
        """
        if name_template is None:
            name_template = default_name_template(file_expansion_keys, ".ass")
        variant_data_iterator = self._molecule_variants_data_iter(
            benchmark, calc_details, file_expansion_keys
        )
        file_name_generator = self._gen_file_names(name_template)
        file_content_generator = self._gen_assignment_content(state_id_key)

        # create a tree representing the folder structure
        # TODO: allow more complex input for more complex folder
        #       structure?
        if flat_structure:
            tree = DummyNode()
        else:
            tree = Node("name")
        folder_structure_generator = self._gen_folder_structure(tree)

        return self._create_files(variant_data_iterator, basepath,
                                  file_name_generator, file_content_generator,
                                  folder_structure_generator)

    def _molecule_variants_data_iter(self, benchmark: MoleculeList[Molecule],
                                     calc_details: dict,
                                     file_expansion_keys: tuple):
        # for each molecule:
        #  find all unique combination of relevant keys, for instance,
        #  perform calculations for different basis sets in independent
        #  files.
        for molecule in benchmark:
            molecule: Molecule
            # list instead of set -> values have not to be hashable
            variants = []
            variant_properties = []
            for property in molecule.state_data.values():
                var = tuple((key, property.get(key, None))
                            for key in file_expansion_keys)
                if any(val is None for _, val in var):
                    continue
                try:  # no new variant -> add the property to the list
                    i = variants.index(var)
                    variant_properties[i].append(property)
                except ValueError:  # new variant
                    variants.append(var)
                    variant_properties.append([property])
            for var, props in zip(variants, variant_properties):
                log.debug(f"Creating file for: {molecule.name} -> {var}.",
                          "Template Constructor")
                # collect all the relevant data
                variant_data: dict = molecule.system_data.copy()
                if "name" in variant_data:
                    log.warning("The key 'name' in the molecules 'system_data'"
                                " is reservedfor the name of the molecule. "
                                "Overwriting existing value "
                                f"{variant_data['name']} with "
                                f"{molecule.name}.", "Template Constructor")
                variant_data["name"] = molecule.name
                # add the relevant subset of state_data for resolving the
                # current variant
                for key, val in var:
                    if key in variant_data:
                        log.warning(f"Found conflicting entry for {key}. "
                                    f"Overwriting existing value "
                                    f"{variant_data[key]} with {val}.",
                                    "TemplateConstructor")
                    variant_data[key] = val
                # add the additional data from the user
                variant_data.update(calc_details)
                # add the data from the config file as backup
                for key, val in config.items():
                    if key not in variant_data:
                        variant_data[key] = val
                # for the assignment files we also need the state data of the
                # states belonging to the current variant
                yield variant_data, props

    def _substitute_template(self, template: str):
        def _substitute(data):
            subvals, _ = data
            return substitute_template(template, subvals)
        return _substitute

    def _gen_file_names(self, name_template: str):
        # generate filenames by resolving the given template using the
        # data for the given data point.
        # if this results in multiple identical file names, a counter
        # starting at 0 is added to the corresponding names.
        def _name_generator(data) -> tuple[str, ...]:
            # resolve the template
            subvals, _ = data
            file_names = substitute_template(name_template, subvals)
            if len(file_names) == 1:
                return file_names
            # more than 1 file name
            # -> if the same file_name is encountered more than once
            #    we add an counter to the name to ensure it is unique:
            #    name.suffix -> name_number.suffix
            name_counter = Counter(file_names)
            f_index = defaultdict(int)
            ret = []
            for fname in file_names:
                if (n_files := name_counter[fname]) > 1:
                    idx = str(f_index[fname]).zfill(len(str(n_files)))
                    f_index[fname] += 1
                    fname = Path(fname)
                    fname = f"{fname.stem}_{idx}{fname.suffix}"
                ret.append(fname)
            return tuple(ret)
        return _name_generator

    def _gen_folder_structure(self, tree: Node):
        # just a wrapper to remove the molecule from the data, which is
        # needed for the assignment file content
        def _gen_folders(data):
            variant_data, _ = data
            return generator(variant_data)
        generator = self._folders_from_tree(tree)
        return _gen_folders

    def _gen_assignment_content(self, state_id_key):
        def _gen_assignment(data: tuple[dict, list]) -> tuple[str]:
            variant_data, properties = data
            # collect all the state id's of properties of the molecule
            # that match the variant data
            state_ids = []
            for prop in properties:
                s_id = prop.get(state_id_key, None)
                if s_id is None:
                    log.warning(f"Property of molecule {variant_data['name']} "
                                f"has no assignment: {prop}",
                                "Template Constructor")
                    continue
                if s_id not in state_ids:  # s_id has not to be hashable
                    state_ids.append(s_id)
            return (new_assignment_file(state_ids),)
        return _gen_assignment


class CompressedTemplateConstructor(TemplateConstructor):

    def __init__(self, template: str):
        super().__init__(template)

    def create_inputs(self, benchmark: MoleculeList[Molecule], basepath: str,
                      calc_details: dict,
                      file_expansion_keys: tuple = ("basis",),
                      flat_structure: bool = False,
                      name_template: str | None = None,
                      reference_path: str = "references.json",
                      compressed_property: str | None = None) -> list:
        # Create compressed benchmark
        # We create a new MoleculeList where each Molecule contains
        # only one geometry.
        # Additionally, we create a dict of references to the individual
        # molecules

        compressed = MoleculeList()
        references = {}

        def _unique(xyz: list, charge: int, mult: int) -> int:
            all_xyzs = [(m.system_data["xyz"],
                         m.system_data["charge"],
                         m.system_data["multiplicity"]) for m in compressed]
            xyz_list = (xyz, charge, mult)
            if xyz_list in all_xyzs:
                return all_xyzs.index(xyz_list)
            return -1

        for mol in benchmark:
            if "xyz_list" not in mol.system_data:
                mol_counter = len(compressed)
                mol.name = f"m{mol_counter:06d}"
                compressed.append(mol)
                references[mol.name] = (mol_counter,)
                mol_counter += 1
                continue

            # Number of Molecules in mol
            n_mols = len(mol.system_data["xyz_list"])
            references[mol.name] = {"molecules": list(),
                                    "factors": list()}
            for i in range(n_mols):
                mol_counter = len(compressed)
                idx = _unique(mol.system_data["xyz_list"][i],
                              mol.system_data["charge_list"][i],
                              mol.system_data["multiplicity_list"][i])

                if idx < 0:
                    # Prepare new Molecule
                    name = f"m{mol_counter:06d}"
                    system_data = dict()
                    for system_dp, system_val in mol.system_data.items():
                        if system_dp.endswith("_list"):
                            sd = system_dp[:-5]  # -5 to cut off "_list"
                            system_data[sd] = system_val[i]
                        else:
                            system_data[system_dp] = system_val
                    newmol = Molecule(name, name, system_data, mol.state_data)

                    compressed.append(newmol)
                    idx = mol_counter

                idx_name = f"m{idx:06d}"
                references[mol.name]["molecules"].append(idx_name)

                if len(mol.state_data) > 1 and compressed_property is None:
                    log.critical("Please specify a property key from which the"
                                 + " stochiometry should be read",
                                 "CompressedTemplateConstructor")
                elif (len(mol.state_data) > 1 and
                      compressed_property not in mol.state_data):
                    log.critical("compressed_property was not found in"
                                 + f" Molecule {mol.name}",
                                 "CompressedTemplateConstructor")
                if compressed_property is None:
                    pkey = list(mol.state_data.keys())[0]
                else:
                    pkey = [k for k, v in mol.state_data.items()
                            if v["type"] == compressed_property][0]

                factor_key = None
                if "stochiometry" in mol.state_data[pkey]:
                    factor_key = "stochiometry"
                elif "factors" in mol.state_data[pkey]:
                    factor_key = "factors"

                if factor_key is None:
                    log.critical(f"Could not find a factor in state {pkey} of "
                                 f"molecule {mol.name}.",
                                 "CompressedTemplateConstructor")
                stoch: list = mol.state_data[pkey][factor_key]
                references[mol.name]["factors"] = stoch

        inputs = super().create_inputs(compressed, basepath, calc_details,
                                       file_expansion_keys, flat_structure,
                                       name_template)

        full_reference_path = Path(basepath) / Path(reference_path)
        with open(full_reference_path, "w") as f:
            json.dump(references, f, ensure_ascii=True, indent=4,
                      sort_keys=True)

        return inputs
