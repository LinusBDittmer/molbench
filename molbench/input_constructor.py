"""Python file for input constructors.

"""

import json
from pathlib import Path
from .assignment import new_assignment_file
from . import logger as log
from .configuration import config
from .functions import substitute_template
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
                      folder_structure_generator: callable,
                      stoch_name_generator: callable):

        basepath: Path = Path(basepath).resolve()
        if not basepath.exists():
            basepath.mkdir(parents=True, exist_ok=True)
        log.debug("Path created successfully.")

        # - data_iterable should return 1 data set for each file that are to be
        #   created
        # - file_name_generator is a callable that takes the data set and
        #   returns the filename as string
        # - file_content_generator is a callable that takes the data set and
        #   returns the content of the file
        # - folder_structure_generator is a callable that produces the
        #   relative path of folders where the file is to be placed:
        #   basepath / folder_path / file_name
        created_files = []
        for data in data_iterable:
            name: str = file_name_generator(data)
            stoch_name: str = stoch_name_generator(data)
            content: str = file_content_generator(data)
            folders: Path = folder_structure_generator(data)
            stochiometry: list = data[1][0].get("stochiometry", None)
            if stochiometry is None:
                stochiometry = [1.0]
            stochiometry_dict: dict = dict()

            # Content is a list of file contents
            for content_idx, (subcontent, subname) in enumerate(zip(content, name)):
                # create the folder and write the file content
                path = basepath / folders 
                if not path.exists():
                    path.mkdir(parents=True, exist_ok=True)
                if len(content) == 1:
                    file = path / subname
                else:
                    fileindex = "_" + str(content_idx).zfill(len(str(len(content))))
                    subname_proper: str = subname[:subname.rfind(".")]
                    subname_proper += fileindex + subname[subname.rfind("."):]
                    file = path / subname_proper
                if file.is_file():
                    log.warning(f"Overwriting existing file {file}.", "Input Constructor")
                stochiometry_dict[str(file)] = stochiometry[content_idx]
                with open(file, "w") as f:
                    f.write(subcontent)
                created_files.append(file)
                        
            if len(stochiometry_dict.keys()) > 1:
                path = basepath / folders
                file = path / stoch_name[0]
                if file.is_file():
                    log.warning(f"Overwriting existing file {file}.", "Input Constructor")
                with open(file, "w") as f:
                    json.dump(stochiometry_dict, f, ensure_ascii=True, 
                              sort_keys=True, indent=2)

        return created_files

    def _folders_from_tree(self, root: Node) -> Path:
        def _folder_structure(variant_data: dict):
            path = Path("")
            for generation in root.traverse_generations():
                folder_name = []
                for node in generation:
                    val = variant_data.get(node.value, None)
                    if val is None:
                        log.critical("Failed to resolve folder path for ",
                                     f"{variant_data}.", "Input Constructor")
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
                      name_template: str = None,
                      stochiometry_template: str = None) -> list:
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
        file_expansion_keys: tuple
            Possibly multiple input files need to be generated for one
            molecule, e.g., for different basis sets and/or multiplicities.
            This parameter defines the keys for which to check the Molecules.
            For each uninque combination of the corresponding values
            a file is generated.
        flat_structure: bool
            True: place all generated files in the same folder
            False: generate a folder for each Molecule
        name_template: str
            Template to define the names of the generated files. By default,
            a template based on the file_expansion_keys is used that ensures
            that all generated files have a unique name.
        """
        if name_template is None:
            name_template = self._default_name_template(file_expansion_keys,
                                                        ".in")
        if stochiometry_template is None:
            stochiometry_template = self._default_name_template(file_expansion_keys,
                                                                ".stoch")
        variant_data_iterator = self._molecule_variants_data_iter(
            benchmark, calc_details, file_expansion_keys
        )
        file_name_generator = self._substitute_template(name_template)
        stoch_filename_generator = self._substitute_template(stochiometry_template)
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
                                  folder_structure_generator, stoch_filename_generator)

    def create_assignments(self, benchmark: MoleculeList[Molecule],
                           basepath: str, calc_details: dict,
                           file_expansion_keys: tuple = ("basis",),
                           flat_structure: bool = False,
                           name_template: str = None,
                           state_id_key="state_id") -> list:
        """
        Create assignment files for the provided set of Molecules.

        Parameters
        ----------
        benchmark: MoleculeList
            The set of Molecules to generate assignment files for.
        basepath: str
            The path where all the generated files are placed.
        calc_details: dict
            Additional information to generate the filenames.
        file_expansion_keys: tuple
            Possibly multiple assignment files need to be generated for one
            molecule, e.g., for different basis sets and/or multiplicities.
            This parameter defines the keys for which to check the Molecules.
            For each uninque combination of the corresponding values
            a file is generated.
        flat_structure: bool
            True: place all generated files in the same folder
            False: generate a folder for each Molecule
        name_template: str
            Template to define the names of the generated files. By default,
            a template based on the file_expansion_keys is used that ensures
            that all generated files have a unique name.
        state_id_key
            The key where the state_ids an be found in the 'state_data' of the
            molecules.
        """
        if name_template is None:
            name_template = self._default_name_template(file_expansion_keys,
                                                        ".ass")
        variant_data_iterator = self._molecule_variants_data_iter(
            benchmark, calc_details, file_expansion_keys
        )
        file_name_generator = self._substitute_template(name_template)
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

    def _default_name_template(self, file_expansion_keys: tuple,
                               file_extension) -> str:
        # construct a name that ensures that each file has a unique name
        name_template = "[[name]]_[[method]]"
        for key in file_expansion_keys:
            if key not in ["name", "method"]:
                name_template += f"_[[{key}]]"
        return name_template + file_extension

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
                log.debug(f"Creating file for: {molecule.name} -> {var}.", "Template Constructor")
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

    def _gen_folder_structure(self, tree: Node):
        # just a wrapper to remove the molecule from the data, which is
        # needed for the assignment file content
        def _gen_folders(data):
            variant_data, _ = data
            return generator(variant_data)
        generator = self._folders_from_tree(tree)
        return _gen_folders

    def _gen_assignment_content(self, state_id_key):
        def _gen_assignment(data: tuple[dict, list]):
            variant_data, properties = data
            # collect all the state id's of properties of the molecule
            # that match the variant data
            state_ids = []
            for prop in properties:
                s_id = prop.get(state_id_key, None)
                if s_id is None:
                    log.warning(f"Property of molecule {variant_data['name']} "
                                f"has no assignment: {prop}", "Template Constructor")
                    continue
                if s_id not in state_ids:  # s_id has not to be hashable
                    state_ids.append(s_id)
            return new_assignment_file(state_ids)
        return _gen_assignment
