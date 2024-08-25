"""Python file for all external parsers.

"""
from . import logger as log
from .molecule import Molecule, MoleculeList
from .assignment import parse_assignment_file
import os.path
from pathlib import Path
import numpy
import json


class ExternalParser:
    """
    Parent class for an Output-API.

    This class is used to load external data from a directory. Any class that
    is supposed to load external data from a directory must inherit from this
    class. The naming convention for subclasses is [SUITE]_[METHOD]_Parser
    (e. g. PySCF_FCI_Parser, QChem_MP2_Parser, ORCA_CCSDT_Parser).

    Methods
    -------
    load(filepath: str, suffix: str = 'out') -> dict
        Load external data from the specified filepath with the given suffix.

    Attributes
    ----------
    None
    """
    _registry: dict[str, 'ExternalParser'] = {}

    def __init_subclass__(cls) -> None:
        # store an instance of each chils class in _registry
        parser = cls.__name__.replace("_Parser", '')
        cls._registry[parser] = cls()

    def parse_file(self, outfile: str) -> dict:
        raise NotImplementedError("Parse file not implemented on the "
                                  "base class")

    def load(self, filepath: str, out_suffix: str = '.out',
             assignment_suffix: str = ".ass", 
             stochiometry_suffix: str = ".stoch") -> MoleculeList:
        """
        Loads the external data in the given folder.

        Parameters
        ----------
        filepath : str
            Path to the external data.
        out_suffix : str, optional
            Files with the given suffix will be treated as output files
            and imported (default: '.out').
        assignment_suffix : str, optional
            Files with the given suffix will be treated as assignment files
            (default: '.ass').
        stochiometry_suffix : str, optional
            Files with the given suffix will be treated as stochiometry files
            (default: '.stoch')
        """
        outfiles = self._fetch_all_outfiles(filepath, out_suffix)
        stochfiles = self._fetch_all_outfiles(filepath, stochiometry_suffix)
        metadata_extractors = OutFile()._registry

        data = MoleculeList()
        for outf in outfiles:
            suite = self._suite_from_outfile(outf)
            extractor = metadata_extractors.get(suite, None)
            if extractor is None:
                log.critical("No extractor available to find a suitable Parser"
                             f" for suite {suite}.", "External Parser")
            # contains all the relevant metadata, suite, method, etc
            # to determine which parser is needed
            parser = extractor.find_parser(outf)
            parser = self._registry.get(parser, None)
            if parser is None:
                log.critical(f"No parser available for parsing {parser}.", "External Parser")
            # load the file and add assignments if available
            mol = self._load_file(
                outfile=outf, out_parser=parser.parse_file,
                assignment_parser=parse_assignment_file,
                assignment_suffix=assignment_suffix
            )
            data.append(mol)
        
        if len(stochfiles) == 0:
            return data
        # Set of files which have stochiometry attached
        corrected_data = MoleculeList()
        correction_set = set()
        for stochfile in stochfiles:
            stoch_dict = json.load(open(stochfile, "r"))
            keylist = [k[:k.rfind(".")] for k in stoch_dict.keys()]
            correction_set.update(set(keylist))
        
        for molecule in data:
            data_id = molecule.data_id[:molecule.data_id.rfind(".")]
            if data_id not in correction_set:
                corrected_data.append(molecule)

        for stochfile in stochfiles:
            stoch_dict = json.load(open(stochfile, "r"))
            global_name: str = None
            system_data = dict()
            state_data: dict = dict()

            current_state_idx: int = 0
            relative_state_data: dict = dict()

            for stochidx, stochkey in enumerate(stoch_dict):
                if global_name is None:
                    global_name = "_".join(stochkey.split("_")[:-1])
                for mol in data:
                    molname = mol.data_id.replace(".out", ".in")
                    if molname == stochkey:
                        active_mol: Molecule = mol
                        break

                if active_mol is None:
                    log.critical(f"Could not find an output file at {stochkey}. Maybe the path was changed?", 
                                 "External Parser")
                
                for key, value in active_mol.system_data.items():
                    if key+"_list" not in system_data:
                        system_data[key+"_list"] = [value]
                    else:
                        system_data[key+"_list"].append(value)

                relative_properties = set([v["type"] for k, v in active_mol.state_data.items()])
                for relprop in relative_properties:
                    if relprop not in relative_state_data:
                        relative_state_data[relprop] = dict()
                for _, active_state in active_mol.state_data.items():
                    current_state_key = f"p{current_state_idx:02d}"
                    current_type = active_state["type"]
                    active_state["component index"] = stochidx
                    relative_state_data[current_type][current_state_key] = stoch_dict[stochkey]
                    state_data[current_state_key] = dict(active_state)
                    current_state_idx += 1
            
            for relprop, stochiometry in relative_state_data.items():
                current_state_key = f"p{current_state_idx:02d}"
                state_data[current_state_key] = {"type": relprop, "stochiometry": stochiometry}
                current_state_idx += 1

            pretty_name = os.path.basename(global_name).split("_")[0]
            mol: Molecule = Molecule(pretty_name, stochkey, system_data, state_data, current_state_idx)
            corrected_data.append(mol)

        return corrected_data

    def _load_file(self, outfile: str, out_parser: callable,
                   assignment_parser: callable,
                   assignment_suffix: str = ".ass") -> Molecule:
        """
        Parses and imports the given outfile, and if possible finds the
        corresponding assignment file and adds the assignments of states.

        Parameters
        ----------
        outfile : str
            The output file to import.
        out_parser : callable
            Parser that imports the outfile to a dictionary.
        assignment_parser : callable
            Parser that imports the assignmentfile to a dictionary.
        assignment_suffix : str, optional
            The suffix of the assignment file (default: '.ass').
            filename.out -> filename.ass
        """
        data: dict = out_parser(outfile)  # read in as dict
        mol = Molecule.from_external(data, outfile)
        ass_file = self._assignmentfile_from_outfile(outfile,
                                                     assignment_suffix)
        if ass_file is not None:  # assignment file exists -> add assignments
            assignments = assignment_parser(ass_file)
            mol.add_assignments(assignments)
        return mol

    def _fetch_all_outfiles(self, path: str, suffix: str = '.out') -> list:
        outfiles = []
        for root, _, files in os.walk(os.path.abspath(path), topdown=True,
                                      followlinks=True):
            for f in files:
                if f.endswith(suffix):
                    fp = os.path.abspath(os.path.join(root, f))
                    outfiles.append(fp)
        return outfiles

    def _assignmentfile_from_outfile(self, outfile: str,
                                     assignment_suffix: str = ".ass"
                                     ) -> str | None:
        """
        Finds the assignmentfile for a given output file. Returns None if no
        assignment file could be found.
        """
        assignmentf = Path(outfile).with_suffix(assignment_suffix)
        return assignmentf if assignmentf.is_file() else None

    def _suite_from_outfile(self, outfile: str) -> str:
        # since we have to open the file multiple times don't use IO here!

        # check JSON
        try:
            json.load(open(outfile, 'r'))
            return "JSON"
        except json.JSONDecodeError:
            pass
        with open(outfile, 'r') as f:
            # check QCHEM
            if any("Welcome to Q-Chem" in line for line in f):
                return "QChem"
        log.critical(f"Could not determine a suite for outfile {outfile}.",
                     "External Parser")


class QChem_MP2_Parser(ExternalParser):

    def parse_file(self, outfile: str) -> dict:
        outname = os.path.basename(outfile)
        # Signature must contain molkey, basis, method
        metadata = {"name": outname.split("_")[0]}
        data = {}
        # Flags
        dip_flag = False
        mulliken_flag = False
        mul_charges = []

        with open(outfile, "r") as outf:
            for line in outf:
                if line.strip() == "":
                    continue
                lsplit = line.strip().split()
                lnw = " ".join(lsplit)
                if "basis" not in metadata:
                    if len(lsplit) >= 2 and lsplit[0].lower() == "basis":
                        b = lsplit[1]
                        if lsplit[1] == "=":
                            b = lsplit[2]
                        data["basis"] = b
                if "method" not in metadata:
                    if len(lsplit) >= 2 and lsplit[0].lower() == "method":
                        m = lsplit[1]
                        if lsplit[1] == "=":
                            m = lsplit[2]
                        data["method"] = m
                if "energy" not in data:
                    if "MP2 total energy =" in lnw:
                        data["energy"] = float(lsplit[-2])
                if "dipole moment" not in data:
                    if not dip_flag and "Dipole Moment (Debye)" in lnw:
                        dip_flag = True
                        continue
                    elif dip_flag:
                        dx, dy, dz = lsplit[1], lsplit[3], lsplit[5]
                        dip_mom = numpy.array(
                            [float(dx), float(dy), float(dz)]
                        )
                        # Debye to au
                        dip_mom /= 0.393430307
                        data["dipole moment"] = tuple(dip_mom)
                        data["total dipole moment"] = (
                            numpy.linalg.norm(dip_mom)
                        )
                if "mulliken charges" not in data:
                    if "Ground-State Mulliken Net Atomic Charges" in lnw \
                            and not mulliken_flag:
                        mulliken_flag = True
                    elif mulliken_flag:
                        if lsplit[0] == "Atom":
                            continue
                        if "------------------------" in lsplit[0] and \
                                len(mul_charges) == 0:
                            continue
                        if "------------------------" in lsplit[0] and \
                                len(mul_charges) > 0:
                            data["mulliken charges"] = tuple(mul_charges)
                        else:
                            mul_charges.append(float(lsplit[-1]))
        metadata["data"] = data
        return {"s0": metadata}
    
class QChem_RIBWS2_Parser(ExternalParser):

    def parse_file(self, outfile: str) -> dict:
        outname = os.path.basename(outfile)
        # Signature must contain molkey, basis, method
        metadata = {"name": outname.split("_")[0]}
        data: dict = dict()
        with open(outfile, "r") as outf:
            for line in outf.readlines():
                lsplit = line.strip().split()
                if "RIBW-S2         total energy" in line:
                    data["energy"] = float(lsplit[-2])
                if "basis" not in metadata:
                    if len(lsplit) >= 2 and lsplit[0].lower() == "basis":
                        b = lsplit[1]
                        if lsplit[1] == "=":
                            b = lsplit[2]
                        metadata["basis"] = b
                if "method" not in metadata:
                    if len(lsplit) >= 2 and lsplit[0].lower() == "method":
                        m = lsplit[1]
                        if lsplit[1] == "=":
                            m = lsplit[2]
                        metadata["method"] = m
                if "scf energy" not in metadata:
                    if "SCF   energy =" in line:
                        data["scf energy"] = float(lsplit[-1])
        metadata["data"] = data
        return {"s0": metadata}

class QChem_RICC2_Parser(ExternalParser):

    def parse_file(self, outfile: str) -> dict:
        outname = os.path.basename(outfile)
        # Signature must contain molkey, basis, method
        metadata = {"name": outname.split("_")[0]}
        data: dict = dict()
        with open(outfile, "r") as outf:
            for line in outf.readlines():
                lsplit = line.strip().split()
                if "RICC2         total energy" in line:
                    data["energy"] = float(lsplit[-2])
                if "basis" not in metadata:
                    if len(lsplit) >= 2 and lsplit[0].lower() == "basis":
                        b = lsplit[1]
                        if lsplit[1] == "=":
                            b = lsplit[2]
                        metadata["basis"] = b
                if "method" not in metadata:
                    if len(lsplit) >= 2 and lsplit[0].lower() == "method":
                        m = lsplit[1]
                        if lsplit[1] == "=":
                            m = lsplit[2]
                        metadata["method"] = m
                if "scf energy" not in metadata:
                    if "SCF   energy =" in line:
                        data["scf energy"] = float(lsplit[-1])
        if metadata["method"] == "ribws2":
            metadata["method"] = "RIBWS-CC2"
        metadata["data"] = data
        return {"s0": metadata} 

class JSON_Parser(ExternalParser):
    def parse_file(self, outfile: str) -> dict:
        return json.load(open(outfile, "r"))

class OutFile:
    """
    Parent class for obtaining metadata from outputfiles.

    This class is used to extract metadate like the used method and basis set
    from an output file. Any class to extract data from a specific output
    file must inherit from this class. The naming convention for subclasses is
    [Suite]OutFile (e.g. QChemOutFile, ORCAOutFile).

    Methods
    -------
    find_parser(outfile: str) -> str
        Gather metadata from the given output file to identify a suitable
        parser for the file.

    Attributes
    ----------
    None
    """

    _registry: dict[str, 'OutFile'] = {}

    def __init_subclass__(cls) -> None:
        # store an instance of each child class in _registry
        suite = cls.__name__.replace("OutFile", "")
        cls._registry[suite] = cls()

    def find_parser(self, outfile: str) -> str:
        pass


class QChemOutFile(OutFile):
    def find_parser(self, outfile: str) -> str:
        parser = None
        with open(outfile, "r") as f:
            for line in f.readlines():
                if "RICC2         total energy" in line:
                    parser = "QChem_RICC2"
                if "RIBW-S2         total energy" in line:
                    parser = "QChem_RIBWS2"
        return parser


class JSONOutFile(OutFile):
    def find_parser(self, outfile: str) -> str:
        return "JSON"
