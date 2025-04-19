"""Parsing class for external data
"""

from . import logger as log
from .molecule import Molecule, MoleculeList
from .assignment import parse_assignment_file
import os.path
import inspect
from pathlib import Path


class ExternalParser:

    def load(self, filepath: str, parser: callable,
             out_suffix: str = ".out",
             assignment_suffix: str = ".ass"):
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
        """
        outfiles = self._fetch_all_outfiles(filepath, out_suffix)

        param_num: int = len(inspect.signature(parser).parameters)

        if param_num not in (1, 2):
            log.critical("Parser callable can only have the following "
                         + "arguments:\n\noutput_filepath : str\n    "
                         + "Path to the output file\nname : str (optional)"
                         + "\n    Name for the molecule object\n\nThe given"
                         + f" callable has {param_num} parameters.",
                         "ExternalParser")

        data = MoleculeList()
        for outf in outfiles:
            # contains all the relevant metadata, suite, method, etc
            # to determine which parser is needed
            if parser is None:
                log.critical(f"No parser available for parsing {parser}.",
                             "External Parser")
            # load the file and add assignments if available
            mol = self._load_file(
                outfile=outf, out_parser=parser,
                assignment_parser=parse_assignment_file,
                assignment_suffix=assignment_suffix,
                requires_name=(param_num == 2)
            )
            data.append(mol)

        return data

    def _load_file(self, outfile: str, out_parser: callable,
                   assignment_parser: callable,
                   assignment_suffix: str = ".ass",
                   requires_name: bool = False) -> Molecule:
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
        # XXX: Why the split?
        # Above in the error message you say that the callable may
        # take an optional name parameter. But you did not mention
        # any suffix separated by '_'
        name: str = Path(outfile).stem.split("_")[0]
        if requires_name:
            data: dict = out_parser(outfile, name)  # read in as dict
        else:
            data: dict = out_parser(outfile)
        mol = Molecule.from_external(data, outfile)
        ass_file = self._assignmentfile_from_outfile(outfile,
                                                     assignment_suffix)
        if ass_file is not None:  # assignment file exists -> add assignments
            assignments = assignment_parser(ass_file)
            mol.add_assignments(assignments)
        return mol

    def _fetch_all_outfiles(self, path: str,
                            suffix: str = '.out') -> list[str]:
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
        return str(assignmentf) if assignmentf.is_file() else None
