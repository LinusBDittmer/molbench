"""Python file for all external parsers.

"""
from . import logger as log
from .molecule import Molecule
import os.path
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
        raise NotImplementedError()

    def load(self, filepath: str, suffix: str = 'out') -> dict:
        outfiles = self._fetch_all_outfiles(filepath, suffix)
        metadata_extractors = OutFile()._registry

        data = []
        for outf in outfiles:
            suite = self._suite_from_outfile(outf)
            extractor = metadata_extractors.get(suite, None)
            if extractor is None:
                log.critical("No extractor available to find a suitable Parser"
                             f" for suite {suite}.")
            # contains all the relevant metadata, suite, method, etc
            # to determine which parser is needed
            parser = extractor.find_parser(outf)
            parser = self._registry.get(parser, None)
            if parser is None:
                log.critical(f"No parser available for parsing {parser}.")
            data.append(Molecule.from_external(parser.parse_file(outf),
                                               outf))
        return data

    def _fetch_all_outfiles(self, path: str, suffix: str = 'out') -> list:
        outfiles = []

        for root, _, files in os.walk(path):
            for f in files:
                if f.endswith(suffix):
                    fp = os.path.abspath(os.path.join(root, f))
                    outfiles.append(fp)
        return outfiles

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
                     self)


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

    def load(self, filepath: str, suffix: str = 'out') -> dict:
        outfiles = self._fetch_all_outfiles(filepath, suffix)

        data = []
        for outfile in outfiles:
            # filter out non qchem MP2 outfiles
            if self._suite_from_outfile(outfile) != "QChem" or \
                    QChemOutFile().find_parser(outfile) != "QChem_MP2":
                continue
            data.append(Molecule.from_external(self.parse_file(outfile),
                                               outfile))
        return data


class JSON_Parser(ExternalParser):
    def parse_file(self, outfile: str) -> dict:
        return json.load(open(outfile, "r"))

    def load(self, filepath: str, suffix: str = 'out') -> dict:
        outfiles = self._fetch_all_outfiles(filepath, suffix)
        return [Molecule.from_external(self.parse_file(outf), outf)
                for outf in outfiles
                if self._suite_from_outfile(outf) == "JSON"]


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
        raise NotImplementedError


class JSONOutFile(OutFile):
    def find_parser(self, outfile: str) -> str:
        return "JSON"
