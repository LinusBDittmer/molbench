"""  # noqa E501
Molbench Initialization File

This file initializes the Molbench Package. Molbench is a Python package used
for benchmarking quantum chemical applications, methods, and suites.

Methods
-------
load_benchmark(benchmark: str) -> dict
    Load a benchmark identified by the provided benchmark string or filepath.

# RETRACTED
load_external(filepath: str, parser: ExternalParser) -> dict
    Load data from an external source, interpreting it using the provided
    parser.

export_xyz(benchmark: dict, filepath: str, flat_structure: bool = False,
           name_template: str = '[[name]].xyz') -> None
    Export xyz files for the specified benchmark to the directory specified at
    filepath.

export_input(benchmark: dict, filepath: str, constructor: InputConstructor,
             flat_structure: bool = False,
             name_template: str = '[[name]]_[[method]]_[[basis]].in') -> None
    Write input files for the given benchmark to the directory located at
    filepath using the provided constructor.

export_comparison(benchmark: dict, external_data: dict, properties: tuple,
                  filepath: str,
                  comparator: Comparator = CsvComparator) -> None
    Compare the given external data with the benchmark and write the
    comparison results to a file.

Classes
-------
ExternalParser
    Parent class for an Output-API. Used to load external data from a
    directory.

InputConstructor
    Parent class for an Inputfile-API. Used to write input files for a
    benchmark.

Comparator
    Parent class for a comparison between external data and a benchmark set.

CsvComparator(Comparator)
    Subclass of Comparator. Creates a csv table for comparison.

Functions
---------
ExternalParser.load(filepath: str, suffix: str = 'out') -> dict
    Load external data from the specified filepath with the given suffix.

InputConstructor.create(benchmark: dict, filepath: str,
                        flat_structure: bool = False,
                        name_template: str = '[[name]]_[[method]]_[[basis]].in') -> None
    Create input files for a benchmark at the specified filepath.

Comparator.compare(benchmark: dict, external_data: dict, properties: tuple) -> str
    Compare benchmark data with external data and return file contents of the comparison.

CsvComparator.compare(benchmark: dict, external_data: dict, properties: tuple) -> str
    Compare benchmark data with external data and return file contents of the
    comparison in CSV format.
"""

from .configuration import config
from .benchmark_parser import JSONBenchmarkParser
from .input_constructor import (InputConstructor, TemplateConstructor,
                                CompressedTemplateConstructor)
from .bash_wrapper import create_bash_files, make_send_script
from .comparison import Comparison
from .statistics import Statistics, register_as_error_measure
from .export import Exporter, LatexExporter
from .external_parser import ExternalParser
from .molecule import Molecule
# from .extrapolate import CBSExtrapolator

__all__ = ["config", "Molecule", "load_benchmark", "InputConstructor",
           "TemplateConstructor", "CompressedTemplateConstructor",
           "create_bash_files", "make_send_script",
           "JSONBenchmarkParser",
           "Comparison", "Statistics", "register_as_error_measure",
           "Exporter", "LatexExporter",
           "ExternalParser"]
__version__ = "0.0.1"
__authors__ = ["Linus Bjarne Dittmer", "Jonas Leitner"]
