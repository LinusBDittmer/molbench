"""
Molbench is a Python package used for benchmarking quantum chemical
applications, methods, and suites. See the class docstrings of
JSONBenchmarkParser, ExternalParser, Comparison, Statistics,
TemplateConstructor, and LatexExporter for the package's API.
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
from .json_encoder import MolbenchJSONEncoder

__all__ = ["config", "Molecule", "InputConstructor",
           "TemplateConstructor", "CompressedTemplateConstructor",
           "create_bash_files", "make_send_script",
           "JSONBenchmarkParser",
           "Comparison", "Statistics", "register_as_error_measure",
           "Exporter", "LatexExporter",
           "ExternalParser", "MolbenchJSONEncoder"]
__version__ = "0.0.1"
__authors__ = ["Linus Bjarne Dittmer", "Jonas Leitner"]
