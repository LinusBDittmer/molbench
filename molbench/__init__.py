"""
Molbench is a Python package used for benchmarking quantum chemical
applications, methods, and suites. See the class docstrings of
JSONBenchmarkParser, ExternalParser, Comparison, Statistics,
TemplateConstructor, and LatexExporter for the package's API.
"""

from .bash_wrapper import create_bash_files, make_send_script
from .benchmark_parser import JSONBenchmarkParser
from .comparison import Comparison
from .configuration import config
from .export import Exporter, LatexExporter
from .external_parser import ExternalParser
from .input_constructor import (
    CompressedTemplateConstructor,
    InputConstructor,
    TemplateConstructor,
)
from .json_encoder import MolbenchJSONEncoder
from .molecule import Molecule
from .statistics import Statistics, register_as_error_measure

__all__ = [
    "Comparison",
    "CompressedTemplateConstructor",
    "Exporter",
    "ExternalParser",
    "InputConstructor",
    "JSONBenchmarkParser",
    "LatexExporter",
    "MolbenchJSONEncoder",
    "Molecule",
    "Statistics",
    "TemplateConstructor",
    "config",
    "create_bash_files",
    "make_send_script",
    "register_as_error_measure",
]
__version__ = "0.0.1"
__authors__ = ["Linus Bjarne Dittmer", "Jonas Leitner"]
