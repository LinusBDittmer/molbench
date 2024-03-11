"""Python script for handling benchmarks.


"""

from . import logger as log
from .molecule import Molecule, MoleculeList
import os
import json


class BenchmarkParser:
    premade_benchmarks = None

    def __init__(self):
        self._collect_premade_benchmarks()

    @classmethod
    def _collect_premade_benchmarks(cls) -> None:
        """
        Collect premade benchmark files from the 'benchmarks' directory.

        Explanation
        -----------
        - `premade_benchmarks`: Global dictionary to store premade benchmark
           files.
        - If `premade_benchmarks` is already populated, exit the function.
        - `wpath`: Get the directory path of the current file.
        - `rpath`: Create the absolute path to the 'benchmarks' directory
           relative to the current file.
        - Initialize `premade_benchmarks` as an empty dictionary.
        - Iterate over files in the 'benchmarks' directory:
            - If the file has a '.json' extension, add it to
              `premade_benchmarks` with the filename
              (without extension) as key and the absolute path as value.
        """
        if cls.premade_benchmarks is not None:
            return
        wpath = os.path.dirname(os.path.abspath(__file__))
        rpath = os.path.join(wpath, "benchmarks")
        cls.premade_benchmarks = {}
        for root, _, files in os.walk(rpath):
            for file in files:
                if file.endswith('.json'):
                    key = os.path.splitext(file)[0]
                    val = os.path.abspath(os.path.join(root, file))
                    cls.premade_benchmarks.update({key: val})

    def load(self, benchmark: str, benchmark_id=None,
             use_local_benchmark: bool = False) -> list[Molecule]:
        if not use_local_benchmark and benchmark in self.premade_benchmarks:
            benchfile = self.premade_benchmarks[benchmark]
        else:
            if not os.path.exists(benchmark):
                log.critical(f"Benchmark file {benchmark} does not exists or "
                             "cannot be seen.", "load_benchmark")
            benchfile = benchmark
        content = self.parse_benchmark(benchfile)
        if benchmark_id is None:
            benchmark_id = benchmark
        return MoleculeList(
            Molecule.from_benchmark(moldata, benchmark_id, molkey)
            for molkey, moldata in content.items()
        )

    def parse_benchmark(self, benchmarkfile: str) -> dict:
        raise NotImplementedError("The parse_benchmark function is not "
                                  "implemented in the BenchmarkParser "
                                  "superclass. Please use a child class "
                                  "instead.")


class JSONBenchmarkParser(BenchmarkParser):
    def parse_benchmark(self, benchmarkfile: str) -> dict:
        """
        Parse a benchmark file.

        Parameters
        ----------
        benchmarkfile : str
            Path to the benchmark file

        Returns
        -------
        dict
            Dictionary containing the loaded benchmark data.
        """
        try:
            return json.load(open(benchmarkfile, "r"))
        except json.JSONDecodeError:
            log.critical(f"Could not read benchmark file {benchmarkfile}.")
