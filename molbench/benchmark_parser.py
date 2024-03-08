"""Python script for handling benchmarks.


"""

import os
import json
import molbench.logger as log

premade_benchmarks = None


def _collect_premade_benchmarks():
    """
    Collect premade benchmark files from the 'benchmarks' directory.

    Explanation
    -----------
    - `premade_benchmarks`: Global dictionary to store premade benchmark files.
    - If `premade_benchmarks` is already populated, exit the function.
    - `wpath`: Get the directory path of the current file.
    - `rpath`: Create the absolute path to the 'benchmarks' directory relative to the current file.
    - Initialize `premade_benchmarks` as an empty dictionary.
    - Iterate over files in the 'benchmarks' directory:
        - If the file has a '.json' extension, add it to `premade_benchmarks` with the filename 
          (without extension) as key and the absolute path as value.

    """
    global premade_benchmarks
    if premade_benchmarks is not None:
        return

    wpath = os.path.dirname(os.path.abspath(__file__))
    rpath = os.path.join(wpath, "benchmarks")
    premade_benchmarks = {}
    for root, _, files in os.walk(rpath):
        for file in files:
            if file.endswith('.json'):
                key = os.path.splitext(file)[0]
                val = os.path.abspath(os.path.join(root, file))
                premade_benchmarks.update({key: val})

class BenchmarkParser:

    def __init__(self):
        _collect_premade_benchmarks()

    def parse(self, benchmark: str, use_local_benchmark: bool = False) -> dict:
        global premade_benchmarks
        if not use_local_benchmark:
            _collect_premade_benchmarks()
        return self._parse_benchmark(benchmark, use_local_benchmark)

    def _load_benchmark(self, benchmark: str, use_local_benchmark: bool = False) -> dict:
        raise NotImplementedError("The load function is not implemented in the BenchmarkParser
                                   superclass. Please use a child class instead.")

class JSONBenchmarkParser(BenchmarkParser):

    def __init__(self):
        super().__init__()

    def parse_benchmark(benchmark: str, use_local_benchmark: bool = False) -> dict:
        """
        Load a benchmark file.

        Parameters
        ----------
        benchmark : str
            Path to the benchmark file or name of a premade benchmark.
        use_local_benchmark : bool, optional
            If True, `benchmark` is interpreted as a local path, otherwise as the name of a premade 
            benchmark. Default is False.

        Returns
        -------
        dict
            Dictionary containing the loaded benchmark data.

        Explanation
        -----------
        - `premade_benchmarks`: Global dictionary storing premade benchmark files.
        - Call `_collect_premade_benchmarks` to populate `premade_benchmarks`.
        - If `premade_benchmarks` is populated and `benchmark` is in `premade_benchmarks`, use the 
          corresponding premade benchmark.
        - Otherwise, if `benchmark` is not a local path or doesn't exist, log a critical error.
        - If `benchmark` doesn't end with '.json', log a warning.
        - Attempt to load the benchmark file as JSON.
        - If an error occurs during JSON decoding, log a critical error.

        """
        global premade_benchmarks

        # If the benchmark is premade, use the premade benchmark
        # Otherwise, interpret as a path
        if not use_local_benchmark and premade_benchmarks is not None and \
                benchmark in premade_benchmarks:
            benchmark = premade_benchmarks[benchmark]
        else:
            if not os.path.exists(benchmark):
                log.critical(f"Benchmark file {benchmark} does not exists or "
                             "cannot be seen.", "load_benchmark")
            if not benchmark.endswith(".json"):
                log.warning(f"Benchmark file {benchmark} is not labelled as JSON. "
                            "Loading may nonetheless be possible.",
                            "load_benchmark")

        try:
            with open(benchmark, "r") as f:
                bm_dict = json.load(f)
        except json.JSONDecodeError:
            log.critical(f"Benchmark file {benchmark} cannot be read.")

        return bm_dict
