"""Python script for handling benchmarks.


"""

from . import logger as log
from .molecule import Molecule
import os
import json

premade_benchmarks = None


def _collect_premade_benchmarks():
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


def load_benchmark(benchmark: str, benchmark_id=None,
                   use_local_benchmark: bool = False) -> dict:
    _collect_premade_benchmarks()
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

    if benchmark_id is None:
        benchmark_id = benchmark
    return [Molecule.from_benchmark(moldata, benchmark_id, molkey)
            for molkey, moldata in bm_dict.items()]
