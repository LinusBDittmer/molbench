# molbench

A Python package for benchmarking quantum chemical applications. molbench loads pre-built benchmark datasets, compares them against your computed results, evaluates statistical error metrics, generates input files for quantum chemistry codes, and exports results as LaTeX tables.

**Authors:** Linus Bjarne Dittmer, Jonas Leitner  
**License:** MIT  
**Python:** 3.12+

---

## Table of Contents

- [Installation](#installation)
- [Included Benchmarks](#included-benchmarks)
- [Core Concepts](#core-concepts)
- [Comparison Workflow](#comparison-workflow)
  - [1. Load a Benchmark](#1-load-a-benchmark)
  - [2. Load Computed Results](#2-load-computed-results)
  - [3. Build a Comparison](#3-build-a-comparison)
  - [4. Compute Statistics](#4-compute-statistics)
  - [5. Export to LaTeX](#5-export-to-latex)
- [Input Generation Workflow](#input-generation-workflow)
  - [Generating Input Files](#generating-input-files)
  - [Generating Assignment Files](#generating-assignment-files)
  - [Relative-Energy Benchmarks (CompressedTemplateConstructor)](#relative-energy-benchmarks-compressedtemplateconstructor)
- [Assignment Files](#assignment-files)
- [HPC: Bash Scripts and Send Scripts](#hpc-bash-scripts-and-send-scripts)
- [Filtering and Slicing MoleculeLists](#filtering-and-slicing-moleculelists)
- [Custom Error Measures](#custom-error-measures)
- [Configuration](#configuration)
- [Data Structures](#data-structures)
  - [Molecule and Datapoint](#molecule-and-datapoint)
  - [Comparison](#comparison)
- [Benchmark JSON Schema](#benchmark-json-schema)
- [Development](#development)

---

## Installation

```bash
# Development install (recommended)
pip install -e .

# or via Poetry
poetry install
```

**Dependencies:** `numpy`, `matplotlib`, `argparse`

---

## Included Benchmarks

| Name | Description | DOI |
|---|---|---|
| `ascdb` | Adiabatic and vertical singlet/triplet energies | 10.1039/C9CP03211H |
| `questdb` | QUEST excitation energy database (TBE/aug-cc-pVTZ) | — |
| `jacquemin18` | Benchmark for excitation energies (Jacquemin 2018) | 10.1021/acs.jctc.8b00406 |
| `jacquemin21` | Extended excitation benchmark (Jacquemin 2021) | 10.1021/acs.jctc.0c01228 |
| `doubles_jacquemin19` | Double excitations (Jacquemin 2019) | 10.1021/acs.jctc.8b01205 |
| `esa_jacquemin25` | Excited state absorption (Jacquemin 2025) | — |

All benchmarks are stored as JSON files inside the package and require no external download.

---

## Core Concepts

molbench organises data around two workflows:

```
Comparison workflow:
  JSONBenchmarkParser.load()   →  MoleculeList (reference)
  ExternalParser.load()        →  MoleculeList (computed)
  Comparison.add(both)         →  nested dict [molecule][basis][method][proptype][data_id]
  Statistics.compare(...)      →  signed errors
  Statistics.evaluate(...)     →  MAE, RMSD, MSE, ...

Input generation workflow:
  JSONBenchmarkParser.load()             →  MoleculeList
  TemplateConstructor.create_inputs()    →  input files (.in)
  TemplateConstructor.create_assignments() →  assignment files (.ass)
  (user runs calculations)
  ExternalParser.load(...)               →  back into comparison workflow
```

---

## Comparison Workflow

### 1. Load a Benchmark

```python
from molbench import JSONBenchmarkParser

# Load a built-in benchmark; benchmark_id is stored as data_id on each molecule
bench = JSONBenchmarkParser().load("questdb", benchmark_id="TBE")
print(f"Loaded {len(bench)} molecules")

# Load a custom JSON file in the same schema
bench = JSONBenchmarkParser().load("/path/to/my_benchmark.json", benchmark_id="ref")
```

### 2. Load Computed Results

You supply a **parser callable** that reads one output file and returns
`(name, system_data, state_data)`. The parser may take 1 or 2 arguments;
when 2 parameters are declared, the second receives the file's stem as a
default name hint.

```python
import json
from molbench import ExternalParser


def my_parser(filepath):
    raw = json.load(open(filepath))
    name = raw["name"]
    system_data = {}
    state_data = {
        "gs": {
            "basis": raw["basis"],
            "method": raw["method"],
            "data": {"energy": {"value": raw["energy"], "unit": "au"}},
        }
    }
    return name, system_data, state_data


computed = ExternalParser().load(
    "/path/to/output/dir",
    parser=my_parser,
    out_suffix=".out",  # scan for files with this extension
    assignment_suffix=".ass",  # companion assignment files (optional)
)
```

For excitation spectra benchmarks you will typically also need
**assignment files** — see [Assignment Files](#assignment-files) below.

### 3. Build a Comparison

```python
from molbench import Comparison

c = Comparison()  # default separators: ("basis", "method")
c.add(bench)
c.add(computed)

# For excitation-energy benchmarks add a transition_id separator:
c = Comparison("basis", "method", "transition_id")
c.add(bench)
c.add(computed)
```

### 4. Compute Statistics

```python
from molbench import Statistics

stats = Statistics(c)

# compare() returns {ref_key: {interest_key: signed_error (Datapoint)}}
errors = stats.compare(
    interest={"method": "my_method"},
    reference={"method": "TBE"},
)

# evaluate() accepts one or more metric names, or "all"
result = stats.evaluate(errors, "mae", "rmsd", "mse", proptype="energy")
# result = {"mae": (value, count), "rmsd": (value, count), ...}

mae, n = result["mae"]
print(f"MAE = {mae:.4f}  (n = {n})")

# Relative errors with optional damping to avoid division by small values
errors_rel = stats.compare(
    interest={"method": "my_method"},
    reference={"method": "TBE"},
    relative=True,
    relative_damping=0.01,  # adds 0.01 to |ref| in denominator
    error_thresh=0.5,  # warn if |error| > threshold
)
```

Available metrics: `"mse"`, `"mae"`, `"rmsd"`, `"sde"`, `"min"`, `"max"`,
`"median_se"`, `"all"`.

### 5. Export to LaTeX

```python
from molbench import LatexExporter
from molbench.formatting import LatexFormatter
from molbench.tree import Node, DummyNode

fmt = LatexFormatter(n_decimals=3)
exporter = LatexExporter(fmt, sort_cols=True, sort_rows=True)

# Build a column tree: one column per method
columns = DummyNode()
Node("my_method", columns)
Node("other_method", columns)

with open("table.tex", "w") as f:
    exporter.export(c, "energy", f, columns=columns)
```

---

## Input Generation Workflow

### Generating Input Files

`TemplateConstructor` reads a template file where placeholders written as
`[[key]]` are substituted from `system_data`, `state_data`, and a
user-supplied `calc_details` dict. Keys from the global
[configuration](#configuration) are used as fallback.

```python
from molbench import JSONBenchmarkParser, TemplateConstructor

bench = JSONBenchmarkParser().load("ascdb", benchmark_id="TBE")

# Use a built-in template (see molbench/templates/)
tc = TemplateConstructor("pyscf_ordmp2")

calc = {
    "method": "hf",
    "basis": "cc-pvdz",
    "verbose": 0,
    "xyz_unit": "A",
    "symmetry": False,
    "spin_unrestricted": False,
    "conv_tol": 1e-9,
    "conv_tol_grad": 1e-6,
    "scf_max_cycle": 100,
}

# Creates one subdirectory per molecule, one .in file per basis/method combo
files = tc.create_inputs(bench, "/path/to/inputs", calc)

# Flat layout: all files in one directory
files = tc.create_inputs(bench, "/path/to/inputs", calc, flat_structure=True)
```

**Built-in templates:**

| Name | Description |
|---|---|
| `pyscf_ordmp2` | PySCF OO-RDMP2 ground state |
| `adcc_mp2` | adcc MP2 |
| `adcc_re2` | adcc RE-ADC(2) excitations |
| `adcc_re-adc` | adcc RE-ADC excitations |
| `mrcc_cc` | MRCC coupled-cluster |
| `qchem_basic` | Q-Chem basic single point |
| `qchem_ribws2` | Q-Chem RI-BWS2 |

**Custom templates** can be plain text files; pass the file path instead of
the built-in name:

```
# my_template.txt
molecule {
  [[xyz]]
  charge [[charge]]
  multiplicity [[multiplicity]]
}

method = [[method]]
basis  = [[basis]]
```

```python
tc = TemplateConstructor("/path/to/my_template.txt")
```

**List-valued keys** (e.g. `xyz_list`, `charge_list` from multi-geometry
molecules) automatically expand the template: one file is generated per
list element. Use `[[key->0]]` to index into a list directly.

### Generating Assignment Files

For benchmarks that contain excited states you need to tell molbench which
computed state corresponds to which reference transition. Generate skeleton
assignment files alongside the input files:

```python
tc.create_assignments(bench, "/path/to/inputs", calc)
```

This creates a `.ass` file next to each `.in` file, pre-filled with the
reference transition IDs on the left. After running your calculations,
fill in the external state IDs on the right — see
[Assignment Files](#assignment-files).

### Relative-Energy Benchmarks (CompressedTemplateConstructor)

For benchmarks like `ascdb` that store relative energies (multi-geometry
molecules), use `CompressedTemplateConstructor`. It deduplicates geometries,
creates one input file per unique geometry, and writes a `references.json`
that maps compressed molecule names back to the original benchmark entries
for post-processing.

```python
from molbench import CompressedTemplateConstructor

tc = CompressedTemplateConstructor("pyscf_ordmp2")
tc.create_inputs(bench, "/path/to/inputs", calc, reference_path="references.json")
```

---

## Assignment Files

When comparing computed excitation spectra against a reference benchmark,
the state ordering from your code may differ from the benchmark. Assignment
files (`.ass`) resolve this mapping.

**Format** — one line per state, `ref_id ==> external_id`:

```
# ref_state_id ==> external_id
s0->s1 ==> null
s0->t1 ==> null
```

After running your calculations, replace `null` with the state label your
code produced (e.g. the root index printed in the output):

```
# ref_state_id ==> external_id
s0->s1 ==> state_001
s0->t1 ==> state_002
```

molbench looks for a `.ass` file with the same stem as the output file and
applies the mapping automatically when loading via `ExternalParser`. States
left as `null` (unassigned) are dropped from the loaded molecule.

The resulting `assigned_transition_id` field is used by `Comparison` in
place of `transition_id`, so the assignment is transparent to the rest of
the workflow.

**Programmatic creation/parsing:**

```python
from molbench.assignment import new_assignment_file, parse_assignment_file

# Create a skeleton file for a list of reference transition IDs
content = new_assignment_file(["s0->s1", "s0->t1"])
with open("mol.ass", "w") as f:
    f.write(content)

# Parse a filled-in file → {external_id: ref_id}
mapping = parse_assignment_file("mol.ass")
# e.g. {"state_001": "s0->s1", "state_002": "s0->t1"}
```

---

## HPC: Bash Scripts and Send Scripts

After generating input files, `create_bash_files` invokes a script-generator
command in each input directory. The command is expected to produce a `.sh`
or `.sbatch` job script. `make_send_script` then writes a single master
script that `cd`s into each directory and submits all jobs.

```python
from molbench import create_bash_files, make_send_script

# Run 'gen_script.py <filename>.in' in each input directory.
# [[key]] placeholders are resolved against the global config.
bash_scripts = create_bash_files(input_files, command="python gen_script.py")

# Write a master submission script
with open("submit_all.sh", "w") as f:
    make_send_script(bash_scripts, send_command="sbatch", sendscript=f)
```

---

## Filtering and Slicing MoleculeLists

`MoleculeList` extends `list` with filtering helpers that inspect both
`system_data` and `state_data` keys:

```python
# Keep only entries where method == "TBE"
tbe = bench.filter("method", "TBE")

# Remove entries with a specific basis
filtered = bench.remove("basis", "sto-3g")

# Keep only entries where n_atoms is between 5 and 20
small = bench.filter_by_range("n_atoms", min=5, max=20)

# Slice like a regular list
subset = bench[:10]

# Combine multi-geometry entries into relative-energy molecules
rel = bench.apply_stochiometry(stochiometry_dict)
```

---

## Custom Error Measures

Register your own statistical measure with the `@register_as_error_measure`
decorator. The function receives `(signed_errors: dict, assign: Callable)`
and must return `(value, count)`:

```python
from molbench import Statistics, register_as_error_measure
from molbench.statistics import _collect_errors
import numpy as np


@register_as_error_measure
def winsorized_mae(signed_errors, assign):
    errors = _collect_errors(signed_errors, assign)
    if not errors:
        return 0.0, 0
    arr = np.clip(np.abs(errors), None, np.percentile(np.abs(errors), 90))
    return float(arr.mean()), len(errors)


result = stats.evaluate(errors, "winsorized_mae", proptype="energy")
```

---

## Configuration

molbench reads a JSON config file at startup. The path is resolved in order:

1. `$MOLBENCH_CONFIG` environment variable (absolute path to a JSON file).
2. `molbench/local_config.json` inside the package directory (default).

**Default `local_config.json`:**

```json
{
    "threads":  1,
    "walltime": "12:00:00",
    "memory":   50000,
    "queue":    "short"
}
```

All keys are available as `[[key]]` placeholders in input templates, so HPC
resource settings (threads, memory, walltime, queue) can be changed in one
place.

**Override at runtime:**

```bash
MOLBENCH_CONFIG=/path/to/hpc_config.json python my_script.py
```

**Verbosity** is controlled by `MOLBENCH_VERBOSE`:
- unset / `0` — warnings and errors only
- `1` — info messages
- `2` — debug output

---

## Data Structures

### Molecule and Datapoint

```
Molecule
├── name          str       — unique molecule identifier
├── data_id       str       — source tag (benchmark_id or file path)
├── system_data   dict      — geometry and system info
│   ├── xyz                   str    coordinates as a single string (Angstrom)
│   ├── charge                int
│   ├── multiplicity          int
│   ├── n_atoms               int
│   └── *_list variants for multi-geometry molecules
└── state_data    dict      — {state_id: state_dict}
    └── state_dict
        ├── method              str
        ├── basis               str
        ├── data                {propname: Datapoint}
        ├── transition_id       str     (excitation benchmarks)
        ├── assigned_transition_id  str (set by add_assignments after loading)
        └── stochiometry        list    (relative-energy benchmarks)

Datapoint(value, unit)
    Supports arithmetic: +, -, *, / (unit-checked for + and -)
    abs(dp) → Datapoint with abs(value)
    dp * 2.0, dp / 3.0 → scale by a dimensionless factor
```

### Comparison

`Comparison` is a nested `dict` keyed as:

```
comparison[molecule_name][basis][method][proptype][data_id] = Datapoint
```

The intermediate separators default to `("basis", "method")` but can be
customised at construction — useful for excitation databases where you also
need to separate by `transition_id`:

```python
c = Comparison("basis", "method", "transition_id")
```

Walk helpers:

```python
# All entries for a specific property across the whole comparison
for path, val_dict in c.walk_by_key("excitation_energy"):
    print(path, val_dict)

# All Datapoint leaf values
for path, dp in c.walk_values():
    print(path, dp.value, dp.unit)
```

---

## Benchmark JSON Schema

Custom benchmark files must follow this schema. For **single-geometry**
molecules:

```json
{
  "molecule_name": {
    "charge": 0,
    "multiplicity": 1,
    "n_atoms": 7,
    "xyz": ["C 0.0 0.0 0.0", "H 1.0 0.0 0.0"],
    "properties": {
      "state_id": {
        "basis":  "aug-cc-pvtz",
        "method": "TBE",
        "data": {
          "excitation_energy": {"value": 4.31, "unit": "eV"},
          "oscillator_strength": {"value": 0.01, "unit": "au"}
        },
        "transition_id": "s0->s1"
      }
    }
  }
}
```

For **multi-geometry** (relative-energy) molecules, replace the
single-valued keys with `_list` variants and add `stochiometry`:

```json
{
  "AE18pE-1": {
    "charge_list":       [0],
    "multiplicity_list": [2],
    "n_atoms_list":      [1],
    "xyz_list":          [["H 0.0 0.0 0.0"]],
    "properties": {
      "p001": {
        "basis": "CBS",
        "method": "TBE",
        "stochiometry": [1.0],
        "data": {"energy": {"value": -0.5, "unit": "au"}}
      }
    }
  }
}
```

---

## Development

```bash
# Install dev dependencies
pip install -e .
pip install pytest flake8

# Lint
flake8

# Run all fast tests
pytest -m "not slow"

# Run pyscf end-to-end tests (requires pyscf)
pip install pyscf
pytest -m slow tests/integration/test_full_pipeline_pyscf.py
```

The test suite is organised as:

```
tests/
├── conftest.py                             shared fixtures
├── unit/                                   one file per module
│   ├── test_datapoint.py
│   ├── test_molecule.py
│   ├── test_molecule_list.py
│   ├── test_comparison.py
│   ├── test_statistics.py
│   ├── test_benchmark_parser.py
│   ├── test_external_parser.py
│   ├── test_input_constructor.py
│   ├── test_assignment.py
│   ├── test_export.py
│   ├── test_formatting.py
│   ├── test_tree.py
│   ├── test_functions.py
│   ├── test_json_encoder.py
│   └── test_bash_wrapper.py
└── integration/
    ├── test_benchmark_to_comparison.py
    ├── test_compare_evaluate.py
    ├── test_input_generation_pipeline.py
    └── test_full_pipeline_pyscf.py         (marked @pytest.mark.slow)
```

CI runs two jobs: fast unit + integration tests on every push to `dev`, and
a separate pyscf end-to-end job. Both jobs use Python 3.12.
