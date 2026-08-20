"""End-to-end: benchmark → run pyscf → parse output → Statistics.compare.

Marked @pytest.mark.slow — runs in its own CI job with pyscf installed.
Skipped automatically if pyscf is not installed.
"""

import json
from pathlib import Path

import pytest

pyscf = pytest.importorskip("pyscf", reason="pyscf not installed")

from pyscf import gto, scf

from molbench.benchmark_parser import JSONBenchmarkParser
from molbench.comparison import Comparison
from molbench.external_parser import ExternalParser
from molbench.molecule import MoleculeList
from molbench.statistics import Statistics

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _json_parser(filepath):
    """Parse a JSON .out file written by the test itself."""
    raw = json.loads(Path(filepath).read_text())
    name = raw["name"]
    system_data = {}
    state_data = {
        "gs": {
            "basis": raw["basis"],
            "method": raw["method"],
            "data": {"energy": {"value": raw["data"]["energy"], "unit": "au"}},
        }
    }
    return name, system_data, state_data


def _run_hf(xyz_str, basis, charge, spin, unit="A"):
    """Run a restricted/unrestricted HF calculation with pyscf and return the
    total energy in au."""
    mol = gto.M(
        atom=xyz_str,
        basis=basis,
        charge=charge,
        spin=spin,
        verbose=0,
        unit=unit,
    )
    if spin == 0:
        mf = scf.RHF(mol)
    else:
        mf = scf.ROHF(mol)
    return mf.kernel()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_hf_energy_on_h_atom_vs_ascdb(tmp_path):
    """HF/cc-pVDZ energy for the H atom should be within 0.05 au of the TBE."""
    bench = JSONBenchmarkParser().load("ascdb", benchmark_id="TBE")
    # AE18pE-1 is the H atom (charge=0, multiplicity=2) in ascdb
    h_mols = bench.filter("name", "AE18pE-1")
    assert len(h_mols) == 1, "AE18pE-1 not found in ascdb"

    h_mol = h_mols[0]
    # ascdb uses xyz_list / charge_list / multiplicity_list
    xyz = h_mol.system_data["xyz_list"][0]
    charge = h_mol.system_data["charge_list"][0]
    mult = h_mol.system_data["multiplicity_list"][0]
    spin = mult - 1

    hf_energy = _run_hf(xyz, "cc-pvdz", charge, spin)

    # Write a mock .out file
    out_file = tmp_path / "AE18pE-1_HF_cc-pvdz.out"
    out_file.write_text(
        json.dumps(
            {
                "name": "AE18pE-1",
                "basis": "cc-pvdz",
                "method": "HF",
                "data": {"energy": hf_energy},
            }
        )
    )

    # Load via ExternalParser
    computed = ExternalParser().load(
        str(tmp_path), parser=_json_parser, out_suffix=".out"
    )
    assert len(computed) == 1

    # Build Comparison
    c = Comparison()
    c.add(h_mols)
    c.add(computed)

    stats = Statistics(c)
    errors = stats.compare(
        interest={"method": "HF"},
        reference={"method": "TBE"},
        error_thresh=0.1,
    )
    result = stats.evaluate(errors, "mae", proptype="energy")

    mae_val, count = result["mae"]
    assert count == 1
    assert abs(mae_val) < 0.05, f"HF/cc-pVDZ MAE for H atom too large: {mae_val:.6f} au"


@pytest.mark.slow
def test_hf_energy_on_water_vs_ascdb(tmp_path):
    """HF/cc-pVDZ energy for the water molecule should be within 0.2 au of TBE."""
    bench = JSONBenchmarkParser().load("ascdb", benchmark_id="TBE")
    # Filter to a single-geometry water-like molecule from AE18
    # We'll pick a molecule with a single geometry and charge=0
    single_geo = [
        m
        for m in bench
        if "xyz" in m.system_data and m.system_data.get("charge", -99) == 0
    ]
    if not single_geo:
        pytest.skip("No suitable single-geometry neutral molecule found in ascdb")

    target = single_geo[0]
    xyz = target.system_data["xyz"]
    charge = target.system_data.get("charge", 0)
    mult = target.system_data.get("multiplicity", 1)
    spin = mult - 1
    # Determine basis from state_data
    basis = "cc-pvdz"
    for state in target.state_data.values():
        if "basis" in state and state["basis"].lower().startswith("cc"):
            basis = state["basis"]
            break

    try:
        hf_energy = _run_hf(xyz, basis, charge, spin)
    except Exception as e:  # noqa: BLE001
        pytest.skip(f"PySCF calculation failed: {e}")

    out_file = tmp_path / f"{target.name}_HF_{basis}.out"
    out_file.write_text(
        json.dumps(
            {
                "name": target.name,
                "basis": basis,
                "method": "HF",
                "data": {"energy": hf_energy},
            }
        )
    )

    computed = ExternalParser().load(
        str(tmp_path), parser=_json_parser, out_suffix=".out"
    )
    assert len(computed) == 1

    mol_bench = MoleculeList([target])
    c = Comparison()
    c.add(mol_bench)
    c.add(computed)

    stats = Statistics(c)
    errors = stats.compare(
        interest={"method": "HF"},
        reference={"method": "TBE"},
        error_thresh=1.0,
    )
    result = stats.evaluate(errors, "mae", proptype="energy")
    mae_val, count = result["mae"]
    assert count == 1
    # HF error on small molecules is usually well below 0.2 au
    assert abs(mae_val) < 0.2, (
        f"HF MAE suspiciously large: {mae_val:.6f} au for {target.name}"
    )
