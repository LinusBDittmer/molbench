"""
Basis set extrapolation to the CBS limit
"""

from .molecule import MoleculeList, Molecule
from .functions import determine_basis_cardinality
from . import logger as log
from math import exp
from math import log as ln

class CBSExtrapolator:

    def __init__(self, basis_cardinality_dict: dict = None):
        self.basis_cardinality: dict = dict()
        if basis_cardinality_dict is not None:
            self.basis_cardinality.update({k.lower(): v for k, v in basis_cardinality_dict.items()})

    def cardinality(self, basis: str):
        if basis.lower() in self.basis_cardinality:
            return self.basis_cardinality[basis.lower()]
        return determine_basis_cardinality(basis)

    def _extrapolate(self, scf_energies: tuple[float], correlation_energies: tuple[float], 
                     cardinalities: tuple[int]) -> float:
        
        # SCF energy: Two-point Heller fit
        # see SI of doi.org/10.1080/1539445X.2020.1714656

        if len(scf_energies) > 3 or len(correlation_energies) > 3:
            log.warning(f"Detected {len(scf_energies)} SCF energies. "
                        + "Only the first three are taken into account.", "CBS Extrapolator")
            scf_energies = scf_energies[:3]
            cardinalities = cardinalities[:3]
            correlation_energies = correlation_energies[:3]

        sorted_cards = list(sorted(cardinalities))
        idx_map = [sorted_cards.index(i) for i in cardinalities]
        scf_energies = [scf_energies[i] for i in idx_map]
        correlation_energies = [correlation_energies[i] for i in idx_map]
        cardinalities = sorted_cards

        ratio: float = (scf_energies[2] - scf_energies[1]) / (scf_energies[1] - scf_energies[0])
        alpha: float = - ln(ratio)
        beta: float = (scf_energies[2] - scf_energies[1]) / (exp(- alpha * cardinalities[1]) * (ratio-1))
        scf_cbs: float = scf_energies[2] - beta * exp(-alpha*cardinalities[1])

        # Correlation energy: Helgaker Polynomial Fit

        n0, n1, n2 = cardinalities
        e0, e1, e2 = correlation_energies

        if (e0**2 + e1**2 + e2**2) < 10**-8:
            return scf_cbs + e2

        def _nr_step(z: float):
            efrac = (e2 - e1) / (e1 - e0) 
            val = efrac * (n0**(-z) - n1**(-z)) - (n1**(-z) - n2**(-z))
            deriv = (efrac * (-ln(n0) * n0**(-z) + ln(n1) * n1**(-z)) 
                     + ln(n1) * n1**(-z) - ln(n2) * n2**(-z))
            return - val / deriv

        # z = max(min(((e1-e2)/(e0-e1))**(1/(n0-n2)), 2), 0.5)
        z = 1.5
        max_stepsize = 10**5
        for i in range(1000):
            zz = 0.5 + i / 100
            local_stepsize = _nr_step(zz)
            if abs(local_stepsize) < max_stepsize:
                max_stepsize = abs(local_stepsize)
                z = zz

        converged = False
        print("Starting new NR iteration")
        for i in range(50):
            step = _nr_step(z)
            print(f"Z: {z}, stepsize: {step}")
            z += step
            if abs(step) < 10**-5:
                converged = True
                break
        
        if not converged:
            log.warning("Numerical solution for correlation energy extrapolation did not converge.", "CBS Extrapolator")
        
        bval = (e1 - e0) / (n1**(-z) - n0**(-z))
        ecorr_cbs = e0 - bval * n0**(-z)
        if z < 0 or abs(z) < 10**-5:
            log.warning("Correlation energies are overly linear, the latest correlation energy is taken", "CBS Extrapolator")
            ecorr_cbs = e2
        
        log.debug(f"Correlation energies: {e0} / {e1} / {e2}", "CBS Extrapolator")
        log.debug(f"Extrapolated SCF energy: {scf_cbs}", "CBS Extrapolator")
        log.debug(f"Extrapolated Correlation energy: {ecorr_cbs}", "CBS Extrapolator")

        return scf_cbs + ecorr_cbs


    def extrapolate(self, mols: MoleculeList, min_basis_set_num: int = 3, 
                    energy_keyword: str = 'energy', scf_energy_keyword: str = 'scf energy',
                    consider_relative_properties: str = False) -> MoleculeList:
        if consider_relative_properties:
            energy_keyword = (energy_keyword, "component " + energy_keyword)
            scf_energy_keyword = (scf_energy_keyword, "component " + scf_energy_keyword)
        else:
            energy_keyword = (energy_keyword,)
            scf_energy_keyword = (scf_energy_keyword,)

        # Molecules sorted by name
        name_dict: dict = dict()
        for mol in mols:
            if mol.name not in name_dict:
                name_dict[mol.name] = [mol]
            else:
                name_dict[mol.name].append(mol)

        def _get_basis_sets(mol: Molecule) -> set:
            bs: set[str] = set()
            for _, state in mol.state_data.items():
                if "basis" in state:
                    bs.add(state["basis"].lower())
            return bs

        # Sorting out all molecules with less than the required number of basis sets
        cbs_dict: dict = dict()
        for name, mollist in name_dict.items():
            if len(mollist) < min_basis_set_num:
                continue
            basis_sets: set = set()
            for mol in mollist:
                basis_sets.update(_get_basis_sets(mol))
            if len(basis_sets) >= min_basis_set_num:
                cbs_dict[name] = mollist
        
        cbs_mollist: MoleculeList = MoleculeList()

        for name, mollist in cbs_dict.items():
            data_id: str = mollist[0].data_id + "-cbs"
            system_data: dict = mollist[0].system_data
            state_data: dict = dict()
            energies: list[float] = []
            scf_energies: list[float] = []
            basis_cardinalities: list[int] = []
            energy_dict: dict = dict()
            scf_energy_dict: dict = dict()
            cardinality_dict: dict = dict()

            for mol in mollist:
                key_offset: int = len(state_data)
                keymap: dict = dict()
                # Transferring data from one Mol into the joint Mol
                for idx, key in enumerate(mol.state_data):
                    i = idx + key_offset
                    newkey = f"p{i:03d}"
                    state_data[newkey] = mol.state_data[key]
                    keymap[key] = newkey
                for _, state in mol.state_data.items():
                    if "stochiometry" in state:
                        state["stochiometry"] = {keymap.get(k, k): v for k, v in state["stochiometry"].items()}
                        continue
                    if state["type"] in energy_keyword and "value" in state:
                        if "component index" in state:
                            if state["component index"] not in energy_dict:
                                energy_dict[state["component index"]] = [state["value"]]
                            else:
                                energy_dict[state["component index"]].append(state["value"])
                            local_c = self.cardinality(state["basis"])
                            if state["component index"] not in cardinality_dict:
                                cardinality_dict[state["component index"]] = [local_c]
                            else:
                                cardinality_dict[state["component index"]].append(local_c)
                        else:
                            energies.append(state["value"])
                            basis_cardinalities.append(self.cardinality(state["basis"]))
                        
                    elif state["type"] in scf_energy_keyword and "value" in state:
                        if "component index" in state:
                            if state["component index"] not in scf_energy_dict:
                                scf_energy_dict[state["component index"]] = [state["value"]]
                            else:
                                scf_energy_dict[state["component index"]].append(state["value"])
                        else:
                            scf_energies.append(state["value"])

            # Split behaviour for relative and absolute energies.
            # Absolute energies

            if len(energy_dict) == 0:
                if len(energies) != len(scf_energies):
                    log.critical("Mismatched reading of energies and SCF energies: "
                                + f"Found {len(energies)} energies and "
                                + f"{len(scf_energies)}", "CBS Extrapolator")

                correlation_energies: list[float] = [e - e_scf for e, e_scf in zip(energies, scf_energies)]
                cbs_energy = self._extrapolate(scf_energies, correlation_energies, basis_cardinalities)
                cbs_key = list(state_data.keys())[0]
                state_data[cbs_key+"_cbs"] = {"type": energy_keyword[0], 
                                            "method": state_data[cbs_key]["method"], 
                                            "basis": "CBS", 
                                            "value": cbs_energy}

                cbs_mollist.append(Molecule(name, data_id, system_data, state_data))
            else:
                for component_index in energy_dict:
                    energies: list[float] = energy_dict[component_index]
                    scf_energies: list[float] = scf_energy_dict[component_index]
                    basis_cardinalities: list[float] = cardinality_dict[component_index]
                    correlation_energies: list[float] = [e - e_scf for e, e_scf in zip(energies, scf_energies)]
                    cbs_energy = self._extrapolate(scf_energies, correlation_energies, basis_cardinalities)
                    etype = energy_keyword[0]
                    method = ""
                    for state, statedata in mol.state_data.items():
                        if "component index" in statedata:
                            if statedata["component index"] == component_index:
                                etype = statedata.get("type", energy_keyword[0])
                                method = statedata.get("method", "")
                                break
                    cbs_index = len(state_data)
                    cbs_key = f"p{cbs_index:03d}"
                    etype = etype.replace("scf", "")
                    etype = " ".join(etype.split())
                    state_data[cbs_key] = {"type": etype,
                                            "method": method,
                                            "basis": "CBS",
                                            "value": cbs_energy,
                                            "component index": component_index}
                # Correcting relative entries
                # We assume duplication and only take the first one
                tmp_state_data = dict()
                for stateidx, statename in enumerate(state_data):
                    statedict = state_data[statename]
                    if not "stochiometry" in statedict:
                        continue
                    stochiometry = statedict["stochiometry"]
                    # keymap = {k: str(int(k[1:])).format("03d") for k in stochiometry}
                    component_indices = {k: state_data[k]["component index"] for k in stochiometry}
                    keymap2 = dict()
                    etype = ""
                    method = ""

                    for finalstate in state_data:
                        if "basis" not in state_data[finalstate]:
                            continue
                        if state_data[finalstate]["basis"].lower() != "cbs":
                            continue
                        if "component index" not in state_data[finalstate]:
                            continue
                        keymap2[state_data[finalstate]["component index"]] = finalstate
                        etype = state_data[finalstate].get("type", "").replace("component", "").strip()
                        method = state_data[finalstate].get("method", "")
                
                    new_stochiometry = {keymap2[component_indices[k]]: v for k, v in stochiometry.items()}
                    cbs_index = len(state_data)
                    cbs_key = f"p{cbs_index:03d}"

                    tmp_state_data[cbs_key] = {"type": etype, "method": method, "stochiometry": new_stochiometry}
                    if stateidx >= len(mollist[0].state_data):
                        break
                state_data.update(tmp_state_data)

                cbs_mollist.append(Molecule(name, data_id, system_data, state_data))

        return cbs_mollist

        

        
