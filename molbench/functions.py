from . import logger as log


def substitute_template(template: str, subvals: dict) -> tuple[str]:
    if "xyz_list" not in subvals:
        return (_substitute_single(template, subvals),)
    num_subs: int = len(subvals["xyz_list"])
    # We remove the "_list" from each key and call _substitute_single
    substitutions: list = []
    for index in range(num_subs):
        subval_dict: dict = dict()
        for key, val in subvals.items():
            newkey: str = key.removesuffix("_list")
            subval_dict[newkey] = val[index] if "_list" in key else val
        substitutions.append(_substitute_single(template, subval_dict))
    return tuple(substitutions)

def _substitute_single(template: str, subvals: dict) -> str:
    while True:
        start = template.find("[[")
        stop = template.find("]]")
        if start == -1 or stop == -1:
            break
        key = template[start+2:stop]
        key_subindex = None
        if "->" in key:
            key, key_subindex = key.split("->")
            key_subindex = int(key_subindex)
        val = subvals.get(key, None)
        if key_subindex is not None:
            val = val[key_subindex]
        if val is None:
            log.error(f"No value for required parameter {key} "
                      f"available. Available are {subvals}.",
                      "Functions: Substitute Template", "KeyError")
        template = template.replace(template[start:stop+2], str(val))
    return template


def walk_dict_by_key(indict: dict, desired_key, prev_keys: tuple = tuple()):
    """Walk an arbitrarily nested dictionary looking for the desired key
       yielding the squence of keys and the corresponding value."""
    for key, val in indict.items():
        if key == desired_key:
            yield prev_keys + (key,), val
        elif isinstance(val, dict):
            for data in walk_dict_by_key(val, desired_key, prev_keys + (key,)):
                yield data


def walk_dict_values(indict: dict, prev_keys: tuple = tuple()):
    """Walk an arbitrarily nested dictionary yielding all non dictionary
       values and the corresponding sequence of keys."""
    for key, val in indict.items():
        if isinstance(val, dict):
            for data in walk_dict_values(val, prev_keys + (key,)):
                yield data
        else:
            yield prev_keys + (key,), val

def determine_basis_cardinality(basis: str):
    # Dunnings basis sets
    def _dunnings(bas: str):
        b: list[str] = bas.split("-")
        cstr: str = b[b.index("cc")+1]
        zetaidx: int = 2
        zetaoffset: int = 0
        if cstr.startswith("pwv") or cstr.startswith("pcv"):
            zetaidx += 1
        if "(" in cstr:
            zetaidx += 1
            zetaoffset += 1
        if cstr[zetaidx] == "d": return 2 + zetaoffset
        if cstr[zetaidx] == "t": return 3 + zetaoffset
        if cstr[zetaidx] == "q": return 4 + zetaoffset
        if cstr[zetaidx].isnumeric(): return int(cstr[zetaidx]) + zetaoffset
        
        log.error(f"Basis set {bas} was interpreted as a Dunning's basis but could not be identified!", 
                  "Functions: Determine Basis Cardinality", "Basis identification error")
        return 0

    # Karlsruhe def2
    def _karlsruhe(bas: str):
        b: list[str] = bas.split("-")
        cstr: str = b[b.index("def2")+1]
        zetaidx: int = 0
        if cstr.startswith("m"):
            zetaidx += 1
        if cstr[zetaidx] == "s": return 1
        if cstr[zetaidx] == "t": return 3
        if cstr[zetaidx] == "q": return 4
        if cstr[zetaidx].isnumeric(): return int(cstr[zetaidx])

        log.error(f"Basis set {bas} was interpreted as a Karlruhe basis but could not be identified!", 
                  "Functions: Determine Basis Cardinality", "Basis identification error")
        return 0

    b: str = basis.lower()
    # Dunnings
    for identifier in ["cc-p", "aug-cc-p", "jun-cc-p", "jul-cc-p", "maug-cc-p"]:
        if b.startswith(identifier):
            return _dunnings(b)
    
    if "def2" in b:
        return _karlsruhe(b)
    
    log.error(f"Unknown basis format for {b}.", "Functions: Determine Basis Cardinality", 
              "Basis identification error")

    return 0

    