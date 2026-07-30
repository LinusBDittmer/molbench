from . import logger as log


def substitute_template(template: str, subvals: dict) -> tuple[str, ...]:
    """
    Substitute placeholders "[[placeholder]]" in the template with the provided
    value in the subvals dictionary {placeholder: value}.
    Furthermore, it is possible to substitute an element of a value list/tuple
    into the template using the syntax "[[placeholder->index]]", where
    "index" refers to the position of the element in the list/tuple.

    For placeholder keys that end with the suffix "_list" and have a value
    of type list/tuple, multiple templates (one for each entry in the
    list/tuple) will be generated. For instance a subvals dictionary
    {charge_list: [0, 1], multiplicity_list: [1, 3], conv_tol: 1e-9}
    will generate two templates with the following values
    - charge = 0; multiplicity = 1, conv_tol = 1e-9
    - charge = 1; multiplicity = 3, conv_tol = 1e-9
    Note that the list "_list" suffix is removed before substitution and
    that an error is thrown if the length of charge_list and
    multiplicity_list are not equal.
    """
    if not any(key.endswith("_list") and isinstance(val, (list, tuple))
               for key, val in subvals.items()):
        return (_substitute_single_template(template, subvals),)
    # split subvals in values to expand and common values
    # and remove the _list suffix
    common = []
    to_expand = []
    for key, val in subvals.items():
        if key.endswith("_list") and isinstance(val, (list, tuple)):
            to_expand.append((key[:-5], val))
        else:
            common.append((key, val))
    # ensure that all expansion lists are of the same length
    n_variants = len(to_expand[0][1])
    if not all(len(val) == n_variants for _, val in to_expand):
        log.critical("List subvals to expand into multiple templates have to "
                     f"be all of the same length. Got\n{to_expand}\n from\n"
                     f"{subvals}", "Functions: Substitute Template")
    if n_variants == 0:
        log.warning("Received (an) empty list(s) to expand for placeholders "
                    f"{[key for key, _ in to_expand]} - no template variants "
                    f"will be generated from\n{subvals}",
                    "Functions: Substitute Template")
    # build subvals dicts for all variants
    variants = [dict(common) for _ in range(n_variants)]
    for key, val_list in to_expand:
        for var, val in zip(variants, val_list):
            if key in var:
                log.error(f"The key {key} generated from {key}_list already "
                          f"exists in the variant\n{var}\nOverwriting existing"
                          " value", "Functions: Substitute Template",
                          "KeyError")
            var[key] = val
    return tuple(
        _substitute_single_template(template, var) for var in variants
    )


def _substitute_single_template(template: str, subvals: dict) -> str:
    """
    Subsitutes placeholders "[[placeholder]]" in the template with the provided
    values in the subvals dictionary {"placeholder": val}.
    Furthermore, it is possible to substitute an element of a value list/tuple
    into the template using the syntax "[[placeholder->index]]", where
    "index" refers to the position of the element in the list/tuple.
    """
    while True:
        start = template.find("[[")
        if start == -1:
            break
        # search for the closing "]]" only after the opening "[[" so a stray
        # "]]" earlier in the template can't be mistaken for the closing
        # bracket of a later placeholder
        stop = template.find("]]", start)
        if stop == -1:
            log.error("Found an unclosed placeholder starting at "
                      f"'{template[start:start+30]}' in the template. "
                      "Leaving it as literal text - check for a missing "
                      "closing ']]'.", "Functions: Substitute Template")
            break
        key = template[start+2:stop]
        # check if we have [[key->number]] and get the key and number
        val_idx = None
        if "->" in key:
            key, val_idx = key.split("->")
            if not val_idx.isnumeric():
                log.critical(f"{val_idx} is not a number. Placeholders of the "
                             f"form {key}->{val_idx} need to be of the form "
                             "'key'->'number'.",
                             "Functions: Substitute Template")
            val_idx = int(val_idx)
        # get the actual value
        val = subvals.get(key, None)
        if val_idx is not None and val is not None:
            if not isinstance(val, (list, tuple)):
                log.critical("The value for a placeholder of the form "
                             "'[[placeholder->number]]' needs to be a 'list' "
                             f"or 'tuple'. Found {val} for key {key}.",
                             "Functions: Substitute Template")
            if val_idx >= len(val):
                log.critical(f"Index {val_idx} is out of range for "
                             f"placeholder '{key}' which only has "
                             f"{len(val)} element(s).",
                             "Functions: Substitute Template")
            val = val[val_idx]
        if val is None:
            log.critical(f"No value available for placeholder '{key}'. "
                         f"Available substitution keys are "
                         f"{list(subvals.keys())}. If '{key}' doesn't look "
                         "like a real placeholder name (e.g. it contains "
                         "spaces or shell syntax), the template likely "
                         "contains literal '[[' / ']]' content - such as "
                         "bash's own '[[ ... ]]' test syntax - that "
                         "collides with the placeholder delimiters and "
                         "needs to be rewritten to avoid double brackets.",
                         "Functions: Substitute Template")
        # update the template
        template = template.replace(template[start:stop+2], str(val))
    return template


def default_name_template(file_expansion_keys: tuple,
                          file_extension: str) -> str:
    # construct a name that ensures that each file has a unique name
    name_template = "[[name]]_[[method]]"
    for key in file_expansion_keys:
        if key not in ["name", "method"]:
            name_template += f"_[[{key}]]"
    return name_template + file_extension


def walk_dict_by_key(indict: dict, desired_key, prev_keys: tuple = tuple()):
    """
    Walk an arbitrarily nested dictionary looking for the desired key
    yielding the squence of keys and the corresponding value.
    """
    for key, val in indict.items():
        if key == desired_key:
            yield prev_keys + (key,), val
        elif isinstance(val, dict):
            for data in walk_dict_by_key(val, desired_key, prev_keys + (key,)):
                yield data


def walk_dict_values(indict: dict, prev_keys: tuple = tuple()):
    """
    Walk an arbitrarily nested dictionary yielding all non dictionary
    values and the corresponding sequence of keys.
    """
    for key, val in indict.items():
        if isinstance(val, dict):
            for data in walk_dict_values(val, prev_keys + (key,)):
                yield data
        else:
            yield prev_keys + (key,), val


def determine_basis_cardinality(basis: str):
    # Dunnings basis sets
    def _dunnings(bas: str):
        try:
            b: list[str] = bas.split("-")
            cstr: str = b[b.index("cc")+1]
            zetaidx: int = 2
            zetaoffset: int = 0
            if cstr.startswith("pwv") or cstr.startswith("pcv"):
                zetaidx += 1
            if "(" in cstr:
                zetaidx += 1
                zetaoffset += 1
            if cstr[zetaidx] == "d":
                return 2 + zetaoffset
            if cstr[zetaidx] == "t":
                return 3 + zetaoffset
            if cstr[zetaidx] == "q":
                return 4 + zetaoffset
            if cstr[zetaidx].isnumeric():
                return int(cstr[zetaidx]) + zetaoffset
        except IndexError:
            pass

        log.error(f"Basis set {bas} was interpreted as a Dunning's basis but "
                  "could not be identified!",
                  "Functions: Determine Basis Cardinality",
                  "Basis identification error")
        return 0

    # Karlsruhe def2
    def _karlsruhe(bas: str):
        try:
            b: list[str] = bas.split("-")
            cstr: str = b[b.index("def2")+1]
            zetaidx: int = 0
            if cstr.startswith("m"):
                zetaidx += 1
            if cstr[zetaidx] == "s":
                return 1
            if cstr[zetaidx] == "t":
                return 3
            if cstr[zetaidx] == "q":
                return 4
            if cstr[zetaidx].isnumeric():
                return int(cstr[zetaidx])
        except IndexError:
            pass

        log.error(f"Basis set {bas} was interpreted as a Karlruhe basis but "
                  "could not be identified!",
                  "Functions: Determine Basis Cardinality",
                  "Basis identification error")
        return 0

    b: str = basis.lower()
    # Dunnings
    dunning_prefs = ["cc-p", "aug-cc-p", "jun-cc-p", "jul-cc-p", "maug-cc-p"]
    if any(b.startswith(pref) for pref in dunning_prefs):
        return _dunnings(b)

    if "def2" in b:
        return _karlsruhe(b)

    log.error(f"Unknown basis format for {b}.",
              "Functions: Determine Basis Cardinality",
              "Basis identification error")
    return 0
