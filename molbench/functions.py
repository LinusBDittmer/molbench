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
                      "substitute_template", KeyError)
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
