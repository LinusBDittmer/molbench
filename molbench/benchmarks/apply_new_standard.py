import json

benchmark_old = "doubles_jacquemin19.json"
benchmark_new = "doubles_jacquemin19_new.json"

# benchmark_old = "jacquemin18.json"
# benchmark_new = "jacquemin18_new.json"

# benchmark_old = "jacquemin21.json"
# benchmark_new = "jacquemin21_new.json"

# benchmark_old = "ascdb.json"
# benchmark_new = "ascdb_new.json"

combine = True

default_units = {
    "excitation energy": "eV",
    "energy": "au",
    "oscillator strength": "au"
}

with open(benchmark_old, "r") as f:
    old = json.load(f)

new_data = dict()
for molecule, moldata in old.items():
    new_data[molecule] = dict(moldata)
    new_data[molecule]["properties"] = dict()

    new_props = dict()
    for prop, propdata in moldata["properties"].items():
        ind_prop = {"data": dict()}
        keytype = None
        keyval = None
        keyunit = None
        for key, val in propdata.items():
            if key in ("method", "basis"):
                ind_prop[key] = val
            elif key == "state_id":
                ind_prop["transition_id"] = "s0->" + val
            elif key == "type":
                keytype = val
            elif key == "value":
                keyval = val
            elif key == "unit":
                keyunit = val
            else:
                ind_prop[key] = val

        if keyunit is None:
            # Defaults
            keyunit = default_units[keytype]

        keytype = keytype.replace(" ", "_")

        ind_prop["data"][keytype] = {"value": keyval, "unit": keyunit}
        new_props[prop] = ind_prop

    if combine:
        combined_props = dict()
        for prop, propdata in new_props.items():
            added = False
            if len(combined_props.keys()) == 0:
                combined_props[prop] = propdata
            else:
                for cpr, cprd in combined_props.items():
                    if cprd["method"] == propdata["method"] and cprd["basis"] == propdata["basis"]:
                        cprd["data"].update(propdata["data"])
                        added = True
                        break
                if not added:
                    combined_props[prop] = propdata
                    
    else:
        combined_props = new_props

    new_data[molecule]["properties"] = combined_props

with open(benchmark_new, "w") as f:
    json.dump(new_data, f, sort_keys=True, indent=4, ensure_ascii=True)

