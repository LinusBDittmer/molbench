from . import logger as log


def new_assignment_file(state_ids: list, comment_token: str = "#",
                        id_separator: str = "->",
                        null_token: str = "null") -> str:
    """
    Creates the content of an assignment file.

    Parameters
    ----------
    state_ids : list
        The state ids of the reference state to write in the new assignment
        file.
    comment_token : str, optional
        Token to mark comments (default: '#').
    id_separator : str, optional
        The separator of reference and external id (default: '->').
    null_token : str, optional
        Identifier for empty, not assigned states.
    """
    content = f"{comment_token} ref_state_id {id_separator} external_id"
    for id in state_ids:
        content += f"\n{id} {id_separator} {null_token}"
    return content


def parse_assignment_file(assignmentfile: str,
                          comment_token: str = "#",
                          id_separator: str = "->",
                          null_token: str = "null",
                          import_external: callable = None,
                          import_ref: callable = None) -> dict:
    """
    Parses an assignment file.

    Parameters
    ----------
    assignmentfile : str
        The file to parse.
    comment_token : str, optional
        Token to mark comments in the assignment file (default: '#').
    id_separator : str, optional
        The token separating reference and external id (default: '->').
    null_token : str, optional
        Identifier for empty, not assigned states (default: 'null').
    import_external: callable, optional
        Imports the external id before inserting the id in the dictionary.
    import_ref : callable, optional
        Imports the reference id before inserting the id in the dictionary.

    Returns
    -------
    dict
        Dictionary of the form {external_id: ref_id}.
    """
    with open(assignmentfile, "r") as f:
        content = f.read()

    state_assignments = {}
    for line in content.splitlines():
        line = line.strip()
        # - remove comments
        if line.startswith(comment_token):
            continue
        assignment = line.split(comment_token, maxsplit=1)[0]
        # - read both state id's and import them if necessary
        assignment = assignment.split(id_separator)
        if len(assignment) != 2:
            log.critical(f"Invalid use of state id separator {id_separator} in"
                         f" line {line} of assignment file {assignmentfile}.",
                         "parse_assignment_file", "Assignment")

        ref = assignment[0].strip()
        external = assignment[1].strip()
        if ref == null_token or external == null_token:  # skip not assigned
            log.warning(f"Unassigned state in file {assignmentfile}.",
                        "Assigment")
            continue

        if import_external is not None:
            external = import_external(external)
        if import_ref is not None:
            ref = import_ref(ref)

        if external in state_assignments:
            log.warning(f"The external state {external} in assignment file "
                        f"{assignmentfile} is assigned twice. Overwriting the "
                        "first assignment.", "Assigment")
        state_assignments[external] = ref
    return state_assignments
