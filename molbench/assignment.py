def new_assignment_file(state_ids: list) -> str:
    content = "# ref_state_id -> external_id"
    for id in state_ids:
        content += f"\n{id} -> null"
    return content
