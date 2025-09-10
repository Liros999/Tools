import json

def parse_groupings(file_path):
    """
    Parses the Subtiwiki groupings file and returns a hierarchical dictionary.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    groupings = {}
    current_path = []

    for line in lines:
        line = line.strip()
        if not line:
            continue

        parts = line.split(' ', 1)
        if len(parts) != 2:
            continue
        dot_id = parts[0].strip('.')
        level = len(dot_id.split('.'))
        name = parts[1].strip()

        while len(current_path) >= level:
            current_path.pop()

        current_dict = groupings
        for key in current_path:
            current_dict = current_dict[key]['children']

        # Preserve SubtiWiki dot-notation for reliable category mapping
        node = {'dot': dot_id, 'children': {}}
        current_dict[name] = node
        current_path.append(name)

    return groupings

def save_groupings_as_json(groupings, json_path):
    """
    Saves the groupings dictionary as a JSON file.
    """
    with open(json_path, 'w') as f:
        json.dump(groupings, f, indent=4)

if __name__ == '__main__':
    groupings = parse_groupings('Subtiwiki_Grouping')
    save_groupings_as_json(groupings, 'static/subtiwiki_groupings.json')
    print("Successfully parsed and saved the Subtiwiki groupings to static/subtiwiki_groupings.json")