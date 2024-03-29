import ast
import json
import os


FOLDER_SOLUTIONS = "src/solutions/"

def save_dictionary(data: dict, file_name: str) -> None:
    """Saves the dict as a json file

    Args:
        data (dict):
    """
    filepath = os.path.join(FOLDER_SOLUTIONS, file_name + ".json")
    
    with open(filepath, 'w') as outfile:
        json.dump(data, outfile, sort_keys=True, indent=4)
        print(f"\n{file_name} saved in {filepath}")

def load_solution(file_name: str) -> dict:
    
    data = dict()
    with open(os.path.join(FOLDER_SOLUTIONS, file_name + ".json")) as f:
        data = json.load(f)
    
    # Remap variables keys to tuple
    converted_data = {key: _convert_to_tuple(value) for key, value in data['variables'].items()}
    data["variables"] = converted_data
    
    return data


# Function to convert string representation of tuple to actual tuple
def _convert_to_tuple(dictionary):
    return {ast.literal_eval(key): value for key, value in dictionary.items()}


if __name__ == "__main__":
    for k,v in load_solution("usaphmp_2023_12_26_18_35_10.json").items():
        print(f"{k}: {v}")
