import ast
import json
import os


FOLDER = "./solutions/"

def save_solution(data: dict, file_name: str) -> None:
    """Saves the dict as a json file

    Args:
        data (dict):
    """
    filepath = os.path.join(FOLDER, file_name + ".json")
    
    with open(filepath, 'w') as outfile:
        json.dump(data, outfile, sort_keys=True, indent=4)
        print(f"\n{file_name} saved in {filepath}")

def load_solution(file_name: str) -> dict:
    
    data = dict()
    with open(os.path.join(FOLDER, file_name + ".json")) as f:
        data = json.load(f)
    
    # Remap variables keys to tuple
    converted_data = {key: __convert_to_tuple__(value) for key, value in data['variables'].items()}
    data["variables"] = converted_data
    
    return data


# Function to convert string representation of tuple to actual tuple
def __convert_to_tuple__(dictionary):
    return {ast.literal_eval(key): value for key, value in dictionary.items()}    



if __name__ == "__main__":
    for k,v in load_solution("usaphmp_2023_12_26_18_35_10.json").items():
        print(f"{k}: {v}")
