import json
import os
from time import strftime

FOLDER = "./solutions/"

def save_solution(data: dict):
    timestamp = strftime("%Y_%m_%d_%H_%M_%S")
    file_name = data['model'] + "_" + timestamp + ".json"
    filepath = os.path.join(FOLDER, file_name) 
    with open(filepath, 'w') as outfile:
        json.dump(data, outfile, sort_keys=False, indent=4)
        print(f"Solution of model {data['model']} saved in {filepath}")
   


# if __name__ == "__main__":
#     save_solution({"model": "m1"})
