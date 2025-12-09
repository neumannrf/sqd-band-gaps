import numpy as np
import os
import pickle
import yaml
import json
import re 

def load_yaml(folder_path,file_name=None):
    if file_name is None:
        with open((folder_path), "r") as f:
            return yaml.safe_load(f)

    else:
        with open(os.path.join(folder_path,file_name), "r") as f:
            return yaml.safe_load(f)
        

def save_pickle(obj,folder_path,file_name=None):
    if file_name is None:
        pickle.dump(obj,open(folder_path,"wb"))
    else:
        pickle.dump(obj,open(os.path.join(folder_path,file_name),"wb"))

def load_pickle(folder_path,file_name=None):
    if file_name is None:
        return pickle.load(open(folder_path,"rb"))
    else:
        return pickle.load(open(os.path.join(folder_path,file_name),"rb"))
    
def save_json(obj,folder_path,file_name=None):
    if file_name is None:
        json.dump(obj,open(folder_path,"w"))
    else:
        json.dump(obj,open(os.path.join(folder_path,file_name),"w"))

def load_json(folder_path,file_name=None):
    if file_name is None:
        return json.load(open(folder_path,"r"))
    else:
        return json.load(open(os.path.join(folder_path,file_name),"r"))

def make_pyscf_slater_det(nelec, num_orbitals):
    a, b = nelec
    n = a + b
    result = [0] * num_orbitals
    i = 0
    while n > 0 and i < num_orbitals:
        add = min(2, n)
        result[i] = add
        n -= add
        i += 1
    return np.array(result)

def threshold_filter(arr, threshold=2.220446049250313e-16):

    if threshold is not None:
    # Create a copy to avoid modifying the original array
        filtered_arr = np.copy(arr)
        
        # Set all elements below the threshold to zero
        filtered_arr[np.abs(filtered_arr) < threshold] = 0

        return filtered_arr
    
    else:
        return arr
    
def filter_hamming_weights(orbitals_dict, reference_down, reference_up):
    filtered_dict = {}

    for binary_string, value in orbitals_dict.items():
        # Split the string into up and down electron halves
        half_length = len(binary_string) // 2
        up_electrons = binary_string[:half_length]
        down_electrons = binary_string[half_length:]

        # Calculate Hamming weights
        up_weight = up_electrons.count('1')
        down_weight = down_electrons.count('1')

        # Check if current string matches the defined reference Hamming weights
        if (up_weight, down_weight) == (reference_up, reference_down):
            filtered_dict[binary_string] = value

    return filtered_dict

def find_max_index(folder_path):
    pattern = re.compile(r"(\d+)")
    max_i = -1

    for filename in os.listdir(folder_path):
        match = pattern.search(filename)  # search instead of match to find number anywhere
        if match:
            i = int(match.group(1))
            if i > max_i:
                max_i = i

    print("Max index: ",max_i)
    return max_i

def load_energy_variance(sqd_results_folder,hci=False):
    # Loop over files in folder
    point_list = []
    for filename in os.listdir(sqd_results_folder):
        if filename.endswith('.json'):
            # match = re.search(r'(\d+)', filename)
            # if match:
            #     index = int(match.group(1))
            #     file_path = os.path.join(folder_path, filename)
            file_path = os.path.join(sqd_results_folder, filename)
            data_dict = load_json(file_path)
            if hci:
                point_list.append([data_dict["variance"],data_dict["hci_energy"],data_dict["subspace_dimension"]])
            else:    
                point_list.append([data_dict["variance"],data_dict["sqd_energy"],data_dict["pyscf_subspace_dimension"]])
        
    return np.array(point_list)

def extract_subspace_energy_pairs(folder_path,key="pyscf_subspace_dimension"):
    result = []
    for filename in os.listdir(folder_path):
        if filename.startswith("results_dict_") and filename.endswith(".json"):
            file_path = os.path.join(folder_path, filename)
            try:
                with open(file_path, 'r') as f:
                    data = json.load(f)
                    subspace_dim = data.get(key)
                    print(subspace_dim)
                    if key=="pyscf_subspace_dimension":
                        sqd_energy = data.get("sqd_energy")
                    else:
                        sqd_energy = data.get("hci_energy")
                    variance = np.abs(data.get("variance"))
                    print("Variance",variance)
                    if subspace_dim is not None and sqd_energy is not None:
                        result.append([subspace_dim, sqd_energy,variance])
            except (json.JSONDecodeError, IOError) as e:
                print(f"Skipping {filename}: {e}")
    # Sort by pyscf_subspace_dimension (the first element of each pair)
    result.sort(key=lambda x: x[0])
    return np.array(result)
