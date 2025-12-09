import pickle
import os
from qiskit_ibm_runtime import QiskitRuntimeService
from collections import defaultdict

from sqd_lattice.util import load_pickle,load_yaml,load_json,filter_hamming_weights

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

def common_keys(dict_a, dict_b):
    """
    Returns the list of keys from dictionary A that are also present in dictionary B.

    Parameters:
    dict_a (dict): The first dictionary.
    dict_b (dict): The second dictionary.

    Returns:
    list: A list of keys from A that are also in B.
    """
    return [key for key in dict_a if key in dict_b]

material = 'zirconia'
ne = "23"

folder_path = f"runs/{material}/{ne}e"

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
results_pre_process = load_pickle(folder_path,"results_pre_process.pkl")
num_orbitals = data_dict['num_orbitals']
ne_a,ne_b = results_pre_process['ne_ab']

## LOADING COUNTS FROM PLATAFORM
# service = QiskitRuntimeService(name='alan-ad')
# job = service.job('czr2tb9qnmvg008tvb5g')
# job_result = job.result()
# merged_dict = defaultdict(int)
# # List of dictionaries to merge
# dicts = [job_result[i].data.c.get_counts() for i in range(len(job_result))]
# # Iterate over each dictionary and sum values for common keys
# for d in dicts:
#     for key, value in d.items():
#         merged_dict[key] += value
# # Convert back to a regular dictionary if needed
# counts = dict(merged_dict)

# LOAD COUNTS FROM FILE
counts = load_json(os.path.join(folder_path,'counts','counts_hardware_torino.json'))

print(f'Number of unique counts: {len(counts)}')

f_counts = filter_hamming_weights(counts,ne_a,ne_b)

print(f"Filtered hardware SQD lenght: {len(f_counts)}")

print(f'Signal to noise ratio: {len(f_counts)/len(counts)}')

from qiskit_addon_sqd.counts import generate_counts_uniform

rand_seed = 42
counts_w = generate_counts_uniform(len(counts), num_orbitals * 2, rand_seed=rand_seed)
filtered_counts_w = filter_hamming_weights(counts_w,ne_a, ne_b)
print(f"White noise {len(filtered_counts_w)/len(counts_w)}")

# SAVE COUNTS TO FILE
# pickle.dump(counts,open(os.path.join(self.data_dir_path,f'counts_hardware_{sum(self.nelec)}e.pkl'),'wb'))

# LOAD NOISELESS COUNTS
# counts_noiseless = load_json(os.path.join(folder_path,'counts','counts_ffsim.json'))

# result = common_keys(counts, counts_noiseless)
# print(f"Noiseless SQD lenght: {len(counts_noiseless)}")
# print(f"Intersection: {len(result)}")