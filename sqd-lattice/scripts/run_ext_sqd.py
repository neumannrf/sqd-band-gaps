import os 
import argparse
import numpy as np
import json
import time
import h5py
import concurrent.futures
import ffsim
import pyscf
from pyscf import ao2mo
#careful with the next few imports. 
# We are exposing methods that are not supposed to be visible by the user
from pyscf.fci.selected_ci import _all_linkstr_index, _as_SCIvector
from pyscf.fci import direct_spin1
from pyscf.fci.selected_ci import enlarge_space

from sqd_lattice.util import load_yaml,load_pickle
from sqd_lattice.sqd import load_completed_batches,load_iteration_from_checkpoint,save_batch_checkpoint,filter_addresses,load_final_results,make_hf_bitstring
from sqd_lattice.ext_sqd import filter_configurations,build_bitstring_matrix_and_probs
from sqd_lattice.ext_sqd import create_all_double_transitions,create_transitions_doubles_and_singles,create_transitions_doubles_across_spin,create_transitions_doubles_within_same_spin,create_transitions_single

try:
    from sqd_lattice.sqd_dice import process_batch_dice as process_batch
    print("DICE solver available!")
    dice = True
except ImportError:
    print("DICE solver not installed, running with PYSCF solver")
    from sqd_lattice.sqd import process_batch
    dice = False

from qiskit_addon_sqd.counts import counts_to_arrays
from qiskit_addon_sqd.configuration_recovery import recover_configurations
from qiskit_addon_sqd.subsampling import postselect_and_subsample
from qiskit_addon_sqd.fermion import bitstring_matrix_to_ci_strs,enlarge_batch_from_transitions

# Loading material name
parser = argparse.ArgumentParser()
parser.add_argument('--material')
parser.add_argument('--ne',help="Total number of electrons to index file")
parser.add_argument('--samples_per_batch',help="samples per batch for SQD")
args = parser.parse_args()

material = args.material
ne = args.ne
samples_per_batch = int(args.samples_per_batch)

folder_path = f"runs/{material}/{ne}e"

print(f"Running Ext-SQD for {material} {ne} electrons and {samples_per_batch} samples per batch",flush=True)

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
results_pre_process = load_pickle(folder_path,"results_pre_process.pkl")
sqd_config = load_yaml(folder_path,"sqd_config.yaml")
ext_sqd_config = load_yaml(folder_path,"sqd_config.yaml")["ext_sqd"]

print("\n",ext_sqd_config,"\n",flush=True)

# Degfining variables from data dicts
if base_config["basis"] != "atomic":
    hcore = results_pre_process["hopping_matrix_rot"] 
    eri = results_pre_process["two_body_tensor_rot"]
else: 
    hcore = data_dict["hopping_matrix"]
    eri = results_pre_process["two_body_tensor"]
num_orbitals = data_dict["num_orbitals"]
num_elec_ab = results_pre_process["ne_ab"]
num_elec_a = num_elec_ab[0]
num_elec_b = num_elec_ab[1]
open_shell = False if num_elec_a==num_elec_b else True
nuclear_repulsion_energy = 0.
iterations = ext_sqd_config["recovery_iterations"]
n_batches = ext_sqd_config["n_batches"]
num_cpus = ext_sqd_config["num_cpus"]
checkpoint_interval = n_batches
max_davidson_cycles = sqd_config["max_davidson_cycles"]
sampling = ext_sqd_config["sampling"]
threshold = ext_sqd_config["amplitude_threshold"]

# Defining results and checkpoints files
results_folder = os.path.join(folder_path,f"results_ext_sqd_{sqd_config['sampling']}")
results_full_folder = os.path.join(results_folder,"results_full")
os.makedirs(results_full_folder,exist_ok=True)
checkpoint_file = os.path.join(results_full_folder,f'results_ext_{samples_per_batch}_checkpoint.h5')
results_file = os.path.join(results_full_folder,f'results_ext_{samples_per_batch}.h5')
old_results_file = os.path.join(folder_path,
                                f"results_sqd_{sqd_config['sampling']}",
                                "results_full",
                                f"results_{ext_sqd_config['samples_per_batch_old']}.h5")


#####################################################################################################

# Load previous SQD data
e_hist, occupancy_hist, final_data = load_final_results(old_results_file)

# Initialize a dictionary to hold merged configurations
merged_configs = {}

# Iterate over all batches
for batch in sorted(final_data.keys()):
    ci_coeffs = final_data[batch]['ci_coeffs']
    addresses_alpha = final_data[batch]['addresses_alpha']
    addresses_beta = final_data[batch]['addresses_beta']

    if ci_coeffs is None or addresses_alpha is None or addresses_beta is None:
        print(f"Skipping batch {batch} due to missing data.")
        continue

    # Ensure that the current batch has the same number of orbitals
    current_max_address = 0
    if len(addresses_alpha) > 0:
        current_max_address = max(current_max_address, addresses_alpha.max())
    if len(addresses_beta) > 0:
        current_max_address = max(current_max_address, addresses_beta.max())
    
    current_norb = int(current_max_address).bit_length() if current_max_address > 0 else 1
    if current_norb != num_orbitals:
        print(f"Batch {batch} has a different number of orbitals ({current_norb}) than expected ({num_orbitals}). Skipping.")
        continue

    # Filter configurations in the current batch
    filtered_configurations_idx, filtered_coeffs = filter_configurations(ci_coeffs, threshold)

    # Iterate over each passing (i, j) pair
    for (i, j), coeff in zip(filtered_configurations_idx, filtered_coeffs):
        alpha_addr = addresses_alpha[i]
        beta_addr  = addresses_beta[j]
        key = (alpha_addr, beta_addr)

        if key in merged_configs:
            merged_configs[key] += coeff
        else:
            merged_configs[key] = coeff

print(f"Total unique configurations after merging all batches: {len(merged_configs)}")

# Build bitstring_matrix and probs_array from merged configurations
bitstring_matrix, probs_array = build_bitstring_matrix_and_probs(merged_configs, num_orbitals)

print("Final bitstring_matrix shape:", bitstring_matrix.shape)
print("Final probs_array shape:", probs_array.shape)
print("Sum of probs_array:", probs_array.sum())

e_hist = np.zeros((iterations, n_batches))
s_hist = np.zeros((iterations, n_batches))
occupancy_hist = np.zeros((iterations, n_batches, 2 * num_orbitals))
occupancies_bitwise = None

identity_transition = np.array([['I'] * (2 * num_orbitals)], dtype='<U1')
if ext_sqd_config["transitions"] == 'singles_and_doubles':
    transitions = create_transitions_doubles_and_singles(num_orbitals)
elif ext_sqd_config["transitions"] == 'all_doubles':
    transitions = create_all_double_transitions(num_orbitals)
elif ext_sqd_config["transitions"] == 'singles_and_doubles_across_spin':
    transitions = create_transitions_doubles_across_spin(num_orbitals)
elif ext_sqd_config["transitions"] == 'singles_and_doubles_within_same_spin':
    transitions = create_transitions_doubles_within_same_spin(num_orbitals)
elif ext_sqd_config["transitions"] == 'singles':
    transitions = create_transitions_single(num_orbitals)
transitions = np.vstack([transitions, identity_transition])

####################################################### MAIN CALCULATION #################################################################################

rand_seed = 469420
hf_string = make_hf_bitstring(num_elec_ab,num_orbitals)

for i in range(iterations):
    print(f"Starting configuration recovery iteration {i}", flush=True)
    ti_total = time.time()

    if occupancies_bitwise is None:
        bs_mat_tmp = bitstring_matrix
        probs_arr_tmp = probs_array
    else:
        bs_mat_tmp, probs_arr_tmp = recover_configurations(
            bitstring_matrix,
            probs_array,
            occupancies_bitwise,
            num_elec_a,
            num_elec_b,
            rand_seed=rand_seed,
        )

    batches = postselect_and_subsample(
        bs_mat_tmp,
        probs_arr_tmp,
        hamming_right=num_elec_a,
        hamming_left=num_elec_b,
        samples_per_batch=samples_per_batch,
        num_batches=n_batches,
        rand_seed=rand_seed,
    )

    completed_batches = load_completed_batches(i,checkpoint_file)
    print(f"Completed batches in iteration {i}: {completed_batches}", flush=True)

    e_tmp, occs_tmp, ci_coeffs_list, addresses_list = load_iteration_from_checkpoint(i, n_batches, num_orbitals,checkpoint_file)
    if ci_coeffs_list[0] is None:
        ci_coeffs_list = [None]*n_batches
    if addresses_list[0] is None:
        addresses_list = [None]*n_batches
    
    for j in range(0, n_batches, checkpoint_interval):
        end_idx = min(j + checkpoint_interval - 1, n_batches - 1)
        batch_indices = [k for k in range(j, end_idx + 1) if k not in completed_batches]
        if not batch_indices:
            print(f"Skipping group of batches {j}-{end_idx} for iteration {i} (already completed)", flush=True)
            continue

        def process_batch_wrapper(k):
            if dice:
                return process_batch(
                    enlarge_batch_from_transitions(batches[k],transitions),
                    hcore,
                    eri,
                    nuclear_repulsion_energy,
                    open_shell,
                    k,
                    num_cpus=num_cpus,
                    hf_string=hf_string
                )
            else:
                return process_batch(
                    enlarge_batch_from_transitions(batches[k],transitions),
                    hcore,
                    eri,
                    nuclear_repulsion_energy,
                    open_shell,
                    k,
                    max_davidson_cycles,
                    hf_string=hf_string
                )
    

        if sqd_config["parallelize_batches"]:
            # Use ProcessPoolExecutor to parallelize the execution
            with concurrent.futures.ProcessPoolExecutor(max_workers=n_batches) as executor:
                futures = []
                for k in batch_indices:
                    time.sleep(2)
                    futures.append(executor.submit(process_batch_wrapper, k))
                results = [future.result() for future in concurrent.futures.as_completed(futures)]
            # Shutdown ensures all processes are completed before the main process exits.
            executor.shutdown(wait=True)

        else:
            results = [process_batch_wrapper(k) for k in batch_indices]

        batch_results = results
        for idx, result in enumerate(batch_results):
            batch_idx = batch_indices[idx]
            e_tmp[batch_idx], occs_tmp[batch_idx, :], ci_coeffs_sci, addresses = result
            ci_coeffs_list[batch_idx] = ci_coeffs_sci
            addresses_list[batch_idx] = addresses

        save_batch_checkpoint(i, j, end_idx, e_tmp, occs_tmp, ci_coeffs_list, addresses_list, checkpoint_file)
        print(f"\n\nSaved checkpoint for iteration {i}, batches {j}-{end_idx}", flush=True)
        
        print(f"Iteration {i} finished! Total computation time: {(time.time() - ti_total)/3600}",flush=True)

    e_hist[i, :] = e_tmp
    occupancy_hist[i,:,:] = occs_tmp

# Save final results including addresses and ci_coeffs
# Since these arrays differ per batch, we store them similarly in a final results file.
with h5py.File(results_file, 'w') as f:
    f.create_dataset('e_hist', data=e_hist)
    f.create_dataset('occupancy_hist', data=occupancy_hist)
    # Also store final ci_coeffs and addresses from last iteration
    # They are lists of length n_batches
    grp = f.create_group('final_iteration_data')
    for b in range(n_batches):
        # Only store if not None
        if ci_coeffs_list[b] is not None:
            grp.create_dataset(f'ci_coeffs_batch_{b}', data=ci_coeffs_list[b])
            grp.create_dataset(f'addresses_alpha_batch_{b}', data=addresses_list[b][0])#, dtype=vlen_str_dtype)
            grp.create_dataset(f'addresses_beta_batch_{b}', data=addresses_list[b][1])#, dtype=vlen_str_dtype)



# ######################################## SAVING ###############################################

# Taking minimum energy from the batches
y1 = []
for energies in e_hist:
    minimal_idx = np.argsort(energies)[0]
    y1.append(energies[minimal_idx])
energy_sqd = float(np.real(y1[-1]))
print(f"\nEnergy through iterations: {y1}",flush=True)
print(f"Energy across batches: {e_hist}")
print(f"Best SQD energy: {energy_sqd}",flush=True)
print(f"CCSD energy: {results_pre_process['ccsd_energy']}",flush=True)
if results_pre_process["fci_energy"] is not None:
    print(f"Error from FCI: {np.abs(results_pre_process['fci_energy'] - energy_sqd)}",flush=True)
print(f"Best batch: {minimal_idx}",flush=True)
print(f"Minimal energy dice subspace dimension: {ci_coeffs_list[minimal_idx].shape[0]*ci_coeffs_list[minimal_idx].shape[1]}",flush=True)

# Saving best occupancy across all batches
best_occs = occupancy_hist[-1,minimal_idx,:]
occs_folder = os.path.join(results_folder,"occs")
os.makedirs(occs_folder,exist_ok=True)
np.save(os.path.join(occs_folder,f'occs_sqd_{samples_per_batch}.npy'),best_occs)

# Getting probabilities from ground state CI Vector
print("Preparing probabilities....",flush=True)
ti = time.time()
probabilities = np.sort(np.square(ci_coeffs_list[minimal_idx]).flatten())[::-1]
idx = np.where(probabilities < 1e-15)[0] # Chop probability array for values below 1e-15
if idx.size > 0:
    arr_chopped = probabilities[:idx[0]]
else:
    arr_chopped = probabilities  
probabilities = arr_chopped
probs_folder = os.path.join(results_folder,"probs")
os.makedirs(probs_folder,exist_ok=True)
np.save(os.path.join(probs_folder,f'probabilities_{samples_per_batch}.npy'), probabilities)
print(f"Probabilities saved in {(time.time() - ti)/3600} hours",flush=True)

## Computing variance
print("Computing variance...",flush=True)
ti = time.time()
myci = pyscf.fci.selected_ci.SelectedCI()

h2e_abs = direct_spin1.absorb_h1e(hcore, eri, num_orbitals, num_elec_ab, .5)
h2e_abs = ao2mo.restore(1, h2e_abs, num_orbitals) #makes sure that the tensor has the right index permutation symms

ci_vector,addresses_alpha,addresses_beta = filter_addresses(ci_coeffs_list[minimal_idx],addresses_list[minimal_idx][0],addresses_list[minimal_idx][1])

print(f"\nFinal subspace dimension: {ci_vector.shape[0]*ci_vector.shape[1]}",flush=True)
print(f"Hilbert space fraction: {ci_vector.shape[0]*ci_vector.shape[1]/ffsim.dim(num_orbitals,tuple(num_elec_ab))}",flush=True)

def hop(c):
    cSCI = _as_SCIvector(c, (addresses_alpha,addresses_beta))
    myci.select_cuoff = None #to trick enlarge space
    myci.ci_coeff_cutoff = -1 #to trick enlarge space
    Cenlarged = enlarge_space(myci, cSCI, h2e_abs, num_orbitals, num_elec_ab)
    Cenlarged = enlarge_space(myci, Cenlarged, h2e_abs, num_orbitals, num_elec_ab) # For safety...

    hc = myci.contract_2e(h2e_abs, Cenlarged, num_orbitals, num_elec_ab)
    return hc, Cenlarged
    
Hc, Cenlarged = hop(ci_vector)
H2 = np.sum(np.conjugate(Hc) * Hc)
e_qsci = np.sum(np.conjugate(Cenlarged) * Hc)

print ('Hamiltonian variance: ' + str((H2-e_qsci**2)/e_qsci**2) + f' Time:{(time.time() - ti)/3600} hours',flush=True)

addresses = bitstring_matrix_to_ci_strs(batches[minimal_idx], open_shell=open_shell)
results_dict = {'sqd_energy':energy_sqd,
                'best_batch':int(minimal_idx),
                'samples_per_batch':samples_per_batch,
                'addresses_alpha': addresses_alpha.shape[0],
                'addresses_beta': addresses_beta.shape[0],
                "dice_subspace_dimension":ci_coeffs_list[minimal_idx].shape[0]*ci_coeffs_list[minimal_idx].shape[1],
                "pyscf_subspace_dimension":ci_vector.shape[0]*ci_vector.shape[1],
                "ground_state_support":probabilities.shape[0],
                "variance": np.float64((H2-e_qsci**2)/e_qsci**2),
                }

if results_pre_process["fci_energy"] is not None:
    results_dict['error_fci'] = np.abs(results_pre_process["fci_energy"]-energy_sqd)

dicts_folder = os.path.join(results_folder,"dicts")
os.makedirs(dicts_folder,exist_ok=True)
json.dump(results_dict,open(os.path.join(dicts_folder,f'results_dict_{samples_per_batch}.json'),'w'))

print('\nFinished!',flush=True)