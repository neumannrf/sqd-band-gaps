import os 
import argparse
import numpy as np
import json
import time
import h5py
import ffsim
import pyscf
from pyscf import ao2mo
#careful with the next few imports. 
# We are exposing methods that are not supposed to be visible by the user
from pyscf.fci.selected_ci import _all_linkstr_index, _as_SCIvector
from pyscf.fci import direct_spin1
from pyscf.fci.selected_ci import enlarge_space

from sqd_lattice.util import load_yaml,load_pickle

from qiskit_addon_sqd.counts import counts_to_arrays
from qiskit_addon_sqd.configuration_recovery import recover_configurations
from qiskit_addon_sqd.subsampling import postselect_and_subsample
from qiskit_addon_sqd.fermion import bitstring_matrix_to_ci_strs

from sqd_lattice.sqd import filter_addresses

from qiskit_addon_dice_solver import solve_hci

# # Loading material name
parser = argparse.ArgumentParser()
parser.add_argument('--material')
parser.add_argument('--ne',help="Total number of electrons to index file")
parser.add_argument('--select_cutoff_idx',help="Select cutoff for HCI, will be iterated over")
args = parser.parse_args()

material = "zirconia"
ne = 24
select_cutoff_idx = int(6)

folder_path = f"runs/{material}/{ne}e"

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
results_pre_process = load_pickle(folder_path,"results_pre_process.pkl")
hci_config = load_yaml(folder_path,"hci_config.yaml")
max_iter = hci_config["max_iter"]

# Degfining variables from data dicts
hcore = np.load("/home/alanduriez/sqd-lattice/notebooks/hcore_rot.npy")
eri = np.load("/home/alanduriez/sqd-lattice/notebooks/eri_pyscf_rot.npy").real

num_orbitals = 20
num_elec_ab = (28,28)
num_elec_a = num_elec_ab[0]
num_elec_b = num_elec_ab[1]
open_shell = False if num_elec_a==num_elec_b else True

print(f"Starting HCI  for {material} {sum(num_elec_ab)} electrons",flush=True)
# print(hci_config)

# Making the list of select cutoffs
sel_cutoff_range = hci_config["select_cutoffs"]
sel_cutoff_step = hci_config["select_cutoff_steps"]
select_cutoff_list = np.linspace(*sel_cutoff_range,sel_cutoff_step)

ti = time.time()
energy_sci, state_sci, avg_occs = solve_hci(
    hcore=hcore,
    eri=eri,
    norb=num_orbitals,
    nelec=num_elec_ab,
    ci_strs=None,
    spin_sq=None,
    energy_tol=1e-12,
    mpirun_options=["-quiet", "-n", f"{1}"],
    select_cutoff=1e-13,
    max_iter=10)

print(f"HCI calculation finished! Time :{time.time()-ti}")
print(f"HCI energy: {energy_sci}")

results_folder = os.path.join(folder_path,"results_hci")

# Saving occupancies
occs_folder = os.path.join(results_folder,"occs")
os.makedirs(occs_folder,exist_ok=True)
np.save(os.path.join(occs_folder,f'occs_hci_{max_iter}.npy'),avg_occs)

# # Saving probs
probabilities = np.sort(np.square(state_sci.amplitudes).flatten())[::-1]
idx = np.where(probabilities < 1e-15)[0] # Chop probability array for values below 1e-15
if idx.size > 0:
    arr_chopped = probabilities[:idx[0]]
else:
    arr_chopped = probabilities  
probabilities = arr_chopped
probs_folder = os.path.join(results_folder,"probs")
os.makedirs(probs_folder,exist_ok=True)
np.save(os.path.join(probs_folder,f'probabilities_{max_iter}.npy'), probabilities)

## Computing variance
print("Computing variance...",flush=True)
ti = time.time()
myci = pyscf.fci.selected_ci.SelectedCI()

h2e_abs = direct_spin1.absorb_h1e(hcore, eri, num_orbitals, num_elec_ab, .5)
h2e_abs = ao2mo.restore(1, h2e_abs, num_orbitals) #makes sure that the tensor has the right index permutation symms

ci_vector,addresses_alpha,addresses_beta = filter_addresses(state_sci.amplitudes,state_sci.ci_strs_a,state_sci.ci_strs_b)

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
variance = (H2-e_qsci**2)/e_qsci**2

print ('Hamiltonian variance: ' + str(variance) + f' Time:{(time.time() - ti)/3600} hours',flush=True)

# # Saving results dict
results_dict = {'hci_energy':energy_sci,
                'max_iter':max_iter,
                'addresses_alpha': addresses_alpha.shape[0],
                'addresses_beta': addresses_beta.shape[0],
                "subspace_dimension":ci_vector.shape[0]*ci_vector.shape[1],
                "ground_state_support":probabilities.shape[0],
                "variance": np.float64(variance),
                }

if results_pre_process["fci_energy"] is not None:
    error_fci = np.abs(results_pre_process["fci_energy"]-energy_sci)
    results_dict['error_fci'] = error_fci
    print(f"Error from FCI: {error_fci}",flush=True)

dicts_folder = os.path.join(results_folder,"dicts")
os.makedirs(dicts_folder,exist_ok=True)
json.dump(results_dict,open(os.path.join(dicts_folder,f'results_dict_{select_cutoff_idx}.json'),'w'))

print("\nFinished!")


    


