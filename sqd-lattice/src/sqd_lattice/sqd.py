import os
import h5py
import numpy as np
import time
from qiskit_addon_sqd.counts import counts_to_arrays

from qiskit_addon_sqd.fermion import solve_fermion,bitstring_matrix_to_ci_strs

def filter_addresses(ci_vector,addresses_alpha,addresses_beta):
    """ Removing padding from DICE SQD output for computing variance """
    entries_alpha = list(np.where(addresses_alpha == 0.0)[0])
    entries_beta = list(np.where(addresses_beta == 0.0)[0])

    addresses_alpha = np.delete(addresses_alpha, entries_alpha)
    addresses_beta = np.delete(addresses_beta, entries_beta)

    ci_vector = np.delete(ci_vector, entries_beta, axis=1)

    # Remove columns
    ci_vector = np.delete(ci_vector, entries_alpha, axis=0)

    return ci_vector,addresses_alpha,addresses_beta

def load_final_results(final_file):
    """
    Load data from the final HDF5 file.
    
    Args:
        final_file (str): Path to the HDF5 file.
    
    Returns:
        tuple: (e_hist, s_hist, occupancy_hist, final_data)
    """
    if not os.path.exists(final_file):
        raise FileNotFoundError("Final results file not found.")

    with h5py.File(final_file, 'r') as f:
        e_hist = f['e_hist'][:]
        # s_hist = f['s_hist'][:]
        occupancy_hist = f['occupancy_hist'][:]

        final_data = {}
        if 'final_iteration_data' in f:
            grp = f['final_iteration_data']
            # Identify batches
            batches = set()
            for key in grp.keys():
                if 'batch_' in key:
                    try:
                        batch_num = int(key.split('_')[-1])
                        batches.add(batch_num)
                    except ValueError:
                        print(f"Unexpected batch key format: {key}")
            batches = sorted(batches)

            for b in batches:
                final_data[b] = {}
                ci_key = f'ci_coeffs_batch_{b}'
                alpha_key = f'addresses_alpha_batch_{b}'
                beta_key = f'addresses_beta_batch_{b}'

                final_data[b]['ci_coeffs'] = grp[ci_key][:] if ci_key in grp else None
                final_data[b]['addresses_alpha'] = grp[alpha_key][:] if alpha_key in grp else None
                final_data[b]['addresses_beta'] = grp[beta_key][:] if beta_key in grp else None

    # return e_hist, s_hist, occupancy_hist, final_data
    return e_hist, occupancy_hist, final_data

def save_batch_checkpoint(iteration, batch_start, batch_end, e_tmp, occs_tmp, ci_coeffs_list, addresses_list, checkpoint_file):
    with h5py.File(checkpoint_file, 'a') as f:
        for batch_num in range(batch_start, batch_end + 1):
            if f.get(f'e_hist_{iteration}_batch_{batch_num}') is None:
                # Save energies, spin, and occupations as normal arrays
                f.create_dataset(f'e_hist_{iteration}_batch_{batch_num}', data=e_tmp[batch_num])
                f.create_dataset(f'occ_hist_{iteration}_batch_{batch_num}', data=occs_tmp[batch_num, :])
                
                # Extract data from lists
                ci_coeffs_sci = ci_coeffs_list[batch_num]  # e.g., a 2D float array
                addresses_alpha = addresses_list[batch_num][0]  # array of ints
                addresses_beta = addresses_list[batch_num][1]   # array of ints

                # Convert to np arrays if needed
                # addresses_alpha = np.array(addresses_alpha, dtype=np.int64)
                # addresses_beta = np.array(addresses_beta, dtype=np.int64)

                # If ci_coeffs_sci is a simple 2D NumPy array (and shape differs batch-to-batch is fine), no vlen needed:
                # Just do:
                f.create_dataset(f'ci_coeffs_{iteration}_batch_{batch_num}', data=ci_coeffs_sci)
                
                # Create datasets for addresses with vlen dtype (not strictly needed if we trust each dataset to have its own shape)
                f.create_dataset(f'addresses_alpha_{iteration}_batch_{batch_num}', data=addresses_alpha)#, dtype=vlen_int)
                f.create_dataset(f'addresses_beta_{iteration}_batch_{batch_num}', data=addresses_beta)#, dtype=vlen_int)

def load_completed_batches(iteration,checkpoint_file):
    completed_batches = set()
    if os.path.exists(checkpoint_file):
        with h5py.File(checkpoint_file, 'r') as f:
            for key in f.keys():
                if key.startswith(f'e_hist_{iteration}_batch_'):
                    # Extract batch number
                    parts = key.split('_')
                    batch_num = int(parts[-1])
                    completed_batches.add(batch_num)
    return completed_batches

def load_iteration_from_checkpoint(iteration, n_batches, num_orbitals,checkpoint_file):
    e_tmp = np.zeros(n_batches)
    occs_tmp = np.zeros((n_batches, 2 * num_orbitals))
    ci_coeffs_list = [None]*n_batches
    addresses_list = [None]*n_batches
    if os.path.exists(checkpoint_file):
        with h5py.File(checkpoint_file, 'r') as f:
            for batch_num in range(n_batches):
                if f.get(f'e_hist_{iteration}_batch_{batch_num}') is not None:
                    e_tmp[batch_num] = f[f'e_hist_{iteration}_batch_{batch_num}'][()]
                    occs_tmp[batch_num, :] = f[f'occ_hist_{iteration}_batch_{batch_num}'][:]

                    # Load ci_coeffs and addresses if present
                    ci_coeff_key = f'ci_coeffs_{iteration}_batch_{batch_num}'
                    addr_alpha_key = f'addresses_alpha_{iteration}_batch_{batch_num}'
                    addr_beta_key = f'addresses_beta_{iteration}_batch_{batch_num}'
                    if (ci_coeff_key in f.keys()) and (addr_alpha_key in f.keys()) and (addr_beta_key in f.keys()):
                        ci_coeffs_list[batch_num] = f[ci_coeff_key][:]
                        alpha_addrs = f[addr_alpha_key][:]
                        beta_addrs = f[addr_beta_key][:]
                        addresses_list[batch_num] = (alpha_addrs, beta_addrs)
    return e_tmp, occs_tmp, ci_coeffs_list, addresses_list

def process_batch(batch, hcore, eri, nuclear_repulsion_energy, open_shell,batch_idx,max_davidson_cycles,hf_string):
    print(f"\nProcessing batch {batch_idx}...", flush=True)
    if not np.any(np.all(batch == hf_string, axis=1)):    
        batch = np.vstack([batch,hf_string])
    addresses = bitstring_matrix_to_ci_strs(batch, open_shell=open_shell)    
    print(f"Subspace dimension: {len(addresses[0]) * len(addresses[1])}", flush=True)
    
    print(f"Open shell:{open_shell}")
    ti = time.time()
    energy_sci, state_sci, avg_occs,spin = solve_fermion(
        batch,
        hcore,
        eri,
        open_shell=open_shell,
        spin_sq=None,
        max_davidson=max_davidson_cycles
    )

    tf = time.time()
    energy_sci += nuclear_repulsion_energy
    occs_tmp = np.append(avg_occs[0], avg_occs[1])

    print(f"Solution alpha strings: {state_sci.ci_strs_a.shape}",flush=True)
    print(f"Solution beta strings: {state_sci.ci_strs_b.shape}",flush=True)
    print(f"Coeffs shape: {state_sci.amplitudes.shape}",flush=True)
    print(f"Dice subspace Dimension: {state_sci.amplitudes.shape[0]*state_sci.amplitudes.shape[1]}",flush=True)
    
    print(f"\nBatch {batch_idx} finished! Computation time: {(tf-ti)/3600} hours",flush=True)
    print(f"SQD Energy: {energy_sci}",flush=True)
    
    return energy_sci, occs_tmp, state_sci.amplitudes, [state_sci.ci_strs_a,state_sci.ci_strs_b]

def make_hf_bitstring(ne_ab, num_orbitals):
    ne_a, na_b = ne_ab
    
    if not (0 <= ne_a <= num_orbitals and 0 <= na_b <= num_orbitals):
        raise ValueError("t1 and t2 must be between 0 and N inclusive.")
    
    left = "0" * (num_orbitals - na_b) + "1" * na_b
    right = "0" * (num_orbitals - ne_a) + "1" * ne_a
    
    bitstring_list = counts_to_arrays({left + right:1})

    return bitstring_list[0][0]
