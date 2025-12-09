import os
import h5py
import numpy as np
import time

from qiskit_addon_sqd.fermion import bitstring_matrix_to_ci_strs
from qiskit_addon_dice_solver import solve_fermion

def process_batch_dice(batch, hcore, eri, nuclear_repulsion_energy, open_shell,batch_idx,num_cpus,hf_string):
    print(f"\nProcessing batch {batch_idx}...", flush=True)
    if not np.any(np.all(batch == hf_string, axis=1)):    
        print("Adding hf string to bacth",flush=True)
        batch = np.vstack([batch,hf_string])
    addresses = bitstring_matrix_to_ci_strs(batch, open_shell=open_shell)
    print(f"Subspace dimension: {len(addresses[0]) * len(addresses[1])}", flush=True)
    print(f"Open shell:{open_shell}")
    ti = time.time()
    energy_sci, state_sci, avg_occs = solve_fermion(
        batch,
        hcore,
        eri,
        open_shell=open_shell,
        mpirun_options=["-quiet","-n", f"{num_cpus}"],
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
