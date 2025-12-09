import ffsim
import numpy as np
from qiskit.circuit import QuantumCircuit,QuantumRegister


def make_lucj_circuit(nelec:tuple,
                      num_orbitals:int,
                      t2:np.ndarray,
                      layers=1,
                      truncated_lucj=False
                      ):
    
    if nelec[0]==nelec[1]:
        alpha_alpha_indices = [(p, p + 1) for p in range(num_orbitals-1)] 
        alpha_beta_indices = [(p, p) for p in range(num_orbitals) if p % 4 == 0]
        interaction_pairs = (alpha_alpha_indices,alpha_beta_indices)
    else:
        alpha_alpha_indices = [(p, p + 1) for p in range(num_orbitals-1)] 
        alpha_beta_indices = [(p, p) for p in range(num_orbitals) if p % 4 == 0] 
        beta_beta_indices =  [(p, p + 1) for p in range(num_orbitals-1)]
        interaction_pairs = (alpha_alpha_indices,alpha_beta_indices,beta_beta_indices)    

    if nelec[0]==nelec[1]:
        if truncated_lucj:
            print("Building truncated LUCJ circut for closed shell",flush=True)
            n_reps = 2
            circ_ffsim = ffsim.UCJOpSpinBalanced.from_t_amplitudes(t2, n_reps=n_reps,interaction_pairs=interaction_pairs)
            circ_ffsim = ffsim.UCJOpSpinBalanced(
                diag_coulomb_mats=circ_ffsim.diag_coulomb_mats[:-1],
                orbital_rotations=circ_ffsim.orbital_rotations[:-1],
                final_orbital_rotation=circ_ffsim.orbital_rotations[-1])
        else:
            circ_ffsim = ffsim.UCJOpSpinBalanced.from_t_amplitudes(t2=t2, n_reps=layers,
                                                            interaction_pairs=interaction_pairs)
    else:
        if truncated_lucj:
            print("Building truncated LUCJ circuit for open shell",flush=True)
            n_reps = 2
            circ_ffsim = ffsim.UCJOpSpinUnbalanced.from_t_amplitudes(t2, n_reps=n_reps,interaction_pairs=interaction_pairs)
            circ_ffsim = ffsim.UCJOpSpinUnbalanced(
                diag_coulomb_mats=circ_ffsim.diag_coulomb_mats[:-1],
                orbital_rotations=circ_ffsim.orbital_rotations[:-1],
                final_orbital_rotation=circ_ffsim.orbital_rotations[-1])            
        
        else:
            circ_ffsim = ffsim.UCJOpSpinUnbalanced.from_t_amplitudes(t2=t2, n_reps=layers,
                                                            interaction_pairs=interaction_pairs)
        
    return circ_ffsim

def make_circuit_qiskit(circ_ffsim,
                        num_orbitals,
                        nelec,
                        basis_rotation,
                        basis='atomic'):
        qubits = QuantumRegister(2*num_orbitals, name="q")
        circ_qiskit = QuantumCircuit(qubits)
        circ_qiskit.append(ffsim.qiskit.PrepareHartreeFockJW(num_orbitals,nelec), qubits)
        
        if nelec[0]==nelec[1]:
            circ_qiskit.append(ffsim.qiskit.UCJOpSpinBalancedJW(circ_ffsim), qubits)
        else:
            circ_qiskit.append(ffsim.qiskit.UCJOpSpinUnbalancedJW(circ_ffsim), qubits)

        if basis=='atomic':
            circ_qiskit.append(ffsim.qiskit.OrbitalRotationJW(num_orbitals,basis_rotation),qubits)
        circ_qiskit.measure_all()

        return circ_qiskit