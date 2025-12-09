import pyscf
import numpy as np
import ffsim 
from sqd_lattice.util import make_pyscf_slater_det
import time
from ffsim import FermionOperator,cre_a,cre_b,des_a,des_b,linear_operator

def make_two_body_tensor(data:dict) -> np.ndarray:
    num_orbitals = data["hopping_matrix"].shape[0]
    two_body_tensor = np.zeros((num_orbitals,num_orbitals,num_orbitals,num_orbitals))
    for k,v_dict in data["hubbard_params"].items():
        if k[0]==k[1]:
            for i in v_dict["index_i"]:
                two_body_tensor[i,i,i,i] = 2*v_dict["V"]
        else:
            for i in v_dict["index_i"]:
                for j in v_dict["index_j"]:
                    two_body_tensor[i,i,j,j] = 2*v_dict["V"]
                    two_body_tensor[j,j,i,i] = 2*v_dict["V"]
    return two_body_tensor

def run_pyscf(data:dict,logger,nelec=None,compute_fci=False) -> dict:

    logger.info(f"Starting calculations for nelec: {nelec}")
    
    mol = pyscf.gto.M(verbose=0)
    mol.nelec = nelec
    mol.incore_anyway = True

    num_orbitals = data["hopping_matrix"].shape[0]

    two_body_tensor = make_two_body_tensor(data)

    ## Hartree fock calculation
    hf_obj = pyscf.scf.RHF(mol)
    hf_obj.get_hcore = lambda *args: np.real(data["hopping_matrix"])#*(1/27.2114079527))
    hf_obj.get_ovlp = lambda *args: np.eye(num_orbitals)
    hf_obj._eri = two_body_tensor
    pyscf_slater_det = make_pyscf_slater_det(nelec,num_orbitals)
    logger.info(pyscf_slater_det)
    hf_obj.get_occ = lambda *args:pyscf_slater_det
    logger.info('\nRunning Hartree Fock...')
    ti = time.time()
    hf_obj.kernel()
    logger.info(f"Hartree fock calculation done! Energy: {hf_obj.e_tot}, Time: {time.time()-ti} seconds")

    ## CCSD calculation
    ccsd_obj = pyscf.cc.CCSD(hf_obj)
    logger.info("\nRunning CCSD...")
    ti = time.time()
    while ccsd_obj.max_cycle > 0:
        try:
            ccsd_obj.kernel() 
            logger.info("CCSD ran successfully.")
            break
        except np.linalg.LinAlgError:
            logger.info(f"LinAlgError encountered with max_cycles = {ccsd_obj.max_cycle}. Retrying...")
            ccsd_obj.max_cycle -= 1
    else:
        logger.info("Failed to run method: max_cycles reduced to 0.")
    logger.info(f"CCSD calculation done! Energy: {ccsd_obj.e_tot}, Correction over HF: {ccsd_obj.e_corr} Time: {time.time() - ti} seconds")
    ccsd_t_correction = ccsd_obj.ccsd_t()
    logger.info(f"CCSD(T) calculation done! Energy: {ccsd_obj.e_tot + ccsd_t_correction}, Correction over CCSD: {ccsd_t_correction}")

    results_pre_process = dict(hf_energy=hf_obj.e_tot,
                        two_body_tensor=two_body_tensor,
                        ne_ab=nelec,
                        mo_coeff=hf_obj.mo_coeff,
                        ccsd_energy=ccsd_obj.e_tot,
                        t1=ccsd_obj.t1,
                        t2=ccsd_obj.t2,
                        ccsdt_energy = ccsd_obj.e_tot + ccsd_t_correction,
                        fci_energy=None)

    if compute_fci:
    ## FCI calculation
        fci_obj = pyscf.fci.FCI(hf_obj)
        logger.info("\nRunning FCI...")
        ti = time.time()
        fci_obj.kernel()
        logger.info(f"FCI calculation done! Energy: {fci_obj.e_tot}    Computation time: {time.time()-ti} seconds") 
        results_pre_process["fci_energy"] = fci_obj.e_tot

    logger.info("\nFinished PYSCF calculations!")

    return results_pre_process


def compute_exp_vals_lucj(nelec,num_orbitals,linear_op_ffsim, circ_ffsim,logger):
    logger.info("Computing Slater det energy...")
    slater_det_state = ffsim.hartree_fock_state(num_orbitals,nelec)
    ti = time.time()
    slater_det_energy = np.real(np.vdot(slater_det_state, linear_op_ffsim @ slater_det_state))
    logger.info(f"Slater determinant energy: {slater_det_energy}, Time: {time.time()-ti} seconds")
    logger.info(f"Computing LUCJ initialization energy...")
    statevector_lucj = ffsim.apply_unitary(slater_det_state,circ_ffsim, norb=num_orbitals, nelec=nelec,copy=True)
    lucj_init_energy = np.real(np.vdot(statevector_lucj, linear_op_ffsim @ statevector_lucj))
    logger.info(f"LUCJ initit energy: {slater_det_energy}, Time: {time.time()-ti} seconds")
    logger.info(f"Finished building circuits!")
    
    return slater_det_energy,lucj_init_energy

def get_ucj_energy(self,circuit_ffsim,operator,reference_state):
    statevector = ffsim.apply_unitary(reference_state,circuit_ffsim, norb=self.num_orbitals, nelec=self.nelec,copy=True)
    ti = time.time()
    self.logger.info('Computing CCSD initialization energy...')
    self.data_dict['ucj_init_energy'] = np.real(np.vdot(statevector, operator @ statevector))

    tf = time.time()
    self.logger.info(f"UCJ initialization energy: {self.data_dict['ucj_init_energy']}  Evaluation time {tf - ti} seconds")
    return statevector

def rotate_tensors(hopping_matrix,two_body_tensor,basis_rotation):
    hopping_matrix_rot = basis_rotation.T @hopping_matrix @ basis_rotation
    two_body_tensor_rot = np.einsum('pqrs,pi,qj,rk,sl->ijkl', two_body_tensor, basis_rotation, basis_rotation, basis_rotation, basis_rotation)
    return hopping_matrix_rot,two_body_tensor_rot

def make_hamiltonian_ffsim(data):
    # Building the Hamiltonian from scratch with fermionic operators
    coeffs: dict[tuple[tuple[bool, bool, int], ...], complex] = defaultdict(float)

    for k,v_dict in data["hubbard_params"].items():
        if k[0]==k[1]:
            for i in v_dict["index_i"]:
                coeffs[cre_a(i), des_a(i), cre_b(i), des_b(i)] = v_dict["V"]
        else:
            for i in v_dict["index_i"]:
                for j in v_dict["index_j"]:
                    coeffs[cre_a(i), des_a(i), cre_a(j), des_a(j)] = v_dict["V"]
                    coeffs[cre_a(i), des_a(i), cre_b(j), des_b(j)] = v_dict["V"]
                    coeffs[cre_b(i), des_b(i), cre_a(j), des_a(j)] = v_dict["V"]
                    coeffs[cre_b(i), des_b(i), cre_b(j), des_b(j)] = v_dict["V"]

    for i in range(data["num_orbitals"]):
        for j in range(data["num_orbitals"]):
            coeffs[cre_a(i),des_a(j)] = data["hopping_matrix"][i,j]
            coeffs[cre_b(i),des_b(j)] = data["hopping_matrix"][i,j]

    return FermionOperator(coeffs)