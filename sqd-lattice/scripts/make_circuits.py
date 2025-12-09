import os 
import pickle
import numpy as np
import matplotlib.pyplot as plt
import importlib.util
import yaml
import logging
import argparse
import ffsim
import time

from qiskit import transpile
from qiskit.transpiler import CouplingMap,generate_preset_pass_manager
from qiskit import qasm3

from sqd_lattice.pre_process import run_pyscf,make_two_body_tensor,compute_exp_vals_lucj,rotate_tensors
from sqd_lattice.circuits import make_lucj_circuit,make_circuit_qiskit
from sqd_lattice.util import save_pickle

# Loading material name
parser = argparse.ArgumentParser()
parser.add_argument('--material')
args = parser.parse_args()

material = args.material

folder_path = f"runs/{material}"

# Loading files for DFT data
orbital_basis = np.load(os.path.join(folder_path,"dft_data", f"basis.npy"), allow_pickle=True) # Orbital ordering
hopping_matrix = np.load(os.path.join(folder_path,"dft_data", f"hopping_matrix.npy"))
num_orbitals = hopping_matrix.shape[0]

# Loading configs input file
with open(os.path.join(folder_path,"base_config.yaml"), "r") as f:
    config = yaml.safe_load(f)

# Load U and V parameters from python script
script_path = os.path.join(folder_path,"dft_data", "hubbard_params.py")
module_name = "hubbard_params"
spec = importlib.util.spec_from_file_location(module_name, script_path)
hubbard_params = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hubbard_params)

# U_hubbard = hubbard_params.U_hubbard
# V_hubbard = hubbard_params.V_hubbard
ne_ab = hubbard_params.ne_ab

hubbard_params = np.load(os.path.join(folder_path,"dft_data","hubbard_params.npy"),allow_pickle=True)

print(f"Starting calculations for Material: {args.material}\n",flush=True)
# print("U Hubbard: ",U_hubbard)
# print("V Hubbard: ",V_hubbard)
print(hubbard_params)

print(config)

data = dict(material=material,
              folder_path=folder_path,
              orbital_basis=orbital_basis,
              num_orbitals=hopping_matrix.shape[0],
              hopping_matrix=hopping_matrix,
            #   U_hubbard=U_hubbard,
            #   V_hubbard=V_hubbard,
              hubbard_params=hubbard_params.item(),
              ne_ab=ne_ab)

data["two_body_tensor"] = make_two_body_tensor(data)

pickle.dump(data,open(os.path.join(folder_path,"data_dict.pkl"),"wb"))

# Running everything for the three electron numbers
ne_ab_list = [(ne_ab[0],ne_ab[0]-1),ne_ab,(ne_ab[0]+1,ne_ab[0])][:]

for ne_ab in ne_ab_list:
    ne_folder_path = os.path.join(folder_path,f"{sum(ne_ab)}e")
    os.makedirs(os.path.join(ne_folder_path,"logs"),exist_ok=True)
    
    # Writing SQD config file
    sqd_config = dict(samples_per_batch=list(config["sqd"]["sqd_samples_per_batch"]),
                      sampling=config["sqd"]["sampling"],
                      n_batches=1,
                      recovery_iterations=1,
                      num_cpus=1,
                      max_davidson_cycles=800,
                      parallelize_batches=False,
                      ext_sqd=dict(transitions="singles",
                                   samples_per_batch_old=int(config["sqd"]["sqd_samples_per_batch"][-1]),
                                   sampling=config["sqd"]["sampling"],
                                   samples_per_batch=[250],
                                   amplitude_threshold=1e-3,
                                   n_batches=1,
                                   recovery_iterations=1,
                                   num_cpus=50))
    with open(os.path.join(ne_folder_path,"sqd_config.yaml"), 'w') as file:
        yaml.dump(sqd_config, file, default_flow_style=False)
    
    # Writing HCI config file
    hci_config = dict(select_cutoffs=list(config["hci"]["select_cutoffs"]),
                      energy_tol=float(config["hci"]["energy_tol"]),
                      max_iter=config['hci']['max_iter'],
                      num_cpus=int(config['hci']['num_cpus']),
                      select_cutoff_steps=int(config['hci']['select_cutoff_steps']))
    with open(os.path.join(ne_folder_path,"hci_config.yaml"), 'w') as file:
        yaml.dump(hci_config, file, default_flow_style=False)
    
    # Making loggers
    log_path = os.path.join(ne_folder_path,"logs", "log_pyscf_ffsim.log")
    logger = logging.getLogger(log_path)  # Unique name for each logger
    logger.setLevel(logging.INFO)

    # Avoid adding multiple handlers to the same logger
    if not logger.handlers:
        fh = logging.FileHandler(log_path,mode='w')
        formatter = logging.Formatter("%(message)s")
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    
    print(f"\nComputing classical benchmarks and making circuits for {sum(ne_ab)} electrons...",flush=True)

    results_pre_process = run_pyscf(data,
                              logger=logger,
                              nelec=ne_ab,
                              compute_fci=config["pyscf"]["compute_fci"])
    
    
    print("Classical benchmarks done! Building circuits and operators...",flush=True)

    print('Changing basis of one and two body tensors...',flush=True)
    ti = time.time()
    hopping_matrix_rot,two_body_tensor_rot = rotate_tensors(hopping_matrix,
                                                            data["two_body_tensor"],
                                                            results_pre_process["mo_coeff"])
    print(f"Change of basis done! Time {time.time() - ti} seconds")
    results_pre_process["hopping_matrix_rot"] = hopping_matrix_rot
    results_pre_process["two_body_tensor_rot"] = two_body_tensor_rot

    # Building Hamiltonian in ffsim
    hamiltonian = ffsim.MolecularHamiltonian(one_body_tensor=hopping_matrix,two_body_tensor=data["two_body_tensor"])
    hamiltonian = hamiltonian.rotated(results_pre_process["mo_coeff"].T) # Changing to molecular

    circ_ffsim = make_lucj_circuit(ne_ab,
                                   data["num_orbitals"],
                                   results_pre_process["t2"],
                                   layers=config["ffsim"]["layers_lucj"],
                                   truncated_lucj=config["ffsim"]["truncated_lucj"])
    
    results_pre_process["circ_ffsim"] = circ_ffsim

    # Storing transpiled qiskit circuit for MPS simulator
    circ_qiskit = make_circuit_qiskit(circ_ffsim,num_orbitals,results_pre_process["ne_ab"],results_pre_process["mo_coeff"],basis=config["basis"])

    pass_manager = generate_preset_pass_manager(optimization_level=3,
                                                basis_gates=['cx', 'rx', 'ry', 'rz'],
                                                approximation_degree=1,
                                                coupling_map=CouplingMap.from_line(circ_qiskit.num_qubits),
                                                seed_transpiler=12345
    )
    
    if config["ffsim"]["compute_exp_vals"]:
        print("Computing exact expectation values of circuits...")
        linear_op_ffsim = ffsim.linear_operator(hamiltonian,norb=num_orbitals,nelec=ne_ab)
        results_pre_process["slater_det_energy"],results_pre_process["lucj_init_energy"] = compute_exp_vals_lucj(ne_ab,
                                                                                                                 num_orbitals,
                                                                                                                 linear_op_ffsim,
                                                                                                                 circ_ffsim,
                                                                                                                 logger)
        print("Expectation values computed!",flush=True)
    
    save_pickle(results_pre_process,ne_folder_path,"results_pre_process.pkl")

    logger.info("Finished!")


print("\nFinished!",flush=True)



