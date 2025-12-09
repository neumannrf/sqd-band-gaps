import os
import numpy as np
import ffsim 
import argparse
import json

from sqd_lattice.util import load_pickle,load_yaml
from sqd_lattice.circuits import make_circuit_qiskit

import time

# Loading material name
parser = argparse.ArgumentParser()
parser.add_argument('--material')
parser.add_argument('--ne',help="Total number of electrons to index file")
args = parser.parse_args()

material = args.material
ne = args.ne

folder_path = f"runs/{material}/{ne}e"

print(f"Sampling circuit (ffsim) for {material} {ne} electrons",flush=True)

base_config = load_yaml(f"runs/{material}","base_config.yaml")
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
results_pre_process = load_pickle(folder_path,"results_pre_process.pkl")

num_orbitals = data_dict["num_orbitals"]

circ_ffsim = results_pre_process["circ_ffsim"]

circ_qiskit = make_circuit_qiskit(circ_ffsim,num_orbitals,results_pre_process["ne_ab"],results_pre_process["mo_coeff"])

print("Sampling circuit...")
sampler = ffsim.qiskit.FfsimSampler(default_shots=base_config["sampler_ffsim"]["shots"], seed=12345)
pub = (circ_qiskit,)
job = sampler.run([pub])
result = job.result()
pub_result = result[0]
counts = pub_result.data.meas.get_counts()

counts_folder = os.path.join(folder_path,"counts")
os.makedirs(counts_folder,exist_ok=True)
json.dump(counts,open(os.path.join(counts_folder,"counts_ffsim.json"),"w"))

print("Finished!",flush=True)
