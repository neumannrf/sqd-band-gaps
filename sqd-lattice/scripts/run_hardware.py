import os
import argparse
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from qiskit_ibm_runtime.options import SamplerOptions,DynamicalDecouplingOptions,EnvironmentOptions
from collections import defaultdict

from qiskit import qpy

from sqd_lattice.util import load_pickle,load_yaml,save_json

# Loading material name
parser = argparse.ArgumentParser()
parser.add_argument('--material')
parser.add_argument('--ne',help="Total number of electrons to index file")
args = parser.parse_args()

material = args.material
ne = args.ne

folder_path = f"runs/{material}/{ne}e"

base_config = load_yaml(f"runs/{material}","base_config.yaml")
hardware_config = load_yaml(f"runs/{material}","hardware_config.yaml")
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
results_pre_process = load_pickle(folder_path,"results_pre_process.pkl")

print(f"Collecting hardware counts for {material} {ne} electrons ",flush=True)
print(hardware_config)

# Loading transpiled circuit
transpiled = qpy.load(open(os.path.join(folder_path,f'transpiled_hardware.qpy'),'rb'))[0]

# Setting options
service = QiskitRuntimeService(name=hardware_config["instance"])
backend = service.backend(hardware_config["processor"],use_fractional_gates=hardware_config["fractional_gates"])

# Setting options
options = SamplerOptions()
if hardware_config["dynamical_decoupling"]:
    dd_options = DynamicalDecouplingOptions(enable=True, sequence_type="XpXm",skip_reset_qubits=True)
    options.dynamical_decoupling = dd_options

job_tags = [backend.name,f'{material}',f'{ne}e',f'{base_config["basis"]}_basis']
if hardware_config["fractional_gates"]:
    job_tags.append('fractional_gates')
if hardware_config["dynamical_decoupling"]:
    job_tags.append('dynamical_decoupling')

options.environment = EnvironmentOptions(job_tags=job_tags)

pub = (transpiled,)
# Sample job execution
sampler = SamplerV2(mode=backend,options=options)
n_pubs = hardware_config["n_pubs"]
pubs = [(transpiled,None,hardware_config["shots_per_pub"]) for _ in range(n_pubs)]
job = sampler.run(pubs)
print(job.job_id(),flush=True)

job_result = job.result()

merged_dict = defaultdict(int)

# List of dictionaries to merge
dicts = [job_result[i].data.c.get_counts() for i in range(n_pubs)]

# Iterate over each dictionary and sum values for common keys
for d in dicts:
    for key, value in d.items():
        merged_dict[key] += value

# Convert back to a regular dictionary if needed
counts = dict(merged_dict)

# Write dictionary to file
name = 'counts_hardware'
if hardware_config["fractional_gates"]:
    name+='_frac'

counts_folder = os.path.join(folder_path,"counts")
os.makedirs(counts_folder,exist_ok=True)

save_json(counts, os.path.join(counts_folder,name+'.json'))

print("Finished!")