import warnings
warnings.filterwarnings("ignore") 
import numpy as np
from PAOFLOW import PAOFLOW


paoflow_basis = PAOFLOW.PAOFLOW(savedir='out/#MATERIAL.save', model=None, outputdir='.', 
                          smearing='gauss',verbose=True)
data_controller = paoflow_basis.data_controller
arry,attr = paoflow_basis.data_controller.data_dicts()
paoflow_basis.projections()
paoflow_basis.projectability(pthr=0.95)
paoflow_basis.pao_hamiltonian(write_binary=True)
data_controller_basis = paoflow_basis.data_controller
arry_basis,attr_basis = paoflow_basis.data_controller.data_dicts()

print(f"Shape of the Hamiltonian tensor: {arry['Hks'].shape}\n")
print(f"Length of the basis vector: {len(arry['basis'])}\n")

np.save('#MATERIAL_basis_#PSEUDO.npy', arry['basis'])
np.save('#MATERIAL_hopping_#PSEUDO.npy', arry['Hks'])
