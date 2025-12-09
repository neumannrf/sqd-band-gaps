# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import warnings
warnings.filterwarnings("ignore") 
from PAOFLOW import PAOFLOW
import sys
     
##
# 1st argument: .save from QE
# 2nd argument: .gnu.dat band structure 
# 3rd argument: Fermi energy

plt.rcParams.update({'font.size': 20})

# %% [markdown]
# ### Load arry files for two MOFs as example

# %%
material = PAOFLOW.PAOFLOW(savedir=sys.argv[1], model=None, outputdir='.', 
                          smearing='gauss',verbose=True)
#paoflow.restart_load('dft-dump')
data_controller = material.data_controller
arry,attr = material.data_controller.data_dicts()
material.projections()
material.projectability(pthr=0.95)
material.pao_hamiltonian(write_binary=True)
# data_controller_basis = paoflow_basis.data_controller
# arry_basis,attr_basis = paoflow_basis.data_controller.data_dicts()
basis = arry['basis']
#paoflow.restart_dump('dft-dump')
#paoflow.finish_execution()

# %% [markdown]
# ### Hooping matrix

# %%
plt.figure(figsize=(10, 10))

# plot only the Gamma point
matrix = arry['Hks'][:,:,0,0,0,0]
factor = 0.5
coeff = True

list_of_orbitals = [f"{item['atom']}+{item['label']}" for item in basis]

absolute_matrix = np.array([[np.abs(item)**2 for item in row] for row in matrix])

# Create a list of labels corresponding to the elements
unique_elements = sorted(set(label.split('+')[0] for label in list_of_orbitals))
element_to_indices = {element: [i for i, label in enumerate(list_of_orbitals) if label.startswith(element)] for element in unique_elements}

# Create a list of orbital labels and their indices
orbital_labels = [label.split('+')[1][1:] for label in list_of_orbitals]
unique_orbitals = sorted(set(orbital_labels))

# Group labels
grouped_labels = []
grouped_index = []
current_label = None

for index,label in enumerate(orbital_labels):
    if label != current_label:
        grouped_labels.append(label)
        grouped_index.append(index)
        current_label = label

orbital_to_indices = {orbital: [i for i, orbital_label in enumerate(orbital_labels) if orbital_label == orbital] for orbital in unique_orbitals}
orbital_positions = {orbital: (min(indices) + max(indices)) / 2 for orbital, indices in orbital_to_indices.items()}

np.fill_diagonal(absolute_matrix,0)
vmin = np.min(abs(absolute_matrix))
vmax = np.max(abs(absolute_matrix))

absolute_matrix = abs(absolute_matrix)/vmax


plt.imshow(absolute_matrix, cmap='plasma', aspect='equal', vmax=factor, vmin=0)
cbar = plt.colorbar(label='$|t_{ij}|^{2}$| (eV$^2$)', shrink=0.5)
cbar.ax.text(1, 1.1, r'$\|t_{ij}\|^{2}_{max}$ = ' + f'{vmax:.2f} eV', transform=cbar.ax.transAxes, ha='left', va='top', fontsize=12)

for element, indices in element_to_indices.items():
    min_index = min(indices)
    max_index = max(indices)
    rect = patches.Rectangle((min_index - 0.5, min_index - 0.5), len(indices), len(indices),
                            linewidth=2, edgecolor='white', facecolor='none')
    plt.axhline(y=min_index - 0.5, color='white', linewidth=1, linestyle='--')
    plt.axhline(y=max_index + 0.5, color='white', linewidth=1, linestyle='--')
    plt.axvline(x=min_index - 0.5, color='white', linewidth=1, linestyle='--')
    plt.axvline(x=max_index + 0.5, color='white', linewidth=1, linestyle='--')
    plt.gca().add_patch(rect)
    plt.text((min_index + max_index) / 2, (min_index + max_index) / 2, element, color='white', 
            fontsize=18, ha='center', va='center', weight='bold')
    

plt.xticks(grouped_index, grouped_labels, rotation=45, ha='center', fontsize=18)
plt.yticks(grouped_index, grouped_labels, rotation=-45, ha='right', fontsize=18)
plt.yticks([])
# plt.title('$\|t_{ij}\|^{2}_{max}$ = ' + f'{vmax:.2f} eV')
plt.ylabel('Orbitals')
plt.xlabel('Orbitals')
plt.tight_layout()
plt.savefig(sys.argv[4]+'_hopping.pdf')
plt.show()


# %%

plt.figure(figsize=(10, 10))


path = 'G-L' 
sym_points = {'G':[0.0, 0.0, 0.0], \
              'L':[0.5, 0.5, 0.5]}
material.bands(ibrav=2, nk =5000 , band_path=path, high_sym_points=sym_points)

data_dft = np.loadtxt(sys.argv[2])

k = np.unique(data_dft[:,0])
kbands = np.reshape(data_dft[:, 1], (-1,len(k)))
k = k/max(k)

for band in range(len(kbands)):
    plt.plot(k, kbands[band,:] -float(sys.argv[3]), '-',linewidth=5, alpha=0.5, color='red', markersize=2)
    # plt.plot(k, kbands[band,:], '-',linewidth=5, alpha=0.5, color='red', markersize=2)


# plot the paoflow bands
eband = arry['E_k']
for ib in range(eband.shape[1]):
    plt.plot(np.linspace(0,1,eband.shape[0]),eband[:,ib],'--',color='black', markersize=0.5)


plt.plot([],[],'-r', label='DFT')
plt.plot([],[],'--k', label='PAOFLOW')
plt.legend()

# High symmetry k-points (check bands_pp.out)
plt.axvline(0, linewidth=0.75, color='k', alpha=0.5)
plt.axvline(1, linewidth=0.75, color='k', alpha=0.5)


plt.xlim(0, 1)
plt.ylabel("Energy (eV)")
# plt.ylim(5,20)

# Annotate the band gap
plt.tight_layout()
plt.savefig(sys.argv[4]+'_bands.pdf')
plt.show()   


