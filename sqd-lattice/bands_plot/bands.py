# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import warnings
warnings.filterwarnings("ignore") 
from PAOFLOW import PAOFLOW
import sys
     
plt.rcParams.update({'font.size': 20})

# %%
data_dft = np.loadtxt(*PATH_BANDS_DFT)

k = np.unique(data_dft[:,0])
kbands = np.reshape(data_dft[:, 1], (-1,len(k)))
k = k/max(k)

for band in range(len(kbands)):
    plt.plot(k, kbands[band,:] - *FERMI ENERGY, '-',linewidth=5, alpha=0.5, color='red', markersize=2)
    # plt.plot(k, kbands[band,:], '-',linewidth=5, alpha=0.5, color='red', markersize=2)

# plot the paoflow bands
eband = np.load(*PATH_BANDS_PAOFLOW)
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
plt.show()   

# %%



