import os 
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
import matplotlib.colors as mcolors
plt.style.use(['science'])

import ffsim

from sqd_lattice.util import load_yaml,load_pickle,load_json,find_max_index,extract_subspace_energy_pairs

# plt.style.use(['science', 'no-latex'])

colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']

# Loading material name
parser = argparse.ArgumentParser()
parser.add_argument('--material')
args = parser.parse_args()

material = args.material

print(f"Plotting subspace dimension for {material}",flush=True)

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
plot_config = load_yaml("plot_config.yaml")["subspace_dimension"]
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
figsize = (2,3)
fontsize = int(plot_config["fontsize"])
ms = plot_config["marker_size"]
linewidth = int(plot_config["linewidth"])

print(ms)

folder_path = f"runs/{material}"

ne_total = sum(data_dict["ne_ab"])
ne_list = [ne_total-1,ne_total,ne_total+1]

results_dict = {}
run_fci = False

for ne in ne_list:
    folder = os.path.join(folder_path,f"{ne}e")
    results_pre_process = load_pickle(os.path.join(folder,"results_pre_process.pkl"))
    results_dict[ne] = {}
    # for method in ['sqd','ext_sqd',"hci"]:

    sqd_results_folder = os.path.join(folder,f"results_hci","dicts")
        
    points = extract_subspace_energy_pairs(sqd_results_folder,key="subspace_dimension")

    print(points.shape)

    fci_dims = ffsim.dim(data_dict["num_orbitals"],results_pre_process["ne_ab"])
    points[:,0] = points[:,0]/fci_dims

    if results_pre_process["fci_energy"] is not None:
        run_fci = True
        points[:,1] = np.abs(points[:,1] - results_pre_process["fci_energy"])
        ylabel = r"$|E_{gs}(\text{HCI})-E_{gs}(\text{FCI})|$ (eV)"
    
    else:
        run_fci = False
        ylabel = "Energy (eV)"

    results_dict[ne]["hci"] = points

    
    # hci_results_folder = os.path.join(folder,"results_hci","dicts")
    # best_hci_result = load_json(os.path.join(hci_results_folder,f"results_dict_{find_max_index(hci_results_folder)}.json"))
    # hci_energy = best_hci_result["hci_energy"]
    # if results_pre_process["fci_energy"] is not None:
    #     hci_energy = np.abs(hci_energy - results_pre_process["fci_energy"])
    # results_dict[ne]["hci"] = np.array([[best_hci_result["subspace_dimension"]/fci_dims,hci_energy]])

figure_folder = os.path.join("plots",f"{material}","subspace_dimension","hci")
os.makedirs(figure_folder,exist_ok=True)

cmap = plt.cm.inferno   # choose a colormap
norm = plt.Normalize()

# idx = 0
# ax.loglog(results_dict[ne_list[idx]]['hci'][:,0],results_dict[ne_list[idx]]['hci'][:,1],marker='o',label=f"HCI $N_e-1$",linestyle="dotted",ms=ms)

# idx = 1
# ax.loglog(results_dict[ne_list[idx]]['hci'][:,0],results_dict[ne_list[idx]]['hci'][:,1],marker='s',label=f"HCI $N_e$",linestyle="dotted",ms=ms)

# idx = 2
# ax.loglog(results_dict[ne_list[idx]]['hci'][:,0],results_dict[ne_list[idx]]['hci'][:,1],marker='>',label=f"HCI $N_e+1$",linestyle="dotted",ms=ms)

# Ne - 1
cmap = plt.cm.inferno

# collect all z values to get global min/max for consistent color scaling
z_all = np.concatenate([results_dict[ne]['hci'][:,2] for ne in ne_list])

print("ZALL",z_all.min())

norm = mcolors.LogNorm(vmin=z_all.min(), vmax=z_all.max())

f, ax = plt.subplots(1, 3, figsize=(8.,2.3))

for idx, marker, label in zip(
    [0,1,2],
    ['o','o','o'],
    [r"$N_e-1$", r"$N_e$", r"HCI $N_e+1$"]
):
    x = results_dict[ne_list[idx]]['hci'][:,0]
    y = results_dict[ne_list[idx]]['hci'][:,1]
    z = results_dict[ne_list[idx]]['hci'][:,2]

    # scatter with log-normalized colormap
    sc = ax[idx].scatter(
        x, y,
        c=z, cmap=cmap, norm=norm,
        marker=marker,
        s=ms**2,
        edgecolors="k", linewidths=0.4
    )

    # log axes
    ax[idx].set_xscale('log')
    if run_fci == True:
        ax[idx].set_yscale('log')

    # labels

    ax[idx].set_ylabel(ylabel, fontsize=fontsize)

    ax[idx].set_xlabel("Fraction of FCI space", fontsize=fontsize)
    ax[idx].tick_params(axis='both', labelsize=fontsize)

    # colorbar on the side of each subplot
    cbar = f.colorbar(sc, ax=ax[idx])
    cbar.set_label("Variance", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

f.tight_layout()

# if material=='hafnium':
# plt.legend(ncol=3,bbox_to_anchor=(0.5,-0.8),loc='lower center',fontsize=fontsize-2)
plt.savefig(os.path.join(figure_folder,"subspace_dimension.png"),format='png')
plt.savefig(os.path.join(figure_folder,"subspace_dimension.pdf"),format='pdf')

print("Finished!")