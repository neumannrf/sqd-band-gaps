import os 
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])

import ffsim

from sqd_lattice.util import load_yaml,load_pickle,load_json,find_max_index,extract_subspace_energy_pairs

# plt.style.use(['science', 'no-latex'])

colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']

# Loading material name
# parser = argparse.ArgumentParser()
# parser.add_argument('--material')
# args = parser.parse_args()

processor_dict={"zirconia_2":"hardware","hafnium_2":"hardware","platinum_V0":"hardware"}

f,ax = plt.subplots(1,2,figsize=(6,3))

material = "hafnium_2"

print(f"Plotting subspace dimension for {material}",flush=True)

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
plot_config = load_yaml("plot_config.yaml")["subspace_dimension"]
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
figsize = tuple(plot_config["figsize"])
fontsize = int(plot_config["fontsize"])
ms = plot_config["marker_size"]
linewidth = int(plot_config["linewidth"])
sampling_sqd = processor_dict[material]
alpha=plot_config["alpha"]

print(ms)

folder_path = f"runs/{material}"

ne_total = sum(data_dict["ne_ab"])
ne_list = [ne_total-1,ne_total,ne_total+1]

results_dict = {}
run_fci = False

hf_energies = []
ccsd_energies = []
ccsdt_energies = []

for ne in ne_list:
    folder = os.path.join(folder_path,f"{ne}e")
    results_pre_process = load_pickle(os.path.join(folder,"results_pre_process.pkl"))
    
    results_dict[ne] = {}
    for method in ['sqd','ext_sqd']:
        sqd_results_folder = os.path.join(folder,f"results_{method}_{sampling_sqd}","dicts")

        points = extract_subspace_energy_pairs(sqd_results_folder,key="pyscf_subspace_dimension")

        fci_dims = ffsim.dim(data_dict["num_orbitals"],results_pre_process["ne_ab"])
        print(points)
        points[:,0] = points[:,0]/fci_dims

        hci_results_folder = os.path.join(folder,"results_hci","dicts")
        best_hci_result = load_json(os.path.join(hci_results_folder,f"results_dict_{find_max_index(hci_results_folder)}.json"))
        hci_energy = best_hci_result["hci_energy"]

        points[:,1] = np.abs(points[:,1] - hci_energy)

        results_dict[ne][method] = points
    
    hf_energies.append(np.abs(hci_energy-results_pre_process["hf_energy"]))
    ccsd_energies.append(np.abs(hci_energy-results_pre_process["ccsd_energy"]))
    ccsdt_energies.append(np.abs(hci_energy-results_pre_process["ccsdt_energy"]))

    
    # hci_results_folder = os.path.join(folder,"results_hci","dicts")
    # best_hci_result = load_json(os.path.join(hci_results_folder,f"results_dict_{find_max_index(hci_results_folder)}.json"))
    # hci_energy = best_hci_result["hci_energy"]
    # if results_pre_process["fci_energy"] is not None:
    #     hci_energy = np.abs(hci_energy - results_pre_process["fci_energy"])
    # results_dict[ne]["hci"] = np.array([[best_hci_result["subspace_dimension"]/fci_dims,hci_energy]])

print("CCSD Energies: ",ccsdt_energies)
print("HF Energies: ",hf_energies)

figure_folder = os.path.join("plots",f"{material}","subspace_dimension",sampling_sqd)
os.makedirs(figure_folder,exist_ok=True)

idx = 0
# ax[0].axhline(y=hf_energies[idx], linestyle="solid",color=colors[2], linewidth=plot_config["linewidth"])
ax[0].axhline(y=ccsdt_energies[idx], linestyle="dotted",color=colors[0], linewidth=plot_config["linewidth"])
ax[0].loglog(results_dict[ne_list[idx]]['sqd'][:,0],results_dict[ne_list[idx]]['sqd'][:,1],color=colors[0],marker='p',linestyle="solid",ms=ms,markerfacecolor='none',alpha=alpha)
ax[0].loglog(results_dict[ne_list[idx]]['ext_sqd'][:,0],results_dict[ne_list[idx]]['ext_sqd'][:,1],color=colors[0],marker='p',linestyle="dashed",ms=ms,markerfacecolor='none',alpha=alpha)

idx = 1
# ax[0].axhline(y=hf_energies[idx], linestyle="solid",color=colors[2], linewidth=plot_config["linewidth"])
ax[0].axhline(y=ccsdt_energies[idx], linestyle="dashed", color=colors[1],linewidth=plot_config["linewidth"])
ax[0].loglog(results_dict[ne_list[idx]]['sqd'][:,0],results_dict[ne_list[idx]]['sqd'][:,1],color=colors[1],marker='X',linestyle="solid",ms=ms,markerfacecolor='none',alpha=alpha)
ax[0].loglog(results_dict[ne_list[idx]]['ext_sqd'][:,0],results_dict[ne_list[idx]]['ext_sqd'][:,1],color=colors[1],marker='X',linestyle="dashed",ms=ms,markerfacecolor='none',alpha=alpha)

idx = 2
# ax[0].axhline(y=hf_energies[idx], linestyle="solid",color=colors[2], linewidth=plot_config["linewidth"])
ax[0].axhline(y=ccsdt_energies[idx]+0.000003, linestyle="dashdot", color=colors[2], linewidth=plot_config["linewidth"])
ax[0].loglog(results_dict[ne_list[idx]]['sqd'][:,0],results_dict[ne_list[idx]]['sqd'][:,1],color=colors[2],marker='^',linestyle="solid",ms=ms,markerfacecolor='none',alpha=alpha)
ax[0].loglog(results_dict[ne_list[idx]]['ext_sqd'][:,0],results_dict[ne_list[idx]]['ext_sqd'][:,1],color=colors[2],marker='^',linestyle="dashed",ms=ms,markerfacecolor='none',alpha=alpha)



ax[0].set_xlabel("Fraction of Hilbert space",fontsize=fontsize-2)
ax[0].set_ylabel(r"$|E_{gs}-E_{gs}(\text{HCI})|$ (eV)",fontsize=fontsize-2)

ax[0].tick_params(axis='both', labelsize=fontsize)

material = "zirconia_2"

print(f"Plotting subspace dimension for {material}",flush=True)

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
plot_config = load_yaml("plot_config.yaml")["subspace_dimension"]
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
figsize = tuple(plot_config["figsize"])
fontsize = int(plot_config["fontsize"])
ms = plot_config["marker_size"]
linewidth = int(plot_config["linewidth"])
sampling_sqd = processor_dict[material]
alpha=plot_config["alpha"]

print(ms)

folder_path = f"runs/{material}"

ne_total = sum(data_dict["ne_ab"])
ne_list = [ne_total-1,ne_total,ne_total+1]

results_dict = {}
run_fci = False

hf_energies = []
ccsd_energies = []
ccsdt_energies = []

for ne in ne_list:
    folder = os.path.join(folder_path,f"{ne}e")
    results_pre_process = load_pickle(os.path.join(folder,"results_pre_process.pkl"))
    results_dict[ne] = {}
    for method in ['sqd','ext_sqd']:
        sqd_results_folder = os.path.join(folder,f"results_{method}_{sampling_sqd}","dicts")

        points = extract_subspace_energy_pairs(sqd_results_folder,key="pyscf_subspace_dimension")

        fci_dims = ffsim.dim(data_dict["num_orbitals"],results_pre_process["ne_ab"])
        print(points)
        points[:,0] = points[:,0]/fci_dims

        hci_results_folder = os.path.join(folder,"results_hci","dicts")
        best_hci_result = load_json(os.path.join(hci_results_folder,f"results_dict_{find_max_index(hci_results_folder)}.json"))
        hci_energy = best_hci_result["hci_energy"]

        points[:,1] = np.abs(points[:,1] - hci_energy)

        results_dict[ne][method] = points

    hf_energies.append(np.abs(hci_energy-results_pre_process["hf_energy"]))
    ccsd_energies.append(np.abs(hci_energy-results_pre_process["ccsd_energy"]))
    ccsdt_energies.append(np.abs(hci_energy-results_pre_process["ccsdt_energy"]))
    
    # hci_results_folder = os.path.join(folder,"results_hci","dicts")
    # best_hci_result = load_json(os.path.join(hci_results_folder,f"results_dict_{find_max_index(hci_results_folder)}.json"))
    # hci_energy = best_hci_result["hci_energy"]
    # if results_pre_process["fci_energy"] is not None:
    #     hci_energy = np.abs(hci_energy - results_pre_process["fci_energy"])
    # results_dict[ne]["hci"] = np.array([[best_hci_result["subspace_dimension"]/fci_dims,hci_energy]])

print("CCSD Energies: ",ccsdt_energies)
print("HF Energies: ",hf_energies)

figure_folder = os.path.join("plots",f"{material}","subspace_dimension",sampling_sqd)
os.makedirs(figure_folder,exist_ok=True)

idx = 0
# ax[1].axhline(y=hf_energies[idx], linestyle="solid",color=colors[0], linewidth=plot_config["linewidth"])
# ax[1].axhline(y=ccsdt_energies[idx], linestyle="dotted",color=colors[0],label=r"CCSD(T) $N_e-1$", linewidth=plot_config["linewidth"])
ax[1].axhline(y=ccsd_energies[idx], linestyle="dotted",color=colors[0],label=r"CCSD $N_e-1$", linewidth=plot_config["linewidth"])
ax[1].loglog(results_dict[ne_list[idx]]['sqd'][:,0],results_dict[ne_list[idx]]['sqd'][:,1],color=colors[0],label=r"QSCI $N_e-1$",marker='p',linestyle="solid",ms=ms,markerfacecolor='none',alpha=alpha)
ax[1].loglog(results_dict[ne_list[idx]]['ext_sqd'][:,0],results_dict[ne_list[idx]]['ext_sqd'][:,1],color=colors[0],marker='p',label=r"Ext-QSCI $N_e-1$",linestyle="dashed",ms=ms,markerfacecolor='none',alpha=alpha)

idx = 1
# ax[1].axhline(y=hf_energies[idx], linestyle="solid",color=colors[1], linewidth=plot_config["linewidth"])
# ax[1].axhline(y=ccsdt_energies[idx], linestyle="dashed", color=colors[1],label=r"CCSD(T) $N_e$",linewidth=plot_config["linewidth"])
ax[1].axhline(y=ccsd_energies[idx], linestyle="dashed", color=colors[1],label=r"CCSD $N_e$",linewidth=plot_config["linewidth"])
ax[1].loglog(results_dict[ne_list[idx]]['sqd'][:,0],results_dict[ne_list[idx]]['sqd'][:,1],color=colors[1],label=r"QSCI $N_e$",marker='X',linestyle="solid",ms=ms,markerfacecolor='none',alpha=alpha)
ax[1].loglog(results_dict[ne_list[idx]]['ext_sqd'][:,0],results_dict[ne_list[idx]]['ext_sqd'][:,1],color=colors[1],label=r"Ext-QSCI $N_e$",marker='X',linestyle="dashed",ms=ms,markerfacecolor='none',alpha=alpha)

idx = 2
# ax[1].axhline(y=hf_energies[idx], linestyle="solid",color=colors[2], linewidth=plot_config["linewidth"])
# ax[1].axhline(y=ccsdt_energies[idx]+0.00000115, linestyle="dashdot", color=colors[2],label=r"CCSD(T) $N_e+1$", linewidth=plot_config["linewidth"])
ax[1].axhline(y=ccsd_energies[idx], linestyle="dashdot", color=colors[2],label=r"CCSD $N_e+1$", linewidth=plot_config["linewidth"])
ax[1].loglog(results_dict[ne_list[idx]]['sqd'][:,0],results_dict[ne_list[idx]]['sqd'][:,1],color=colors[2],label=r"QSCI $N_e+1$",marker='^',linestyle="solid",ms=ms,markerfacecolor='none',alpha=alpha)
ax[1].loglog(results_dict[ne_list[idx]]['ext_sqd'][:,0],results_dict[ne_list[idx]]['ext_sqd'][:,1],color=colors[2],label=r"Ext-QSCI $N_e+1$",marker='^',linestyle="dashed",ms=ms,markerfacecolor='none',alpha=alpha)



ax[1].set_xlabel("Fraction of Hilbert space",fontsize=fontsize-2)
ax[1].set_ylabel(r"$|E_{gs}-E_{gs}(\text{HCI})|$ (eV)",fontsize=fontsize-2)

ax[1].tick_params(axis='both', labelsize=fontsize)

f.legend(ncol=3,bbox_to_anchor=(0.57,-0.15),loc='lower center',fontsize=fontsize-2)

plt.tight_layout()

plt.savefig(os.path.join(figure_folder,"subspace_dimension.png"),format='png')
plt.savefig(os.path.join(figure_folder,"subspace_dimension.pdf"),format='pdf', bbox_inches="tight", pad_inches=0.01)

print("Finished!")