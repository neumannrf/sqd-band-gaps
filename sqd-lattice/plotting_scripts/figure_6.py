import os 
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])
from sqd_lattice.util import save_json,load_yaml,load_pickle,load_json,find_max_index

processor_dict={"zirconia_2":"hardware","hafnium_2":"hardware","platinum_V0":"hardware"}
ylims={"zirconia_2":0,"hafnium_2":2.5,"platinum_V0":1.7}
ylims2={"zirconia_2":2,"hafnium_2":2,"platinum_V0":0.5}
colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']
linestyles = 2*["solid","dotted","dashdot","dashed"]

def compute_band_gap(energy_list):
    return energy_list[0]+energy_list[2]-2*(energy_list[1])

# Loading material name
# parser = argparse.ArgumentParser()
# parser.add_argument('--material')
# args = parser.parse_args()
# Create horizontal bar chart
f,ax = plt.subplots(1,2,figsize=(6, 3))

material = "hafnium_2"

print(f"Plotting band gap comparison for {material}",flush=True)

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
plot_config = load_yaml("plot_config.yaml")["figure_6"]
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
bandgap_dict = load_json(f"runs/{material}","bandgap_dict.json")
figsize = tuple(plot_config["figsize"])
fontsize = int(plot_config["fontsize"])

folder_path = f"runs/{material}"

ne_total = sum(data_dict["ne_ab"])

hf_energies = []
ccsd_energies = []
ccsdt_energies = []
hci_energies = []
fci_energies = []
sqd_energies = []
ext_sqd_energies = []

for ne in [ne_total-1,ne_total,ne_total+1]:
    folder = os.path.join(folder_path,f"{ne}e")
    results_pre_process = load_pickle(os.path.join(folder,"results_pre_process.pkl"))

    hf_energies.append(results_pre_process["hf_energy"])
    ccsd_energies.append(results_pre_process["ccsd_energy"])
    ccsdt_energies.append(results_pre_process["ccsdt_energy"])

    hci_results_folder = os.path.join(folder,"results_hci","dicts")
    best_hci_result = load_json(os.path.join(hci_results_folder,f"results_dict_{find_max_index(hci_results_folder)}.json"))
    hci_energies.append(best_hci_result["hci_energy"])

    sampling_sqd = "sqd_"+processor_dict[material]
    print(f"Loaded SQD results for {sampling_sqd} sampling")
    sqd_results_folder = os.path.join(folder,f"results_{sampling_sqd}","dicts")
    best_sqd_result = load_json(os.path.join(sqd_results_folder,f"results_dict_{find_max_index(sqd_results_folder)}.json"))
    sqd_energies.append(best_sqd_result["sqd_energy"])
    print(sqd_energies)

    sampling_sqd = "ext_sqd_"+processor_dict[material]
    print(f"Loaded SQD results for {sampling_sqd} sampling")
    sqd_results_folder = os.path.join(folder,f"results_{sampling_sqd}","dicts")
    best_sqd_result = load_json(os.path.join(sqd_results_folder,f"results_dict_{find_max_index(sqd_results_folder)}.json"))
    ext_sqd_energies.append(best_sqd_result["sqd_energy"])

    if results_pre_process["fci_energy"] is not None:
        fci_energies.append(results_pre_process["fci_energy"])


bandgap_hf = compute_band_gap(hf_energies)
bandgap_ccsd = compute_band_gap(ccsd_energies)
bandgap_ccsdt = compute_band_gap(ccsdt_energies)
bandgap_hci = compute_band_gap(hci_energies)

bandgap_dict_all = bandgap_dict.copy()

bandgap_dict_all["hf"] = bandgap_hf
bandgap_dict_all["ccsd"] = bandgap_ccsd
bandgap_dict_all["ccsdt"] = bandgap_ccsdt
bandgap_dict_all["hci"] = bandgap_hci

if len(fci_energies)==3:
    bandgap_ci = compute_band_gap(fci_energies)
    bandgap_dict_all["fci"] = bandgap_ci
    
    ci_label = 'FCI'
else:
    bandgap_ci = compute_band_gap(hci_energies)
    ci_label = 'HCI'

bandgap_sqd = compute_band_gap(sqd_energies)
bandgap_ext_sqd = compute_band_gap(ext_sqd_energies)

print("Error in bandgap: ", np.abs(bandgap_ext_sqd-bandgap_hci))

bandgap_dict_all["sqd"] = bandgap_sqd
bandgap_dict_all["ext_sqd"] = bandgap_ext_sqd

# save_json(bandgap_dict_all,f"runs/{material}","bandgap_dict.json")

# Example data
labels = ['CCSD(T)', "HCI",'SQD',"Ext-SQD"] 
values = [bandgap_ccsdt, bandgap_ci, bandgap_sqd,bandgap_ext_sqd]

ax[0].axhspan(bandgap_dict["Expt"][0], bandgap_dict["Expt"][1], color='grey', alpha=0.3)  # grey region between x_start and x_end    
ax[0].axhline(y=bandgap_dict["Expt"][0], color='grey',label="Experimental", linestyle='dashed', linewidth=1.2,zorder=1)
ax[0].axhline(y=bandgap_dict["Expt"][1], color='grey', linestyle='dashed', linewidth=1.2,zorder=1)
# ax.text(0.28,0.86, 'Experiment', color='black', ha='center', transform=ax.transAxes,fontsize=fontsize-2)
# ax.text(1.1,0.8, 'GW', color='orange', ha='center', transform=ax.transAxes,fontsize=fontsize-2)
ax[0].axhline(y=bandgap_dict["GW"], color="orange",label="GW", linestyle="dotted", linewidth=plot_config["linewidth"])


ax[0].bar(labels, values)


# Add labels and title
ax[0].set_ylabel(r'Band gap at $\Gamma$ (eV)',fontsize=fontsize-2)

ax[0].tick_params(axis='both',labelsize=fontsize-2,labelrotation=90.)

for tick in ax[0].xaxis.get_major_ticks():
    tick.tick1line.set_visible(False)
    tick.tick2line.set_visible(False)
for tick in ax[0].yaxis.get_major_ticks():
    tick.tick1line.set_visible(False)
    tick.tick2line.set_visible(False)

material = "zirconia_2"

print(f"Plotting band gap comparison for {material}",flush=True)

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
plot_config = load_yaml("plot_config.yaml")["figure_6"]
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
bandgap_dict = load_json(f"runs/{material}","bandgap_dict.json")
figsize = tuple(plot_config["figsize"])
fontsize = int(plot_config["fontsize"])

folder_path = f"runs/{material}"

ne_total = sum(data_dict["ne_ab"])

hf_energies = []
ccsd_energies = []
ccsdt_energies = []
hci_energies = []
fci_energies = []
sqd_energies = []
ext_sqd_energies = []

for ne in [ne_total-1,ne_total,ne_total+1]:
    folder = os.path.join(folder_path,f"{ne}e")
    results_pre_process = load_pickle(os.path.join(folder,"results_pre_process.pkl"))

    hf_energies.append(results_pre_process["hf_energy"])
    ccsd_energies.append(results_pre_process["ccsd_energy"])
    ccsdt_energies.append(results_pre_process["ccsdt_energy"])

    hci_results_folder = os.path.join(folder,"results_hci","dicts")
    best_hci_result = load_json(os.path.join(hci_results_folder,f"results_dict_{find_max_index(hci_results_folder)}.json"))
    hci_energies.append(best_hci_result["hci_energy"])

    sampling_sqd = "sqd_"+processor_dict[material]
    print(f"Loaded SQD results for {sampling_sqd} sampling")
    sqd_results_folder = os.path.join(folder,f"results_{sampling_sqd}","dicts")
    best_sqd_result = load_json(os.path.join(sqd_results_folder,f"results_dict_{find_max_index(sqd_results_folder)}.json"))
    sqd_energies.append(best_sqd_result["sqd_energy"])
    print(sqd_energies)

    sampling_sqd = "ext_sqd_"+processor_dict[material]
    print(f"Loaded SQD results for {sampling_sqd} sampling")
    sqd_results_folder = os.path.join(folder,f"results_{sampling_sqd}","dicts")
    best_sqd_result = load_json(os.path.join(sqd_results_folder,f"results_dict_{find_max_index(sqd_results_folder)}.json"))
    ext_sqd_energies.append(best_sqd_result["sqd_energy"])

    if results_pre_process["fci_energy"] is not None:
        fci_energies.append(results_pre_process["fci_energy"])


bandgap_hf = compute_band_gap(hf_energies)
bandgap_ccsd = compute_band_gap(ccsd_energies)
bandgap_ccsdt = compute_band_gap(ccsdt_energies)
bandgap_hci = compute_band_gap(hci_energies)

bandgap_dict_all = bandgap_dict.copy()

bandgap_dict_all["hf"] = bandgap_hf
bandgap_dict_all["ccsd"] = bandgap_ccsd
bandgap_dict_all["ccsdt"] = bandgap_ccsdt
bandgap_dict_all["hci"] = bandgap_hci

if len(fci_energies)==3:
    bandgap_ci = compute_band_gap(fci_energies)
    bandgap_dict_all["fci"] = bandgap_ci
    
    ci_label = 'FCI'
else:
    bandgap_ci = compute_band_gap(hci_energies)
    ci_label = 'HCI'

bandgap_sqd = compute_band_gap(sqd_energies)
bandgap_ext_sqd = compute_band_gap(ext_sqd_energies)


print("Error in bandgap: ", np.abs(bandgap_ext_sqd-bandgap_hci))

bandgap_dict_all["sqd"] = bandgap_sqd
bandgap_dict_all["ext_sqd"] = bandgap_ext_sqd

# save_json(bandgap_dict_all,f"runs/{material}","bandgap_dict.json")

# Example data
labels = ['CCSD(T)', "HCI",'SQD',"Ext-SQD"] 
values = [bandgap_ccsdt, bandgap_ci, bandgap_sqd,bandgap_ext_sqd]

ax[1].axhspan(bandgap_dict["Expt"][0], bandgap_dict["Expt"][1], color='grey', alpha=0.3)  # grey region between x_start and x_end    
ax[1].axhline(y=bandgap_dict["Expt"][0], color='grey', linestyle='dashed', linewidth=1.2,zorder=1)
ax[1].axhline(y=bandgap_dict["Expt"][1], color='grey', linestyle='dashed', linewidth=1.2,zorder=1)
# ax.text(0.28,0.86, 'Experiment', color='black', ha='center', transform=ax.transAxes,fontsize=fontsize-2)
# ax.text(1.1,0.8, 'GW', color='orange', ha='center', transform=ax.transAxes,fontsize=fontsize-2)
ax[1].axhline(y=bandgap_dict["GW"], color="orange", linestyle="dotted", linewidth=plot_config["linewidth"])


ax[1].bar(labels, values)


# Add labels and title
ax[1].set_ylabel(r'Band gap at $\Gamma$ (eV)',fontsize=fontsize-2)

ax[1].tick_params(axis='both',labelsize=fontsize-2,labelrotation=90.)

for tick in ax[1].xaxis.get_major_ticks():
    tick.tick1line.set_visible(False)
    tick.tick2line.set_visible(False)
for tick in ax[1].yaxis.get_major_ticks():
    tick.tick1line.set_visible(False)
    tick.tick2line.set_visible(False)

# Show the plot
f.legend(loc='lower center',bbox_to_anchor=((0.52,-0.075)),ncols=2)
f.tight_layout()

# Make folder
figure_folder = os.path.join("plots",f"{material}","figure_6")
os.makedirs(figure_folder,exist_ok=True)
plt.savefig(os.path.join(figure_folder,"figure_6.png"),format='png')
plt.savefig(os.path.join(figure_folder,"figure_6.pdf"),format='pdf')

print("Finished!",flush=True)






    


