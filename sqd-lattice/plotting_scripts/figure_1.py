import os 
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])
from sqd_lattice.util import save_json,load_yaml,load_pickle,load_json,find_max_index

colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']
linestyles = ["solid","dotted","dashdot","dashed"]

# Loading material name
parser = argparse.ArgumentParser()
parser.add_argument('--material')
args = parser.parse_args()

f, ax = plt.subplots(1,2,figsize=(6,3))

material = "hafnium_2"

print(f"Plotting figure 1 for {material}",flush=True)

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
plot_config = load_yaml("plot_config.yaml")["figure_1"]
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
bandgap_dict = load_json(f"runs/{material}","bandgap_dict.json")

fontsize = int(plot_config["fontsize"])

folder_path = f"runs/{material}"

print(bandgap_dict)

if material=='zirconia':
    ci_label = "hci"
else:
    ci_label = "fci"

cases = ['TB', '$U$', '$V$', '$U+V$']

bar_width = 0.3
x = np.arange(len(cases))     # Create positions for all cases

method_A,method_B = [[bandgap_dict["DFT"],bandgap_dict["DFT+U"],bandgap_dict["DFT+V"],bandgap_dict["DFT+U+V"]],[bandgap_dict["tb"],bandgap_dict["u"],bandgap_dict["v"],bandgap_dict["hci"]]]

print(method_A)
print(method_B)

# ax[0].axhspan(bandgap_dict["Expt"][0], bandgap_dict["Expt"][1], color='grey', alpha=0.3)  # grey region between x_start and x_end    
# ax[0].axhline(y=bandgap_dict["Expt"][0], color='grey', linestyle='dashed', linewidth=1.2,zorder=1)
# ax[0].axhline(y=bandgap_dict["Expt"][1], color='grey', linestyle='dashed', linewidth=1.2,zorder=1)
# ax[0].axhline(y=bandgap_dict["GW"], color="orange", linestyle="dotted", linewidth=plot_config["linewidth"])
# ax[0].text(0.2,0.72, 'Experiment', color='black', ha='center', transform=ax[0].transAxes,fontsize=fontsize-1)
# ax[0].text(0.12,0.85, 'GW', color='orange', ha='center', transform=ax[0].transAxes,fontsize=fontsize-1)

# Plot methods A and B for first 4 cases
ax[0].bar(x - bar_width/2, method_A, width=bar_width, label='DFT', color=colors[0],zorder=2,edgecolor='black'
)

ax[0].bar(x + bar_width/2, method_B, width=bar_width, label='Lattice', color=colors[1],zorder=2,edgecolor='black')
ax[0].tick_params(axis='both',labelsize=fontsize-2,labelrotation=90.)

# Customize plot
ax[0].set_xticks(x, cases)
# plt.xlabel('Cases')
ax[0].set_ylabel(r'Band gap at $\Gamma$ (eV)',fontsize=fontsize-2)

# ax[0].set_ylim(0,7.5)
ax[0].set_ylim(0,6)

material = "zirconia_2"

print(f"Plotting figure 1 for {material}",flush=True)

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
plot_config = load_yaml("plot_config.yaml")["figure_1"]
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
bandgap_dict = load_json(f"runs/{material}","bandgap_dict.json")

fontsize = int(plot_config["fontsize"])

folder_path = f"runs/{material}"

print(bandgap_dict)

if material=='zirconia':
    ci_label = "hci"
else:
    ci_label = "fci"

cases = ['TB', '$U$', '$V$', '$U+V$']

bar_width = 0.3
x = np.arange(len(cases))     # Create positions for all cases

method_A,method_B = [[bandgap_dict["DFT"],bandgap_dict["DFT+U"],bandgap_dict["DFT+V"],bandgap_dict["DFT+U+V"]],[bandgap_dict["tb"],bandgap_dict["u"],bandgap_dict["v"],bandgap_dict["hci"]]]

print(method_A)
print(method_B)

# ax[1].axhspan(bandgap_dict["Expt"][0], bandgap_dict["Expt"][1], color='grey', alpha=0.3)  # grey region between x_start and x_end    
# ax[1].axhline(y=bandgap_dict["Expt"][0], color='grey', linestyle='dashed', linewidth=1.2,zorder=1)
# ax[1].axhline(y=bandgap_dict["Expt"][1], color='grey', linestyle='dashed', linewidth=1.2,zorder=1)
# ax[1].axhline(y=bandgap_dict["GW"], color="orange", linestyle="dotted", linewidth=plot_config["linewidth"])
# ax[1].text(0.2,0.86, 'Experiment', color='black', ha='center', transform=ax[1].transAxes,fontsize=fontsize-1)
# ax[1].text(0.12,0.75, 'GW', color='orange', ha='center', transform=ax[1].transAxes,fontsize=fontsize-1)

# Plot methods A and B for first 4 cases
ax[1].bar(x - bar_width/2, method_A, width=bar_width, color=colors[0],zorder=2,edgecolor='black')

ax[1].bar(x + bar_width/2, method_B, width=bar_width, color=colors[1],zorder=2,edgecolor='black')
ax[1].tick_params(axis='both',labelsize=fontsize-2,labelrotation=90.)

# Customize plot
ax[1].set_xticks(x, cases)
# plt.xlabel('Cases')
ax[1].set_ylabel(r'Band gap at $\Gamma$ (eV)',fontsize=fontsize-2)

ax[1].set_ylim(0,7.5)

f.legend(ncol=2,fontsize=fontsize-2,bbox_to_anchor=(0.55,-0.06),loc='lower center')
# Make folder

figure_folder = os.path.join("plots",f"{material}","figure_1")

os.makedirs(figure_folder,exist_ok=True)
plt.tight_layout()
# plt.savefig(os.path.join(figure_folder,f"figure_1.pdf"),format='pdf', bbox_inches="tight", pad_inches=0.01)
plt.savefig(os.path.join(figure_folder,f"figure_1.pdf"),format='pdf')
plt.savefig(os.path.join(figure_folder,f"figure_1.png"))
plt.close()

print("Finished!")

