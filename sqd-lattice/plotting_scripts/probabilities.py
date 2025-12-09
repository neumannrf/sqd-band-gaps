import os 
import argparse
import numpy as np
import matplotlib.pyplot as plt

import ffsim

from sqd_lattice.util import load_yaml,load_pickle,load_json,find_max_index

def format_sci_notation_latex(n):
    exponent = len(str(n)) - 1
    base = n / (10 ** exponent)
    base = round(base, 1)
    return rf"${base} \times 10^{{{exponent}}}$"

colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']

# Loading material name
parser = argparse.ArgumentParser()
parser.add_argument('--material')
args = parser.parse_args()

material = args.material

print(f"Plotting probabilities for {material}",flush=True)

# Loading data
base_config = load_yaml(f"runs/{material}","base_config.yaml")
plot_config = load_yaml("plot_config.yaml")["probabilities"]
data_dict = load_pickle(f"runs/{material}","data_dict.pkl")
figsize = tuple(plot_config["figsize"])
fontsize = int(plot_config["fontsize"])
linewidth = int(plot_config["linewidth"])

folder_path = f"runs/{material}"

ne_total = sum(data_dict["ne_ab"])

probabilities_sqd = []
ne_list = []
fci_dimensions = []

for ne in [ne_total-1,ne_total,ne_total+1]:
    folder = os.path.join(folder_path,f"{ne}e")
    results_pre_process = load_pickle(os.path.join(folder,"results_pre_process.pkl"))

    sampling_sqd = plot_config["folder"]
    print(f"Loaded SQD results for {sampling_sqd}",flush=True)
    sqd_results_folder = os.path.join(folder,f"results_{sampling_sqd}","probs")
    best_sqd_probs = np.load(os.path.join(sqd_results_folder,f"probabilities_{find_max_index(sqd_results_folder)}.npy"))
    probabilities_sqd.append(best_sqd_probs)

    fci_dimensions.append(ffsim.dim(data_dict["num_orbitals"],results_pre_process["ne_ab"]))
    ne_list.append(ne)

def plot_probabilities(prob_list,file_name):
    # # Create horizontal bar chart
    f,ax = plt.subplots(1,1,figsize=figsize)
    ax.loglog(prob_list[0],label=f"$N_e={ne_list[0]}$, "+r"$D_{\text{full}} \simeq$ "+f"{format_sci_notation_latex(int(fci_dimensions[0]))}",color=colors[0],linestyle='dashed',linewidth=linewidth)
    ax.loglog(prob_list[1],label=f"$N_e={ne_list[1]}$, "+r"$D_{\text{full}} \simeq$ "+f"{format_sci_notation_latex(int(fci_dimensions[1]))}",color=colors[1],linestyle='dashdot',linewidth=linewidth)
    ax.loglog(prob_list[2],label=f"$N_e={ne_list[2]}$, "+r"$D_{\text{full}} \simeq$ "+f"{format_sci_notation_latex(int(fci_dimensions[2]))}",color=colors[2],linestyle='dotted',linewidth=linewidth)

    ax.set_xlabel('$\mathbf{x}$',fontsize=fontsize)
    ax.set_ylabel('$|c_{\mathbf{x}}|^2$',fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize-2)

    f.legend(loc='lower center',fontsize=fontsize-2,bbox_to_anchor=(0.57,-0.01))    

    figure_folder = os.path.join("plots",f"{material}","probs",sampling_sqd)

    plt.tight_layout(rect=[0, 0.25, 1, 1])

    os.makedirs(figure_folder,exist_ok=True)
    plt.savefig(os.path.join(figure_folder,f"{file_name}.png"),format='png')
    plt.savefig(os.path.join(figure_folder,f"{file_name}.pdf"),format='pdf')

plot_probabilities(probabilities_sqd,"probs_sqd")

print("Finished!")

    


