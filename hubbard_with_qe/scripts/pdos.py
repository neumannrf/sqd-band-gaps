import re
import pandas as pd
import sys

##########################################
#Functions


def read_proj_outfile(name):
    with open(name, 'r') as file:
        content = file.readlines()
    return content

def get_parameters(content):
    parameters = {}

    for i_line, line in enumerate(content):
        if 'Problem Sizes' in line:
            for i in content[i_line+1:i_line+5]:
                parameters[i.split()[0]] = i.split()[-1]
            break
    return parameters

def get_info_states(content, parameters):

    header = [['state', 'atom', 'wfc', 'l', 'm']]

    for i_line, line in enumerate(content):
        if 'Atomic states used for projection' in line:
            for i in content[i_line+3:i_line+3+int(parameters['natomwfc'])]:
                info_si = [int(re.sub(r'\D', '',i.split()[j])) for j in [2,4,8,9]]
                header.append(info_si)
            break
    return header


def get_fermi_energy(data_scf):
    for line in data_scf:
        if 'highest' in line:
            return(float(line.split()[-2]))

def get_kpoints(content):

    k_points = []

    for i, line in enumerate(content):
        if 'number of k points=' in line:
            nk = int(line.split('=')[-1])
            for ni in range(nk):
                k_points.append([ni+1, float(content[i+ni+2].split()[-1])])
            break
#            k_points.append([float(k) for k in i.split()[2:]])
    return k_points

def get_occupations_wfc(content,parameters):

    header = [['k', 'band', 'ao', 'occ', 'energy', 'wfc']]
    count_k = 1

    for i_line, line in enumerate(content):
            if 'k = ' in line:
                while count_k <= int(parameters['nkstot']):
                    if 'eV' in content[i_line]:
                        band = int(re.sub(r'\D', '',content[i_line].split()[2]))
                        energy = float(content[i_line].split()[4])

                    if 'psi =' in content[i_line]:
                        occ = []
                        ao = []
                        while '|psi|' not in content[i_line]:
                            for i in content[i_line].split('+'):
                                if i.strip():
                                    if '=' in i:
                                        i = i.split('=')[-1]
                                    occ.append(float(i.split('*')[0]))
                                    ao.append(int(re.sub(r'\D', '',i.split('*')[1])))
                            i_line += 1
                    if '|psi|' in content[i_line]:
                        wfc = float(content[i_line].split()[-1])
                        for i,j in enumerate(occ):
                            # print(count_k, band, ao[i], j, energy, wfc) 
                            header.append([count_k, band, ao[i], j, energy, wfc])
                        if band == int(parameters['nbnd']):
                            count_k += 1
                    i_line += 1
                    
    return header

################################################################
# MAIN
# Provide SCF and PROJWFC.X output files, respectively.

################################################################
# Read files and build dataframe

content = read_proj_outfile(sys.argv[2])
data_scf = read_proj_outfile(sys.argv[1])
k_weights = get_kpoints(data_scf)
parameters = get_parameters(content)
full_data = get_occupations_wfc(content, parameters)
fermi = get_fermi_energy(data_scf)

df_weights = pd.DataFrame(k_weights, columns=['k', 'k_weight'])
df_data = pd.DataFrame(full_data[1:], columns=full_data[0])
df_full = pd.merge(df_data, df_weights, on='k')

# Save parsing in pdos_info.txt and pdos_info.pkl
df_full.to_pickle('pdos_info.pkl')
df_full.to_csv('pdos_info.txt', sep=' ', index=False)
################################################################

# Compute occupation for every k point
occ_ao={}

df_occupied = df_full[df_full['energy'] <= fermi]

for ao in range(1,int(parameters['natomwfc'])+1):
    bkup_ao = df_occupied[df_occupied['ao'] == ao] 
    sum_ao = 0

    for index, line in bkup_ao.iterrows():
        sum_ao += line['occ']*line['k_weight']

    occ_ao[ao] = sum_ao

occupation_all = pd.DataFrame(list(occ_ao.items()), columns=['orbital', 'occupation'])

# Save occupation for every k point in pdos_all

occupation_all.to_pickle('pdos_all.pkl')

total = occupation_all['occupation'].sum()

with open('pdos_all.txt', 'w') as f:
    f.write(occupation_all.to_string(index=False))  
    f.write(f"\nTotal sum: {total}\n")

################################################################

# Compute occupation for every Gamma point
occ_gamma={}

df_occupied = df_full[df_full['energy'] <= fermi]
df_occupied_gamma = df_occupied[df_occupied['k']==1]

for ao in range(1,int(parameters['natomwfc'])+1):
    bkup_ao = df_occupied_gamma[df_occupied_gamma['ao'] == ao] 
    sum_ao = 0

    for index, line in bkup_ao.iterrows():
        sum_ao += line['occ']*line['k_weight']/(line['k_weight']/2)

    occ_gamma[ao] = sum_ao

occupation_gamma = pd.DataFrame(list(occ_gamma.items()), columns=['orbital', 'occupation'])

# Save occupation for every Gamma point in pdos_all

occupation_gamma.to_pickle('pdos_gamma.pkl')

total = occupation_gamma['occupation'].sum()

with open('pdos_gamma.txt', 'w') as f:
    f.write(occupation_gamma.to_string(index=False))  
    f.write(f"\nTotal sum: {total}\n")