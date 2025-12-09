# %%
import warnings
warnings.filterwarnings("ignore") 
import numpy as np
from PAOFLOW import PAOFLOW
import sys
import itertools

material = sys.argv[1]

paoflow = PAOFLOW.PAOFLOW(savedir=f'out/{material}.save', model=None, outputdir='.', 
                          smearing='gauss', acbn0=False, verbose=False)
paoflow.projections()
paoflow.projectability(pthr=0.95)
paoflow.pao_hamiltonian(write_binary=True)
data_controller = paoflow.data_controller
arry,attr = paoflow.data_controller.data_dicts()
basis = arry['basis']
coordinates = arry['tau']


def read_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    lines = [line.strip('\n').split() for line in lines]
    return lines

def extract_U_block(full_v):

    control_flag = False 

    for index, i in enumerate(full_v):
        if 'site' in i:
            ini_ind = index
            control_flag = True

        if not i and control_flag:
            final_ind = index   
            break

    format_block = [[int(i[0]), int(i[1]), i[2], int(i[3]), int(i[4]), i[5], i[6], float(i[7])] for i in full_v[ini_ind+1:final_ind]]
    
    return format_block

def extract_V_block(full_v):
    control_flag = False

    for index, i in enumerate(full_v):
        if "Atom" in full_v[index-2]:
            ini_ind = index
            control_flag = True
        if i and '=' in i[0] and control_flag:
            final_ind = index   
            break

    filtered_list = [sublist for sublist in full_v[ini_ind:final_ind] if sublist] 

    format_block = [[int(i[0]), i[1], int(i[2]), i[3], float(i[4]), float(i[5])] for i in filtered_list]
    
    return format_block

def map_index(vblock, n):

    for index,i in enumerate(vblock):
        if i[2] <= n:
            continue
        elif i[2] % n == 0:
            new_number = n
            vblock[index][2] = new_number
        else:
            new_number = i[2] % n
            vblock[index][2] = new_number

    return vblock

def compute_paoflox_index(basis, coordinates, index_qe, orbital):


    indexes_paoflow = []

    for index_paoflow,atom_orb in enumerate(basis):
        if np.allclose(atom_orb['tau'],coordinates[index_qe - 1]) and np.char.upper(atom_orb['label']) == np.char.upper(orbital):
            indexes_paoflow.append(index_paoflow)


    return indexes_paoflow



full_v=read_file(f'../hubbard_computation/{material}.Hubbard_parameters.dat')

n = attr['natoms']

vblock = extract_V_block(full_v)

ublock = extract_U_block(full_v)

vblock_mapped = map_index(vblock, n)

combinations = list(itertools.combinations(range(1, n+1), 2))

combinations_full = [(i,i) for i in range(1, n+1)]

for i in combinations:
    combinations_full.append(i)

interactions_v = []

for i in vblock_mapped:
    
    if (i[0], i[2]) in combinations:
        interactions_v.append(i)
        combinations.remove((i[0], i[2]))

hubb_dict = {}

for i in combinations_full:
    hubb_dict[i] = {'species_i': None, 'orbital_i': None, 
                    'species_j': None, 'orbital_j': None, 
                    'V': None, 'distance': None, 
                    'index_i': None, 'index_j': None}

for i in ublock:
    key = (i[0], i[0])
    if key in hubb_dict:
        hubb_dict[key]['species_i'] = i[2]
        hubb_dict[key]['orbital_i'] = i[6]
        hubb_dict[key]['species_j'] = i[2]
        hubb_dict[key]['orbital_j'] = i[6]
        hubb_dict[key]['V'] = i[7]
        hubb_dict[key]['distance'] = 0
        hubb_dict[key]['index_i'] = compute_paoflox_index(basis, coordinates, i[0], i[6])
        hubb_dict[key]['index_j'] = compute_paoflox_index(basis, coordinates, i[0], i[6])

for i in interactions_v:
    key = (i[0], i[2])

    for j in hubb_dict.keys():
        if (i[0],i[0]) == j:
            orb_i = hubb_dict[j]['orbital_i']
        elif (i[2],i[2]) == j:
            orb_j = hubb_dict[j]['orbital_j']

    if key in hubb_dict:
        hubb_dict[key]['species_i'] = i[1]
        hubb_dict[key]['orbital_i'] = orb_i
        hubb_dict[key]['species_j'] = i[3]
        hubb_dict[key]['orbital_j'] = orb_j
        hubb_dict[key]['V'] = i[5]
        hubb_dict[key]['distance'] = i[4]
        hubb_dict[key]['index_i'] = compute_paoflox_index(basis, coordinates, i[0], orb_i)
        hubb_dict[key]['index_j'] = compute_paoflox_index(basis, coordinates, i[2], orb_j)

np.save(f'hubbard_info_{material}.npy', hubb_dict)

# %%



