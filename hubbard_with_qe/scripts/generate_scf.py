# %%
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import sys 
import numpy as np
import pymatgen.core as pm
from pymatgen.io.pwscf import PWInput
import json

# %%
def write_file(filename, content):
    with open(filename, "w") as f:
        f.write(content)

def find_lattice(a, b, c, alpha, beta, gamma):
    """
    Identify the Bravais lattice from lattice parameters.
    """
    thr = 1e-4
    reticolo = ""

    if abs(alpha - 90.0) < thr and abs(gamma - 90.0) < thr:
        if abs(beta - 90.0) < thr:
            if abs(a - b) < thr and abs(a - c) < thr:
                reticolo = "cubic"
            elif abs(a - b) < thr:
                reticolo = "tetragonal"
            else:
                reticolo = "orthorhombic"
        else:
            reticolo = "monoclinic"
    elif abs(alpha - 90.0) < thr and abs(beta - 90.0) < thr and abs(gamma - 120.0) < thr:
        reticolo = "hexagonal"
    elif abs(alpha - beta) < thr and abs(alpha - gamma) < thr and abs(a - b) < thr and abs(a - c) < thr:
        reticolo = "rhombohedral"
    else:
        reticolo = "triclinic"

    return reticolo

def find_ibrav(spacegroup_symbol, reticolo):
    """
    Map Bravais lattice and space group centering to QE's ibrav value.
    """
    primitive = spacegroup_symbol.startswith("P")
    bodycentered = spacegroup_symbol.startswith("I")
    facecentered = spacegroup_symbol.startswith("F")
    basecentered = spacegroup_symbol.startswith("C")

    match reticolo:
        case "cubic":
            if primitive:
                ibrav = 1
            if facecentered:
                ibrav = 2
            if bodycentered:
                ibrav = 3
        case "tetragonal":
            if primitive:
                ibrav = 6
            if bodycentered:
                ibrav = 7
        case "orthorhombic":
            if primitive:
                ibrav = 8
            if basecentered:
                ibrav = 9
            if facecentered:
                ibrav = 10
            if bodycentered:
                ibrav = 11
        case "monoclinic":
            if primitive:
                ibrav = 12
            if basecentered:
                ibrav = 13
        case "triclinic":
            ibrav = 14
        case "hexagonal":
            ibrav = 4
        case "rhombohedral":
            if primitive:
                ibrav = 4
            else:
                ibrav = 5
        case _:
            ibrav = 0

    return ibrav

def priority(element):
    ''''
    Objective: attributes a value to the element according to the classification:
    has_d_electrons, nearest-neighbors of d_element, others

    Input: element (string)
    Output: value 2, 1, 0 (integer)
    '''

    if pm.Element(element).is_transition_metal:
        value = 1
    else:
        value = 0

    return value


# %%
material=sys.argv[1]
# Load structure from a file (CIF, POSCAR, etc.)
structure = Structure.from_file(material+'.cif')
primitive_structure = structure.get_primitive_structure()

# %%
# Analyze space group and lattice type
spa = SpacegroupAnalyzer(primitive_structure)

# Get the space group number and lattice system
spacegroup_number = spa.get_space_group_number()
spacegroup_symbol = spa.get_space_group_symbol()  # Corrected!
lattice_type = spa.get_lattice_type()
crystal_system = spa.get_crystal_system()

# Get conventional cell (important for QE)
conventional_structure = spa.get_conventional_standard_structure()

# Get lattice parameters
lattice = conventional_structure.lattice
a, b, c = lattice.a, lattice.b, lattice.c
alpha, beta, gamma = lattice.alpha, lattice.beta, lattice.gamma

# Convert lengths to Bohr (QE uses Bohr)
bohr_conversion = 1.88973  # 1 Å = 1.88973 Bohr
a_bohr, b_bohr, c_bohr = a * bohr_conversion, b * bohr_conversion, c * bohr_conversion

cell = find_lattice(a, b, c, alpha, beta, gamma)
ibrav = find_ibrav(spacegroup_symbol, cell)

# Compute celldm parameters
celldm1 = a_bohr
celldm2 = b / a
celldm3 = c / a
celldm4 = np.cos(alpha)
celldm5 = np.cos(beta)
celldm6 = np.cos(gamma)

celldms = [celldm1, celldm2, celldm3, celldm4, celldm5, celldm6]

# %%


with open(f'structure_info_{material}.dat', 'w') as f:
    # Print results
    print(f"Space Group Number: {spacegroup_number}", file=f)
    print(f"Space Group Symbol: {spacegroup_symbol}", file=f)
    print(f"Lattice Type: {lattice_type}", file=f)
    print(f"Crystal System: {crystal_system}", file=f)
    print(f"Lattice Parameters: a={a:.4f}, b={b:.4f}, c={c:.4f}", file=f)
    print(f"Angles: α={alpha:.2f}, β={beta:.2f}, γ={gamma:.2f}", file=f)

    # Print results
    print(f"Detected Space Group: {spacegroup_symbol} ({spacegroup_number})", file=f)
    print(f"Crystal System: {crystal_system}", file=f)
    print(f"Lattice Type: {lattice_type}", file=f)
    print(f"Suggested ibrav: {ibrav}", file=f)

    print("\nQE celldm parameters:", file=f)
    print(f"celldm(1) = {celldm1:.6f}  ! a (in Bohr)", file=f)
    if ibrav in [4, 6, 7, 8, 9, 10, 11, 12, 13, 14]: print(f"celldm(2) = {celldm2:.6f}  ! b/a", file=f)
    if ibrav in [4, 6, 7, 8, 9, 10, 11, 12, 13, 14]: print(f"celldm(3) = {celldm3:.6f}  ! c/a", file=f)
    if ibrav in [12, 13, 14]: print(f"celldm(4) = {celldm4:.6f}  ! cos(α)", file=f)
    if ibrav in [14]: print(f"celldm(5) = {celldm5:.6f}  ! cos(β)", file=f)
    if ibrav in [14]: print(f"celldm(6) = {celldm6:.6f}  ! cos(γ)\n", file=f)


    print("Atoms in the primitive cell: \n", file=f)
    for i in primitive_structure:
        print(i.species, i.frac_coords, file=f)

# symmetrized_structure = spa.get_symmetrized_structure()
# for i, group in enumerate(symmetrized_structure.equivalent_sites):
#     print(f"Inequivalent site {i+1}: {group[0].species}, Fractional Coordinates: {group[0].frac_coords}, Real Coordinates: {group[0].coords}")

# coords = [symmetrized_structure.equivalent_sites[0][0].frac_coords, symmetrized_structure.equivalent_sites[1][0].frac_coords, symmetrized_structure.equivalent_sites[1][1].frac_coords] 
# species = [symmetrized_structure.equivalent_sites[0][0].species, symmetrized_structure.equivalent_sites[1][0].species, symmetrized_structure.equivalent_sites[1][1].species] 


# %%
def make_input(material, structure, ibrav, celldms, pseudo):
    '''
    Objective: make input file with PWInput

    Input: 
    Output: 
    '''

    name_output = material + '.scf.in'

    system={"ibrav":ibrav,
        "ecutwfc": '#ECUTWFC',
        "ecutrho": '#ECUTRHO',
        "occupations":"fixed",
        "nspin":1,
        "celldm(1)": celldms[0],
  	"!nbnd": 'nbnd' }

    if ibrav in [4, 6, 7, 8, 9, 10, 11, 12, 13, 14]: system["celldm(2)"] = celldms[1]
    if ibrav in [4, 6, 7, 8, 9, 10, 11, 12, 13, 14]: system["celldm(3)"] = celldms[2]
    if ibrav in [12, 13, 14]: system["celldm(4)"] = celldms[3]
    if ibrav in [14]: system["celldm(5)"] = celldms[4]
    if ibrav in [14]: system["celldm(6)"] = celldms[5]

    primitive_structure = structure.get_primitive_structure()
    species = []
    coords = []
    for i in primitive_structure.sites:
        species.append(i.specie.name)
        coords.append(i.frac_coords)

    lattice = structure.lattice

    reduced_structure = Structure(lattice, species, coords)

    reduced_structure.sort(key= lambda x: priority(x.specie.name), reverse=True)

    # pseudo = {}

    # for element in reduced_structure.composition.elements:
    #     pseudo[element.name] = f"PSEUDO_{element.name}"
    

    pwi = PWInput(reduced_structure, pseudo=pseudo ,  
            control={"calculation":"scf", 
                        "restart_mode":"from_scratch",
                        "pseudo_dir":'#PATH',
                        "outdir":"./out/",
                        "prefix":material},
            system=system,
            electrons={"diagonalization":"david",
                        "conv_thr":1.0e-8,
                        "mixing_beta":0.1,
                        "electron_maxstep":200},
            ions={"ion_dynamics": "bfgs"},
            kpoints_grid="#K#K#K",
            kpoints_shift=[0, 0, 0]                     
            )
    pwi.write_file(name_output)

    # write the Hubbard card

    if int(sys.argv[3]) == 1:
     
        # initialize input
        input_hubbard = "\nHUBBARD {ortho-atomic}\n"

        for element in reduced_structure.composition.elements:
                if element.is_transition_metal:
                    
                    metal = element

                    #finds the index of the first atom from element species
                    for index, atom in enumerate(reduced_structure):
                        if element == atom.specie:
                            index_metal = index+1
                            break

                    #write the U line for the d-element
                    for orbital in element.electronic_structure.split('.'):
                        if 'd'in orbital:
                            orb_element = orbital[:2]
                            input_hubbard = input_hubbard + f"V {element.name}-{orb_element} {element.name}-{orb_element} {index_metal} {index_metal} 1.0d-8\n"

        for neighbor in reduced_structure.composition.elements:
                
                if neighbor is not metal:

                    valence_nn = pm.Element(neighbor).electronic_structure.split('.')

                    #finds the index of the first atom from this element species
                    for index, atom in enumerate(reduced_structure):
                        if str(neighbor) == str(atom.specie.name):
                            index_nn = index+1
                            break

                    orbital_nn = valence_nn[-1]
                    input_hubbard = input_hubbard + f"V {metal.name}-{orb_element} {neighbor}-{orbital_nn[:2]} {index_metal} {index_nn} 1.0d-8\n"

        #append the Hubbard card in file
        with open(name_output, 'a') as f:
            f.write(input_hubbard)
        f.close()

    with open(name_output, 'r') as f:
            data = f.readlines()

    for index, line in enumerate(data):
            if "CELL_PARAMETERS" in line:
                    save_index = index
                    break
    del data[save_index:save_index+4]

    with open(name_output, 'w') as f:
        for i in data:
            f.write(i)
    f.close()


pseudo_string = sys.argv[2]

pseudo = json.loads(pseudo_string)


# %%
make_input(material, structure, ibrav, celldms, pseudo)

# %%



