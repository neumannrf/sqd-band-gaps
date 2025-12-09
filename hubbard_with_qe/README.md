# Hubbard Parameters with Quantum ESPRESSO

A set of scripts for calculating Hubbard parameters by running DFT+U+V calculations and computing their respective band gaps using Quantum ESPRESSO. This repo is targeted towards metal oxides of the type MOX, where M is the metal, O is the oxygen and X is the number of O atoms. 

## Overview

This repository contains scripts to:

1. Generate QE input files from CIF structures
2. Run DFT calculations and compute band gaps
3. Calculate Hubbard parameters 
4. Perform DFT+U+V, DFT+U and DFT+V calculations
5. Perform Tight-Binding (TB) projection using PAOFLOW
6. Parse results
7. Extract Hubbard information and copy files for QC workflow in the outputs_for_QC/ folder

## Prerequisites

- Quantum ESPRESSO
- Python with PAOFLOW package
- MPI environment

## Directory Structure

```
.
├── input                   # Main input configuration file
├── main.sh                 # Main execution script
├── scripts/
│   ├── generate_scf.py                  # Generate Quantum Espresso input for provided CIF file      
│   ├── run_dft.sh                       # DFT calculations and geometry optimization
│       ├── convergence_ecut.sh          # Convergence of energy cut-offs and k-grid
│   ├── run_hubbard_computation.sh       # Hubbard parameter calculation
│   ├── run_hubbard.sh                   # DFT+U, DFT+V and DFT+U+V calculations 
│   ├── run_paoflow.sh                   # PAOFLOW projections
│       ├── make_hopping.py              # Compute TB projection with PAOFLOW and save basis and matrix in npy file
│       ├── plot_hopping_bands.py        # Band structure and TB matrix plotting
│       ├── pdos.py                      # Compute atomic occupation per orbital at the Gamma point
│   ├── parse.sh                         # Parse final results to csv file
│   ├── extract_hubbard.sh               # Create npy file with hubbard information to be used un QC workflow       
│   └── submit_hp.sh                     # HP convergence calculations of the q-grid
├── pseudos/                             # Pseudopotential files
 

```

## Usage

1. Configure calculation parameters in the `input` file
2. Place pseudopotential files in `pseudos/` directory 
3. Place the CIF file of the desired material in the folder
4. Run the main script:

```bash
./main.sh
```

Alternatively, one can run the script in blocks by leaving only the desired block uncommented. Notice that this code allows for one to test several pseudopotential files for the metal. It is not adapted to do the same for oxygen. 

## Input File Structure

The `input` file contains sections for:

- Material name
- Pseudopotential paths and specifications
- k-points for band structure
- Parameters for bands calculation
- HP calculation parameters
- Projection parameters
- Number of processors

## Output

The calculations produce:
- Band structure data
- Hubbard parameters
- Projected DOS
- Final results in CSV format

## Installation

### Environment Setup

Intall the dependencies using the `hubbard-with-qe/Pipfile`. 

## To-be-done

- Adapt the geometry optimization to other ibrav values (currently working on 1,2,3,4 and 6)
- Compute the band gap at specific K points for number_Kohn_Sham_states >= 20
- Change the K-grid to be consistent with crystal axis (currently kx=ky=kz)
- Organize the code with comments. 

### Additional Requirements

- Quantum ESPRESSO (version 7.2 or higher - except 7.4 which is not compatible with PAOFLOW)
- MPI implementation (e.g., OpenMPI)

To install Quantum ESPRESSO, please follow the instructions at:
https://www.quantum-espresso.org/installation/

Note: Replace this section with specific QE version requirements for your scripts.
