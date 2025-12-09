# sqd-lattice
My implementation of the Sample-based Quantum diagonalization and its variants applied to lattice models.

## Installation

To install the source code of the repository as a package and dependencies run `pip install -e .`. Next, install the qiskit-addon-sqd following the instructions in https://github.com/Qiskit/qiskit-addon-sqd. 

## Overview

The scripts are executed with the `main.zsh` script. After selecting the material with the variable `material`, set the flags of each calculation to `true` or `false`. To execute all calculations, turn all flags to `true`, and the script will automatically wait for the preceding calculation to start the next one. 

The flags of the `main.zsh`script perform the following actions

1. `initialize`: Perform the initial Hartree fock and CCSD calculations, and also build the corresponding circuits for sampling configurations. This script will create files for the electron-number subspaces chosen, and store the results in those folders.
2. `sample_qiskit`: Sample noiseless counts using the ffsim simulator
3. `run_hci`: Performs the HCI calculation for benchmarking (REQUIRES DICE)
4. `run_sqd`: Performs the SQD calculation
5. `run_ext_sqd`: Performs the Extended-SQD calculation

The following flags plot the results. 
