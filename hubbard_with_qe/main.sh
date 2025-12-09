#!/bin/bash

if test -f out; then rm out; fi

path_scripts="$PWD/scripts"
file='input'
material=$(awk 'f && !/END MATERIAL/ {print} /BEGIN MATERIAL/ {f=1} /END MATERIAL/ {f=0}' ${file})
metal="Zn"

######################################################################################

# Block 1: Create QE input from CIF using python script

echo "Generating QE input file from CIF" >> out

python ${path_scripts}/generate_scf.py ${material} "{\"${metal}\": \"NAME1\", \"O\": \"NAME2\"}" 1

if [ ! -f ${material}.scf.in ]; then
    echo "Error: File ${material}.scf.in not found."
    exit 1
fi

######################################################################################

# Block 2: Run DFT calculations and get band gaps

echo "Running DFT calculations and computing band gaps" >> out

${path_scripts}/run_dft.sh

######################################################################################

# Block 3: Run Hubbard calculations

echo "Computing Hubbard parameters" >> out

${path_scripts}/run_hubbard_computation.sh

######################################################################################

# Block 4: Run DFT+U+V, DFT+U and DFT+V

echo "Running DFT+U+V, DFT+U and DFT+V" >> out

${path_scripts}/run_hubbard.sh

######################################################################################

# Block 5: PAOFLOW projection and compute occupations

echo "Computing PAOFLOW projections and occupations" >> out

${path_scripts}/run_paoflow.sh

######################################################################################

# Block 6: Parse final results

echo "Parsing results and dumping in final_results.csv" >> out

${path_scripts}/parse.sh

# ######################################################################################

# Block 7: Final files for QC

selected_pseudo='Zn.pbe-d-hgh'

echo "Copying files for selected pseudo in outputs_for_QC/" >> out

cd ${selected_pseudo}/paoflow

python ${path_scripts}/extract_hubbard.py ${material}

cd ../../

mkdir outputs_for_QC
cp ${selected_pseudo}/paoflow/*.npy outputs_for_QC/

echo "Done" >> out
