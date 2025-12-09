#!/bin/zsh

# Give the name of the material (same as folder)
# material="hafnium"
# material="hafnium_2"
# material="hafnium_mol"
# material="zirconia"
material="zirconia_2"
# material="platinum_2"
# material="platinum_V0"
# material="platinum_V0"
# material="silicon"

# Select calculations to be 
initialize=false
sample_qiskit=false
run_hci=false
run_sqd=false
run_ext_sqd=false

# Selecting which electron number to consider for calculations
# For all electron numbers do selected_indexes=(1 2 3)
selected_indexes=(1 2 3) 

# Select plots to be done
energy=false
subspace_dimension=true
subspace_dimension_hci=false

# Figures for paper
figure_1=false
figure_6=false

##############################################################################################################

base_dir="runs/$material"

echo "Running for $material"
echo "Initilize $initialize"
echo "Sample qiskit $sample_qiskit"
echo "Run SQD $run_sqd"

# Run pyscf and make circuits
if [[ $initialize == "true" ]]; then
  echo "Running initialization..."
  nohup python scripts/make_circuits.py --material $material > "$base_dir/log_make_circuits.log" 2>&1 &
  PID_make_circuits=$!
  wait $PID_make_circuits
  echo "Pre processing done!"
fi

# Selecting the folders of electron numbers
folders=($base_dir/*e(/))
selected_folders=()
for i in $selected_indexes; do
  selected_folders+="$folders[$i]"
done
###############################################################################################################

# Run HCI
if [[ $run_hci == "true" ]]; then
  pids=()
  for dir in $selected_folders; do
    folder_name=${dir:t}      
    number=${folder_name%[a-zA-Z]}  # remove the trailing 'e'
    echo "Running HCI for $number electrons in folder: $dir"

    # Get all iterations values from hci_config.yaml
    # iterations_list=("${(@f)$(yq '.max_iter[]' "$dir/hci_config.yaml")}")

    i=$(yq '.select_cutoff_steps' "$dir/hci_config.yaml")

    for j in {0..$((i - 1))}; do
      echo "  → Running with select_cutoff = $j"
      nohup python scripts/run_hci.py \
        --material "$material" \
        --ne "$number" \
        --select_cutoff_idx "$j" \
        > "$dir/logs/log_hci_${j}.log" 2>&1 &
      sleep 1
      pids+=($!)
    done
  done

  for pid in $pids; do
    wait $pid
  done
  echo "All HCI jobs finished!"
fi

##############################################################################################################

# Sample with Qiskit sampler
if [[ $sample_qiskit == "true" ]]; then
  pids=()  # Array to store process IDs
  for dir in $selected_folders; do
    folder_name=${dir:t}
    number=${folder_name%[a-zA-Z]}  # remove the trailing 'e'
    echo "Sampling qiskit for $number electrons..."
    nohup python scripts/sample_qiskit.py --material $material --ne $number > "$dir/logs/log_sample_qiskit.log" 2>&1 &
    pids+=($!)  # Store the PID of the background process
  done
  # Wait for all background processes to finish
  for pid in $pids; do
    wait $pid
  done
  echo "All Qiskit sampling scripts finished!"
fi

##############################################################################################################

# Generate MPS approximations of the circuit
if [[ $generate_mps == "true" ]]; then
  pids=()  # Array to store process IDs
  for dir in $selected_folders; do
    folder_name=${dir:t}
    number=${folder_name%[a-zA-Z]}  # remove the trailing 'e'
    echo "Generating MPS for $number electrons..."
    nohup python scripts/generate_mps.py --material $material --ne $number > "$dir/logs/log_generate_mps.log" 2>&1 &
    pids+=($!)  # Store the PID of the background process
  done
  # Wait for all background processes to finish
  for pid in $pids; do
    wait $pid
  done
  echo "All MPS generation scripts finished!"
fi

##############################################################################################################

# Generate MPS approximations of the circuit
if [[ $sample_mps == "true" ]]; then
  pids=()  # Array to store process IDs
  for dir in $selected_folders; do
    folder_name=${dir:t}
    number=${folder_name%[a-zA-Z]}  # remove the trailing 'e'
    echo "Sampling MPS for $number electrons..."
    nohup python scripts/sample_mps.py --material $material --ne $number > "$dir/logs/log_sample_mps.log" 2>&1 &
    pids+=($!)  # Store the PID of the background process
  done
  # Wait for all background processes to finish
  for pid in $pids; do
    wait $pid
  done
  echo "All MPS sampling scripts finished!"
fi

##############################################################################################################

# Run SQD
if [[ $run_sqd == "true" ]]; then
  pids=()
  for dir in $selected_folders; do
    folder_name=${dir:t}      
    number=${folder_name%[a-zA-Z]}  # remove the trailing 'e'
    echo "Processing $number electrons in folder: $dir"

    # Get all samples_per_batch values from sqd_config.yaml
    sample_list=("${(@f)$(yq '.samples_per_batch[]' "$dir/sqd_config.yaml")}")

    for spb in $sample_list; do
      echo "  → Running with samples_per_batch = $spb"
      nohup python scripts/run_sqd.py \
        --material "$material" \
        --ne "$number" \
        --samples_per_batch "$spb" \
        > "$dir/logs/log_sqd_${spb}.log" 2>&1 &
      sleep 5
      pids+=($!)
    done
  done

  for pid in $pids; do
    wait $pid
  done
  echo "All SQD jobs finished!"
fi

echo "Done!"

##############################################################################################################

# Create extended counts
if [[ $build_excited_counts == "true" ]]; then
  pids=()  # Array to store process IDs
  for dir in $selected_folders; do
    folder_name=${dir:t}
    number=${folder_name%[a-zA-Z]}  # remove the trailing 'e'
    echo "Building excited counts for $number electrons..."
    nohup python scripts/build_excited_counts.py --material $material --ne $number > "$dir/logs/log_build_excited_counts.log" 2>&1 &
    pids+=($!)  # Store the PID of the background process
  done
  # Wait for all background processes to finish
  for pid in $pids; do
    wait $pid
  done
  echo "All MPS sampling scripts finished!"
fi

##############################################################################################################

# Run Ext-SQD
if [[ $run_ext_sqd == "true" ]]; then
  pids=()
  for dir in $selected_folders; do
    folder_name=${dir:t}      
    number=${folder_name%[a-zA-Z]}  # remove the trailing 'e'
    echo "Running Ext-SQD for $number electrons in folder: $dir"

    # Get all samples_per_batch values from sqd_config.yaml
    sample_list=("${(@f)$(yq '.ext_sqd.samples_per_batch[]' "$dir/sqd_config.yaml")}")

    for spb in $sample_list; do
      echo "  → Running Ext-SQD with samples_per_batch = $spb"
      nohup python scripts/run_ext_sqd.py \
        --material "$material" \
        --ne "$number" \
        --samples_per_batch "$spb" \
        > "$dir/logs/log_ext_sqd_${spb}.log" 2>&1 &
      sleep 5
      pids+=($!)
    done
  done

  for pid in $pids; do
    wait $pid
  done
  echo "All Ext-SQD jobs finished!"
fi

echo "Done!"


########################################################################################################
########################################################################################################
########################################################################################################

if [[ $bandap_comparison == "true" ]]; then
  echo "Plotting band gap comparison..."
  nohup python plotting_scripts/bandgap_comparison.py --material $material > "$base_dir/log_plot_bandgap_comparison.log" 2>&1 &
  echo "Band gap comparison done!"
fi

########################################################################################################

if [[ $energy == "true" ]]; then
  echo "Plotting energies..."
  nohup python plotting_scripts/energy.py --material $material > "$base_dir/log_plot_energies.log" 2>&1 &
  echo "Energy plot done!"
fi

########################################################################################################

if [[ $probabilities == "true" ]]; then
  echo "Plotting band gap comparison..."
  nohup python plotting_scripts/probabilities.py --material $material > "$base_dir/log_plot_probabilities.log" 2>&1 &
  echo "Band gap comparison done!"
fi

##########################################################################################################

if [[ $energy_variance == "true" ]]; then
  for dir in $selected_folders; do
    folder_name=${dir:t}
    number=${folder_name%[a-zA-Z]}  # remove the trailing 'e'
    echo "Plotting energy variance for $number electrons..."
    nohup python plotting_scripts/energy_variance.py --material $material --ne $number > "$base_dir/log_plot_energy_variance_$number.log" 2>&1 &
    pids+=($!)  # Store the PID of the background process
  done
  echo "All energy variance plots finished!"
fi

echo "Finished!" 

##############################################################################################################

if [[ $subspace_dimension_hci == "true" ]]; then
  echo "Plotting subspace dimension..."
  nohup python plotting_scripts/subspace_dimension_hci.py --material $material > "$base_dir/log_plot_subpsace_dimension_hci.log" 2>&1 &
  echo "Subspace dimension HCI done!"
fi

##############################################################################################################

if [[ $subspace_dimension == "true" ]]; then
  echo "Plotting subspace dimension..."
  nohup python plotting_scripts/subspace_dimension.py --material $material > "$base_dir/log_plot_subpsace_dimension.log" 2>&1 &
  echo "Subspace dimension done!"
fi

##############################################################################################################

if [[ $subspace_comparison == "true" ]]; then
  echo "Plotting subspace comparison..."
  nohup python plotting_scripts/subspace_comparison.py --material $material > "$base_dir/log_plot_subpsace_comparison.log" 2>&1 &
  echo "Subspace dimension done!"
fi

########################################################################################################
########################################################################################################
########################################################################################################

if [[ $figure_1 == "true" ]]; then
  echo "Plotting figure 1..."
  nohup python plotting_scripts/figure_1.py --material $material > "$base_dir/log_plot_figure_1.log" 2>&1 &
  echo "Figure 1 done!"
fi

########################################################################################################

if [[ $figure_6 == "true" ]]; then
  echo "Plotting figure 6..."
  nohup python plotting_scripts/figure_6.py --material $material > "$base_dir/log_plot_figure_6.log" 2>&1 &
  echo "Figure 6 done!"
fi