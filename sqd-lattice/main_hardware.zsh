#!/bin/zsh

# Give the name of the material (same as folder)
# material="zirconia"
# material="platinum_V0"
# material="hafnium"
# material="hafnium_custom"
material="zirconia_molecular"
# material="test"
# material="silicon_test"
# material="silicon"
# material="silicon_custom"
# material="silicon_complex"

# Selecting which electron number to consider for calculations
# For all electron numbers do selected_indexes=(1 2 3)
selected_indexes=(1 2 3) 

base_dir="runs/$material"

echo "Collecting hardware counts for material for $material"
##############################################################################################################
# Selecting the folders of electron numbers
folders=($base_dir/*e(/))
selected_folders=()
for i in $selected_indexes; do
  selected_folders+="$folders[$i]"
done

echo "Running for electron numbers"
echo $selected_folders

###############################################################################################################

# Running hardware scripts
for dir in $selected_folders; do
  folder_name=${dir:t}
  number=${folder_name%[a-zA-Z]}  # remove the trailing 'e'
  echo "Collecting hardware counts for $number electrons..."
  nohup python scripts/run_hardware.py --material $material --ne $number > "$dir/logs/log_hardware.log" 2>&1 &
done

echo "All hardware jobs finished!"

echo "Finished!" 
