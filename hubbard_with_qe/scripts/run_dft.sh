#!/bin/bash

path_scripts="$PWD/scripts"


path_init=$PWD

# File containing the input information (change name here)
file="input"

# Creation of variables with the material information read from the input
material=$(awk 'f && !/END MATERIAL/ {print} /BEGIN MATERIAL/ {f=1} /END MATERIAL/ {f=0}' ${file})
list_pseudos=$(awk 'f && !/END PSEUDO/ {print} /BEGIN PSEUDO/ {f=1} /END PSEUDO/ {f=0}' ${file} | awk '{print $1}') 
bands_kpoints=$(awk 'f && !/END KPOINTS_BANDS/ {print} /BEGIN KPOINTS_BANDS/ {f=1} /END KPOINTS_BANDS/ {f=0}' ${file}) 
path_pseudo=$(awk 'f && !/END PATH/ {print} /BEGIN PATH/ {f=1} /END PATH/ {f=0}' ${file})
path_qe=$(awk 'f && !/END QE_PATH/ {print} /BEGIN QE_PATH/ {f=1} /END QE_PATH/ {f=0}' ${file})
first_kpoint=$(awk 'f && !/END KPOINTS_BANDS/ {print} /BEGIN KPOINTS_BANDS/ {f=1} /END KPOINTS_BANDS/ {f=0}' ${file} | tail -n 2 | head -n 1 | awk '{print $NF}' | tr -d '!')
second_kpoint=$(awk 'f && !/END KPOINTS_BANDS/ {print} /BEGIN KPOINTS_BANDS/ {f=1} /END KPOINTS_BANDS/ {f=0}' ${file} | tail -n 1 | awk '{print $NF}' | tr -d '!')

hpin=$(awk 'f && !/END HP/ {print} /BEGIN HP/ {f=1} /END HP/ {f=0}' ${file})
bandsin=$(awk 'f && !/END BANDS/ {print} /BEGIN BANDS/ {f=1} /END BANDS/ {f=0}' ${file})
procs=$(awk 'f && !/END PROCS/ {print} /BEGIN PROCS/ {f=1} /END PROCS/ {f=0}' ${file})


# Remove files to not be overwritten
if test -f compare_pseudos; then rm compare_pseudos; fi
if test -f output_dft; then rm output_dft; fi

# Make header of file with results
echo "pseudo 	opt_a 	opt_c/a 	DFT_indirect	DFT_gap${second_kpoint} 	DFT_gap${first_kpoint}" > aux
awk '{printf "%30s%30s%30s%30s%30s%30s\n", $1, $2, $3, $4, $5, $6}' aux >> compare_pseudos
rm aux

echo "Reading information from input..." >> ${path_init}/output_dft

echo "Material name: ${material}" >> ${path_init}/output_dft
echo "Pseudopotentials for the transition metal:" >> ${path_init}/output_dft
printf "%s\n" "${list_pseudos}" >> ${path_init}/output_dft
echo "Path to pseudopotential files: ${path_pseudo}" >> ${path_init}/output_dft
echo "Kpoints to compute bands:" >> ${path_init}/output_dft
printf "%s\n" "${bands_kpoints}" | tail -n +3 >> ${path_init}/output_dft

printf "Starting the calculations.. \n" >> ${path_init}/output_dft

for p in ${list_pseudos}; do
	echo "Running calculations for pseudopotential ${p}..." >> ${path_init}/output_dft

	folder=$(basename "${p}" .UPF | sed 's/\.upf$//')
	echo $folder
	mkdir $folder; cd $folder
	rm *.in
	cp ../${material}.scf.in .

	# Replace pseudopotential name and path in the input file for QE
	sed -i "s|NAME2|O.pbe-n-kjpaw_psl.0.1.UPF|g" ${material}.scf.in
	sed -i "s|NAME1|${p}|g" ${material}.scf.in
	sed -i "s|#PATH|${path_pseudo}|g" ${material}.scf.in

	# CONVERGENCE TESTING ECUT
	printf "Running convergence tests for ECUT...\n" >> ${path_init}/output_dft

	mkdir testing_ecut; cd testing_ecut 

	head ../${material}.scf.in -n -3 > scf.in
	cp ${path_scripts}/convergence_ecut.sh .
	factor_ecut=$(grep ${p} ${path_init}/${file} | awk '{print $2}')
	echo "Considering the ratio: ecutrho = ${factor_ecut} times ecutwfc" >> ${path_init}/output_dft
	sed -i "s|#FACTOR_ECUT|${factor_ecut}|g" convergence_ecut.sh
	sed -i "s|#PROCS|${procs}|g" convergence_ecut.sh
	sed -i "s|#PATH_QE|${path_qe}|g" convergence_ecut.sh
	./convergence_ecut.sh

	cd ../

	# STRUCTURAL RELAXATION

	echo "Running structural relaxation..." >> ${path_init}/output_dft
	cp testing_ecut/converged.scf.in ${material}.vc-relax.in
	sed -i 's|scf|vc-relax|g' ${material}.vc-relax.in
	mpirun -np ${procs} ${path_qe}/pw.x < ${material}.vc-relax.in > ${material}.vc-relax.out

	nbnd=$(grep Sham ${material}.vc-relax.out | head -n 1 | awk '{print $NF}')
	
	echo "Number of Kohn-Sham states: ${nbnd}" >> ${path_init}/output_dft

	cp ${material}.vc-relax.in ${material}.scf.in
	sed -i 's|!nbnd|nbnd|g' ${material}.scf.in
	sed -i "s|nbnd *= *'nbnd'|nbnd = $((${nbnd}+5))|" ${material}.scf.in
	sed -i 's|vc-relax|scf|g' ${material}.scf.in

	ibrav=$(grep ibrav ${material}.vc-relax.in | awk '{print $NF}' | tr -d ',')

	grep -A30 "Begin final coordinates" ${material}.vc-relax.out | grep -B30 "End final coordinates" > final_coordinates

	v1x=$(grep -A1 "CELL_PARAMETERS" final_coordinates | tail -n 1 | awk '{print $1}'|  tr -d '-')
	v1y=$(grep -A1 "CELL_PARAMETERS" final_coordinates | tail -n 1 | awk '{print $2}'|  tr -d '-')
	v1z=$(grep -A1 "CELL_PARAMETERS" final_coordinates | tail -n 1 | awk '{print $3}'|  tr -d '-')

	v2x=$(grep -A2 "CELL_PARAMETERS" final_coordinates | tail -n 1 | awk '{print $1}'|  tr -d '-')
	v2y=$(grep -A2 "CELL_PARAMETERS" final_coordinates | tail -n 1 | awk '{print $2}'|  tr -d '-')
	v2z=$(grep -A2 "CELL_PARAMETERS" final_coordinates | tail -n 1 | awk '{print $3}'|  tr -d '-')	

	v3x=$(grep -A3 "CELL_PARAMETERS" final_coordinates | tail -n 1 | awk '{print $1}'|  tr -d '-')
	v3y=$(grep -A3 "CELL_PARAMETERS" final_coordinates | tail -n 1 | awk '{print $2}'|  tr -d '-')
	v3z=$(grep -A3 "CELL_PARAMETERS" final_coordinates | tail -n 1 | awk '{print $3}'|  tr -d '-')	

	cell=$(grep "CELL_PARAMETERS" final_coordinates | awk -F'=' '{print $2}' | tr -d ')')

	if [ $ibrav -eq 1 ]; then
		celldm1=$(echo "scale=10; $v1x * $cell" | bc)
		sed -i "s|^ *celldm(1) *= *[0-9.eE+-]*,|  celldm(1) = ${celldm1},|" ${material}.scf.in
		celldm3="1"

	elif [ $ibrav -eq 2 ] || [ $ibrav -eq 3 ]; then
		celldm1=$(echo "scale=10; ($v1x / 0.5) * $cell" | bc)
		sed -i "s|^ *celldm(1) *= *[0-9.eE+-]*,|  celldm(1) = ${celldm1},|" ${material}.scf.in
		celldm3="1"

	elif [ $ibrav -eq 4 ] || [ $ibrav -eq 6 ]; then
		celldm1=$(echo "scale=10; $v1x * $cell" | bc)
		celldm3=${v3z}
		sed -i "s|^ *celldm(1) *= *[0-9.eE+-]*,|  celldm(1) = ${celldm1},|" ${material}.scf.in
		sed -i "s|^ *celldm(3) *= *[0-9.eE+-]*,|  celldm(3) = ${celldm3},|" ${material}.scf.in
	else
		echo "ibrav ${ibrav} not supported. Please check the input file." >> ${path_init}/output_dft
		exit 1
	fi

	awk '/ATOMIC_POSITION/{flag=1; next} /End final coordinates/{flag=0} flag' final_coordinates > atomic_positions

	grep -B1000 ATOMIC_POSITIONS ${material}.scf.in > new_scf.in
	cat atomic_positions >> new_scf.in
	grep -A1000 K_POINTS ${material}.scf.in >> new_scf.in

	mv new_scf.in ${material}.scf.in

	# SCF AND BANDS CALCULATION
	

	echo "Running bands calculations..." >> ${path_init}/output_dft
	head ${material}.scf.in -n -2 > ${material}.bands.in
	printf "%s\n" "${bands_kpoints}" >> ${material}.bands.in
	sed -i 's|scf|bands|g' ${material}.bands.in
	
	printf "%s\n" "${bandsin}" > bands.in 
	sed -i "s|material|${material}|g" bands.in

	mpirun -np ${procs} ${path_qe}/pw.x < ${material}.scf.in > ${material}.scf.out
	mpirun -np ${procs} ${path_qe}/pw.x < ${material}.bands.in > ${material}.bands.out
	mpirun -np ${procs} ${path_qe}/bands.x < bands.in > bands.out

	vb=$(grep highest ${material}.scf.out | awk '{print $NF}')
	cb=$(grep highest ${material}.scf.out | awk '{print $(NF-1)}')
	indirect_gap=$(echo "${vb} - ${cb}" | bc -l)

	echo "Computing band gaps..." >> ${path_init}/output_dft

	if [ "${nbnd}" -gt 10 ]; then
		nx=1 
		ng=4
		col=$((${nbnd}-10))

		linex=$(tail -n ${nx} ${material}_bands.dat | head -n 1)
		lineg=$(sed -n "${ng}p" ${material}_bands.dat)


		value1=$(echo "$linex" | awk -v k="$(($col+1))" '{print $k}')
		value2=$(echo "$linex" | awk -v k="$(($col))" '{print $k}')

		differencex=$(echo "$value1 - $value2" | bc)

		value1=$(echo "$lineg" | awk -v k="$((${col}+1))" '{print $k}')
		value2=$(echo "$lineg" | awk -v k="${col}" '{print $k}')

		differenceg=$(echo "$value1 - $value2" | bc)
	elif [ "${nbnd}" -lt 10 ]; then
		nx=2 
		ng=3
		col=${nbnd}

		linex=$(tail -n ${nx} ${material}_bands.dat | head -n 1)
		lineg=$(sed -n "${ng}p" ${material}_bands.dat)


		value1=$(echo "$linex" | awk -v k="$(($col+1))" '{print $k}')
		value2=$(echo "$linex" | awk -v k="$(($col))" '{print $k}')

		differencex=$(echo "$value1 - $value2" | bc)

		value1=$(echo "$lineg" | awk -v k="$((${col}+1))" '{print $k}')
		value2=$(echo "$lineg" | awk -v k="${col}" '{print $k}')

		differenceg=$(echo "$value1 - $value2" | bc)
	elif [ "${nbnd}" -eq 10 ]; then

		value1=$(tail -n 1 ${material}_bands.dat | head -n 1 | awk '{print $1}')
		value2=$(tail -n 2 ${material}_bands.dat | head -n 1 | awk '{print $NF}')

		differencex=$(echo "$value1 - $value2" | bc)

		value1=$(sed -n "4p" ${material}_bands.dat | head -n 1 | awk '{print $1}')
		value2=$(sed -n "3p" ${material}_bands.dat | head -n 1 | awk '{print $NF}')

		differenceg=$(echo "$value1 - $value2" | bc)
	fi
	echo "${p}	${celldm1}	${celldm3}	${indirect_gap} 	${differencex}		${differenceg}" > aux
	awk '{printf "%30s%30.4f%30.4f%30.4f%30.4f%30.4f\n", $1, $2, $3, $4, $5, $6}' aux >> ../compare_pseudos
	rm aux	
	echo "Complete!" >> ${path_init}/output_dft
	cd ../
done
