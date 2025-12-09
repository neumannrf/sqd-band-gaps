#!/bin/bash

if test -f compare_gaps; then rm compare_gaps; fi
if test -f output_gaps_hubbard; then rm output_gaps_hubbard; fi


path_init=$PWD

file="input"

material=$(awk 'f && !/END MATERIAL/ {print} /BEGIN MATERIAL/ {f=1} /END MATERIAL/ {f=0}' ${file})
bandsin=$(awk 'f && !/END BANDS/ {print} /BEGIN BANDS/ {f=1} /END BANDS/ {f=0}' ${file})
list_pseudos=$(awk 'f && !/END PSEUDO/ {print} /BEGIN PSEUDO/ {f=1} /END PSEUDO/ {f=0}' ${file} | awk '{print $1}') 
procs=$(awk 'f && !/END PROCS/ {print} /BEGIN PROCS/ {f=1} /END PROCS/ {f=0}' ${file})
path_qe=$(awk 'f && !/END QE_PATH/ {print} /BEGIN QE_PATH/ {f=1} /END QE_PATH/ {f=0}' ${file})

echo "Starting calculations..." >> ${path_init}/output_gaps_hubbard

for p in ${list_pseudos}; do
	echo ${p}

	echo "Pseudopotential ${p}..." >> ${path_init}/output_gaps_hubbard

	folder=$(basename "${p}" .UPF | sed 's/\.upf$//')
	cd $folder

	header=$(head hubbard_computation/HUBBARD.dat -n 2| tail -n 1)
	tail hubbard_computation/HUBBARD.dat -n +3 > file_dftuv
	head hubbard_computation/HUBBARD.dat -n 3| tail -n 1 > file_dftu
	head hubbard_computation/HUBBARD.dat -n 4| tail -n 1 > file_dftv
	nbnd=$(grep Sham ${material}.vc-relax.out | head -n 1 | awk '{print $NF}')

	echo ${p} >> ../compare_gaps

	for run in dftuv dftu dftv; do

		echo "Calculation type: ${run}" >> ${path_init}/output_gaps_hubbard

		mkdir ${run}; cd ${run}
		rm *.in
		cp ../${material}.scf.in .
		cp ../${material}.bands.in .
		cp ../bands.in .


		printf "%s\n" "${header}" >> ${material}.scf.in
		cat ../file_${run} >> ${material}.scf.in

		printf "%s\n" "${header}" >> ${material}.bands.in
		cat ../file_${run} >> ${material}.bands.in


		echo "Running scf and bands calculations..." >> ${path_init}/output_gaps_hubbard

		mpirun -np ${procs} ${path_qe}/pw.x < ${material}.scf.in > ${material}.scf.out
		mpirun -np ${procs} ${path_qe}/pw.x < ${material}.bands.in > ${material}.bands.out
		mpirun -np ${procs} ${path_qe}/bands.x < bands.in > bands.out

		echo "Computing the bandgaps..." >> ${path_init}/output_gaps_hubbard

		vb=$(grep highest ${material}.scf.out | awk '{print $NF}')
		cb=$(grep highest ${material}.scf.out | awk '{print $(NF-1)}')
		indirect_gap=$(echo "${vb} - ${cb}" | bc -l)

		# Extract the second line and third line containing the values

		if [ "${nbnd}" -gt 10 ]; then
			nx=1 
			ng=4
			col=$((${nbnd}-10))

			linex=$(tail -n ${nx} ${material}_bands.dat | head -n 1)
			lineg=$(sed -n "${ng}p" ${material}_bands.dat)


			# Extract specific values using awk
			value1=$(echo "$linex" | awk -v k="$(($col+1))" '{print $k}')
			value2=$(echo "$linex" | awk -v k="$(($col))" '{print $k}')

			# Compute the difference
			differencex=$(echo "$value1 - $value2" | bc)

			# Extract specific values using awk
			value1=$(echo "$lineg" | awk -v k="$((${col}+1))" '{print $k}')
			value2=$(echo "$lineg" | awk -v k="${col}" '{print $k}')

			# Compute the difference
			differenceg=$(echo "$value1 - $value2" | bc)

		elif [ "${nbnd}" -lt 10 ]; then
			nx=2 
			ng=3
			col=${nbnd}

			linex=$(tail -n ${nx} ${material}_bands.dat | head -n 1)
			lineg=$(sed -n "${ng}p" ${material}_bands.dat)


			# Extract specific values using awk
			value1=$(echo "$linex" | awk -v k="$(($col+1))" '{print $k}')
			value2=$(echo "$linex" | awk -v k="$(($col))" '{print $k}')

			# Compute the difference
			differencex=$(echo "$value1 - $value2" | bc)

			# Extract specific values using awk
			value1=$(echo "$lineg" | awk -v k="$((${col}+1))" '{print $k}')
			value2=$(echo "$lineg" | awk -v k="${col}" '{print $k}')

			# Compute the difference
			differenceg=$(echo "$value1 - $value2" | bc)
		elif [ "${nbnd}" -eq 10 ]; then

			# Extract specific values using awk
			value1=$(tail -n 1 ${material}_bands.dat | head -n 1 | awk '{print $1}')
			value2=$(tail -n 2 ${material}_bands.dat | head -n 1 | awk '{print $NF}')

			# Compute the difference
			differencex=$(echo "$value1 - $value2" | bc)

			# Extract specific values using awk
			value1=$(sed -n "4p" ${material}_bands.dat | head -n 1 | awk '{print $1}')
			value2=$(sed -n "3p" ${material}_bands.dat | head -n 1 | awk '{print $NF}')

			# Compute the difference
			differenceg=$(echo "$value1 - $value2" | bc)
		fi

		# Print the result
		echo "${run} ${indirect_gap} 	${differencex}		${differenceg}" >> ../../compare_gaps

		echo "Complete!" >> ${path_init}/output_gaps_hubbard
		
		cd ../
	done
	cd ../
done
