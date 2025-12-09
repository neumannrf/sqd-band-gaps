#!/bin/bash

if test -f output_projection; then rm output_projection; fi

# source /opt/share/intel/composerxe/bin/compilervars.sh intel64
# export PATH=$PATH:/opt/share/openmpi-4.0.1/x86_64/bin:/dccstor/quant4clim/qe-7.2/bin:/dccstor/quant4clim/qe-7.2/PW/tools
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/share/openmpi-4.0.1/x86_64/lib

path_init=$PWD
path_scripts="$PWD/scripts"
file="input"

material=$(awk 'f && !/END MATERIAL/ {print} /BEGIN MATERIAL/ {f=1} /END MATERIAL/ {f=0}' ${file})
list_pseudos=$(awk 'f && !/END PSEUDO/ {print} /BEGIN PSEUDO/ {f=1} /END PSEUDO/ {f=0}' ${file} | awk '{print $1}') 
procs=$(awk 'f && !/END PROCS/ {print} /BEGIN PROCS/ {f=1} /END PROCS/ {f=0}' ${file})
projin=$(awk 'f && !/END PROJ/ {print} /BEGIN PROJ/ {f=1} /END PROJ/ {f=0}' ${file})
path_qe=$(awk 'f && !/END QE_PATH/ {print} /BEGIN QE_PATH/ {f=1} /END QE_PATH/ {f=0}' ${file})

echo "Starting calculations..." >> ${path_init}/output_projection

for p in ${list_pseudos}; do
	echo ${p}

	echo "Pseudopotential ${p}..." >> ${path_init}/output_projection

	folder=$(basename "${p}" .UPF | sed 's/\.upf$//')
	cd $folder; mkdir paoflow; cd paoflow
	rm *.in
	cp ../${material}.scf.in .

	# If one wants to increase the k-point grid to perform PAOFLOW projection, uncomment the next two lines:
	# head ../${material}.scf.in -n -1 >  ${material}.scf.in
	# echo "14 14 14 0 0 0" >>  ${material}.scf.in

	echo "Running scf calculation..." >> ${path_init}/output_projection

	mpirun -np ${procs} ${path_qe}/pw.x < ${material}.scf.in > ${material}.scf.out

	cp ${path_scripts}/make_hopping.py .
	cp ${path_scripts}/plot_hopping_bands.py .
	sed -i "s/#PSEUDO/${folder}/g" make_hopping.py
	sed -i "s/#MATERIAL/${material}/g" make_hopping.py

	echo "Performing PAOFLOW projection..." >> ${path_init}/output_projection

	python make_hopping.py > ${material}_${folder}_paoflow_output
	efermi=$(grep highest ${material}.scf.out | awk '{print $(NF-1)}')

	echo "Plotting hopping matrix and bands..." >> ${path_init}/output_projection

	python plot_hopping_bands.py out/${material}.save/ ../${material}_bands.dat.gnu ${efermi} ${material}_${folder}

	cd ../

	echo "Computing the atomic occupations..." >> ${path_init}/output_projection

	for sub in paoflow dftu dftv dftuv; do
		echo "For ${sub}..." >> ${path_init}/output_projection

		cd $sub; mkdir proj; cd proj
		rm *.in

		printf "%s\n" "${projin}" > proj.in
		sed -i "s|material|${material}|g" proj.in

		cp ../${material}.scf.in .
		sed -i "/restart_mode/a   verbosity = 'high'"  ${material}.scf.in
		cp ${path_scripts}/pdos.py .

		mpirun -np ${procs} ${path_qe}/pw.x < ${material}.scf.in > ${material}.scf.out
		mpirun -np ${procs} ${path_qe}/projwfc.x < proj.in > proj.out
		python pdos.py ${material}.scf.out proj.out

		mv pdos_all.txt ${material}_pdos_all_${folder}_${sub}.txt
		mv pdos_all.pkl ${material}_pdos_all_${folder}_${sub}.pkl

		mv pdos_gamma.txt ${material}_pdos_gamma_${folder}_${sub}.txt
		mv pdos_gamma.pkl ${material}_pdos_gamma_${folder}_${sub}.pkl

		mv pdos_info.txt ${material}_pdos_info_${folder}_${sub}.txt
		mv pdos_info.pkl ${material}_pdos_info_${folder}_${sub}.pkl

		cd ../../
	done
	cd ../
	echo "Complete!" >> ${path_init}/output_projection
done
