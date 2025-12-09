#!/bin/bash

if test -f output_hubbard_comp; then rm output_hubbard_comp; fi


file="input"

material=$(awk 'f && !/END MATERIAL/ {print} /BEGIN MATERIAL/ {f=1} /END MATERIAL/ {f=0}' ${file})
hpin=$(awk 'f && !/END HP/ {print} /BEGIN HP/ {f=1} /END HP/ {f=0}' ${file})
procs=$(awk 'f && !/END PROCS/ {print} /BEGIN PROCS/ {f=1} /END PROCS/ {f=0}' ${file})
list_pseudos=$(awk 'f && !/END PSEUDO/ {print} /BEGIN PSEUDO/ {f=1} /END PSEUDO/ {f=0}' ${file} | awk '{print $1}')
path_qe=$(awk 'f && !/END QE_PATH/ {print} /BEGIN QE_PATH/ {f=1} /END QE_PATH/ {f=0}' ${file})

path_init=$PWD

echo "Starting calculations..." >> ${path_init}/output_hubbard_comp

for p in ${list_pseudos}; do
	echo ${p}

	echo "Pseudopotential ${p}..." >> ${path_init}/output_hubbard_comp

	folder=$(basename "${p}" .UPF | sed 's/\.upf$//')
	cd $folder; mkdir hubbard_computation; cd hubbard_computation
	rm *.in
	head ../${material}.scf.in -n -1 >  ${material}.scf.in

	echo "Selecting k-mesh ..." >> ${path_init}/output_hubbard_comp

	grep ! ../testing_ecut/scf.k*.out | sort -t'k' -k2,2n | awk '{printf "%5.3f\n", $5}' > aux
	nl=$(awk 'NR==1 { ref=$1; next } $1 != ref { print NR; found=1; exit } END { if (!found) print 1 }' aux)
	k=$(grep ! ../testing_ecut/scf.k*.out | grep ! ../testing_ecut/scf.k*.out | sort -t'k' -k2,2n |  sed -n "${nl}p" | awk '{print $1}' | cut -d 'k' -f 2 | cut -d '.' -f 1)
	rm aux

	echo "$k $k $k 0 0 0" >>  ${material}.scf.in
    tail ${path_init}/${material}.scf.in -n 3 >>  ${material}.scf.in
    sed -i '/ntyp/a   nosym = .true.'  ${material}.scf.in

	printf "%s\n" "${hpin}" > hp.in
	sed -i "s|material|${material}|g" hp.in

	echo "Running scf calculation with k-mesh ${k}x${k}x${k}" >> ${path_init}/output_hubbard_comp
	mpirun -np ${procs} ${path_qe}/pw.x < ${material}.scf.in > ${material}.scf.out

	q=$(($k/2))
	sed -i "s/#q/${q}/g" hp.in

	echo "Computing Hubbard parameters with q-grid ${q}x${q}x${q}" >> ${path_init}/output_hubbard_comp
	mpirun -np ${procs} ${path_qe}/hp.x < hp.in > hp.out
	cd ${path_init}

	echo "Complete!" >> ${path_init}/output_hubbard_comp
done
