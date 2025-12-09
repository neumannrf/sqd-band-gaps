#!/bin/bash

path_init=$PWD

if test -f output_hp_convergence; then rm output_hp_convergence; fi
if test -f convergence_vs; then rm convergence_vs; fi

# File containing the input information (change name here)
file="input"

# Creation of variables with the material information read from the input
material=$(awk 'f && !/END MATERIAL/ {print} /BEGIN MATERIAL/ {f=1} /END MATERIAL/ {f=0}' ${file})
list_pseudos=$(awk 'f && !/END PSEUDO/ {print} /BEGIN PSEUDO/ {f=1} /END PSEUDO/ {f=0}' ${file} | awk '{print $1}') 
hpin=$(awk 'f && !/END HP/ {print} /BEGIN HP/ {f=1} /END HP/ {f=0}' ${file})
procs=$(awk 'f && !/END PROCS/ {print} /BEGIN PROCS/ {f=1} /END PROCS/ {f=0}' ${file})

echo "Running convergence tests for HP" >> ${path_init}/output_hp_convergence

for p in ${list_pseudos}; do

    echo "Pseudopotential ${p}..."  >> ${path_init}/output_hp_convergence
	folder=$(basename "${p}" .UPF | sed 's/\.upf$//')
	cd $folder

    mkdir testing_hp; cd testing_hp

    printf "%s\n" "${hpin}" > hp.in 
	sed -i "s|material|${material}|g" hp.in
    head ../${material}.scf.in -n -1 > scf.in

    echo "Selecting k-mesh for scf calculation..." >> ${path_init}/output_hp_convergence

    grep ! ../testing_ecut/scf.k*.out | sort -t'k' -k2,2n | awk '{printf "%5.3f\n", $5}' > aux
    nl=$(awk 'NR==1 { ref=$1; next } $1 != ref { print NR; found=1; exit } END { if (!found) print 1 }' aux)
    k=$(grep ! ../testing_ecut/scf.k*.out | sort -t'k' -k2,2n |  sed -n "${nl}p" | awk '{print $1}' | cut -d 'k' -f 2 | cut -d '.' -f 1)
	rm aux

    echo "$k $k $k 0 0 0" >> scf.in
    tail ${path_init}/${material}.scf.in -n 3 >> scf.in
    sed -i '/ntyp/a   nosym = .true.' scf.in

    echo "Running scf with ${k}x${k}x${k} mesh..." >> ${path_init}/output_hp_convergence
    mpirun -np ${procs} pw.x < scf.in > scf.out

    for q in $(($k/2 - 1)) $(($k/2)) $(($k/2 + 1)); do
        echo "Running hp.x with q-grid ${q}x${q}x${q}" >> ${path_init}/output_hp_convergence
        cp hp.in hp.q$q.in
        sed -i "s/#q/$q/g" hp.q$q.in
        mpirun -np ${procs} hp.x < hp.q$q.in > hp.q$q.out
        mv  *.Hubbard_parameters.dat Hubbard_parameters_q$q.dat
        mv HUBBARD.dat HUBBARD_q$q.dat

        v=$(sed -n '4p' HUBBARD_q$q.dat| awk '{print $6}')
        echo "$q $v" >> convergence_vs
    done

    nl=$(awk 'NR==1 {prev=$2; next} $2 != prev {print NR, $2; exit}' convergence_vs | awk '{print $1}')
    q_final=$(sed -n "${nl}p" convergence_vs | awk '{print $1}')

    cp hp.q${q_final}.in hp_correct.in
    cp hp.q${q_final}.out hp_correct.out
    cp HUBBARD.q${q_final}.dat HUBBARD_correct.dat

    echo "Complete"
    cd ${path_init}
done