#!/bin/bash

path_qe=#PATH_QE
factor_ecut=#FACTOR_ECUT
procs=#PROCS


if test -f scf.ecut.in; then rm scf.ecut*.in; fi
sed "s/# K # K # K/ 8 8 8/g" scf.in > scf.ecut.in


for ecut in 80 100 120 140 160 180 200; do
      cp scf.ecut.in scf.ecut${ecut}.in
      sed -i "s/ecutwfc *= *'#ECUTWFC'/ecutwfc = ${ecut}/" scf.ecut${ecut}.in
		sed -i "s/ecutrho *= *'#ECUTRHO'/ecutrho = $((${factor_ecut}*${ecut}))/g" scf.ecut${ecut}.in
     	mpirun -np ${procs} ${path_qe}/pw.x < scf.ecut${ecut}.in > scf.ecut${ecut}.out
done

grep '!' scf.ecut*.out | sed -E 's/.*scf\.ecut([0-9]+)\.out.*/\1 &/' | sort -n | cut -d' ' -f2- | awk '{printf "%5.3f\n", $5}'> aux
nl=$(awk 'NR==1 { ref=$1; next } $1 != ref { print NR; found=1; exit } END { if (!found) print 1 }' aux)
ecut_final=$(grep '!' scf.ecut*.out | sed -E 's/.*scf\.ecut([0-9]+)\.out.*/\1 &/' | sort -n | cut -d' ' -f2- | sed -n "${nl}p" | awk '{print $1}' | cut -d '.' -f 2 | cut -d 't' -f 2)
rm aux

if test -f scf.k.in; then rm scf.k*.in; fi
cp scf.ecut${ecut_final}.in scf.k.in

for k in 4 6 8 10 12 14 16; do
	cp scf.k.in scf.k$k.in
	sed -i "s/8 8 8/$k $k $k/g" scf.k$k.in
	mpirun -np ${procs} ${path_qe}/pw.x < scf.k$k.in > scf.k$k.out
done

grep ! scf.k*.out | sort -t'k' -k2,2n | awk '{printf "%5.3f\n", $5}' > aux
nl=$(awk 'NR==1 { ref=$1; next } $1 != ref { print NR; found=1; exit } END { if (!found) print 1 }' aux)
k=$(grep ! scf.k*.out | sort -t'k' -k2,2n |  sed -n "${nl}p" | awk '{print $1}' | cut -d 'k' -f 2 | cut -d '.' -f 1)

sed "s/ecutwfc *= *'#ECUTWFC'/ecutwfc = ${ecut_final}/" scf.in > converged.scf.in
sed -i "s/ecutrho *= *'#ECUTRHO'/ecutrho = $((${factor_ecut}*${ecut_final}))/" converged.scf.in

# PAOFLOW needs a lot of K-points to do the projection correctly (here defined as a 14x14x14 mesh)
sed -i "s/# K/14/g" converged.scf.in
# You can also set the K-points to the value you found in the convergence test
# sed -i "s/# K/${k}/g" converged.scf.in

grep ! scf.k*.out > k_output
grep ! scf.ecut*.out > ecut_output
