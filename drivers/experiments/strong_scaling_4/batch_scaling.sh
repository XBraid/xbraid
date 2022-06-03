#!/bin/bash

##### These lines are for SLURM
#SBATCH -N 19
#SBATCH -p pbatch
#SBATCH -A paratime
#SBATCH -t 5
#SBATCH -o out.%j
#SBATCH -e err.%j
##### These are shell commands

echo -n 'This machine is '; hostname
echo -n 'My jobid is '; echo $SLURM_JOBID
echo -n 'Timestamp START: ';date

# number of cores
ncores="1 2 4 8 16 32 64 128 256 512 1024"
# ncores="1 2 4 8"

# levels
mlevels="4"

# coarsening factors
cfactors="4"

# fixed arguments
fargs="-tf 4 -nt 4096 -theta -nu 1 -niters 2 -tol 1e-6"

# path to executable
ex="../../drive-lorenz-Delta"

# output directory
outd="."

# output fname
outn="lorenz_theta"

# serial run
# mpirun -n 1 ${ex} ${fargs} -ml 1 > ${outd}/${outn}_ml1
echo "srun -N 1 -n 1 -o ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1"
srun -N 1 -n 1 -o ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1

for nc in $ncores; do
   for ml in $mlevels; do
      for cf in $cfactors; do
         echo "srun -N 19 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_cf${cf}_ml${ml} ${ex} ${fargs} -cf ${cf} -ml ${ml} -Delta"
         echo "srun -N 19 -n ${nc} -o ${outd}/${outn}_nc${nc}_cf${cf}_ml${ml}       ${ex} ${fargs} -cf ${cf} -ml ${ml}"
         srun -N 19 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_cf${cf}_ml${ml} ${ex} ${fargs} -cf ${cf} -ml ${ml} -Delta
         srun -N 19 -n ${nc} -o ${outd}/${outn}_nc${nc}_cf${cf}_ml${ml}       ${ex} ${fargs} -cf ${cf} -ml ${ml}
         # mpirun -n ${nc} ${ex} ${fargs} -ml ${ml} -cf ${cf} -Delta > ${outd}/${outn}_Delta_nc${nc}_cf${cf}_ml${ml}
         # mpirun -n ${nc} ${ex} ${fargs} -ml ${ml} -cf ${cf}        > ${outd}/${outn}_nc${nc}_cf${cf}_ml${ml}
      done
   done
done

echo 'Done'
echo -n 'Timestamp END: ';date
