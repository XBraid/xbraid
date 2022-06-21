#!/bin/bash

##### These lines are for SLURM
#SBATCH -N 19
#SBATCH -p pbatch
#SBATCH -A paratime
#SBATCH -t 30
#SBATCH -o out.%j
#SBATCH -e err.%j
##### These are shell commands

echo -n 'This machine is '; hostname
echo -n 'My jobid is '; echo $SLURM_JOBID
echo -n 'Timestamp START: ';date

# number of cores
ncores="16 32 64 128 256 512 1024"

# levels
mlevels="4"

# coarsening factors
cfactors="4"

# Delta ranks
ranks="2 4 8 16 32 64"

# fixed arguments
fargs="-tf 4 -nt 2048 -nx 64 -nu 1 -nu0 1 -tol 1e-6 -theta -Deltalvl 1"

# path to executable
ex="../../../drive-ks"

# output directory
outd="."

# output fname
outn="ks"

# serial run
# mpirun -n 1 ${ex} ${fargs} -ml 1 > ${outd}/${outn}_ml1
echo "srun -N 1 -n 1 -o ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1"
srun -N 1 -n 1 -o ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1

for nc in $ncores; do
   for ml in $mlevels; do
      for cf in $cfactors; do
         for rank in $ranks; do
	    echo "srun -N 19 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_cf${cf}_ml${ml}     ${ex} ${fargs} -cf ${cf} -ml ${ml} -Delta -rank ${rank}"
            srun -N 19 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_cf${cf}_ml${ml}     ${ex} ${fargs} -cf ${cf} -ml ${ml} -Delta -rank ${rank}
         done
         echo "srun -N 19 -n ${nc} -o ${outd}/${outn}_nc${nc}_cf${cf}_ml${ml}           ${ex} ${fargs} -cf ${cf} -ml ${ml}"
         srun -N 19 -n ${nc} -o ${outd}/${outn}_nc${nc}_cf${cf}_ml${ml}           ${ex} ${fargs} -cf ${cf} -ml ${ml}
         # mpirun -n ${nc} ${ex} ${fargs} -ml ${ml} -cf ${cf} -Delta > ${outd}/${outn}_Delta_nc${nc}_cf${cf}_ml${ml}
         # mpirun -n ${nc} ${ex} ${fargs} -ml ${ml} -cf ${cf}        > ${outd}/${outn}_nc${nc}_cf${cf}_ml${ml}
      done
   done
done

echo 'Done'
echo -n 'Timestamp END: ';date
