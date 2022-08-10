#!/bin/bash

##### These lines are for SLURM
#SBATCH -N 19
#SBATCH -p pbatch
#SBATCH -A paratime
#SBATCH -t 60
#SBATCH -o out.%j
#SBATCH -e err.%j
##### These are shell commands

echo -n 'This machine is '; hostname
echo -n 'My jobid is '; echo $SLURM_JOBID
echo -n 'Timestamp START: ';date

# number of cores
ncores="16 32 64 128 256 512 1024"

# Delta ranks
ranks="2 4 8 9 16 32 64"

# fixed arguments
fargs="-tf 4 -nt 8192 -nx 512 -cf 4 -cf0 16 -nu 1 -nu0 1 -tol 1e-8 -niters 2 -mi 25"
dargs="-theta -Delta -Deltalvl 1"

# path to executable
ex="../../../drive-ks"

# output directory
outd="."

# output fname
outn="ks"

# serial run
# mpirun -n 1 ${ex} ${fargs} -ml 1 > ${outd}/${outn}_ml1
echo "srun -N 1 -n 1 -o ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1"
srun -N 1 -n 1 -o       ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1

for nc in $ncores; do
   #	echo "srun -N 19 -n ${nc} -o ${outd}/${outn}_nc${nc}       ${ex} ${fargs} -ml 2"
   #	srun  -N 19 -n ${nc} -o ${outd}/${outn}_nc${nc}            ${ex} ${fargs} -ml 2

   #	echo "srun -N 19 -n ${nc} -o ${outd}/${outn}_theta_nc${nc} ${ex} ${fargs} -ml 3 -theta"
   #	srun -N 19 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}       ${ex} ${fargs} -ml 3 -theta

   for rank in $ranks; do
      echo "srun -N 19 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc} ${ex} ${fargs} -ml 3 ${dargs} -rank ${rank}"
      srun -N 19 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}       ${ex} ${fargs} -ml 3 ${dargs} -rank ${rank} 
   done
   for rank in $ranks; do
      echo "srun -N 19 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_cf16 ${ex} ${fargs} -ml 3 ${dargs} -niters 3 -cf 16 -rank ${rank}"
      srun -N 19 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_cf16       ${ex} ${fargs} -ml 3 ${dargs} -niters 3 -cf 16 -rank ${rank}
   done
done

echo 'Done'
echo -n 'Timestamp END: ';date
