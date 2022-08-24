#!/bin/bash

##### These lines are for SLURM
#SBATCH -N 40
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
ncores="32 64 128 256 512 1024 2048"

# Delta ranks
# ranks="4 8 9 16 32"
ranks="8 9 16"

# fixed arguments
fargs="-tf 4 -nt 8192 -nx 256 -tol 1e-8"
dargs="-theta -Delta -cglv"

# path to executable
ex="../../../drive-ks"

# output directory
outd="."

# output fname
outn="ks"

# serial run
# mpirun -n 1 ${ex} ${fargs} -ml 1 > ${outd}/${outn}_ml1
# echo "srun -N 1 -n 1 -o ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1"
# srun -N 1 -n 1 -o       ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1

for nc in $ncores; do
   	# echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_theta_fmg_nc${nc} ${ex} ${fargs} -ml 5 -theta -fmg"
   	# srun -N 40 -n ${nc} -o ${outd}/${outn}_theta_fmg_nc${nc}       ${ex} ${fargs} -ml 5 -theta -fmg

   	# echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_theta_rfmg_nc${nc} ${ex} ${fargs} -ml 5 -theta -refine -fmg"
   	# srun -N 40 -n ${nc} -o ${outd}/${outn}_theta_rfmg_nc${nc}       ${ex} ${fargs} -ml 5 -theta -refine -fmg

   for rank in $ranks; do
      # echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_rf_nc${nc} ${ex} ${fargs} -ml 5 ${dargs} -rank ${rank} -refine"
      # srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_rf_nc${nc}       ${ex} ${fargs} -ml 5 ${dargs} -rank ${rank} -refine

      echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_fmg_nc${nc} ${ex} ${fargs} -ml 5 ${dargs} -rank ${rank} -fmg"
      srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_fmg_nc${nc}       ${ex} ${fargs} -ml 5 ${dargs} -rank ${rank} -fmg

      echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_rfmg_nc${nc} ${ex} ${fargs} -ml 5 ${dargs} -rank ${rank} -fmg -refine"
      srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_rfmg_nc${nc}       ${ex} ${fargs} -ml 5 ${dargs} -rank ${rank} -fmg -refine
   done
done

echo 'Done'
echo -n 'Timestamp END: ';date
