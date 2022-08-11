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
ncores="16 32 64 128 256 512 1024 2048"
# ncores="2048 4096"

# Delta ranks
ranks="9"

# fixed arguments (coarse grid prop of LVs is important for convergence here)
fargs="-tf 8 -nt 16384 -nx 256 -cf 4 -cf0 16 -tol 1e-8 -mi 25 -cglv"

# path to executable
ex="../../../drive-ks"

# output directory
outd="."

# output fname
outn="ks"

# serial run
# echo "srun -N 1 -n 1 -o ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1"
# srun -N 1 -n 1 -o       ${outd}/${outn}_ml1 ${ex} ${fargs} -ml 1

for nc in $ncores; do
   # echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_nc${nc}       ${ex} ${fargs} -ml 2"
   # srun  -N 40 -n ${nc} -o ${outd}/${outn}_nc${nc}            ${ex} ${fargs} -ml 2

   # echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_theta_nc${nc} ${ex} ${fargs} -ml 2 -theta"
   # srun -N 40 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}       ${ex} ${fargs} -ml 2 -theta

   for rank in $ranks; do
      # echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_ml3 ${ex} ${fargs} -ml 3 -theta -Delta -rank ${rank}"
      # srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_ml3       ${ex} ${fargs} -ml 3 -theta -Delta -rank ${rank} 

      # echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_ml4 ${ex} ${fargs} -ml 4 -theta -Delta -rank ${rank}"
      # srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_ml4       ${ex} ${fargs} -ml 4 -theta -Delta -rank ${rank} 

      echo "srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_ml4_fmg ${ex} ${fargs} -ml 4 -theta -Delta -rank ${rank} -fmg"
      srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta${rank}_nc${nc}_ml4_fmg       ${ex} ${fargs} -ml 4 -theta -Delta -rank ${rank} -fmg
   done
done

echo 'Done'
echo -n 'Timestamp END: ';date
