#!/bin/bash

##### These lines are for SLURM
#SBATCH -N 19
#SBATCH -p pbatch
#SBATCH -A paratime
#SBATCH -t 4
#SBATCH -o out.%j
#SBATCH -e err.%j
##### These are shell commands

echo -n 'This machine is '; hostname
echo -n 'My jobid is '; echo $SLURM_JOBID
echo -n 'Timestamp START: ';date

# number of cores
ncores="1 2 4 8 16 32 64 128 512 1024"
# ncores="1 2 4 8"

# levels
mlevels="3 4"

# fixed arguments
fargs="-tf 4 -nt 4098 -cf 4 -theta -nu 1"

# path to executable
ex="../../drive-lorenz-Delta"

# output directory
outd="."

# output fname
outn="lorenz_theta"

# serial run
# mpirun -n 1 drive-lorenz-Delta ${fargs} -ml 1 > ${outd}/${outn}_Delta_nc${nc}_ml1
srun -N 1 -n 1 -o ${outd}/${outn}_Delta_nc${1}_ml1 ${ex} ${fargs} -ml 1

for nc in $ncores; do
   for ml in $mlevels; do
      # echo "srun -N 2 -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} -n ${nc} ${ex} ${fargs} -ml ${ml} -Delta"
      srun -N 19 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -ml ${ml} -Delta
      srun -N 19 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -ml ${ml}
      # mpirun -n ${nc} ${ex} ${fargs} -ml ${ml} -Delta > ${outd}/${outn}_Delta_nc${nc}_ml${ml}
      # mpirun -n ${nc} ${ex} ${fargs} -ml ${ml}        > ${outd}/${outn}_nc${nc}_ml${ml}
   done
done

echo 'Done'
echo -n 'Timestamp END: ';date
