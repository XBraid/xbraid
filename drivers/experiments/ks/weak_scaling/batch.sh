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
echo -n 'Timestamp START: '; date

# number of cores
ncores="8 32 128 512"

# timegrid size
ntimes="128 512 2048 8192"

# levels
mlevels="2 3 4 4"

# fixed arguments
fargs="-tf 4 -nx 128 -cf 4"
dargs="-theta -Delta -rank 10"

# path to executable
ex="../../../drive-ks"

# output directory
outd="."

# output fname
outn="ks"

# 8 procs
nt="128"; nc="8"; ml="2"
srun -N 1 -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nt ${nt} -ml 1
srun -N 1 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} -theta
srun -N 1 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} ${dargs} -Deltalvl 0

# 32 procs
nt="512"; nc="32"; ml="3"
srun -N 1 -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nt ${nt} -ml 1
srun -N 1 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nt ${nt} -ml 2 -cf0 8
srun -N 1 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} -theta
srun -N 1 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} ${dargs}

# 128 procs
nt="2048"; nc="128"; ml="4"
srun -N 1 -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nt ${nt} -ml 1
srun -N 3 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nt ${nt} -ml 3 -cf0 8
srun -N 3 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} -theta
srun -N 3 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} ${dargs}

# 512 procs
nt="8192"; nc="512"; ml="4"
srun -N 1  -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nt ${nt} -ml 1
srun -N 10 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nt ${nt} -ml 4     -cf0 8
srun -N 10 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} -cf0 16 -theta
srun -N 10 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} -cf0 16 ${dargs} 


echo 'Done'
echo -n 'Timestamp END: ';date
