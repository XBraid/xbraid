#!/bin/bash

##### These lines are for SLURM
#SBATCH -N 40
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
ncores="8 32 128 512 2048"

# timegrid size
ntimes="128 512 2048 8192 32768"

# levels
mlevels="2 2 3 4 5"

# fixed arguments
fargs="-tf 4 -nx 512 -cf 4 -cf0 8 -mi 50"
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
nt="512"; nc="32"; ml="2"
srun -N 1 -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nt ${nt} -ml 1
srun -N 1 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nt ${nt} -ml ${ml} 
srun -N 1 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} -theta
srun -N 1 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} ${dargs}

# 128 procs
nt="2048"; nc="128"; ml="3"
srun -N 1 -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nt ${nt} -ml 1
srun -N 3 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nt ${nt} -ml ${ml} 
srun -N 3 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} -theta
srun -N 3 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} ${dargs}

# 512 procs
nt="8192"; nc="512"; ml="4"
srun -N 1  -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nt ${nt} -ml 1
srun -N 10 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nt ${nt} -ml ${ml} 
srun -N 10 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} -theta
srun -N 10 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} ${dargs} 

# 2048 procs
nt="32768"; nc="2048"; ml="5"
srun -N 1  -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nt ${nt} -ml 1
srun -N 40 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nt ${nt} -ml ${ml} 
srun -N 40 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} -theta
srun -N 40 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nt ${nt} -ml ${ml} ${dargs} 

echo 'Done'
echo -n 'Timestamp END: ';date
