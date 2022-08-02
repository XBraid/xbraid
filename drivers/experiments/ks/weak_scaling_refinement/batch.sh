#!/bin/bash

##### These lines are for SLURM
#SBATCH -N 76
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
ncores="8 64 512 4096"

# timegrid size
ntimes="128 512 2048 8192"

# levels
mlevels="2 3 4 5"

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
# nx="128"; nt="128"; nc="8"; ml="2"
# srun -N 1 -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml 1
# srun -N 1 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml ${ml} -cf0 2
# srun -N 1 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml ${ml} -theta
# srun -N 1 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml ${ml} ${dargs} -Deltalvl 0

# 64 procs
# nx="256"; nt="512"; nc="64"; ml="3"
# srun -N 1 -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml 1
# srun -N 2 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml 2 -cf0 8
# srun -N 2 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml ${ml} -theta
# srun -N 2 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml ${ml} ${dargs}

# 512 procs
# nx="512"; nt="2048"; nc="512"; ml="4"
# srun -N 1  -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml 1
# srun -N 10 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml 3 -cf0 8
# srun -N 10 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml ${ml} -theta
# srun -N 10 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml ${ml} ${dargs}

# 4096 procs
nx="1024"; nt="8192"; nc="4096"; ml="5"
srun -N 1  -n 1     -o ${outd}/${outn}_nc${nc}_ml1           ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml 1
srun -N 76 -n ${nc} -o ${outd}/${outn}_nc${nc}_ml${ml}       ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml 4 -cf0 8
srun -N 76 -n ${nc} -o ${outd}/${outn}_theta_nc${nc}_ml${ml} ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml ${ml} -theta
srun -N 76 -n ${nc} -o ${outd}/${outn}_Delta_nc${nc}_ml${ml} ${ex} ${fargs} -nx ${nx} -nt ${nt} -ml ${ml} ${dargs}

echo 'Done'
echo -n 'Timestamp END: ';date
