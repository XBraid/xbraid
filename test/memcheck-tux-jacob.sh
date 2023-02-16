#!/bin/bash
#BHEADER**********************************************************************
#
# Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. Written by 
# Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
# Dobrev, et al. LLNL-CODE-660355. All rights reserved.
# 
# This file is part of XBraid. For support, post issues to the XBraid Github page.
# 
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2.1 dated February 1999.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with this program; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#EHEADER**********************************************************************

# scriptname holds the script name, with the .sh removed
scriptname=`basename $0 .sh`

# Echo usage information
case $1 in
   -h|-help)
      cat <<EOF

   $0 [-h|-help] 

   where: -h|-help   prints this usage information and exits

   This script runs basic memory tests of Braid for a constant coefficient 2D heat 
   equation. The basic strategy is to grep the output of valgrind for any instances of 
   *tw*. The output is written to $scriptname.out, $scriptname.err and 
   $scriptname.dir. This test passes if $scriptname.err is empty.

   Example usage: ./test.sh $0 

EOF
      exit
      ;;
esac

# Determine csplit and mpirun command for this machine 
OS=`uname`
case $OS in
   Linux*) 
      MACHINES_FILE="hostname"
      if [ ! -f $MACHINES_FILE ] ; then
         hostname > $MACHINES_FILE
      fi
      RunString="LD_PRELOAD=/home/jbschroder/.local/lib/valgrind/libmpiwrap-amd64-linux.so  mpirun -machinefile $MACHINES_FILE $*"
      csplitcommand="csplit"
      ;;
   Darwin*)
      csplitcommand="gcsplit"
      RunString="mpirun --hostfile ~/.machinefile_mac"
      ;;
   *)
      RunString="mpirun"
      csplitcommand="csplit"
      ;;
esac


# Setup
example_dir="../examples"
driver_dir="../drivers"
test_dir=`pwd`
output_dir=`pwd`/$scriptname.dir
rm -fr $output_dir
mkdir -p $output_dir


# compile the regression test drivers against valgrind 
echo "Compiling regression test drivers"
cd $example_dir
make clean
make valgrind=yes
cd $driver_dir
make clean
make valgrind=yes
cd $test_dir

#store -2, -1, 0, 1

# Run the following regression tests
valgrind_opts="--log-fd=1 --leak-check=full"
TESTS=( "MPIWRAP_DEBUG=quiet $RunString -np 1 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 1 1 1 -nt 30 -ml 15 -skip 0 -store -2" \
        "MPIWRAP_DEBUG=quiet $RunString -np 4 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 1 -nt 30 -ml 15 -skip 1 -store -1" \
        "MPIWRAP_DEBUG=quiet $RunString -np 4 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 1 -nt 30 -ml 15 -skip 1 -store 0" \
        "MPIWRAP_DEBUG=quiet $RunString -np 4 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 1 -nt 30 -ml 15 -skip 0 -store 1" \
        "MPIWRAP_DEBUG=quiet $RunString -np 8 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 2 -nt 30 -ml 15 -skip 0 -cf0 1 -store -2" \
        "MPIWRAP_DEBUG=quiet $RunString -np 1 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 1 1 1 -nt 32 -ml 15 -scoarsen 1 -skip 0 -store -2" \
        "MPIWRAP_DEBUG=quiet $RunString -np 4 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 1 -nt 32 -ml 15 -scoarsen 1 -skip 1 -store -1" \
        "MPIWRAP_DEBUG=quiet $RunString -np 4 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 1 -nt 32 -ml 15 -scoarsen 1 -skip 1 -store 0" \
        "MPIWRAP_DEBUG=quiet $RunString -np 4 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 1 -nt 32 -ml 15 -scoarsen 1 -skip 1 -store 1" \
        "MPIWRAP_DEBUG=quiet $RunString -np 2 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 1 1 2 -nt 32 -cf 4 -ml 15 -fmg 1 -skip 0 -store -2" \
        "MPIWRAP_DEBUG=quiet $RunString -np 8 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 2 -nt 33 -cf 4 -ml 15 -fmg 1 -skip 1 -store -1" \
        "MPIWRAP_DEBUG=quiet $RunString -np 8 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 2 -nt 33 -cf 4 -ml 15 -fmg 1 -skip 0 -store 0" \
        "MPIWRAP_DEBUG=quiet $RunString -np 8 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 2 -nt 33 -cf 4 -ml 15 -fmg 1 -skip 0 -store 1" \
        "MPIWRAP_DEBUG=quiet $RunString -np 2 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 1 1 2 -nt 33 -cf 4 -ml 15 -fmg 1 -scoarsen 2 -skip 0 -store -2" \
        "MPIWRAP_DEBUG=quiet $RunString -np 8 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 2 -nt 33 -cf 4 -ml 15 -fmg 1 -scoarsen 2 -skip 1 -store -1" \
        "MPIWRAP_DEBUG=quiet $RunString -np 8 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 2 -nt 33 -cf 4 -ml 15 -fmg 1 -scoarsen 2 -skip 0 -store 0" \
        "MPIWRAP_DEBUG=quiet $RunString -np 8 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 2 2 2 -nt 33 -cf 4 -ml 15 -fmg 1 -scoarsen 2 -skip 1 -store 1" \
        "MPIWRAP_DEBUG=quiet $RunString -np 2 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 1 1 2 -nt 33 -cf 4 -ml 15 -fmg 1 -access_level 2 -skip 0" )
        #"MPIWRAP_DEBUG=quiet $RunString -np 4 valgrind $valgrind_opts $driver_dir/drive-diffusion-2D -pgrid 1 1 4 -nt 32 -cf 4 -ml 15 -fmg 1 " )

# The below commands will then dump each of the tests to the output files 
#   $output_dir/unfiltered.std.out.0, 
#   $output_dir/std.out.0, 
#   $output_dir/std.err.0,
#    
#   $output_dir/unfiltered.std.out.1,
#   $output_dir/std.out.1, 
#   $output_dir/std.err.1,
#   ...
#
# The unfiltered output is the direct output of the script, whereas std.out.*
# is filtered by a grep for the lines that are to be checked.  
#
lines_to_check=" _tw.* | tw.* |^  time steps.*|^  number of levels.*|^  iterations.*|^ Fine level spatial problem size.*"
#
# Then, each std.out.num is compared against stored correct output in 
# $scriptname.saved.num, which is generated by splitting $scriptname.saved
#
TestDelimiter='# Begin Test'
$csplitcommand -n 1 --silent --prefix $output_dir/$scriptname.saved. $scriptname.saved "%$TestDelimiter%" "/$TestDelimiter.*/" {*}
#
# The result of that diff is appended to std.err.num. 

# Run regression tests
counter=0
for test in "${TESTS[@]}"
do
   echo "Running Test $counter"
   eval "$test" 1>> $output_dir/unfiltered.std.out.$counter  2>> $output_dir/std.out.$counter
   cd $output_dir
   egrep -o "$lines_to_check" unfiltered.std.out.$counter > std.out.$counter
   diff -U3 -B -bI"$TestDelimiter" $scriptname.saved.$counter std.out.$counter >> std.err.$counter
   cd $test_dir
   counter=$(( $counter + 1 ))
done 

# Echo to stderr all nonempty error files in $output_dir.  test.sh
# collects these file names and puts them in the error report
for errfile in $( find $output_dir ! -size 0 -name "*.err.*" )
do
   echo $errfile >&2
done


# remove machinefile, if created
if [ -n $MACHINES_FILE ] ; then
   rm $MACHINES_FILE 2> /dev/null
fi

rm braid.out.cycle 2> /dev/null
rm braid_timings.* 2> /dev/null
