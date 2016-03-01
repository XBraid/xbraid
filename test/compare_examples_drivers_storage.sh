#!/bin/sh
#BHEADER**********************************************************************
#
# Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. Written by 
# Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
# Dobrev, et al. LLNL-CODE-660355. All rights reserved.
# 
# This file is part of XBraid. Email xbraid-support@llnl.gov for support.
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

   This script runs basic tests of Braid for a constant coefficient 2D heat 
   equation. The output is written to $scriptname.out, $scriptname.err and 
   $scriptname.dir. This test passes if $scriptname.err is empty.

   Example usage: ./test.sh $0 

EOF
      exit
      ;;
esac

# Determine mpi command and mpi setup for this machine
HOST=`hostname`
case $HOST in
   tux*) 
      MACHINES_FILE="hostname"
      if [ ! -f $MACHINES_FILE ] ; then
         hostname > $MACHINES_FILE
      fi
      RunString="mpirun -mca btl ^openib -machinefile $MACHINES_FILE $*"
      ;;
      *) 
         RunString="mpirun"
         ;;
esac


# Setup
example_dir="../examples"
driver_dir="../drivers"
test_dir=`pwd`
output_dir=`pwd`/$scriptname.dir
rm -fr $output_dir
mkdir -p $output_dir


# compile the regression test drivers 
echo "Compiling regression test drivers"
cd $example_dir
make clean
make 
cd $driver_dir
make clean
make 
cd $test_dir

# Run the following regression tests -- here we run a comparison of ex-02 and
# drive-02 The first three tests compare drive-02, ex-02 and ex-02-serial,
# making sure that all three versions have the same discretization error when
# run to a very small tolerance.  The subsequent tests check the discretization
# error after only 2 iterations, so that we can test various XBraid options.
# This way the printed discretization error is printed before convergence and
# hence `unique' to this particular option combination. 
TESTS=( "$RunString -np 8 $driver_dir/drive-02 -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage -1 -skip 0" \
        "$RunString -np 8 $example_dir/ex-02   -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage -1 -skip 0" 
        "$RunString -np 8 $driver_dir/drive-02 -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage -1 -skip 1" \
        "$RunString -np 8 $example_dir/ex-02   -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage -1 -skip 1" \
        "$RunString -np 8 $driver_dir/drive-02 -pgrid 1 1 8 -nt 128 -nx 17 17 -fmg 1 -mi 4 -storage -1 -skip 0" \
        "$RunString -np 8 $example_dir/ex-02   -pgrid 1 1 8 -nt 128 -nx 17 17 -fmg   -mi 4 -storage -1 -skip 0" \
        "$RunString -np 8 $driver_dir/drive-02 -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage 0 -skip 0" \
        "$RunString -np 8 $example_dir/ex-02   -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage 0 -skip 0" \
        "$RunString -np 8 $driver_dir/drive-02 -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage 0 -skip 1" \
        "$RunString -np 8 $example_dir/ex-02   -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage 0 -skip 1" \
        "$RunString -np 8 $driver_dir/drive-02 -pgrid 1 1 8 -nt 128 -nx 17 17 -fmg 1 -mi 4 -storage 0 -skip 0" \
        "$RunString -np 8 $example_dir/ex-02   -pgrid 1 1 8 -nt 128 -nx 17 17 -fmg   -mi 4 -storage 0 -skip 0" \
        "$RunString -np 8 $driver_dir/drive-02 -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage 1 -skip 0" \
        "$RunString -np 8 $example_dir/ex-02   -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage 1 -skip 0" \
        "$RunString -np 8 $driver_dir/drive-02 -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage 1 -skip 1" \
        "$RunString -np 8 $example_dir/ex-02   -pgrid 1 1 8 -nt 128 -nx 17 17 -mi 4 -storage 1 -skip 1" \
        "$RunString -np 8 $driver_dir/drive-02 -pgrid 1 1 8 -nt 128 -nx 17 17 -fmg 1 -mi 4 -storage 1 -skip 0" \
        "$RunString -np 8 $example_dir/ex-02   -pgrid 1 1 8 -nt 128 -nx 17 17 -fmg   -mi 4 -storage 1 -skip 0" )

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
lines_to_check="^  Discr. error at final time.*|.*discr. error at final time.*"
#
# Then, each std.out.num is compared against stored correct output in 
# $scriptname.saved.num, which is generated by splitting $scriptname.saved
#
TestDelimiter='# Begin Test'
csplit -n 1 --silent --prefix $output_dir/$scriptname.saved. $scriptname.saved "%$TestDelimiter%" "/$TestDelimiter.*/" {*}
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


# Additional tests can go here comparing the output from individual tests,
# e.g., two different std.out.* files from identical runs with different
# processor layouts could be identical ...


# Echo to stderr all nonempty error files in $output_dir.  test.sh
# collects these file names and puts them in the error report
for errfile in $( find $output_dir ! -size 0 -name "*.err.*" )
do
   echo $errfile >&2
done


# remove machinefile, if created
if [ -n $MACHINES_FILE ] ; then
   rm $MACHINES_FILE
fi
