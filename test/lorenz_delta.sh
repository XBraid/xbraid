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

   This script runs basic tests of Braid for a constant coefficient 2D heat 
   equation. The output is written to $scriptname.out, $scriptname.err and 
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
      RunString="mpirun -machinefile $MACHINES_FILE $*"
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


# compile the regression test drivers 
echo "Compiling regression test drivers"
cd $example_dir
make clean
make ex-07
cd $test_dir

fixed_args="-tstop 0.5 -ntime 256 -tol 1e-4 "
# Run the following regression tests 
TESTS=( "$RunString -np 1 $example_dir/ex-07 -test" \
        "$RunString -np 2 $example_dir/ex-07 $fixed_args -ml 2 -cf 2 -rank 0" \
        "$RunString -np 2 $example_dir/ex-07 $fixed_args -ml 2 -cf 2 -rank 1" \
        "$RunString -np 2 $example_dir/ex-07 $fixed_args -ml 2 -cf 2 -rank 2" \
        "$RunString -np 2 $example_dir/ex-07 $fixed_args -ml 2 -cf 2 -rank 3" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 2 -cf 4 -rank 1" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 2 -cf 4 -rank 2" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 2 -cf 4 -rank 3" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 0" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 1" \
        "$RunString -np 2 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 2" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 2" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 3" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 1 -defer-lvl 1" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 2 -defer-lvl 1" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 3 -defer-lvl 1" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 1 -defer-iter 1 -mi 3" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 2 -defer-iter 1 -mi 3" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 3 -defer-iter 1 -mi 3" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 1 -nu 1" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 2 -nu 1 -fcfLyap -noLyap" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 3 -noLyap" \
        "$RunString -np 4 $example_dir/ex-07 $fixed_args -ml 3 -cf 2 -rank 3 -noLyap -fcfLyap" )

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
lines_to_check=".*residual norm.*|.*Passed.*|.*passed.*|.*Average Lyapunov Exponents.*"
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
   rm $MACHINES_FILE 2> /dev/null
fi

rm braid.out.cycle 2> /dev/null
rm ex-02*.out.* 2> /dev/null
