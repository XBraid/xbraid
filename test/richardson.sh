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
rm -fr $output_dir 2> /dev/null
mkdir -p $output_dir


# compile the regression test drivers 
# note that there are a lot of unavoidable Fortran warnings about 
# unused app structures, hence we ignore those for ex-01b-f
echo "Compiling regression test drivers"
cd $example_dir
make clean
make ex-02 
make ex-06 
cd $test_dir


# Run the following regression tests 
TESTS=( "$RunString -np 1 $example_dir/ex-06 -nt 19" \
        "$RunString -np 3 $example_dir/ex-06 -nt 19" \
        "$RunString -np 1 $example_dir/ex-06 -nt 19 -time_discr 2 -mi 2" \
        "$RunString -np 3 $example_dir/ex-06 -nt 19 -time_discr 2 -mi 2" \
        "$RunString -np 1 $example_dir/ex-06 -nt 19 -richardson" \
        "$RunString -np 3 $example_dir/ex-06 -nt 19 -richardson" \
        "$RunString -np 1 $example_dir/ex-06 -nt 19 -time_discr 2 -richardson -mi 2" \
        "$RunString -np 3 $example_dir/ex-06 -nt 19 -time_discr 2 -richardson -mi 2" \
        "$RunString -np 1 $example_dir/ex-06 -nt 19 -refinet -max_tpts 400 -mi 9" \
        "$RunString -np 3 $example_dir/ex-06 -nt 19 -refinet -max_tpts 400 -mi 9" \
        "$RunString -np 1 $example_dir/ex-06 -nt 19 -time_discr 2 -refinet -max_tpts 100 -mi 4" \
        "$RunString -np 3 $example_dir/ex-06 -nt 19 -time_discr 2 -refinet -max_tpts 100 -mi 4" \
        "$RunString -np 1 $example_dir/ex-06 -nt 19 -refinet -max_tpts 100 -richardson -mi 6" \
        "$RunString -np 3 $example_dir/ex-06 -nt 19 -refinet -max_tpts 100 -richardson -mi 6" \
        "$RunString -np 1 $example_dir/ex-06 -nt 19 -time_discr 2 -refinet -max_tpts 50 -richardson -mi 5" \
        "$RunString -np 3 $example_dir/ex-06 -nt 19 -time_discr 2 -refinet -max_tpts 50 -richardson -mi 5" \
        "$RunString -np 1 $example_dir/ex-02 -ntime 17 -richardson -mi 4" \
        "$RunString -np 3 $example_dir/ex-02 -ntime 17 -richardson -mi 4" \
        "$RunString -np 1 $example_dir/ex-02 -ntime 17 -refinet 1e-2" \
        "$RunString -np 3 $example_dir/ex-02 -ntime 17 -richardson -CWt 1.2" \
        "$RunString -np 1 $example_dir/ex-06 -nt 37 -richardson -periodic -mi 5" \
        "$RunString -np 3 $example_dir/ex-06 -nt 37 -richardson -periodic -mi 5" )

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
lines_to_check="^  time steps.*|^  number of levels.*|^  max number of levels.*|^  iterations.*|^  residual norm.*|^Finished braid_TestAll: no fails detected, however some results must be|.*Braid: Temporal refinement occurred.*|^  num_syncs.*|^  periodic.*"
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
   rm ex-01*.out.* timegrid.* 2> /dev/null
   echo "Running Test $counter"
   eval "$test" 1>> $output_dir/unfiltered.std.out.$counter  2>> $output_dir/std.out.$counter
   cd $output_dir
   egrep -o "$lines_to_check" unfiltered.std.out.$counter > std.out.$counter
   diff -U3 -B -bI"$TestDelimiter" $scriptname.saved.$counter std.out.$counter >> std.err.$counter
   # check whether test was successfull (portable on UNIX systems)
   if [ `du std.err.$counter | cut -f1` -gt 0 ]; then
      echo "...did not pass."
   fi
   cd $test_dir
   cat ex-01*.out.* > $output_dir/solutionvector.out.$counter 2> /dev/null
   cat timegrid.* > $output_dir/timegrid.$counter 2> /dev/null
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
rm braid_timings.* 2> /dev/null
rm ex-06*.out.* 2> /dev/null
rm ex-02*.out.* 2> /dev/null
rm timegrid.* 2> /dev/null
