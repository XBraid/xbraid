#!/bin/sh
#BHEADER**********************************************************************
# Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# This file is part of WARP.  See file COPYRIGHT for details.
#                                                                          
# WARP is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.                    
#
#EHEADER**********************************************************************


# scriptname holds the script name, with the .sh removed
scriptname=`basename $0 .sh`
example_dir="../examples"
test_to_run="drive-02"


# Echo usage information
case $1 in
   -h|-help)
      cat <<EOF

   $0 [-h|-help] 

   where: -h|-help   prints this usage information and exits

   This script runs basic tests of warp for a constant coefficient 2D heat equation. 
   The output is written to diffusion2D.out, diffusion2D.err and diffusion2D.dir.
   This test passes if diffusion2D.err is empty.

   Example usage: test.sh $0 

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
output_dir=`pwd`/$scriptname.dir
rm -fr $output_dir
mkdir -p $output_dir


# Split $scriptname.saved into $scriptname.saved.0, $scriptname.saved.1, ...,
# where each file is slice from $scriptname.saved representing a single test.
# Each test is delimited by the line "# Begin Test .."  
TestDelimiter='# Begin Test'
csplit -n 1 --silent --prefix $output_dir/$scriptname.saved. $scriptname.saved "%$TestDelimiter%" "/$TestDelimiter.*/" {*}


# Run the regression tests and move all output to the output_dir
#
# The tests are dumped to $output_dir/std.out.0, $output_dir/std.out.1, ...
# Each std.out.num file is diff-ed with the corresponding .saved.num file with 
# the result dumped to std.err.num. 
lines_to_check="^  time steps|^  number of levels|^  coarsening factor|^  num F-C relaxations|^  num rels on level 0|^  iterations|^spatial problem size|  expl"

# Test 0
echo "Running Test 0"
$RunString -np 4 $example_dir/drive-02 -pgrid 1 1 4 -nt 256 -ml 15 1>> $output_dir/std.out.0  2>> $output_dir/std.err.0
egrep "$lines_to_check" $output_dir/std.out.0 > temp.out
mv temp.out $output_dir/std.out.0
diff -U3 -B -bI"$TestDelimiter" $output_dir/$scriptname.saved.0 $output_dir/std.out.0 >> $output_dir/std.err.0
#
# Test 1
echo "Running Test 1"
$RunString -np 4 $example_dir/drive-02 -pgrid 1 1 4 -nt 256 -ml 15 -fmg 1>> $output_dir/std.out.1  2>> $output_dir/std.err.1
egrep "$lines_to_check" $output_dir/std.out.1 > temp.out
mv temp.out $output_dir/std.out.1
diff -U3 -B -bI"$TestDelimiter" $output_dir/$scriptname.saved.1 $output_dir/std.out.1 >> $output_dir/std.err.1
#
# Test 2
echo "Running Test 2"
$RunString -np 8 $example_dir/drive-05  -pgrid 1 1 8 -ml 15 -nt 128 -nx 33 33 -mi 100 -expl -scoarsen >> $output_dir/std.out.2  2>> $output_dir/std.err.2
egrep "$lines_to_check" $output_dir/std.out.2 > temp.out
mv temp.out $output_dir/std.out.2
diff -U3 -B -bI"$TestDelimiter" $output_dir/$scriptname.saved.2 $output_dir/std.out.2 >> $output_dir/std.err.2
#
# Test ...


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
