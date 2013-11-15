#!/bin/sh
#BHEADER**********************************************************************
# Copyright (c) 2008,  Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# This file is part of HYPRE.  See file COPYRIGHT for details.
#
# HYPRE is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.
#
# $Revision: 1.3 $
#EHEADER**********************************************************************


# scriptname holds the script name, with the .sh removed
scriptname=`basename $0 .sh`
example_dir="../examples"
test_to_run="drive-02"

# Echo usage information
case $1 in
   -h|-help)
      cat <<EOF

   $0 [-h|-help] {src_dir}

   where: {src_dir}  is the hypre source directory
          -h|-help   prints this usage information and exits

   This script runs basic tests of warp for a constant coefficient 2D Laplacian

   Example usage: test.sh $0 

EOF
      exit
      ;;
esac

# mpi command and mpi setup for this machine
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

# Run make in the example directory
cd $example_dir
make clean
make $test_to_run

# Run the regression tests and move all output to the output_dir
# The diff command below is for comparing generated output to saved output
# In particular, we need to ignore the timing statements, i.e. -I"time",
# between the stored output files and the generated files
cd ../test
# Test 1
echo "Running Test 1"
$RunString -np 4 $example_dir/drive-02 -pgrid 1 1 4 -nt 256 -ml 15 1>> $output_dir/std.out 2>> $output_dir/std.err
diff -U3 -bI"time"  $scriptname.saved $output_dir/std.out >> $output_dir/std.err
# Test 2 ...

# Echo to stderr all nonempty error files in $output_dir.  This
# is critical, because test.sh assumes this.
for errfile in $( find $output_dir ! -size 0 -name "*.err" )
do
   cat $errfile >&2
done

# remove machinefile, if created
if [ -n $MACHINES_FILE ] ; then
   rm $MACHINES_FILE
fi
