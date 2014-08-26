#!/bin/sh
#BHEADER**********************************************************************
# Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# This file is part of XBraid.  See file COPYRIGHT for details.
#                                                                          
# XBraid is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.                    
#
#EHEADER**********************************************************************

testname=`basename $0 .sh`

# Echo usage information
case $1 in
   -h|-help)
      cat <<EOF

   **** Only run this script on one of the tux machines. ****

   $0 [-h|-help] 

   where: -h|-help   prints this usage information and exits

   This script runs a number of Braid tests suitable for the tux machines.

   Example usage: $0 ..


EOF
      exit
      ;;
esac

# Setup
test_dir=`pwd`
output_dir=`pwd`/$testname.dir
rm -fr $output_dir
mkdir -p $output_dir

# Run the following regression tests 
TESTS=( "diffusion2D.sh " \
        "diffusion2D_scaling.sh " \
        "diffusion2D_check_rnorm.sh " \
        "machine-tux-checkout-compile.sh " \
        "docs.sh" \
        "memcheck-tux-jacob.sh ")

# Run regression tests
for test in "${TESTS[@]}"
do
   thistest=`basename $test .sh`
   eval "./test.sh $test"
   mv -f $thistest.dir $output_dir
   mv -f $thistest.out $output_dir
   mv -f $thistest.err $output_dir
done 

# Echo to stderr all nonempty error files in $output_dir
for errfile in $( find $output_dir ! -size 0 -name "*.err" )
do
   echo $errfile >&2
done

