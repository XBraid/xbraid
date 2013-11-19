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

testname=`basename $0 .sh`

# Echo usage information
case $1 in
   -h|-help)
      cat <<EOF

   **** Only run this script on one of the tux machines. ****

   $0 [-h|-help] 

   where: -h|-help   prints this usage information and exits

   This script runs a number of warp tests suitable for the tux machines.

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

# Test 1
./test.sh diffusion2D.sh 
mv -f diffusion2D.dir $output_dir
mv -f diffusion2D.out $output_dir
mv -f diffusion2D.err $output_dir

# Test 2
# ...

# Echo to stderr all nonempty error files in $output_dir
for errfile in $( find $output_dir ! -size 0 -name "*.err" )
do
   echo $errfile >&2
done

