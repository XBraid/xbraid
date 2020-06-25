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
TESTS=( "diffusion1D.sh " \
        "diffusion1D_scaling.sh" \
        "diffusion1D_check_rnorm.sh" \
        "diffusion2D.sh " \
        "diffusion2D_storage.sh " \
        "diffusion2D_scaling.sh " \
        "diffusion2D_scaling_storage.sh " \
        "diffusion2D_check_rnorm.sh " \
        "diffusion2D_check_rnorm_storage.sh " \
        "compare_examples_drivers.sh " \
        "compare_examples_drivers_storage.sh " \
        "mfem.sh" \
        "test-checkout-compile.sh " \
        "adjoint.sh " \
        "shellvector_bdf2.sh "\
        "memcheck-tux-jacob.sh "\
        "ode1D.sh"\
        "ode1D_cython.sh")
#       Need to fix the issues with refinement = 2 
#        "ode1D-refine-periodic.sh"\

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

