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

##
# Set these variables to what you want
this_version="2.0.0"
destination_dir="/usr/casc/hypre/braid/share/braid_tarballs"
docs_destination_dir="/usr/casc/hypre/braid/share/braid_manuals/"
braid_version_to_checkout='HEAD'
temp_dir="$HOME/braid_temp"
date="`date +%Y.%m.%d` at `date +%H.%M`"

##
# Create archive name
hash=$(git rev-parse $braid_version_to_checkout)
hash_short=${hash:0:8}
archive_name="braid_`date +%Y-%m-%d`_$hash_short.tar.gz"
touch $destination_dir/$archive_name

(
   ##
   # Create export directory
   rm -rf $temp_dir
   mkdir $temp_dir
   cd $temp_dir
   git clone ssh://git@mystash.llnl.gov:7999/xbraid/xbraid.git braid
   cd braid
   git reset --hard $braid_version_to_checkout

   ##
   # Dump a version file into the new directory
   head -n 22 COPYRIGHT > VERSION
   echo "XBraid Version: " >> VERSION
   echo $this_version >> VERSION
   echo " " >> VERSION
   echo "XBraid Git Hash: " >> VERSION
   git rev-parse $braid_version_to_checkout >> VERSION
   echo " " >> VERSION
   echo "Date: " >> VERSION
   echo $date >> VERSION

   ## 
   #clean up the export directory
   rm -rf .git
   make clean

   ##
   # clean up drivers
   cd drivers
   rm -rf hyperdrive
   cd ../

   ##
   # clean up test
   cd test
   rm autotest.err
   rm autotest.out
   rm -rf autotest_finished
   cd ../

   ##
   # clean up docs
   cd docs
   make clean
   make user_manual
   make developer_manual
   rm -rf developer_manual
   rm -rf user_manual
   rm *.tex~
   rm -rf img
   rm *.conf
   rm *.tex
   rm *.sty
   rm Makefile
   cp user_manual.pdf $docs_destination_dir
   cp developer_manual.pdf $docs_destination_dir
   rm user_manual.pdf
   rm developer_manual.pdf
   cd ../

   ##
   # Remove all other tar balls for this hash, but leave an empty one for todays tarball
   cd $destination_dir
   rm *"$hash_short"*
   touch $archive_name
   # Remove all but the 10 most recent tarballs and go back to the directory where this code block started
   rm `ls -t *.tar.gz | awk 'NR>10'`
   cd $temp_dir/braid

   ##
   # zip it up! 
   cd ../ 
   tar -cf braid.tar braid
   gzip braid.tar
   mv braid.tar.gz $destination_dir/$archive_name
   chgrp hypre $destination_dir/$archive_name
   chmod g+rw $destination_dir/$archive_name
   rm -rf $temp_dir
) 1>> /dev/null 2>> /dev/null

# Check size of archive greater than 5MB
ArchiveSize=$(du -k $destination_dir/$archive_name | cut -f 1)
if [ $ArchiveSize -le 200 ] ; then 
    echo "" 1>&2
    echo "Tried creating Braid source archive, but " 1>&2
    echo "$destination_dir/$archive_name is too small" 1>&2
else
    echo ""
    echo "Successfully created Braid source archive, "
    echo "$destination_dir/$archive_name is large enough"
fi 
