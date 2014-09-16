#!/bin/sh

##
# To test that the archive works, do something like
# $ mkdir braid_temp
# $ tar -xvf braid-2014.09.16.tar.gz -C ./braid_temp
# $ cd braid_temp/braid
# $ gvim examples/Makefile
# 
# Change,
#  MFEM_DIR = ../../../mfem
#  METIS_DIR = ../../../metis-4.0

#  HYPRE_DIR = ../../../linear_solvers/hypre
#
# Then, make and do autotests

##
# Set these variables to what you want
this_version="1.0 beta"
destination_dir="/usr/casc/hypre/braid/share/braid_tarballs"
braid_version_to_checkout='HEAD'
temp_dir="$HOME/braid_temp"
date=`date +%Y.%m.%d`
archive_name="braid-`date +%Y.%m.%d`.tar.gz"

touch $destination_dir/$archive_name

(
   ##
   # Create export directory
   rm -rf $temp_dir
   mkdir $temp_dir
   cd $temp_dir
   git clone /usr/casc/hypre/braid/git/braid braid
   cd braid
   git reset --hard $braid_version_to_checkout

   ##
   # Dump a version file into the new directory
   cd braid
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
   # clean up examples
   cd examples
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
   cd ../

   ##
   # Remove all source tarballs older than 10 days
   find $destination_dir/braid-*.tar.gz -mtime +10 -exec rm {} \;

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

# Check size of archive greater than 10MB
ArchiveSize=$(du -k $destination_dir/$archive_name | cut -f 1)
if [ $ArchiveSize -le 10000 ] ; then 
    echo "" 1>&2
    echo "Tried creating Braid source archive, but " 1>&2
    echo "$destination_dir/$archive_name is too small" 1>&2
else
    echo ""
    echo "Successfully created Braid source archive, "
    echo "$destination_dir/$archive_name is large enough"
fi 
