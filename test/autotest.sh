#!/bin/sh
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


# Turn on Shell Script Debugging
#set -xv

# Setup
test_dir=`pwd`
source_dir=`cd ..; pwd`
finished_dir="autotest_finished"
#remote_dir="/usr/casc/hypre/braid/testing/"
remote_dir="/home/jbschroder/joint_repos/braid_test_history/"
remote_subdir="AUTOTEST-`date +%Y.%m.%d-%a`"

# Email Setup
summary_file="Braid_SUMMARY.html"
summary_subject="Braid Autotest Summary `date +%Y-%m-%d`"
email_list="jbschroder@unm.edu, rfalgout@llnl.gov"
#email_list="jbschroder@unm.edu"
timebraid_logo="______           _     _ 
| ___ \         (_)   | |
| |_/ /_ __ __ _ _  __| |
| ___ \ '__/ _\` | |/ _\` |
| |_/ / | | (_| | | (_| |
\____/|_|  \__,_|_|\__,_| "

# Print Help Message
case $1 in
   -h|-help)
      cat <<EOF

   $0  [-init | -{test} ... | -remote-copy | -summary-email | -create-tarball]

   where:

      -init               Initializes for autotests and does a fresh pull and 
                          recompile of the repository. Should be called before 
                          running tests.
      -{test}             Run a single indicated test associated with a specific
                          machine name (e.g., -tux343, -vulcan).
      -remote-copy        Copy the autotest results to the remote archive in
                          /home/jbschroder/joint_repos/braid_test_history".  
                          A second argument can be passed in this instance 
                          that tells autotest to tunnel through a machine during 
                          the copy.
      -summary-email      Sends out a summary email to the developers for all the 
                          tests run today.  This command cannot be run over ssh, 
                          and must be run from a machine with direct access to 
                          /home/jbschroder/joint_repos/braid_test_history
      -create-tarball     Creates a nice share-able tarball of braid, with 
                          unnecessary files removed, fresh documentation and a
                          VERSION file to identify this checkout

   with options:

      -h|-help       prints this usage information and exits
      -t|-trace      echo each command

   The main purpose of this script is to organize the automatic testing process
   and to ensure that all related files have the appropriate permissions.

   Example usage: $0 -init
                  $0 -tux149
                  $0 -remote-copy 'tux343'
                  $0 -summary-email
                  $0 -create-tarball

EOF
      exit
      ;;
   -t|-trace)
      set -xv
      shift
      ;;
   *)
      break
      ;;
esac


# Begin main switch statement over $1
case $1 in
   
   -init)
      # Initialize autotest error and out files
      if [ -e autotest.err ]; then
         rm autotest.err
         touch autotest.err
      fi
      if [ -e autotest.out ]; then
         rm autotest.out
         touch autotest.out
      fi
      
      # Update this copy of the repo
      (
         rm -rf $finished_dir
         mkdir $finished_dir
         cd $source_dir
         git pull
         make clean
         make all
      ) 1>> $test_dir/autotest.out 2>> $test_dir/autotest.err

      break
      ;;

   # Generate a summary based on the remote directory
   -summary-email)
      
      # move all but the last 10 autotest results into yearly subdirectories
      cd $remote_dir
      files=`echo AUTOTEST-2*.*`
      count=`echo $files | wc | awk '{print $2}'`
      for i in $files
      do
         if [ $count -le 10 ]; then
            break;
         fi
         dir=`echo $i | awk -F '.' '{print $1}'`
         if [ ! -d $dir ]; then
            mkdir $dir
            chmod -fR a+rX,ug+w,o-w $dir
         fi
         mv $i $dir/$i
         count=`expr $count - 1`
      done
      
      # Create Summary File
      cd $remote_subdir
      echo "To: $email_list"           >  $summary_file
      echo "From: Jacob Schroder <jbschroder@unm.edu>" >>  $summary_file
      echo "Subject: $summary_subject" >> $summary_file
      echo "Content-Type: text/html"   >> $summary_file
      echo "MIME-Version: 1.0"         >> $summary_file
      echo " "                         >> $summary_file
      echo "<html>"                    >> $summary_file
      echo "<head> </head>"            >> $summary_file
      echo "<PRE>"                     >> $summary_file
      
      # echo logo
      echo "$timebraid_logo" >> $summary_file
      echo  " "     >> $summary_file
      echo $summary_subject            >> $summary_file
      echo ""         >> $summary_file

      # all top-level tests with empty error files are reported as "passed",
      # not including the cron autotest logs
      echo "[PASSED]" >> $summary_file
      for test in $( find . -maxdepth 1 -size 0 -name "*.err" ! -name "*cron*" )
      do
         testname=`basename $test .err`
         echo "-${testname#machine-}" >> $summary_file
      done

      # all top-level tests with non-empty error files are reported as "failed",
      # including the cron autotest logs
      echo " "     >> $summary_file
      echo "[FAILED]" >> $summary_file
      for test in $( find . -maxdepth 1 ! -size 0 -name "*.err" )
      do
         testname=`basename $test .err`
         for prefix in "machine-" "autotest-";
         do
            testname="${testname#$prefix}"
         done
         echo "-$testname" >> $summary_file
      done

      # keep a time stamp of last runs and report if more than 10 days
      echo " "     >> $summary_file
      echo "[DETAILS]" >> $summary_file
      egrep "began at|ended at" *.out >> $summary_file

      # list all non-empty error files in today's output directory
      echo "  "     >> $summary_file
      echo "[ERROR FILES]" >> $summary_file
      output_dir="$remote_dir$remote_subdir" 
      for test in $( find $output_dir ! -size 0 -name "*.err" | sort -r )
      do
         echo "<a href=\"file://$test\">$test</a>" >> $summary_file
      done

      echo "</PRE>"  >> $summary_file
      echo "</html>" >> $summary_file

      # send the email
      /usr/sbin/sendmail -t < $summary_file

      break
      ;;

    -remote-copy)
       
       # Copy all finished tests to the remote directory, in a subdirectory 
       # named by this year, month and day
       (
          # Create new autotest directory, if necessary over ssh
          if [ -n "$2" ]; then
            ssh $2 "mkdir -p $remote_dir$remote_subdir"
          else
            mkdir -p $remote_dir$remote_subdir
          fi
          
          # Copy results, if necessary over ssh
          cd $finished_dir
          if [ -n "$2" ]; then
             scp -q -r *  "$2:$remote_dir$remote_subdir"
             ssh $2 chmod -fR a+rX,ug+w,o-w "$remote_dir$remote_subdir"
          else
             scp -q -r *  $remote_dir$remote_subdir
             chmod -fR a+rX,ug+w,o-w "$remote_dir$remote_subdir"
          fi
       ) 1>> $test_dir/autotest.out 2>> $test_dir/autotest.err

       break
       ;;
   
    -create-tarball)
        # Create a nice share-able tarball
        (
            create_release.sh
        ) 1>> $test_dir/autotest.out 2>> $test_dir/autotest.err

        break
        ;;
 
   *)
      (
         # A machine name is the command line parameter
         case $1 in
            -tux[0-9]*)
               name="tux"
               ;;

            -mac)
               name="mac"
               ;;

               *)
               echo "Unknown machine option: $1" >&2
               exit 1
               ;;
         esac
         
         # Run machine test 
         echo "Test [machine-$name] began at  `date +%T` on `date +%D`"
         echo "Test [machine-$name] began at  `date +%T` on `date +%D`" 1> machine-$name.out
         ./machine-$name.sh   1>> machine-$name.out   2>> machine-$name.err
         echo "Test [machine-$name] ended at  `date +%T` on `date +%D`"
         echo "Test [machine-$name] ended at  `date +%T` on `date +%D`"  1>> machine-$name.out

         # Move all results to $finished_dir
         mv machine-$name.dir $finished_dir
         mv machine-$name.err $finished_dir
         mv machine-$name.out $finished_dir
      ) 1>> $test_dir/autotest.out 2>> $test_dir/autotest.err

      ;;
esac

