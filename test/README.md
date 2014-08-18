## Regression Testing
<!--
 - Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 - Produced at the Lawrence Livermore National Laboratory.
 - This file is part of XBraid.  See file COPYRIGHT for details.
 -
 - XBraid is free software; you can redistribute it and/or modify it under the
 - terms of the GNU Lesser General Public License (as published by the Free
 - Software Foundation) version 2.1 dated February 1999.
 -->

### Overview

-  There are three levels in the testing framework.  At each level, the fine-grain output
   from a `testscript.sh` is dumped into a directory `testscript.dir`, with the 
   standard out and error stored in `testscript.out` and `testscript.err`.  The test
   `testscript.sh` passes if `testscript.err` is empty (nothing is written to standard
   error).

-  Basic instructions:  run a test with a command like
   
         $ ./test.sh diffusion2D.sh

   Then, see if `diffusion2D.err` is of size 0.  If it is not, look at it's contents
   to see which test failed.

-  To add a new regression test, create a new lowest level script like
   `diffusion2D.sh` and then call it from a machine script at level 2.

-  Regression tests should be run before pushing code.  It is recommended to run the 
   basic (lowest level) tests like `diffusion2d.sh` or machine test like `machine-tux.sh`


### Lowest Level Test Scripts 

As an example, here we look at one of the lowest level tests, the diffusion2d test.\n

Files used:
-  `test.sh`
-  `diffusion2D.sh`
-  `diffusion2D.saved`

Output:
-  `diffusion2D.dir`
-  `diffusion2D.err`
-  `diffusion2D.out`

At this level, we execute 
   
      $ ./test.sh diffusion2D.sh

or just
   
      $ ./diffusion2D.sh

The script `diffusion2D.sh` must create `diffusion2D.dir` and place all
fine-grain test output in this directory.  `test.sh` captures the standard out
and error in `diffusion2D.out` and `diffusion2D.err`.  The test `diffusion2D.sh` passes
if `diffusion2D.err` is empty (nothing is written to standard error).

The strategy for low level scripts like `diffusion2D.sh` is to run a sequence of 
test drivers such as
      
      $ mpirun -np 1 ../examples/drive-02 -pgrid 1 1 1 -nt 256
      $ mpirun -np 4 ../examples/drive-02 -pgrid 1 1 4 -nt 256

The output from the first mpirun test must then be written to files named

      diffusion2D.dir/unfiltered.std.out.0
      diffusion2D.dir/std.out.0
      diffusion2D.dir/std.err.0

and the second mpirun test similarly writes the files

      diffusion2D.dir/unfiltered.std.out.1
      diffusion2D.dir/std.out.1
      diffusion2D.dir/std.err.1

Subsequent tests are written to higher numbered files.  The 
`unfiltered.std.out.num` file contains all of the standard out for the 
test, while `std.out.num` contains filtered output (usually from a `grep`
command) and could contain the output lines such as iteration numbers
and number of levels.  The file `std.err.num` contains the standard error 
output.

To see if a test ran correctly, `std.out.num` is compared to saved output in
`diffusion2D.saved`.  The file `diffusion2D.saved` contains the concatenated
output from all the tests that `diffusion2D.sh` will run.  For the above
example, this file could look like 

      # Begin Test 1
      number of levels     = 6
      iterations           = 16
      # Begin Test 2
      number of levels     = 4
      iterations           = 8

This saved output is split into an individual file for each test (using `#
Begin Test` as a delimiter) and these new files are placed in
`diffusion2D.dir`.  So, after running these two regression tests,
`diffusion2D.dir` will contain

      diffusion2D.saved.0
      diffusion2D.saved.1
      unfiltered.std.out.0
      std.out.0
      std.err.0
      unfiltered.std.out.1
      std.out.1
      std.err.1

An individual test has passed if `std.err.num` is empty.  The file
`std.err.num` contains a diff between `diffusion2D.save.num` and `std.out.num`
(the diff ignores whitespace and the delimiter `# Begin Test`).

Last in the directy where you ran `./test.sh diffusion2d.sh`, the files 
      
      diffusion2D.err
      diffusion2D.out

will be created.  If all the tests passed then `diffusion2D.err` will be empty.
Otherwise, it will contain the filenames of the `std.err.num` files that are 
non-empty, representing failed tests. 


### Level 2 Scripts

As an example, here we look at one of the Level 2 tests, the machine-tux test that Jacob
runs. \n

Files used:
-  `machine-tux.sh`

Output:
-  `machine-tux.dir`
-  `machine-tux.err` (only generated if `autotest.sh` is used to run `machine-tux.sh`)
-  `machine-tux.out` (only generated if `autotest.sh` is used to run `machine-tux.sh`)


At this level, we execute 
   
      ./machine-tux.sh

The autotest framework (`autotest.sh`) calls machine scripts in this way.  Each
machine script should be short and call lower-level scripts like
`diffusion2D.sh`.  The output from lower-level scripts must be moved to
`machine-tux.dir` like this: 

      $ ./test.sh diffusion2D.sh 
      $ mv -f diffusion2D.dir machine-tux.dir
      $ mv -f diffusion2D.out machine-tux.dir
      $ mv -f diffusion2D.err machine-tux.dir

All error files from `diffusion2D.sh` will be placed in `machine-tux.dir`, so
if `machine-tux.dir` has all zero `*.err` files, then the `machine-tux` test
has passed.

To begin testing on a new machine, like vulcan, add a new machine script
similar to `machine-tux.sh` and change `autotest.sh` to recognize and run the new
machine test.  To then use `autotest.sh` with the machine script, you'll have to
set up a passwordless connection from the new machine to

      /usr/casc/hypre/braid/testing


### Level 3 Script

Here we look at the highest level, where `autotest.sh` runs all of the level 2 
machine tests and emails out the results. \n 

Files used:
-  `autotest.sh`

Output:
-  `test/autotest_finished`
-  `/usr/casc/hypre/braid/testing/AUTOTEST-20**.**.**-Day` 
- Email to recipients listed in `autotest.sh`

At the highest level sits `autotest.sh` and is called automatically as a cronjob.
If you just want to check to see if you've broken anything with a commit, just use lower
level scripts.

There are four steps to running autotest.

-  Step 1

         $ ./autotesh.sh -init

   will do a pull from master for the current working repository and recompile 
   Braid.  The autotest output files (`autotest.err` and `autotest.out`) and the output 
   directory (autotest_finished) are initialized.

-  Step 2
      
         $ ./autotest.sh -tux343

   will run the autotests on tux343.  This command will look for a 
   `machine-tux.sh`, and execute it, moving the resulting 
   
         machine-tux.dir
         machine-tux.err
         machine-tux.out
      
   into `test/autotest_finished`.

-  Step 3

         $ ./autotest.sh -remote-copy

   will copy `/test/autotest_finished/*` to a time-stamped directory
   such as\n  `/usr/casc/hypre/braid/testing/AUTOTEST-2013.11.18-Mon` \n

   Alternatively,

         $ ./autotesh.sh -remote-copy tux343
   
   will ssh through tux343 to copy to `/usr/casc`.
   Multiple machines may independently be running regression tests and then
   copy to `AUTOTEST-2013.11.18-Mon`.

-  Step 4

         $ ./autotest.sh -summary-email

   will email everyone listed in the `$email_list` (an `autotest.sh` variable) 



### Cronfile 

To add entries to your crontab, First, put your new cronjob lines into 
`cronfile`.  Then see what you already have in your crontab file with

      $ crontab -l

Next, append to cronfile whatever you already have 

      $ crontab -l >> cronfile

Finally, tell crontab to use your cronfile

      $ crontab cronfile

Then make sure it took affect with
      
      $ crontab -l

Crontab entry format uses '*' to mean "every" and '*/m' to 
mean "every m-th". The first five entries on each line 
correspond respectively to:
-   minute (0-56)
-   hour (0-23)
-   day of month (1-31)
-   month (1-12)
-   day of week (0-6)(0=Sunday)

Jacob's crontab (on tux343):

      05 14 * * * source /etc/profile; source $HOME/.bashrc; cd $HOME/joint_repos/braid/test; ./autotest.sh -init
      05 14 * * * source /etc/profile; source $HOME/.bashrc; cd $HOME/joint_repos/braid/test; ./autotest.sh -tux343
      05 14 * * * source /etc/profile; source $HOME/.bashrc; cd $HOME/joint_repos/braid/test; ./autotest.sh -remote-copy
      05 14 * * * source /etc/profile; source $HOME/.bashrc; cd $HOME/joint_repos/braid/test; ./autotest.sh -summary-email

