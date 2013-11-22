## Regression Testing

-  There are three levels in the testing framework.  At each level, the output
   from a `testscript.sh` is dumped into a directory `testscript.dir` with the
   standard out and error stored in `testscript.out` and `testscript.err`.
   `testscript.sh` passes if `testscript.err` is empty (nothing is written to standard
   error).

-  Basic instructions:
   
         $ ./test.sh diffusion2D.sh

   Then, see if `diffusion2D.err` is of size 0.

-  To add a new regression test:
   Create a new low level script like `diffusion2D.sh` and then call it from a machine
   script at level 2.

-  Regression tests should be run before pushing code.


### Lowest Level Test Scripts 

Files used:
-  `test.sh`
-  `diffusion2D.sh`
-  `diffusion2D.saved`
-  `diffusion2D.filters`

Output:
-  `diffusion2D.dir`
-  `diffusion2D.err`
-  `diffusion2D.out`

At this level, we execute 
   
      $ ./test.sh diffusion2D.sh

or just
   
      $ ./diffusion2D.sh

The script `diffusion2D.sh` must create `diffusion2D.dir` and place all fine-grain
test output in this directory.  `test.sh` captures the standard out and error in
`diffusion2D.out` and `diffusion2D.err`.  If you expect any known error messages in
`diffusion2D.err`, you can filter these with regular expressions in
`diffusion2D.filters` (`test.sh` processes the filters file). `diffusion2D.sh` 
passes if `diffusion2D.err` is empty (nothing is written to standard error).

The strategy for low level scripts like diffusion2D.sh is to run commands like

      $ mpirun -np 4 ../examples/drive-02 -pgrid 1 1 4 -nt 256 -ml 15 
   
and then compare the output to saved output in `diffusion2D.saved`.  The file
`diffusion2D.saved` contains the concatenated output from all the tests that
`diffusion2D.sh` will run, i.e., the file looks like

      # Begin Test 1
      ... actual test 1 output ...
      # Begin Test 2
      ... actual test 2 output ...
      # Begin Test 3
      ... 

The directory `diffusion2D.dir` will contain files like

      diffusion2D.saved.0
      diffusion2D.saved.1
      std.err.0
      std.err.1
      std.out.0
      std.out.1

where `diffusion2D.saved.num` is just the expected output for test num and is a
slice from `diffusion2D.saved`.  `std.out.num` is the output from the first test,
and `std.err.num` is the diff between `std.out.num` and `diffusion2D.saved.num`.


### Level 2 Scripts

Files used:
-  `machine-tux.sh`

Output:
-  `machine-tux.dir`
-  `machine-tux.err`
-  `machine-tux.out`

The autotest framework (`autotest.sh`) calls machine scripts, like
`machine-tux.sh`.  Each machine script should be short and call low level
scripts like `diffusion2D.sh`.  The output must then be moved to
`machine-tux.dir`: 

      $ ./test.sh diffusion2D.sh 
      $ mv -f diffusion2D.dir machine-tux.dir
      $ mv -f diffusion2D.out machine-tux.dir
      $ mv -f diffusion2D.err machine-tux.dir

Any nonzero error files from `diffusion2D.sh` will be listed in `machine-tux.err`,
with the detailed output from `diffusion2D.sh` stored in `machine-tux.dir`.  If 
`machine-tux.err` is zero then the test has passed.

To begin testing on a new machine, like vulcan, add a new machine script
similar to `machine-tux.sh` and change `autotest.sh` to recognize and run the new
machine.  To then use `autotest.sh` with the machine script, you'll have to
set up a passwordless connection from the new machine to

      /usr/casc/hypre/warp/testing


### Level 3 Script

Files used:
-  `autotest.sh`

Output:
-  `test/autotest_finished`
-  `/usr/casc/hypre/warp/testing/AUTOTEST-20**.**.**-Day` 

At the highest level sits `autotest.sh` and is called automatically as a cronjob.
If you just want to check to see if you've broken anything with a commit, just use lower
level scripts.

There are four steps to running autotest.

-  Step 1

         $ ./autotesh.sh -init

   will do a pull from master for the current working repository and recompile 
   warp.  The autotest output files (`autotest.err` and `autotest.out`) and the output 
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
   such as  `/usr/casc/hypre/warp/testing/AUTOTEST-2013.11.18-Mon`

         $ ./autotesh.sh -remote-copy tux343
   
   will ssh through tux343 to copy to `/usr/casc`.
   Multiple machines may independently be running regression tests and then
   copy to `AUTOTEST-2013.11.18-Mon`.

-  Step 4

         $ ./autotest.sh -summary-email

   will email everyone listed in in `$email_list` in `autotesh.sh`.



### Cronfile 

To add entries to your crontab, put your new cronjob lines in cronfile. 
Then see what you already have in your crontab file with

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

      05 14 * * * source /etc/profile; source $HOME/.bashrc; cd $HOME/joint_repos/warp/test; ./autotest.sh -init
      05 14 * * * source /etc/profile; source $HOME/.bashrc; cd $HOME/joint_repos/warp/test; ./autotest.sh -tux343
      05 14 * * * source /etc/profile; source $HOME/.bashrc; cd $HOME/joint_repos/warp/test; ./autotest.sh -remote-copy
      05 14 * * * source /etc/profile; source $HOME/.bashrc; cd $HOME/joint_repos/warp/test; ./autotest.sh -summary-email

