from __future__ import print_function
from scipy import *
from matplotlib import pyplot as mpl
from os import sys

##
# Run like 
#  $ python viz-ex-06.py serial|braid
#
#  Vizualize the files from ex-06.c 
#       - ex-04.out.***.###    [ output from ex-04, *** is the step number, and ### is the proc number ]
# 
#  The files ex-04.out.***.###  This file contains six lines, 
#
#  <float> time-value      
#  <int>   time index
#  <int>   global num time points
#  <float> numerical solution value
#  <float> Richardson-based error estimate 
#  <float> True error
#
##

if len(sys.argv) != 1:
    print(" Please read this script for information on how to run ") 
    sys.exit(1)

# Set the braid iteration number and number of steps
file_stem = 'ex-06.out.'
current_rank = 0

# Load the initial solution file and extract the problem size
tstart = 0
tstop = 0.5
step = 0
fname = file_stem + "%04d"%step + '.' + "%03d"%current_rank
data = loadtxt(fname)
nsteps = int(data[2])

# Allocate arrays to plot
soln = zeros((nsteps,))
tpts = zeros((nsteps,))
err_est = zeros((nsteps,))
err = zeros((nsteps,))


# Load data, noting that we don't know the MPI ranks ahead of
# time that generated the files, so we guess :-)
for step in range(nsteps):
    try:
        fname = file_stem + "%04d"%step + '.' + "%03d"%current_rank
        data = loadtxt(fname)
    except:
        current_rank = current_rank + 1
        fname = file_stem + "%04d"%step + '.' + "%03d"%current_rank
        data = loadtxt(fname)
    ##
    tpts[step] = data[0]
    soln[step] = data[3]
    err_est[step] = data[4]
    err[step] = data[5]


# Visualize results
mpl.figure(0)
mpl.plot(tpts, soln, '-k+') 
mpl.ylabel('Numerical Solution')
mpl.xlabel('time')
mpl.title('Example 6')

mpl.figure(1)
mpl.semilogy(tpts, err, '-k+')
if(err_est[0] != -1.0):
    mpl.semilogy(tpts, err_est, '-r+') 
    mpl.legend(['True Error', 'Richardson Error Estimate'])
else:
    mpl.legend(['True Error', 'Richardson Error Estimate'])

mpl.ylabel('Error')
mpl.xlabel('time')
mpl.title('Example 6 Error')

mpl.show()

