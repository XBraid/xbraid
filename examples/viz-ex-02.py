from __future__ import division
from scipy import *
from matplotlib import pyplot as mpl
from os import sys

##
# Run like 
#  $ python viz.py
#
#  Vizualize the files from ex-02.c.  The output is assumed to have 
#  the following format.  The filenames are 
#       file_stem + '.' step_number + '.' + rank
#  So, the filename
#       ex-02.out.0000350.00000
#  is the output from step 350 from processor rank 0.  Right now,
#  there is no parallelism, so it's always rank 0.
#
#  Each output file has the format
#     
#     ntime_steps
#     tstart
#     tstop
#     current_time
#     nspace_points
#     xstart
#     xstop
#     x[0]
#     x[1]
#       .
#       .
#       .
#     x[k]
##

# Set the braid iteration number and number of steps
#file_stem = 'ex-02-serial.out.'
file_stem = 'ex-02.out.'
current_rank = 0

# Load the initial solution file and extract the problem size in space and time and the grid spacings
step = 0
fname = file_stem + "%07d"%step + '.' + "%05d"%current_rank
data = loadtxt(fname)
nsteps = int(data[0])
tstart = float(data[1])
tstop = float(data[2])
tvalues = zeros((nsteps,))
nspace = int(data[4])
xstart = float(data[5])
xstop = float(data[6])
mesh = linspace(xstart, xstop, nspace)
data = zeros((nsteps,data.shape[0]-7))

# Load space-time solution, noting that we don't know the MPI ranks ahead of
# time that generated the files, so we guess :-)
for step in range(nsteps):
    try:
        fname = file_stem + "%07d"%step + '.' + "%05d"%current_rank
        data[step,:] = (loadtxt(fname))[7:]
        tvalues[step] = (loadtxt(fname))[3]
    except:
        current_rank = current_rank + 1
        fname = file_stem + "%07d"%step + '.' + "%05d"%current_rank
        data[step,:] = (loadtxt(fname))[7:]
        tvalues[step+1] = (loadtxt(fname))[3]


# IMSHOW not as useful for adaptive grids, so it's turned off by default
#mpl.figure(0)
#mpl.imshow(data,origin='lower',extent=(xstart, xstop, tstart, tstop))
#mpl.colorbar()
#mpl.ylabel('time')
#mpl.xlabel('space')


for j, idx in enumerate([0, nsteps//5, 2*(nsteps//5), 3*(nsteps//5), 4*(nsteps//5), nsteps-1]):
    mpl.figure(j+1)
    mpl.plot(mesh, data[idx,:], '-o')
    mpl.ylabel('u')
    mpl.xlabel('space')
    mpl.title('Step=%d,  Time=%1.4e'%(idx, tvalues[idx]) )
    mpl.ylim(data[idx,:].min()-1.0, data[idx,:].max()+1.0)
    mpl.xlim(xstart, xstop) 

mpl.figure(0)
mpl.plot(tvalues, zeros_like(tvalues), 'k.')
mpl.title("Time Grid")

mpl.show()

