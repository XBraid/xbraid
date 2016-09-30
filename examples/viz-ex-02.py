from scipy import *
from matplotlib import pyplot as mpl
from os import sys

##
# Run like 
#  $ python viz.py
#
#  Vizualize the files from ex-burgers.c.  The output is assumed to have 
#  the following format.  The filenames are 
#       file_stem + '.' step_number + '.' + rank
#  So, the filename
#       ex-burgers.out.0000350.00000
#  is the output from step 350 from processor rank 0.  Right now,
#  there is no parallelism, so it's always rank 0.
#
#  Each output file has the format
#     
#     ntime_steps
#     tstart
#     tstop
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
rank = 0 #int(sys.argv[1])
file_stem = 'ex-burgers.out.'

# Find out size of problem in space and time, and the 
# grid spacings
step = 0
fname = file_stem + "%07d"%step + '.' + "%05d"%rank
data = loadtxt(fname)
nsteps = int(data[0])
tstart = float(data[1])
tstop = float(data[2])
nspace = int(data[3])
xstart = float(data[4])
xstop = float(data[5])
mesh = linspace(xstart, xstop, nspace)
data = zeros((nsteps,data.shape[0]-6))

# Load space-time solution
for step in range(nsteps):
    fname = file_stem + "%07d"%step + '.' + "%05d"%rank
    data[step,:] = (loadtxt(fname))[6:]

mpl.figure(0)
mpl.imshow(data,origin='lower',extent=(xstart, xstop, tstart, tstop))
mpl.colorbar()
mpl.ylabel('time')
mpl.xlabel('space')

mpl.figure(1)
mpl.plot(mesh, data[0,:], '-o')
mpl.ylabel('u')
mpl.xlabel('space')
mpl.title('Step 0')
mpl.ylim(data[0,:].min()-1.0, data[0,:].max()+1.0)

mpl.figure(2)
mpl.plot(mesh, data[nsteps/2,:], '-o')
mpl.ylabel('u')
mpl.xlabel('space')
mpl.title('Step %d'%(nsteps/2))
mpl.ylim(data[nsteps/2,:].min()-1.0, data[nsteps/2,:].max()+1.0)

mpl.figure(3)
mpl.plot(mesh, data[nsteps-1,:], '-o')
mpl.ylabel('u')
mpl.xlabel('space')
mpl.title('Step %d'%(nsteps-1))
mpl.ylim(data[nsteps-1,:].min()-1.0, data[nsteps-1,:].max()+1.0)


mpl.show()

