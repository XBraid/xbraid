from scipy import *
from matplotlib import pyplot as mpl
from os import sys

##
# Run like 
#  $ python viz.py  <it>  
#
#  where <it> is the XBraid iteration that you want to visualize
#
#  The output files are assumed to all begin with the value given 
#  by file_stem.  And each outputted state vector is assumed to have
#  an integer representing the total number of time steps on line 0
#  followed by a double on each line representing the 1D spatial grid
#
##

# Set the braid iteration number and number of steps
it = int(sys.argv[1])
file_stem = 'ex-burgers.out.'

# Find out size of problem in space and time
step = 0
fname = file_stem + "%07d"%step + '.' + "%05d"%it
data = loadtxt(fname)
nsteps = int(data[0])
data = zeros((nsteps+1,data.shape[0]-1))
# Load space-time solution
for step in range(nsteps+1):
    fname = file_stem + "%07d"%step + '.' + "%05d"%it
    data[step,:] = (loadtxt(fname))[1:]

mpl.figure(0)
mpl.imshow(data,origin='lower')
mpl.colorbar()
mpl.ylabel('time')
mpl.xlabel('space')

mpl.figure(1)
mpl.plot(data[0,:])
mpl.ylabel('u')
mpl.xlabel('space')
mpl.title('Step 0')
mpl.ylim(data[0,:].min()-1.0, data[0,:].max()+1.0)

mpl.figure(2)
mpl.plot(data[nsteps/2,:])
mpl.ylabel('u')
mpl.xlabel('space')
mpl.title('Step %d'%(nsteps/2))
mpl.ylim(data[nsteps/2,:].min()-1.0, data[nsteps/2,:].max()+1.0)

mpl.figure(3)
mpl.plot(data[nsteps,:])
mpl.ylabel('u')
mpl.xlabel('space')
mpl.title('Step %d'%nsteps)
mpl.ylim(data[nsteps,:].min()-1.0, data[nsteps,:].max()+1.0)


mpl.show()

