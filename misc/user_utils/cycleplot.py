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

from scipy import loadtxt, array, zeros, sqrt, arange, log10, setdiff1d
import matplotlib as mpl
import matplotlib.pyplot as plt
from sys import argv, exit
import time

##
# Helper routine for setting the plotting parameters just the way we want them
def return_rcparams(fig_width=5, fig_height=5, fontsize=28, fontfamily='sans-serif'):

    params = {'backend': 'ps',
              'font.family'       : fontfamily,    # can specify many more font and math font
              'font.serif'        :  'cm',          #properties, see ye ole info superhighway  
              'font.sans-serif'   : 'arial',
              'axes.labelsize'    : fontsize,
              'font.size'         : fontsize,
              'axes.titlesize'    : fontsize,
              'legend.fontsize'   : fontsize-6,
              'xtick.labelsize'   : fontsize,
              'ytick.labelsize'   : fontsize,
              'text.usetex'       : True,
              #'figure.figsize'    : [fig_width,fig_height],
              'lines.linewidth'   : 2,
              'lines.markersize'  : 10,
              'xtick.major.size'  :  12, 
              'xtick.minor.size'  : 0,
              'ytick.major.size'  :  12, 
              'ytick.minor.size'  : 0}
    return params

##
# Plot the data
if __name__ == "__main__":
   
   if( (len(argv) == 2)  and (argv[1] == '-help' or argv[1] == '--help' or argv[1] == 'help')):
         print ''' 
               Cycle plotting visualization for XBraid

               Usage 1
               -------
               $$ python cycleplot.py help

               Prints this message
               
               
               Usage 2
               -------
               $$ mpirun -np # braid_driver
               $$ python cycleplot.py
               
               Simultaneously runs your braid_driver and plots how braid
               traverses the levels in your cycle. 
               
               
               Usage 3
               -------
               $$ mpirun -np # braid_driver
               $$ python cycleplot.py <interactive>
               
               Same as Usage 2, only if interactive is > 0, then the plot will refresh every
               so <interactive> number of seconds.


               Usage 4
               -------
               $$ mpirun -np # braid_driver
               $$ python cycleplot.py <interactive> <niter>  <nlevel>
               
               Same as Usage 3, only the plot is scaled to accommodate niter iterations
               and nlevel levels.  This is nice so that your plot isn't always resized. 
               By default, niter=12, and nlevels=5.


               Notes
               ------
               If you need to change the cycling output file, then edit the script 
               to change "fname".
               
               '''
         exit()
   elif( len(argv) == 2 ):
      interactive = float(argv[1])
      niter_max = 12
      nlevel_max = 5
   elif( len(argv) == 4 ):
      interactive = float(argv[1])
      niter_max = int(argv[2])
      nlevel_max = int(argv[3])
   else:
      interactive = 0
      niter_max = 12
      nlevel_max = 5

   ##
   # Outer loop that halts once a residual is read that meets the tolerance
   rnorm = 1; tol = 0; counter = -1; current_lvl = -1
   while(rnorm > tol or current_lvl != 0 ):
      counter += 1

      ##
      #                             0       1       2      3      4     5
      # The data has the columns: level, nrefine, iter, gupper, ||r||, tol
      try:
         fname = "braid.out.cycle"
         data = loadtxt(fname)
         a = data[:,4].min()
      except:
         raise ValueError("Your cycling output file " + fname + " has not been created yet!")
      
      ##
      # For halting purposes, find the current rnorm, tol and current_lvl 
      if interactive > 0:
         rnorm = data[:,4]
         rnorm = rnorm[rnorm >= 0].min()
         tol = data[:,5]
         tol = tol[tol >= 0].min()
         current_lvl = data[-1,0]
      else:
         current_lvl = 0
         rnorm = -1

      ##
      # Set plotting parameters and figure size based on anticipated number of iterations and levels
      # So long as the 
      #     number of iterations <= niter
      # and
      #     number of levels     <= nlevels
      # then the plot window shouldn't change size.  Thus, you may want to adjust the minimum sizes
      # of niter and nlevels below for your simulation
      if interactive > 0:
         niter   = max(niter_max, data[:,2].max())
         nlevels = max(nlevel_max, data[:,0].max() - data[:,0].min())
      else:
         niter   = data[:,2].max()
         nlevels = data[:,0].max() - data[:,0].min()

      fontsize = 22
      fig_size = [niter*1.7, nlevels*1.0]
      if counter == 0:
         fig, ax = plt.subplots(1,1, figsize=fig_size)
      ## For some reason, these functions don't work
      #fig.set_figwidth(fig_size[0])
      #fig.set_figheight(fig_size[1])
      #fig.tight_layout()
      params = return_rcparams(fig_width=fig_size[0], fig_height=fig_size[1], fontfamily='serif', fontsize=fontsize)
      for key,entry in params.items():
          mpl.rcParams[key] = entry


      ##
      # The 'true' level must account for any shifts in level number due to nrefine
      level = data[:,0] - data[:,1]
      
      ##
      # Plot levels, noting that level 0 is the finest and level k is a coarse level.
      # So, plot (-levels) and adjust the y-ticks accordingly
      level = -level
      maxl = max(level)
      minl = maxl - nlevels
      ax.set_ylim([minl-1, maxl+1])
      
      ##
      # Print xticks only when the iteration changes
      new_iterations = (level == maxl).nonzero()[0]
      new_iterations = new_iterations[ ((new_iterations[1:] - new_iterations[:-1]) != 1).nonzero()[0] ]
      xticks = arange( level.shape[0] )[new_iterations]
      ax.set_xticks(xticks)
      ax.set_xticklabels(["%d"%(k+1) for k in range(xticks.shape[0])])
      ax.set_xlim(0, max((2*nlevels+2)*niter, len(level)) )
      ax.plot(level, '-ko')

      ##
      # yticks (on the left) are just the level
      yticks = range(int(min(level)), int(max(level)+1))
      ax.set_yticks(yticks, ["%d"%(-k) for k in yticks])
      ax.set_ylabel('Level', color='k')
      ax.set_xlabel('Iteration')
      ax.figure.canvas.draw()
      
      ##
      # Do a second plot of the residual norms, if we have any residuals available 
      r = data[:,4]
      r[r == -1] = 0
      if r.max() > 0:
         dr = r[:-1] - r[1:]
         r_indices = setdiff1d((dr != 0).nonzero()[0], (r==0).nonzero()[0]-1)  # x-locations to print
         r_to_print = r[r_indices+1]
         r_level = level[r_indices]         # corresponding level numbers
         ax2 = ax.twinx()
         ax2.set_xlim(0, max((2*nlevels+2)*niter, len(level)) )
         ax2.semilogy(r_indices, r_to_print,'-bo')
         ax2.set_ylabel('$\| r_k \|$', color='b')
         
         ##
         # Plot 4 y-ticks on the right for ||r||
         tols = data[:,5]
         mi = min( tols[tols>0].min()/500., r_to_print.min())
         ma = r_to_print.max()
         yticks = [ma, 10**((log10(mi)+log10(ma))*1./3.), 10**((log10(mi)+log10(ma))*2./3.), mi ]
         ax2.set_yticks(yticks, ["%1.1e"%tick for tick in yticks])
         for tl in ax2.get_yticklabels():
            tl.set_color('b')
         ##
         ax2.figure.canvas.draw()


      if counter == 0:
         plt.title('XBraid Cycling')
         plt.show(block=False)
      
      time.sleep(interactive)
   
   ##
   # Force the user to close the plotting window 
   plt.show(block=True)
