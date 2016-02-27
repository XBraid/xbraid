from scipy import loadtxt, array, zeros, sqrt, arange, log10
import matplotlib as mpl
import matplotlib.pyplot as plt
from sys import argv, exit
import time

##
# For some reason, Python doesn't want to do the plots interactively...stupid Python
# It may be the version of Matplotlib on JBS tuxbox
interactive = False

##
# Helper routine for setting the plotting parameters just the way we want them
def return_rcparams(fig_width_pt=700.0, fig_height_pt=700.0, fontsize=28, fontfamily='sans-serif'):
    inches_per_pt = 1.0/72.27                  # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt     # width in inches
    fig_height = fig_height_pt*inches_per_pt   # height in inches
    fig_size =  [fig_width,fig_height]

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
              'figure.figsize'    : fig_size,
              'lines.linewidth'   : 2,
              'lines.markersize'  : 14,
              'xtick.major.size'  :  14, 
              'xtick.minor.size'  : 0,
              'ytick.major.size'  :  14, 
              'ytick.minor.size'  : 0}
    return params

##
# Plot the data
if __name__ == "__main__":
   
   if( len(argv) == 2 ):
      if argv[1] == '-help' or argv[1] == '--help' or argv[1] == 'help':
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
               $$ python cycleplot.py <niter>  <nlevel>
               
               Same as Usage 2, only the plot is scaled to accommodate niter iteration
               and nlevel levels.  This is nice so that your plot isn't always resized. 
               By default, niter=12, and nlevels=5.


               Notes
               ------
               If you need to change the cycling output file, then edit the script 
               to change "fname".
               
               '''
         exit()
   elif( len(argv) == 3 ):
      niter_max = int(argv[1])
      nlevel_max = int(argv[2])
   else:
      niter_max = 12
      nlevel_max = 5

   ##
   # Outer loop that halts once a residual is read that meets the tolerance
   rnorm = 1; tol = 0;
   while(rnorm > tol):
   
      ##
      #                             0       1       2      3      4     5
      # The data has the columns: level, nrefine, iter, gupper, ||r||, tol
      try:
         fname = "braid.out.cycle"
         data = loadtxt(fname)
      except:
         raise ValueError("Your cycling output file " + fname + " has not been created yet!")
      
      ##
      # For halting purposes, find the current rnorm and tol
      if interactive == True:
         rnorm = data[:,4]
         rnorm = rnorm[rnorm >= 0].min()
         tol = data[:,5]
         tol = tol[tol >= 0].min()
      else:
         rnorm = -1

      ##
      # Set plotting parameters and figure size based on anticipated number of iterations and levels
      # So long as the 
      #     number of iterations <= niter
      # and
      #     number of levels     <= nlevels
      # then the plot window shouldn't change size.  Thus, you may want to adjust the minimum sizes
      # of niter and nlevels below for your simulation
      if interactive == True:
         niter   = max(niter_max, data[:,2].max())
         nlevels = max(nlevel_max, data[:,0].max() - data[:,0].min())
      else:
         niter   = data[:,2].max()
         nlevels = data[:,0].max() - data[:,0].min()

      fontsize = 28
      params = return_rcparams(fig_width_pt=niter*110, fig_height_pt=nlevels*90., fontfamily='serif', fontsize=fontsize)
      for key,entry in params.items():
          mpl.rcParams[key] = entry


      ##
      # The 'true' level must account for any shifts in level number due to nrefine
      level = data[:,0] - data[:,1]
      
      ##
      # Plot levels, noting that level 0 is the finest and level k is a coarse level.
      # So, plot (-levels) and adjust the y-ticks accordingly
      level = -level
      plt.plot(level, '-ko')
      maxl = max(level)
      minl = maxl - nlevels 
      plt.ylim([minl-1, maxl+1])

      ##
      # yticks (on the left) are just the level
      yticks = range(int(min(level)), int(max(level)+1))
      plt.yticks(yticks, ["%d"%(-k) for k in yticks])
      plt.ylabel('Level', color='k')
      plt.xlabel('Iteration')
      
      ##
      # Do a second plot of the residual norms, if we have any residuals available 
      r = data[:,4]
      r[r == -1] = 0
      if r.max() > 0:
         dr = r[:-1] - r[1:]
         r_indices = (dr < 0).nonzero()[0]  # x-locations to print
         r_to_print = -dr[r_indices]
         r_level = level[r_indices]         # corresponding level numbers
         import pdb; pdb.set_trace()
         ax2 = plt.twinx()
         ax2.semilogy(r_indices, r_to_print,'-bo')
         plt.ylabel('$\| r_k \|$', color='b')
         
         ##
         # Plot 4 y-ticks on the right for ||r||
         tols = data[:,5]
         mi = min( tols[tols>0].min()/500., r_to_print.min())
         ma = r_to_print.max()
         yticks = [ma, 10**((log10(mi)+log10(ma))*1./3.), 10**((log10(mi)+log10(ma))*2./3.), mi ]
         plt.yticks(yticks, ["%1.1e"%tick for tick in yticks])
         for tl in ax2.get_yticklabels():
            tl.set_color('b')


      ##
      # Print xticks only when the iteration changes
      dlevel = level[:-1] - level[1:]
      dlevel[dlevel == -1] = 0
      ddlevel = dlevel[:-1] - dlevel[1:]
      new_iterations = (ddlevel == -1).nonzero()[0] + 1
      xticks = arange( level.shape[0] )[new_iterations]
      plt.xticks(xticks, ["%d"%(k+1) for k in range(xticks.shape[0])])
      plt.xlim(0, (2*nlevels+1)*niter)


      plt.title('XBraid Cycling')
      plt.gcf().subplots_adjust(bottom=0.15)
      plt.gcf().subplots_adjust(right=0.82)

      plt.show(block=False)
      time.sleep(1)
   
   ##
   # Force the user to close the plotting window 
   plt.show(block=True)
