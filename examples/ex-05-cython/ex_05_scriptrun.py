# Script for MPI runs, 
#
#  $ mpirun -np K  python3 ex_05_run.py
#

import ex_05
from viz_ex_05 import viz_ex_05

# Can choose to print help message
print_help = False


time_discr = ['BE', 'SDIRK3']
advect_discr = ['upwind', 'central', 'fourth', 'fourth_diss', 'fourth_diss_sq']
eps = [0, 1.0]
a = [0, 1.0]

for td in time_discr:
    for ad in advect_discr:
        for ee in eps:
            for aa in a:
                print("\n\n Time, space, eps, a:  " + td + "   " + ad + "   " + str(ee) + "   " + str(aa) )
                core, app = ex_05.InitCoreApp(print_help=print_help, ml=15, nu=1, nu0=1,
                                              CWt=1.0, skip=0, nx=128, ntime=256, eps=ee, 
                                              a=aa, tol=1e-6, cf=2, mi=30, sc=0, fmg=0, 
                                              advect_discr=ad, time_discr=td)

                if print_help == False:
                    ex_05.run_Braid(core, app)
            
                viz_ex_05()


