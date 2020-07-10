# Script for MPI runs, 
#
#  $ mpirun -np K  python3 ex_05_run.py
#

import ex_05
from viz_ex_05 import viz_ex_05

# Can choose to print help message
print_help = False

core, app = ex_05.InitCoreApp(print_help=print_help, ml=15, nu=1, nu0=1,
                              CWt=1.0, skip=0, nx=256, ntime=512, eps=0.0, 
                              a=1.0, tol=1e-6, cf=2, mi=30, sc=0, fmg=0, 
                              advect_discr='fourth_diss_sq', diff_discr='second_order',
                              time_discr='SDIRK3')

if print_help == False:
    ex_05.run_Braid(core, app)

viz_ex_05()


