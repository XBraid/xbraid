# Script for MPI runs, 
#
#  $ mpirun -np K  python3 ex_05_run.py
#

import ex_05
from viz_ex_05 import viz_ex_05

core, app = ex_05.InitCoreApp(a=0.0)
ex_05.run_Braid(core, app)
viz_ex_05()

