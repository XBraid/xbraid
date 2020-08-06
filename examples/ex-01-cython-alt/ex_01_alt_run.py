# Script for MPI runs, 
#
#  $ mpirun -np K  python3 ex_01_alt_run.py
#


import ex_01_alt
core, app = ex_01_alt.InitCoreApp()
ex_01_alt.run_Braid(core, app)

