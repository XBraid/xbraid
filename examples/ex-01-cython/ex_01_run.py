# Run with 
#  $$ mpirun -np K  python3 ex_01_run.py
#

import ex_01
core, app = ex_01.InitCoreApp()
ex_01.run_Braid(core, app)

