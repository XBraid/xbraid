## Compiling the examples

- Type

   drive-0* -help

for instructions on how to run any driver

-  To run the examples
   
    mpirun -np 4 drive-*  [args]


- drive-01 is the simplest example.  It implements a scalar ODE and can be compiled and run
with no outside dependencies.

- drive-02 implements the 2D heat equation on a regular grid.  You must have
[hypre](https://computation-rnd.llnl.gov/linear_solvers/software.php)
installed and these variables in the Makefile set correctly
    
    HYPRE_DIR = ../../linear_solvers/hypre
    HYPRE_FLAGS = -I$(HYPRE_DIR)/include
    HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE


- drive-03 implements the 3D heat equation on a regular grid, and assumes 
hypre is installed just like drive-02.

- drive-05 implements the 2D heat equation on a regular grid, but it
uses with spatial coarsening.  This allows you to use explicit time stepping
on each Warp level, regardless of time step size.  It assumes 
hypre is installed just like drive-02.

- drive-04 is sophisticated test bed for various PDEs, mostly parabolic.  It relies
on the 
[mfem](https://code.google.com/p/mfem/)
package to create general finite element discretizations for the spatial problem.
Other packages must be installed in this order.
  + Unpack and install [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
  + Unpack and install hypre
  + Unpack and install mfem.  Make the serial version of mfem first by only typing ``make``.
    Then make sure to set these variables correctly in the Makefile
       
       USE_METIS_5 = YES
       HYPRE_DIR  = where_ever_linear_solvers_is/hypre 
   
  + Make [GLVIS](https://code.google.com/p/glvis/), which needs serial mfem.
  + Go back to mfem 
  
         make clean
         make parallel

  + Got to warp/examples and set these Makefile variables, and them ``make drive-04``.

         METIS_DIR = ../../metis-5.1.0/lib
         MFEM_DIR = ../../mfem
         MFEM_FLAGS = -I$(MFEM_DIR)
         MFEM_LIB = -L$(MFEM_DIR) -lmfem -L$(METIS_DIR) -lmetis

- To run drive-04 and glvis, open two windows.  In one, start a glvis session
      
         ./glvis
  
  Then, in the other window, run drive-04
      
         mpirun -np ... drive-04 [args]



