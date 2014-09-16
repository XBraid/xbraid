## Compiling and running the examples
<!--
  - Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
  - Produced at the Lawrence Livermore National Laboratory. Written by 
  - Jacob Schroder schroder2@llnl.gov, Rob Falgout falgout2@llnl.gov,
  - Tzanio Kolev kolev1@llnl.gov, Ulrike Yang yang11@llnl.gov, 
  - Veselin Dobrev dobrev1@llnl.gov, et al. 
  - LLNL-CODE-660355. All rights reserved.
  - 
  - This file is part of XBraid. Email schroder2@llnl.gov on how to download. 
  - 
  - This program is free software; you can redistribute it and/or modify it under
  - the terms of the GNU General Public License (as published by the Free Software
  - Foundation) version 2.1 dated February 1999.
  - 
  - This program is distributed in the hope that it will be useful, but WITHOUT ANY
  - WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
  - PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
  - License for more details.
  - 
  - You should have received a copy of the GNU Lesser General Public License along
  - with this program; if not, write to the Free Software Foundation, Inc., 59
  - Temple Place, Suite 330, Boston, MA 02111-1307 USA
 -->

Type

      drive-0* -help

for instructions on how to run any driver.

To run the examples, type
   
      mpirun -np 4 drive-*  [args]


1. drive-01 is the simplest example.  It implements a scalar ODE and can be
  compiled and run with no outside dependencies.

2. drive-02 implements the 2D heat equation on a regular grid.  You must have
   [hypre](https://computation-rnd.llnl.gov/linear_solvers/software.php)
   installed and these variables in examples/Makefile set correctly
    
          HYPRE_DIR = ../../linear_solvers/hypre
          HYPRE_FLAGS = -I$(HYPRE_DIR)/include
          HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE


3. drive-03 implements the 3D heat equation on a regular grid, and assumes 
   [hypre](https://computation-rnd.llnl.gov/linear_solvers/software.php) 
   is installed just like drive-02.

4. drive-05 implements the 2D heat equation on a regular grid, but it
   uses spatial coarsening.  This allows you to use explicit time stepping
   on each Braid level, regardless of time step size.  It assumes 
   [hypre](https://computation-rnd.llnl.gov/linear_solvers/software.php) 
   is installed just like drive-02.

5. drive-04 is a sophisticated test bed for various PDEs, mostly parabolic.  It relies
   on the [mfem](https://code.google.com/p/mfem/)
   package to create general finite element discretizations for the spatial problem.
   Other packages must be installed in this order.
     + Unpack and install [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
     + Unpack and install 
       [hypre](https://computation-rnd.llnl.gov/linear_solvers/software.php) 
     + Unpack and install 
       [mfem.](https://code.google.com/p/mfem/)
       Make the serial version of mfem first by only typing ``make``.
       Then make sure to set these variables correctly in the mfem Makefile:
       
             USE_METIS_5 = YES
             HYPRE_DIR  = where_ever_linear_solvers_is/hypre 
   
     + Make [GLVIS](https://code.google.com/p/glvis/), which needs serial mfem.
       Set these variables in the glvis makefile
            
            MFEM_DIR   = mfem_location
            MFEM_LIB   = -L$(MFEM_DIR) -lmfem

     + Go back to the mfem directory and type
     
            make clean
            make parallel

     + Go to braid/examples and set these Makefile variables, 

            METIS_DIR = ../../metis-5.1.0/lib
            MFEM_DIR = ../../mfem
            MFEM_FLAGS = -I$(MFEM_DIR)
            MFEM_LIB = -L$(MFEM_DIR) -lmfem -L$(METIS_DIR) -lmetis
       
       then type
            
            make drive-04

     + To run drive-04 and glvis, open two windows.  In one, start a glvis session
      
               ./glvis
  
         Then, in the other window, run drive-04
      
               mpirun -np ... drive-04 [args]
         
         Glvis will listen on a port to which drive-04 will dump visualization information.


