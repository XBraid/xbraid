## Drivers: compiling and running
<!--
  - Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
  - Produced at the Lawrence Livermore National Laboratory. Written by 
  - Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
  - Dobrev, et al. LLNL-CODE-660355. All rights reserved.
  - 
  - This file is part of XBraid. Email xbraid-support@llnl.gov for support.
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


1. drive-02 implements the 2D heat equation on a regular grid.  You must have
   [hypre](https://computation.llnl.gov/project/linear_solvers/software.php)
   installed and these variables in examples/Makefile set correctly
    
          HYPRE_DIR = ../../linear_solvers/hypre
          HYPRE_FLAGS = -I$(HYPRE_DIR)/include
          HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE
   
   This driver also support spatial coarsening and explicit time stepping.
   This allows you to use explicit time stepping on each Braid level, 
   regardless of time step size.  

2. drive-03 implements the 3D heat equation on a regular grid, and assumes 
   [hypre](https://computation.llnl.gov/project/linear_solvers/software.php)
   is installed just like drive-02.  This driver does not support spatial 
   coarsening, and thus if explicit time stepping is used, the time stepping 
   switchs to implicit on coarse XBraid grids when the CFL condition is violated.

3. drive-04 is a sophisticated test bed for various PDEs, mostly parabolic.  It relies
   on the [mfem](http://mfem.org)
   package to create general finite element discretizations for the spatial problem.
   Other packages must be installed in this order.
     + Unpack and install [Metis](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
     + Unpack and install 
       [hypre](https://computation.llnl.gov/project/linear_solvers/software.php)
     + Unpack [mfem.](http://mfem.org)
       Then make sure to set these variables correctly in the mfem Makefile:
       
             USE_METIS_5 = YES
             HYPRE_DIR   = where_ever_linear_solvers_is/hypre 
   
     + Make the parallel version of mfem first by typing

            make parallel

     + Make [GLVIS](http://glvis.org).
       Set these variables in the glvis makefile
            
            MFEM_DIR   = mfem_location
            MFEM_LIB   = -L$(MFEM_DIR) -lmfem

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


