## Drivers: compiling and running
<!--
  - Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
  - Produced at the Lawrence Livermore National Laboratory. Written by 
  - Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
  - Dobrev, et al. LLNL-CODE-660355. All rights reserved.
  - 
  - This file is part of XBraid. For support, post issues to the XBraid Github page.
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

      drive-* -help

for instructions on how to run any driver.

To run the examples, type
   
      mpirun -np 4 drive-*  [args]


1. *drive-diffusion-2D* implements the 2D heat equation on a regular grid.  You must have
   [hypre](https://computation.llnl.gov/project/linear_solvers/software.php)
   installed and these variables in examples/Makefile set correctly
    
          HYPRE_DIR = ../../linear_solvers/hypre
          HYPRE_FLAGS = -I$(HYPRE_DIR)/include
          HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE
   
   This driver also support spatial coarsening and explicit time stepping.
   This allows you to use explicit time stepping on each Braid level, 
   regardless of time step size.  

2. *drive-burgers-1D* implements Burger's equation (and also linear advection) in 1D
   using forward or backward Euler in time and Lax-Friedrichs in space.  Spatial
   coarsening is supported, allowing for stable time stepping on coarse time-grids.

   See also *viz-burgers.py* for visualizing the output.

3. *drive-diffusion* is a sophisticated test bed for finite element discretizations of the 
   heat equation. It relies on the [mfem](http://mfem.org)
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
            
            make drive-diffusion

     + To run drive-diffusion and glvis, open two windows.  In one, start a glvis session
      
               ./glvis
  
         Then, in the other window, run drive-diffusion
      
               mpirun -np ... drive-diffusion [args]
         
         Glvis will listen on a port to which drive-diffusion will dump visualization 
         information.

4. The other drive-.cpp files use MFEM to implement other PDEs

     + *drive-adv-diff-DG*:  implements advection(-diffusion) with a discontinuous Galerkin
       discretization.  This driver is under developement.

     + *drive-diffusion-1D-moving-mesh*:  implements the 1D heat equation, but with a
       moving mesh that adapts to the forcing function so that the mesh 
       equidistributes the arc-length of the solution.

     + *drive-diffusion-1D-moving-mesh-serial*:  implements a serial time-stepping version
       of the above problem.

     + *drive-pLaplacian*:  implements the 2D the \f$p\f$-Laplacian (nonlinear diffusion).

     + *drive-diffusion-ben*:  implements the 2D/3D diffusion equation with time-dependent
       coefficients. This is essentially equivalent to drive-diffusion, and could be 
       removed, but we're keeping it around because it implements linear diffusion in the 
       same way that the p-Laplacian driver implemented nonlinear diffusion.  This makes it
       suitable for head-to-head timings.

     + *drive-lin-elasticity*:  implements time-dependent linearized elasticity and is under
       development.

     + *drive-nonlin-elasticity*:  implements time-dependent nonlinear elasticity and is under
       development.

5. Directory drive-adv-diff-1D-Cython/  solves a simple 1D advection-diffussion equation 
   using the Cython interface and numerous spatial and temporal discretizations
   
                drivers/drive-adv-diff-1D-Cython/drive_adv_diff_1D.pyx
      
      and
      
                drivers/drive-adv-diff-1D-Cython/drive_adv_diff_1D-setup.py

6. Directory drive-Lorenz-Delta/ implements the chaotic Lorenz system, with its trademark
   butterfly shaped attractor. The driver uses the Delta correction feature and Lyapunov
   estimation to solve for the backward Lyapunov vectors of the system and to accelerate
   XBraid convergence. Visualize the solution and the Lyapunov vectors with *vis_lorenz_LRDelta.py*
   Also see example 7 (examples/ex-07.c). *This driver is in a broken state, and needs*
   *updating for compatibility with new Delta correction implementation.*

7. Directory drive-KS-Delta/ solves the chaotic Kuramoto-Sivashinsky equation in 1D, using 
   fourth order finite differencing in space and the Lobatto IIIC fully implicit RK method 
   in time. The driver also uses Delta correction and Lyapunov estimation to accelerate
   convergence and to generate estimates to the unstable Lyapunov vectors for the system.
