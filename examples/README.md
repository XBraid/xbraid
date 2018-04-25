## Examples: compiling and running
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

      ex-* -help

for instructions on how to run any example.

To run the examples, type
   
      mpirun -np 4 ex-*  [args]


1. ex-01 is the simplest example.  It implements a scalar ODE and can be
  compiled and run with no outside dependencies.  See Section (@ref exampleone)
  for more discussion of this example.  There are seven versions of this example,
  
    + *ex-01.c*:  simplest possible implementation, start reading this example first
    
    + *ex-01-expanded.c*:  same as ex-01.c but adds more XBraid features
    
    + *ex-01-expanded-bdf2.c*:  same as ex-01-expanded.c, but uses BDF2 instead 
      of backward Euler
    
    + *ex-01-expanded-f.f90*:  same as ex-01-expanded.c, but implemented in f90

    + *ex-01-refinement.c*: same as ex-01.c, but adds the refinement feature

    + *ex-01-adjoint.c*: adds adjoint-based gradient computation to ex-01.c

    + *ex-01-optimization.c*: gradient-based optimization cycle for ex-01-c

2. ex-02 implements the 1D heat equation on a regular grid, using a very simple
   implementation.  This is the next example to read after the various ex-01
   cases.

3. ex-03 implements the 2D heat equation on a regular grid.  You must have
   [hypre](https://computation.llnl.gov/project/linear_solvers/software.php)
   installed and these variables in examples/Makefile set correctly
    
          HYPRE_DIR = ../../linear_solvers/hypre
          HYPRE_FLAGS = -I$(HYPRE_DIR)/include
          HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE

   Only implicit time stepping (backward Euler) is supported.  See Section
   (@ref examplethree) for more discussion of this example.  The driver

          drivers/drive-diffusion
   
   is a more sophisticated version of this simple example that supports
   explicit time stepping and spatial coarsening.

4. ex-04 solves an optimal control problem with time-dependent design variable using a simple steepest-descent optimization iteration. 
   

