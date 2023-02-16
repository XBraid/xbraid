## Examples: compiling and running
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

For C/C++/Fortran examples, type

      ex-* -help

for instructions on how to run.  To run the C/C++/Fortran examples, type
   
      mpirun -np 4 ex-*  [args]

For the Cython examples, see the corresponding *.pyx file. 



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

    + *ex-01-cython/*: is a  directory containing an example using the Braid-Cython
       interface defined in braid.pyx ( braid/braid.pyx ).  It solves the same
       scalar ODE equation as the ex-01 series described above.  This example uses
       a Python-like syntax, in contrast to the ex-01-cython-alt example, which uses
       a C-style syntax.
    
       For instructions on running and compiling, see 
       
                examples/ex-01-cython/ex_01.pyx
       
       and
       
                examples/ex-01-cython/ex_01-setup.py
   
    + *ex-01-cython-alt/*: is a directory containing another example using the
      Braid-Cython interface defined in braid.pyx ( braid/braid.pyx ).  It solves
      the same scalar ODE equation as the ex-01 series described above.  This example
      uses a lower-level C-like syntax for most of it's code, in contrast to the 
      ex-01-cython example, which uses a Python-style syntax.  
      
      For instructions on running and compiling, see
      
                examples/ex-01-cython-alt/ex_01_alt.pyx
      
      and
      
                examples/ex-01-cython-alt/ex_01_alt-setup.py
   
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

4. ex-04 solves a simple optimal control problem with time-dependent design variable
   using a simple steepest-descent optimization iteration.  
   
5. Directory ex-05-cython/  solves a simple 1D heat equation using the Cython interface 
   
                examples/ex-05-cython/ex_05.pyx
      
      and
      
                examples/ex-05-cython/ex_05-setup.py

6. ex-06 solves a simple scalar ODE, but allows for use of the built-in
   Richardson-based error estimator and accuracy improving extrapolation. With
   the "-refinet" option, the error estimator allows for adaptive refinement in
   time, and with the "-richardson" option, Richardson extrapolation is used
   improve the solution at fine-level C-points. 

   The viz script, 
               
               examples/viz-ex-06.py
   
   allows you to visualize the solution, error, and error estimate.  The use of
   "-richardson" notably improves the accuracy of the solution.
    
   The Richardson-based error estimates and/or extrapolation are only available
   after the first Braid iteration, in that the coarse level solution must be
   available to compute the error estimate and extrapolation.  Thus, after an
   adaptive refinement (and new hierarchy is constructed), another iteration is
   again required for the error estimate to be available.  If the error
   estimate isn't available, Braid returns a value of -1.  See this example and
   the comments therein for more details.

7. ex-07 solves the chaotic Lorenz system, utilizing the Delta correction
   feature to accelerate Braid convergence while estimating the Lyapunov vectors
   and Lyapunov exponents.

   The viz script,

               examples/viz-ex-07.py
   
   Plots the solution trajectory in 3D along with the estimated Lyapunov basis
   vectors computed by Braid. The Lyapunov vectors define a basis for the stable,
   neutral, and unstable manifolds of the system, and the Lyapunov exponents give
   qualitative information about the dynamics of the system.

   The command line argument "-rank" controls the number of Lyapunov vectors
   which are tracked, with "-rank 0" turning Delta correction off, and "-rank 3"
   giving a full-rank Delta correction, since the Lorenz system is 3-dimensional.
   The "-defer-lvl" and "-defer-iter" arguments control whether the Delta
   correction is deferred to a coarse level, or later iteration, respectively.
   For more information about these options, use "$ examples/ex-07 -help".
