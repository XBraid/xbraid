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

![](docs/img/logo_with_subtext_2_inch.png)

## What is XBraid?

XBraid is a parallel-in-time software package.  It implements an
optimal-scaling multigrid solver for the (non)linear systems that arise from
the discretization of problems with evolutionary behavior. 

This code and associated algorithms are developed at Lawrence Livermore
National Laboratory 
(see [website](https://computation.llnl.gov/projects/parallel-time-integration-multigrid/),
and at collaborating [academic institutions](https://github.com/XBraid/xbraid/wiki/Team). 

For our publication list, please go [here](https://github.com/XBraid/xbraid/wiki/Project-Publications).

### About XBraid

Typically, solution algorithms for evolution equations are based on a
time-marching approach, solving sequentially for one time step after the other.
Parallelism in these traditional time-integration techniques is limited to
spatial parallelism.  However, current trends in computer architectures are
leading towards systems with more, but not faster, processors, i.e., clock
speeds are stagnate.  Therefore, faster overall runtimes must come from greater
parallelism. Our approach to achieve such parallelism in time is with multigrid.

In this software, we implement a non-intrusive, optimal-scaling time-parallel
method based on multigrid reduction techniques (multigrid-reduction-in-time or
MGRIT).  A few important points about XBraid are as follows.

- The algorithm enables a scalable parallel-in-time approach by applying multigrid to the time dimension.

- It is designed to be nonintrusive. That is, users apply their existing
  sequential time-stepping code according to our interface, and then XBraid
  does the rest. Users have spent years, sometimes decades, developing the
  right time-stepping scheme for their problem. XBraid allows users to keep
  their schemes, but enjoy parallelism in the time dimension.

- XBraid solves exactly the same problem that the existing sequential
  time-stepping scheme does.

- XBraid is flexible, allowing for a variety of time stepping, relaxation, and
  temporal and spatial coarsening options.

- The full approximation scheme multigrid approach is used to accommodate
  nonlinear problems.

- XBraid written in MPI/C with C++ and Fortran 90 interfaces.

- XBraid is released under LGPL 2.1.


## Documentation

- For examples of using XBraid, see the
  [examples/](https://github.com/XBraid/xbraid/tree/master/examples) and
  [drivers/](https://github.com/XBraid/xbraid/tree/master/drivers) directories,
  and in particular examples/ex-01-*.c

- See the [release](https://github.com/XBraid/xbraid/releases) page for links
  to precompiled documentation PDFs that go through, step-by-step how to use
  XBraid 

- For tutorials, see the bottom of our publications 
[page](https://github.com/XBraid/xbraid/wiki/Project-Publications).

## Building XBraid

-  To specify the compilers, flags and options for your machine, edit
   makefile.inc.  For now, we keep it simple and avoid using configure or
   cmake.

-  To make the library, libbraid.a,
   
         $ make

-  To make the examples
   
         $ make all

-  The makefile lets you pass some parameters like debug with 
   
         $ make debug=yes
   
   or
   
         $ make all debug=yes
   
   It would also be easy to add additional parameters, e.g., to compile with
   insure.  


- To set compilers and library locations, look in makefile.inc
  where you can set up an option for your machine to define simple
  stuff like

       CC = mpicc
       MPICC = mpicc
       MPICXX = mpiCC
       LFLAGS = -lm


## License

This project is released under the LGPL v2.1 license. See files COPYRIGHT and
LICENSE file for full details.

LLNL Release Number: LLNL-CODE-660355

