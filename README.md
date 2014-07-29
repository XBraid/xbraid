## Building Warp
<!--
 - Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
 - Produced at the Lawrence Livermore National Laboratory.
 - This file is part of WARP.  See file COPYRIGHT for details.
 -
 - WARP is free software; you can redistribute it and/or modify it under the
 - terms of the GNU Lesser General Public License (as published by the Free
 - Software Foundation) version 2.1 dated February 1999.
 -->
 
-  To specify the compilers, flags and options for your machine, edit
   makefile.inc.  For now, we keep it simple and avoid using configure or
   cmake.

-  To make the library, libwarp.a,
   
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
