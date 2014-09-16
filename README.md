## Building XBraid
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

- Copyright information and licensing restrictions can be found in the files 
  COPYRIGHT and LICENSE.

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
