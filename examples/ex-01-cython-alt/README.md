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


This is an example using the Braid-Cython interface defined in braid.pyx (
braid/braid.pyx ).  This example uses a lower-level C-style syntax than the
other basic Cython example in `examples/ex-01-cython`

It solves the same scalar ODE equation 
 
 u' = lambda u, 
 with lambda=-1 and y(0) = 1

as the ex-01 series in the `examples` directory. 

For instructions on running and compiling, see 

   examples/ex-01-cython-alt/ex_01_alt-setup.py

and 

   examples/ex-01-cython-alt/ex_01_alt.pyx

