#BHEADER**********************************************************************
#
# Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. Written by 
# Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
# Dobrev, et al. LLNL-CODE-660355. All rights reserved.
# 
# This file is part of XBraid. Email xbraid-support@llnl.gov for support.
# 
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2.1 dated February 1999.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with this program; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#EHEADER**********************************************************************

# Was DEBUG specified? 
ifeq ($(debug),no)
   optlevel = OPTIMIZE
else
   optlevel = DEBUG
endif

# Was VALGRIND specified? 
ifeq ($(valgrind),yes)
   lvalgrind = YES
else
   lvalgrind = NO
endif

# Default compiler options for different platforms and known machines
# Add an else statement for your machine, if the compile options need adjustment

HOSTNAME := $(shell hostname)

ifeq ($(findstring tux343,$(HOSTNAME)),tux343)
   CC = mpicc
   MPICC = mpicc
   MPICXX = mpiCC
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -lm
      CXXFLAGS = -g -Wall -Wno-reorder -lm
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall -lm
      CXXFLAGS = -O -Wall -Wno-reorder -lm
      FORTFLAGS = -O1 -Wall
   endif
   ifeq ($(lvalgrind),YES)
      LFLAGS += /usr/lib64/valgrind/libmpiwrap-amd64-linux.so 
   endif
else ifeq ($(findstring fourier,$(HOSTNAME)),fourier)
   # for Anders' laptop
   CC = gcc
   MPICC = mpicc
   MPICXX = mpicxx
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
else ifeq ($(findstring sierra,$(HOSTNAME)),sierra)
   # for Sierra
   MPICC = mpiicc
   MPICXX = mpiicpc
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
else ifeq ($(findstring cab,$(HOSTNAME)),cab)
   # for Cab
   MPICC = mpiicc
   MPICXX = mpiicpc
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif   
else ifeq ($(findstring vulcan,$(HOSTNAME)),vulcan)
   # for Vulcan
   MPICC = mpixlc
   MPICXX = mpixlcxx
   MPIF90 = mpixlf90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
else ifeq ($(shell uname -s),Linux)
   MPICC = mpicc
   MPICXX = mpiCC
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -g -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
else
   MPICC = mpicc
   MPICXX = mpiCC
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -g -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
endif

# Rob's flags for running against a memory checker

# MPICC    = insure -g
# MPICXX   = insure -g
# CFLAGS   = -I/home/falgout/codes/mpich2-1.4.1-install/include
# CXXFLAGS = -I/home/falgout/codes/mpich2-1.4.1-install/include
# LFLAGS   = -L/home/falgout/codes/mpich2-1.4.1-install/lib -Wl,-rpath,/home/falgout/codes/mpich2-1.4.1-install/lib -lmpich -lopa -lmpl -lrt -lpthread  -lstdc++ -lm
