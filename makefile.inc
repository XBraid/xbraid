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

# Three compile time options
# make debug=yes|no
# make valgrind=yes|no
# make sequential=yes|no

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

# Was MPI use specified?
ifeq ($(sequential),yes)
   sequential = YES
else
   sequential = NO
endif

# Default compiler options for different platforms and known machines
# Add an else statement for your machine, if the compile options need adjustment

HOSTNAME := $(shell hostname)

ifeq ($(findstring tux408,$(HOSTNAME)),tux408)
   # Jacobs old LLNL tux-boxes
   CC = mpicc
   MPICC = mpicc
   MPICXX = mpicxx
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -lm -fPIC
      CXXFLAGS = -g -Wall -Wno-reorder -lm -fPIC
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall -lm -fPIC
      CXXFLAGS = -O -Wall -Wno-reorder -lm -fPIC
      FORTFLAGS = -O1 -Wall
   endif
   ifeq ($(lvalgrind),YES)
      #only needed for openmpi: LFLAGS += /usr/lib64/valgrind/libmpiwrap-amd64-linux.so 
   endif
else ifeq ($(shell uname -s), Darwin)
   # for Jacob's MacBook running Homebrew 
   # Need to specifically include lstdc++ (!!)
   CC = gcc
   MPICC = mpicc
   MPICXX = mpicxx
   MPIF90 = mpif90
   LFLAGS = -lm -lstdc++
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -fPIC
      CXXFLAGS = -g -Wall -fPIC
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O3 -Wall -fPIC
      CXXFLAGS = -O3 -Wall -fPIC
      FORTFLAGS = -O1 -Wall
   endif
else ifeq ($(findstring cab,$(HOSTNAME)),cab)
   # for Cab
   MPICC = mpiicc
   MPICXX = mpiicpc
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -fPIC
      CXXFLAGS = -g -Wall -fPIC
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall -fPIC
      CXXFLAGS = -O -Wall -fPIC
      FORTFLAGS = -O1 -Wall
   endif   
else ifeq ($(findstring vulcan,$(HOSTNAME)),vulcan)
   # for Vulcan
   MPICC = mpixlc
   MPICXX = mpixlcxx
   MPIF90 = mpixlf90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -fPIC
      CXXFLAGS = -g -Wall -fPIC
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall -fPIC
      CXXFLAGS = -O -Wall -fPIC
      FORTFLAGS = -O1 -Wall
   endif
else ifeq ($(shell uname -s),Linux)
   MPICC = mpicc
   MPICXX = mpiCC
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -fPIC
      CXXFLAGS = -g -Wall -fPIC
      FORTFLAGS = -g -Wall
   else
      CFLAGS = -O3 -Wall -fPIC
      CXXFLAGS = -O3 -Wall -fPIC
      FORTFLAGS = -O1 -Wall
   endif
##
# Add your machine name here
##
else
   MPICC = mpicc
   MPICXX = mpiCC
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -fPIC
      CXXFLAGS = -g -Wall -fPIC
      FORTFLAGS = -g -Wall
   else
      CFLAGS = -O3 -Wall -fPIC
      CXXFLAGS = -O3 -Wall -fPIC
      FORTFLAGS = -O1 -Wall
   endif
endif

# Compiler options when compiling without MPI
ifeq ($(sequential),YES)
   MPICC = gcc
   MPICXX = g++
   MPIF90 = gfortran
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -D braid_SEQUENTIAL
      CXXFLAGS = -g -Wall -D braid_SEQUENTIAL
      FORTFLAGS = -g -Wall -D braid_SEQUENTIAL
   else
      CFLAGS = -O -Wall -D braid_SEQUENTIAL
      CXXFLAGS = -O -Wall -D braid_SEQUENTIAL
      FORTFLAGS = -O1 -Wall -D braid_SEQUENTIAL
   endif
endif

