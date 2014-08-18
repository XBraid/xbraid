#BHEADER**********************************************************************
# Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# This file is part of XBraid.  See file COPYRIGHT for details.
#
# XBraid is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.
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
	LFLAGS = -lm
	ifeq ($(optlevel),DEBUG)
		CFLAGS = -g -Wall -lm
		CXXFLAGS = -g -Wall -Wno-reorder -lm
	else
		CFLAGS = -O -Wall -lm
		CXXFLAGS = -O -Wall -Wno-reorder -lm
	endif
	ifeq ($(lvalgrind),YES)
		LFLAGS += /usr/lib64/valgrind/libmpiwrap-amd64-linux.so 
   endif
else ifeq ($(findstring yorkville,$(HOSTNAME)),yorkville)
  # for Anders' laptop
	MPICC = openmpicc
	MPICXX = openmpicxx
	LFLAGS = -lm
	ifeq ($(optlevel),DEBUG)
		CFLAGS = -g -Wall
		CXXFLAGS = -g -Wall
	else
		CFLAGS = -O -Wall
		CXXFLAGS = -O -Wall
	endif
else ifeq ($(findstring sierra,$(HOSTNAME)),sierra)
  # for Sierra
	MPICC = mpiicc
	MPICXX = mpiicpc
	LFLAGS = -lm
	ifeq ($(optlevel),DEBUG)
		CFLAGS = -g -Wall
		CXXFLAGS = -g -Wall
	else
		CFLAGS = -O -Wall
		CXXFLAGS = -O -Wall
	endif
else ifeq ($(shell uname -s),Linux)
	MPICC = mpicc
	MPICXX = mpiCC
	LFLAGS = -lm
	ifeq ($(optlevel),DEBUG)
		CFLAGS = -g -Wall
		CXXFLAGS = -g -Wall
	else
		CFLAGS = -O -Wall
		CXXFLAGS = -O -Wall
	endif
else
	MPICC = mpicc
	MPICXX = mpiCC
	LFLAGS = -lm
	ifeq ($(optlevel),DEBUG)
		CFLAGS = -g -Wall
		CXXFLAGS = -g -Wall
	else
		CFLAGS = -O -Wall
		CXXFLAGS = -O -Wall
	endif
endif

# Rob's flags for running against a memory checker
#
# MPICC = insure -g
# CFLAGS = -I/home/falgout/codes/mpich2-1.4.1-install/include
#
# MPICC = insure -g
# CFLAGS = -I/home/falgout/codes/mpich2-1.4.1-install/include
# LFLAGS = -L/home/falgout/codes/mpich2-1.4.1-install/lib -Wl,-rpath,/home/falgout/codes/mpich2-1.4.1-install/lib -lmpich -lopa -lmpl -lrt -lpthread  -lstdc++ -lm

