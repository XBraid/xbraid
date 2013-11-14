#BHEADER**********************************************************************
# Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# This file is part of WARP.  See file COPYRIGHT for details.
#
# WARP is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.
#
#EHEADER**********************************************************************

##################################################################
# Import machine specific compilers, options, flags, etc.. 
##################################################################

include ../makefile.inc


##################################################################
# Build exmaples 
##################################################################

HYPRE_DIR = ../../linear_solvers/hypre
HYPRE_FLAGS = -I$(HYPRE_DIR)/include
HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE

MFEM_DIR = ../../mfem
METIS_DIR = ../../metis-4.0
MFEM_FLAGS = -I$(MFEM_DIR)
MFEM_LIB = -L$(MFEM_DIR) -lmfem -L$(METIS_DIR) -lmetis

WARP_FLAGS = -I..
WARP_LIB_FILE = ../libwarp.a

C_NOHYPRE = drive-01
C_EXAMPLES = drive-02 drive-03
# Note: .cpp examples will be linked with mfem
CXX_EXAMPLES = drive-04

.SUFFIXES:
.SUFFIXES: .c .cpp

# put this rule first so it becomes the default
all: $(C_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES)

# Rule for building drive-01
drive-01: drive-01.c $(WARP_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(WARP_FLAGS) $(@).c -o $@ $(WARP_LIB_FILE) $(LFLAGS)

# Rule for compiling .c files
%: %.c $(WARP_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(WARP_FLAGS) $(HYPRE_FLAGS) $(@).c -o $@\
 $(WARP_LIB_FILE) $(HYPRE_LIB) $(LFLAGS)

# Rule for compiling .cpp files; links with mfem
%: %.cpp $(WARP_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICXX) $(CFLAGS) $(WARP_FLAGS) $(HYPRE_FLAGS) $(MFEM_FLAGS)\
 $(@).cpp -o $@ $(MFEM_LIB) $(HYPRE_LIB) $(WARP_LIB_FILE) $(LFLAGS)

clean: cleanout
	rm -f *.o $(C_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES)

cleanout:
	rm -f drive*.out.*
