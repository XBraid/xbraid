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

##################################################################
# Import machine specific compilers, options, flags, etc.. 
##################################################################

BRAID_DIR=..
include $(BRAID_DIR)/makefile.inc


##################################################################
# Build exmaples 
##################################################################

HYPRE_DIR = $(BRAID_DIR)/../linear_solvers/src/hypre
HYPRE_FLAGS = -I$(HYPRE_DIR)/include
HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE

MFEM_DIR = $(BRAID_DIR)/../mfem
MFEM_CONFIG_MK = $(MFEM_DIR)/config/config.mk
MFEM_LIB_FILE = mfem_is_not_built
-include $(MFEM_CONFIG_MK)

BRAID_FLAGS = -I$(BRAID_DIR)
BRAID_LIB_FILE = $(BRAID_DIR)/libbraid.a

C_NOHYPRE = ex-01 ex-burgers ex-lorenz
F_NOHYPRE = ex-01f ex-fburgers
C_EXAMPLES = ex-02 ex-02-serial
# Note: .cpp examples will be linked with mfem
#CXX_EXAMPLES = ex-04

.PHONY: all clean cleanout

.SUFFIXES:
.SUFFIXES: .c .cpp

# put this rule first so it becomes the default
all: $(C_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES)

# Rule for building ex-burgers
ex-burgers: ex-burgers.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

# Rule for building ex-lorenz
ex-lorenz: ex-lorenz.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

# Rule for building ex-01
ex-01: ex-01.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(@).c -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

# Rule for building ex-01f
ex-01f: ex-01f.f90 $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPIF90) $(FORTFLAGS) $(@).f90 -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

# Rule for building ex-fburgers
ex-fburgers: ex-fburgers.f90 $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPIF90) $(FORTFLAGS) $(@).f90 -o $@ $(BRAID_LIB_FILE) $(LFLAGS)

# Rule for compiling .c files
%: %.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(HYPRE_FLAGS) $(@).c -o $@\
 $(BRAID_LIB_FILE) $(HYPRE_LIB) $(LFLAGS)

# Rule for compiling .cpp files; links with mfem
%: %.cpp $(BRAID_LIB_FILE) $(MFEM_LIB_FILE) $(MFEM_CONFIG_MK)
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) $(BRAID_FLAGS) $(MFEM_FLAGS) \
	$< -o $@ $(MFEM_LIBS) $(BRAID_LIB_FILE) $(LFLAGS)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: cleanout
	rm -f *.o $(C_NOHYPRE) $(F_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES) $(F_EXAMPLES) *ror_norm* *_err_* *_mesh* *_sol_*
	rm -rf *.dSYM

cleanout:
	rm -f ex*.out.*
