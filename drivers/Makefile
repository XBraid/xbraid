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

BRAID_FLAGS = -I..
BRAID_LIB_FILE = ../libbraid.a

C_EXAMPLES = drive-02 #drive-03
# Note: .cpp examples will be linked with mfem
CXX_EXAMPLES = drive-04 #drive-05 drive-06 drive-07

.SUFFIXES:
.SUFFIXES: .c .cpp

# put this rule first so it becomes the default
all: $(C_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES)

# Rule for compiling .c files
%: %.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(HYPRE_FLAGS) $(@).c -o $@\
 $(BRAID_LIB_FILE) $(HYPRE_LIB) $(LFLAGS)

# Rule for compiling .cpp files; links with mfem
%: %.cpp $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) $(BRAID_FLAGS) $(HYPRE_FLAGS) $(MFEM_FLAGS)\
 $(@).cpp -o $@ $(MFEM_LIB) $(HYPRE_LIB) $(BRAID_LIB_FILE) $(LFLAGS)

clean: cleanout
	rm -f *.o $(C_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES) *ror_norm* *_err_* *_sol_*
	rm -rf *.dSYM

cleanout:
	rm -f drive*.out.*
