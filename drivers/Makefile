#BHEADER**********************************************************************
#
# Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. Written by 
# Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
# Dobrev, et al. LLNL-CODE-660355. All rights reserved.
# 
# This file is part of XBraid. For support, post issues to the XBraid Github page.
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

BRAID_DIR=../braid
include ../makefile.inc


##################################################################
# Build exmaples 
##################################################################

HYPRE_DIR = ../../hypre/src/hypre
HYPRE_FLAGS = -I$(HYPRE_DIR)/include
HYPRE_LIB = -L$(HYPRE_DIR)/lib -lHYPRE

MFEM_INTERFACE = ./braid_mfem.hpp ./braid_mfem_block.hpp ../braid/braid.hpp
MFEM_DIR = ../../mfem
MFEM_CONFIG_MK = $(MFEM_DIR)/config/config.mk
MFEM_LIB_FILE = mfem_is_not_built
-include $(MFEM_CONFIG_MK)

BRAID_FLAGS = -I$(BRAID_DIR)
BRAID_LIB_FILE = $(BRAID_DIR)/libbraid.a

C_EXAMPLES = drive-burgers-1D drive-diffusion-2D drive-lorenz drive-solve-adjoint-with-xbraid
# Note: .cpp examples will be linked with mfem
CXX_EXAMPLES = drive-diffusion-1D-moving-mesh drive-diffusion-1D-moving-mesh-serial \
					drive-diffusion-ben drive-diffusion drive-lin-elasticity \
					drive-pLaplacian drive-adv-diff-DG drive-nonlin-elasticity

# clean all possible
ALL_EXAMPLES := $(patsubst %.cpp,%,$(wildcard *.c *.cpp))
ALL_EXAMPLES := $(patsubst %.c,%,$(ALL_EXAMPLES))

.PHONY: all clean cleanout

.SUFFIXES:
.SUFFIXES: .c .cpp

# put this rule first so it becomes the default
all: $(C_NOHYPRE) $(C_EXAMPLES) $(CXX_EXAMPLES)

# Rule for compiling .c files
%: %.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) $(HYPRE_FLAGS) $< -o $@\
 $(BRAID_LIB_FILE) $(HYPRE_LIB) $(LFLAGS)

# Rule for compiling .cpp files; links with mfem
%: %.cpp $(BRAID_LIB_FILE) $(MFEM_LIB_FILE) $(MFEM_CONFIG_MK) $(MFEM_INTERFACE)
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) $(BRAID_FLAGS) $(MFEM_FLAGS) \
 $< -o $@ $(MFEM_LIBS) $(BRAID_LIB_FILE) $(LFLAGS)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: cleanout
	rm -f *.o $(ALL_EXAMPLES) *ror_norm* *_err_* *_sol_*
	rm -rf *.dSYM

cleanout:
	rm -f drive*.out.*
