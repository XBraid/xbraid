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

BRAID_DIR=../../braid
include ../../makefile.inc


##################################################################
# Build
##################################################################

SUITE_SPARSE = /usr/include/suitesparse/

MFEM_INTERFACE = ./braid_mfem.hpp ./braid_mfem_block.hpp ../braid/braid.hpp
MFEM_DIR = ../../mfem
MFEM_CONFIG_MK = $(MFEM_DIR)/config/config.mk
MFEM_LIB_FILE = mfem_is_not_built
-include $(MFEM_CONFIG_MK)

BRAID_FLAGS = -I$(BRAID_DIR) -Wno-deprecated-declarations
BRAID_LIB_FILE = $(BRAID_DIR)/libbraid.a

# clean all possible
ALL_EXAMPLES := $(patsubst %.cpp,%,$(wildcard *.c *.cpp))
ALL_EXAMPLES := $(patsubst %.c,%,$(ALL_EXAMPLES))

.PHONY: all clean cleanout

.SUFFIXES:
.SUFFIXES: .c .cpp

drive-lorenz-Delta: drive-lorenz-Delta.cpp $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) $(BRAID_FLAGS) -Idrive-lorenz-Delta-lib -IEigen \
 $< -o $@ $(BRAID_LIB_FILE) lorenz_lib.cpp $(LFLAGS)

drive-lorenz-LRDelta: drive-lorenz-LRDelta.cpp $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) $(BRAID_FLAGS) -Idrive-lorenz-Delta-lib -IEigen \
 $< -o $@ $(BRAID_LIB_FILE) lorenz_lib.cpp $(LFLAGS)

clean: cleanout
	rm -f *.o $(ALL_EXAMPLES) *ror_norm* *_err_* *_sol_*
	rm -rf *.dSYM

cleanout:
	rm -f drive*.out.*
