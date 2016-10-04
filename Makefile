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

BRAID_DIR=.
include makefile.inc


##################################################################
# Targets
##################################################################

BRAID_HEADERS = _braid.h braid.h _util.h braid_test.h braid_status.h braid_defs.h

BRAID_FILES = _util.c braid.c _braid.c braid_test.c _braid_status.c braid_F90_iface.c

ifeq ($(sequential),yes)
	BRAID_HEADERS += mpistubs.h
	BRAID_FILES += mpistubs.c
	SEQFLAGS = -Dbraid_SEQUENTIAL
else
	SEQFLAGS = 
endif

BRAID_OBJ = $(BRAID_FILES:.c=.o)

.PHONY: all examples drivers clean info
.SUFFIXES:
.SUFFIXES: .c .o

# Rule for compiling .c files
%.o: %.c $(BRAID_HEADERS)
	$(MPICC) $(SEQFLAGS) $(CFLAGS) -c $< -o $@

libbraid.a: $(BRAID_HEADERS) $(BRAID_OBJ)
	@echo "Building" $@ "..."
	ar cruv libbraid.a $(BRAID_OBJ)
	ranlib libbraid.a

all: libbraid.a examples drivers

examples: libbraid.a
	cd examples; $(MAKE)

drivers: libbraid.a
	cd drivers; $(MAKE)

clean:
	rm -f *.o libbraid.a
	cd examples; $(MAKE) clean
	cd drivers; $(MAKE) clean

info:
	@echo "MPICC     = `which $(MPICC)`"
	@echo "MPICXX    = `which $(MPICXX)`"
	@echo "MPIF90    = `which $(MPIF90)`"
	@echo "CFLAGS    = $(CFLAGS)"
	@echo "CXXFLAGS  = $(CXXFLAGS)"
	@echo "FORTFLAGS = $(FORTFLAGS)"
	@echo "LFLAGS    = $(LFLAGS)"
