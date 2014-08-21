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

##################################################################
# Import machine specific compilers, options, flags, etc.. 
##################################################################

include makefile.inc


##################################################################
# Targets
##################################################################

BRAID_HEADERS = _braid.h braid.h util.h braid_test.h braid_status.h braid_defs.h

BRAID_FILES = util.c braid.c _braid.c braid_test.c braid_status.c

BRAID_OBJ = $(BRAID_FILES:.c=.o)

.PHONY: examples
.SUFFIXES:
.SUFFIXES: .c .o

# Rule for compiling .c files
%.o: %.c $(BRAID_HEADERS)
	$(MPICC) $(CFLAGS) -c $< -o $@

libbraid.a: $(BRAID_HEADERS) $(BRAID_OBJ)
	@echo "Building" $@ "..."
	ar cruv libbraid.a $(BRAID_OBJ)
	ranlib libbraid.a

all: libbraid.a examples

examples: libbraid.a
	cd examples; $(MAKE)

clean:
	rm -f *.o libbraid.a
	cd examples; $(MAKE) clean
