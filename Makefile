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

MPICC ?= mpicc
CFLAGS = -g -Wall

# MPICC = insure -g
# CFLAGS = -I/home/falgout/codes/mpich2-1.4.1-install/include

##################################################################
# Targets
##################################################################

WARP_HEADERS = _warp.h warp.h util.h

WARP_FILES = util.c warp.c _warp.c

WARP_OBJ = $(WARP_FILES:.c=.o)

.PHONY: examples
.SUFFIXES:
.SUFFIXES: .c .o

# Rule for compiling .c files
%.o: %.c
	$(MPICC) $(CFLAGS) -c $< -o $@

libwarp.a: $(WARP_HEADERS) $(WARP_OBJ)
	@echo "Building" $@ "..."
	ar cruv libwarp.a $(WARP_OBJ)
	ranlib libwarp.a

all: libwarp.a examples

examples: libwarp.a
	cd examples; $(MAKE)

clean:
	rm -f *.o libwarp.a
	cd examples; $(MAKE) clean
