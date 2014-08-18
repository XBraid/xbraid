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

include ../../makefile.inc

.SUFFIXES:
.SUFFIXES: .c .o 


BRAID_FLAGS = -I../..
BRAID_LIB = ../../libbraid.a

OBJ5 = advect.o c_array.o bop6g.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o assign_gp.o \
       dwdx.o dwdt.o advect_data.o explicit_rk4_stepper.o

OBJT = test-gf.o c_array.o bop6g.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o assign_gp.o \
       dwdx.o dwdt.o advect_data.o

OBJR = restrict.o c_array.o bop6g.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o assign_gp.o \
       dwdx.o dwdt.o advect_data.o explicit_rk4_stepper.o

OBJ1 = hyperdrive-01.o c_array.o bop6g.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o assign_gp.o \
       dwdx.o dwdt.o c-twforce.o advect_data.o explicit_rk4_stepper.o save_grid_fcn.o

.c.o:
	$(CC) -c $(CFLAGS) $(BRAID_FLAGS) $<

all: advect hyperdrive-01 test-gf restrict

advect: $(OBJ5) 
	$(CC) -o $@ $(OBJ5) $(LFLAGS)

hyperdrive-01: $(OBJ1) 
	$(CC) -o $@ $(OBJ1) $(LFLAGS) $(BRAID_LIB)

restrict: $(OBJR) 
	$(CC) -o $@ $(OBJR) $(LFLAGS) $(BRAID_LIB)

test-gf: $(OBJT) 
	$(CC) -o $@ $(OBJT) $(LFLAGS) $(BRAID_LIB)

clean:
	rm -f *.o *.a
