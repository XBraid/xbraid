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
	$(MPICC) -c $(CFLAGS) $(BRAID_FLAGS) $<

all: advect hyperdrive-01 test-gf restrict

advect: $(OBJ5) 
	$(MPICC) -o $@ $(OBJ5) $(LFLAGS) $(BRAID_LIB)

hyperdrive-01: $(OBJ1) 
	$(MPICC) -o $@ $(OBJ1) $(LFLAGS) $(BRAID_LIB)

restrict: $(OBJR) 
	$(MPICC) -o $@ $(OBJR) $(LFLAGS) $(BRAID_LIB)

test-gf: $(OBJT) 
	$(MPICC) -o $@ $(OBJT) $(LFLAGS) $(BRAID_LIB)

clean:
	rm -f *.o *.a
