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

CC = mpicc -g -Wall
LFLAGS = -lm

# CC = insure -g
# LFLAGS = -I/home/falgout/codes/mpich2-1.4.1-install/include -L/home/falgout/codes/mpich2-1.4.1-install/lib -Wl,-rpath,/home/falgout/codes/mpich2-1.4.1-install/lib -lmpich -lopa -lmpl -lrt -lpthread  -lstdc++ -lm

HYPRE_DIR = ../warp-hypre/linear_solvers/hypre
HYPRE_FLAGS = -I$(HYPRE_DIR)/include -L$(HYPRE_DIR)/lib -lHYPRE

##################################################################
# Targets
##################################################################

WARP_HEADERS = _warp.h warp.h

WARP_FILES = util.c warp.c

all: drive-01 drive-02

clean: cleanout
	rm -f *.o drive-01 drive-02
cleanout:
	rm -f drive*.out.*

##################################################################
# Rules
##################################################################

drive-01: drive-01.c ${WARP_FILES}
	@echo  "Building" $@ "... "
	${CC} -o $@ $@.c ${WARP_FILES} ${LFLAGS}

drive-02: drive-02.c ${WARP_FILES}
	@echo  "Building" $@ "... "
	${CC} -o $@ $@.c ${WARP_FILES} ${HYPRE_FLAGS} ${LFLAGS}
