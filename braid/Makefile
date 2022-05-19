#BHEADER**********************************************************************
# Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory.
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

include ../makefile.inc

##################################################################
# Targets
##################################################################

BRAID_HEADERS =\
 _braid.h\
 base.h\
 status.h\
 tape.h\
 util.h\
 braid.h\
 braid_status.h\
 braid_test.h

BRAID_FILES =\
 access.c\
 adjoint.c\
 base.c\
 braid.c\
 braid_status.c\
 braid_F90_iface.c\
 braid_test.c\
 communication.c\
 distribution.c\
 drive.c\
 grid.c\
 hierarchy.c\
 interp.c\
 mpistubs.c\
 norm.c\
 refine.c\
 relax.c\
 residual.c\
 restrict.c\
 space.c\
 step.c\
 tape.c\
 util.c\
 uvector.c

ifeq ($(sequential),yes)
	BRAID_HEADERS += mpistubs.h
	BRAID_FILES += mpistubs.c
	SEQFLAGS = -Dbraid_SEQUENTIAL
else
	SEQFLAGS = 
endif

BRAID_OBJ = $(BRAID_FILES:.c=.o)

.PHONY: all clean
.SUFFIXES:
.SUFFIXES: .c .o

# Rule for compiling .c files
%.o: %.c $(BRAID_HEADERS)
	$(MPICC) $(SEQFLAGS) $(CFLAGS) -c $< -o $@

libbraid.a: $(BRAID_HEADERS) $(BRAID_OBJ)
	@echo "Building" $@ "..."
	ar cruv libbraid.a $(BRAID_OBJ)
	ranlib libbraid.a

all: libbraid.a 

clean:
	rm -f *.o libbraid.a

