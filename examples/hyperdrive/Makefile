##################################################################
# Import machine specific compilers, options, flags, etc.. 
##################################################################

include ../../makefile.inc

.SUFFIXES:
.SUFFIXES: .c .o 


WARP_FLAGS = -I../..
WARP_LIB = ../../libwarp.a

OBJ5 = advect.o c_array.o bop6g.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o assign_gp.o \
       dwdx.o dwdt.o advect_data.o explicit_rk4_stepper.o

OBJT = test-gf.o c_array.o bop6g.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o assign_gp.o \
       dwdx.o dwdt.o advect_data.o

OBJ1 = hyperdrive-01.o c_array.o bop6g.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o assign_gp.o \
       dwdx.o dwdt.o c-twforce.o advect_data.o explicit_rk4_stepper.o save_grid_fcn.o

.c.o:
	$(CC) -c $(CFLAGS) $(WARP_FLAGS) $<

all: advect hyperdrive-01 test-gf

advect: $(OBJ5) 
	$(CC) -o $@ $(OBJ5) $(LFLAGS)

hyperdrive-01: $(OBJ1) 
	$(CC) -o $@ $(OBJ1) $(LFLAGS) $(WARP_LIB)

test-gf: $(OBJT) 
	$(CC) -o $@ $(OBJT) $(LFLAGS) $(WARP_LIB)

clean:
	rm -f *.o *.a
