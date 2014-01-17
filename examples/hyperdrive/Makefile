##################################################################
# Import machine specific compilers, options, flags, etc.. 
##################################################################

include ../../makefile.inc

.SUFFIXES:
.SUFFIXES: .c .o 

WARP_FLAGS = -I../..
WARP_LIB = ../../libwarp.a

OBJ5 = kreiss.o c_array.o bop6g.o sbpghost.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o c-bckreiss.o \
       c-dwdtkreiss.o c-twforce.o kreiss_data.o explicit_rk4_stepper.o

OBJ1 = hyperdrive-01.o c_array.o bop6g.o sbpghost.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o c-bckreiss.o \
       c-dwdtkreiss.o c-twforce.o kreiss_data.o explicit_rk4_stepper.o

.c.o:
	$(CC) -c $(CFLAGS) $(WARP_FLAGS) $<

kreiss: $(OBJ5) 
	$(CC) -o $@ $(OBJ5) $(LFLAGS)

hyperdrive-01: $(OBJ1) 
	$(CC) -o $@ $(OBJ1) $(LFLAGS) $(WARP_LIB)

clean:
	rm -f *.o *.a
