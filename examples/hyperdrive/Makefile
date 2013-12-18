.SUFFIXES:
.SUFFIXES: .c .o 

CC = gcc
CFLAGS = -g

.c.o:
	$(CC) -c $(CFLAGS) $<

OBJ5 = hyperdrive-01.o c_array.o bop6g.o sbpghost.o c-exact.o c-bdata.o c-twbndry.o c-evalerr.o c-bckreiss.o \
       c-dwdtkreiss.o c-twforce.o kreiss_data.o explicit_rk4_stepper.o

hyperdrive-01: $(OBJ5) 
	$(CC) -o $@ $(OBJ5) 

clean:
	rm -f *.o *.a
