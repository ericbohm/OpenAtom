
all : libhpm.a

libhpm.a : hpm.o timebase.o
	ar rv libhpm.a hpm.o timebase.o

hpm.o : hpm.c 
	mpicc -c -g -O  -I/bgsys/drivers/ppcfloor/arch/include hpm.c

timebase.o : timebase.c
	mpicc -c -g -O  -I/bgsys/drivers/ppcfloor/arch/include timebase.c

clean :
	rm -f *.a *.o
