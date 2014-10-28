
all : libhpm.a

libhpm.a : hpm.o timebase.o
ar rv libhpm.a hpm.o timebase.o

hpm.o : hpm.c 
mpixlc -c -g -O  -I/bgsys/drivers/ppcfloor/arch/include hpm.c

timebase.o : timebase.c
mpixlc -c -g -O  -I/bgsys/drivers/ppcfloor/arch/include timebase.c

clean :
rm -f *.a *.o
