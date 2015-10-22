# Makefile  

# July. 2015  minjung.kim@yale.edu

# last update: Sep. 2015 by Minjung

# compilers
CC = g++
CFLAGS = -g -O3
FFTWPATH = /sw/lib
FFTLIB = -L$(FFTWPATH) -lfftw3 -lm
FFTINC = -I/sw/include
LAPACKLIB = -framework Accelerate

MYFLAGS = -DUSE_LAPACK -DDEBUG


.PHONY: all clean 
.SUFFIXES: .C .o

# source codes
SRCF= matrix.C mylapack.C\
      iter_invmtrx.C coulomb.C calc_Eps.C gspace.C kqindex.C\
      fft_routines.C my_fftw.C calc_P.C do_fft_states.C\
      get_fftsize.C read_states.C read_usrinput.C\
      mtrxop_3x3.C read_sysinfo.C util.C epsilon.C 




OBJF= $(SRCF:%.C=%.o)

# target
epsilon: $(OBJF)
	$(CC) $(OBJF) -o $@ $(FFTINC) $(FFTLIB) $(LAPACKLIB)

.C.o: 
	$(CC) $(CFLAGS) $(MYFLAGS) -c $< -o $@

clean:
	rm *.o epsilon 

depend:
	$(CC) -MM $(SRCF) 



## dependency

