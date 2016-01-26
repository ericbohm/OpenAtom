# Makefile  

# July. 2015  minjung.kim@yale.edu

# last update: Sep. 2015 by Minjung

# compilers
CC = g++
CFLAGS = -g -O3 
FFTWPATH = /Users/Eric/Programming/ppl/fftw3
FFTLIB = -L$(FFTWPATH)/lib -lfftw3 -lm
FFTINC = -I$(FFTWPATH)/include
LAPACKLIB = -framework Accelerate

MYFLAGS = $(FFTINC) -DUSE_LAPACK -DDEBUG #-DBGW


.PHONY: all clean 
.SUFFIXES: .C .o

# source codes
SRCF= print_util.C matrix.C mylapack.C\
      iter_invmtrx.C coulomb.C calc_Eps.C gspace.C kqindex.C\
      fft_routines.C my_fftw.C calc_P.C do_fft_states.C\
      get_fftsize.C read_states.C\
      check_inputs.C mtrxop_3x3.C sysinfo.C usrinput.C util.C\
      epsilon.C # main.C 




OBJF= $(SRCF:%.C=%.o)

# target
epsilon: $(OBJF)
	$(CC) $(OBJF) -o $@ $(FFTLIB) $(LAPACKLIB)

.C.o: 
	$(CC) $(CFLAGS) $(MYFLAGS) -c $< -o $@

clean:
	rm *.o epsilon

depend:
	$(CC) -MM $(SRCF) 



## dependency

