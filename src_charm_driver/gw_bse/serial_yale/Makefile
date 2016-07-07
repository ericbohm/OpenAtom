# Makefile

# July. 2015  minjung.kim@yale.edu

# last update: Jan. 2016 by Minjung

# compilers
CC = g++
CFLAGS = -g -O3 
FFTWPATH = /sw/lib
FFTLIB = -L$(FFTWPATH) -lfftw3 -lm
FFTINC = -I/sw/include
LAPACKLIB = -framework Accelerate

MYFLAGS = -DUSE_LAPACK #-DUSE_P_INTERPOLATION #-DGPP 
DEBFLAGS = -DDEBUG -DVERBOSE #-DBGW


.PHONY: all clean 
.SUFFIXES: .C .o

# source codes
SRCF= debug_gsize.C\
      vcoulb.C\
      print_util.C\
      gpp.C matrix.C mylapack.C\
      iter_invmtrx.C coulomb.C calc_Eps.C gspace.C kqindex.C\
      interpolator.C\
      fft_routines.C my_fftw.C calc_P.C do_fft_states.C\
      get_fftsize.C read_states.C\
      check_inputs.C mtrxop_3x3.C sysinfo.C usrinput.C util.C\
      main.C




OBJF= $(SRCF:%.C=%.o)

# target
gw.x: $(OBJF)
	$(CC) $(OBJF) -o $@ $(FFTINC) $(FFTLIB) $(LAPACKLIB)

.C.o: 
	$(CC) $(CFLAGS) $(MYFLAGS) $(DEBFLAGS) -c $< -o $@

clean:
	rm *.o gw.x

depend:
	$(CC) -MM $(SRCF) 



## dependency

