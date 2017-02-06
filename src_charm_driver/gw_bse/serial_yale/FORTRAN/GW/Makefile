# makefile
# Sep. 2014  minjung.kim@yale.edu
# 
# serial code
# tested on MacOSX

# last updated: July. 2015

# Executables
EXECUTABLES= epsilon gpp 

# Compile options
OPTS = -O3 -ffree-form 
# compiler (I'm using GNU compiler)
FC = gfortran $(OPTS)

# libraries
FFTWPATH = /sw/lib
FFTLIB = -L$(FFTWPATH) -lfftw -lm 
#LAPACKPATH = /sw2/lib/lapack-3.5.0/
#LAPACKLIB = -L$(LAPACKPATH) -lblas -llapack
LAPACKLIB = -framework Accelerate # Yosemite


# source files 
# I guess the order of source files is important??
SRCF_EPSILON = set_FFTsize.f90 utility.f90 readwfn.f90 \
	mtrx_op.f90 \
	kqinfo.f90 create_mtrx.f90 fft_subroutines.f90 do_fft_all.f90 \
	shifted_grid.f90 calc_P.f90 gspace.f90 calc_eps.f90 inverse.f90\
	usrinput.f90
SRCF_GPP = gpp_subroutines.f90
MOD = constant.o structures.o usrinput.o

# target
epsilon: $(MOD)  
	$(FC) $(SRCF_EPSILON) epsilon.f90 $(FFTLIB) $(LAPACKLIB) -o epsilon

gpp: $(MOD)	
	$(FC) $(SRCF_EPSILON) $(SRCF_GPP) gpp_main.f90 $(FFTLIB) $(LAPACKLIB) -o gpp


# modules
usrinput.o: constant.o
	$(FC) -c usrinput.f90
structures.o: constant.o 
	$(FC) -c structures.f90
constant.o: 
	$(FC) -c constant.f90


# Utility targets  
# .PHONY means targets are not associated with any source files
.PHONY: clean
clean:
	rm -f *.o
	rm -f *.mod
	rm -f ${EXECUTABLES} 
