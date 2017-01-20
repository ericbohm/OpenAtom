# makefile
# Sep. 2014  minjung.kim@yale.edu
# 
# serial code
# tested on MacOSX


# Executables
EXECUTABLES= sigma ##epsilon

# Compile options
OPTS = #-O3 -ffree-form 
VERSION = -fsf-4.9
#OPTS = -ffree-form -ffree-line-length-none -fno-second-underscore
##OPTS = -g -ffree-form
# compiler (I'm using GNU compiler)
FC = gfortran$(VERSION) $(OPTS)

# libraries
FFTWPATH = /sw/lib
FFTLIB = -L$(FFTWPATH) -lfftw -lm 
#BLASPATH = /sw2/lib/lapack-3.5.0/BLAS
#BLASLIB = -L$(BLASPATH) -lblas
LAPACKPATH = /Users/subhasishmandal//Downloads/lapack-3.5.0 ###   /sw2/lib/lapack-3.5.0/
#LAPACKLIB = -L$(LAPACKPATH) -lblas -llapack
LAPACKLIB = -framework Accelerate


%.o : %.f90
	$(FC) -c $< -o $@

# source files 
# I guess the order of source files is important??
SRCF = set_FFTsize.f90 utility.f90 readwfn.f90 \
	mtrx_op.f90 \
	kqinfo.f90 create_mtrx.f90 fft_subroutines.f90 do_fft_all.f90 \
	shifted_grid.f90 umkalpp.f90 calc_v.f90 gspace.f90 calc_sig.f90 inverse.f90\
	usrinput.f90 structures.f90 constant.f90 sigma.f90
OBJS=$(patsubst %.f90, %.o, $(SRCF))
MOD = constant.o structures.o usrinput.o

# target
sigma: $(MOD) $(OBJS)
	$(FC) $(OBJS) $(FFTLIB) $(LAPACKLIB) -o sigma


# Utility targets  
# .PHONY means targets are not associated with any source files
.PHONY: clean
clean:
	rm -f *.o
	rm -f *.mod
	rm -f ${EXECUTABLES} 
