####### Machine specific settings for the build #######

# Location of the charm installation
CHARMBASE     = $(HOME)/charm
# Location of the FFTW library installation
FFT_HOME	  = $(HOME)/fftw


# Extra flags for compiling and linking the whole code. 
# Note that most necessary flags are specified in other makefiles. Only 
# mention any machine speciic stuff here.
# If you realllly want, you can clear and reset the flags here or on the command line
#---------------------------------------------------------------
CPPFLAGS     += -DFORTRANUNDERSCORE 
CXXFLAGS     += -mdynamic-no-pic -mpowerpc-gpopt -mtune=G5 -mcpu=G5 -mpowerpc64
CFLAGS       +=
FFLAGS       += -qstrict -qextname -funroll-all-loops -fsched-interblock \
                -falign-loops=16 -falign-jumps=16 -falign-functions=16 \
				-falign-jumps-max-skip=15 -falign-loops-max-skip=15 \
				-force_cpusubtype_ALL -ffast-math
LDFLAGS      += -memory gnu
LDLIBS       += -lm


# Options related to math routines
#---------------------------------------------------------------
# Should we use dual fft or not
DUAL_FFTW     = -DDUAL_FFTW_OFF
# Extra math libraries to be linked in if DUAL_FFTW is off
MATH_LIB      =
# Where the linker should find these extra math libraries
MATH_LIB_PATH = -L$(MKL_HOME)/lib/em64t
# Which math library sources (that OpenAtom lugs around) need to be compiled
src_math      = $(src_blas) $(src_lapack) $(src_eispack) 
# Special optimization options to used for compiling fastadd.C
fastadd.o:        CXXFLAGS  +=
# Should we pass -D_IBM_ESSL_ as a preprocessor flag when building ibm_essl_dummy.o
ibm_essl_dummy.o: CPPFLAGS  +=


# Options related to the physics module
#---------------------------------------------------------------
# Where should we look for standard_include.h
STANDARD_INC  = $(physics)/include/pentium_par
# What flags do we use when compiling the fragile portions of piny
$(libphysics): OPT_CARE      = -O2
# Flags for compiling the fortran portions of the physics code
$(libphysics): FFLAGS       += -fno-second-underscore

