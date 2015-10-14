#===============================================================================
#------------------------- Path to Charm++ and FFTW ----------------------------
  # Location of the charm installation
  CHARMBASE	= $(CHARM_DIR)
  # Location of the FFTW library installation
  FFT_HOME	= /dcsdata/home/bohm/fftw64

#===============================================================================
# Flags, include paths, libraries etc. on a per-target basis

# CPPFLAGS - Flags used for all preprocessing, compilation and linking
# FFLAGS   - Flags used for compiling and linking fortran code
# CFLAGS   - Flags used for compiling and linking C code
# CXXFLAGS - Flags used for compiling and linking C++ code
# LDFLAGS  - Flags used only for the link stage
# LDLIBS   - Extra libraries to be linked in

#-------------------------------------------------------------------------------
#------------------------- Flags for the whole code ----------------------------
  # Optimization level and debug (Dont add other flags to OPT)
  OPT       = -O3 -optimize -production
  # What flags do we use when compiling the fragile portions of piny
  OPT_CARE  = -O2
  CPPFLAGS += $(DUAL_FFTW) -DFORTRANUNDERSCORE -I$(FFT_HOME)/include
  FFLAGS   += $(OPT)
  CFLAGS   += $(OPT)
  CXXFLAGS += $(OPT)

#-------------------------------------------------------------------------------
#------------------------------ Flags for linking ------------------------------
  LDFLAGS  += -L$(FFT_HOME)/lib -memory gnu
  LDLIBS   += -module CkMulticast -lz -lconv-util -lm

#-------------------------------------------------------------------------------
#----------------- Flags and settings just for the driver code -----------------
$(libdriver):  CPPFLAGS += -I. -I$(driver) -I$(base) -I$(base)/include \
			   -I$(STANDARD_INC)
$(libdriver):  FFLAGS   +=
$(libdriver):  CFLAGS   +=
$(libdriver):  CXXFLAGS +=

#-------------------------------------------------------------------------------
#----------------- Flags and settings just for the physics code ----------------
$(libphysics): CPPFLAGS += -I$(STANDARD_INC)
$(libphysics): FFLAGS   += -fno-second-underscore
$(libphysics): CFLAGS   += 
$(libphysics): CXXFLAGS += 

# Where should we look for standard_include.h
STANDARD_INC             = $(physics)/include/pentium_par

#-------------------------------------------------------------------------------
#------------------ Flags and settings just for the math libs ------------------
$(libmath):    CPPFLAGS +=
$(libmath):    FFLAGS   +=
$(libmath):    CFLAGS   += -seq
$(libmath):    CXXFLAGS += -seq

#===============================================================================
#-------------------------------------------------------------------------------
  # Should we use dual fft or not
  DUAL_FFTW	= -DDUAL_FFTW_OFF
  # Which math library sources (that OpenAtom lugs around) need to be compiled
  src_math	= $(src_blas) $(src_lapack) $(src_eispack)
  # Special optimization options to used for compiling fastadd.C
  fastadd.o:     CXXFLAGS +=
  # Should we pass -D_IBM_ESSL_ as a preprocessor flag when building 
  # ibm_essl_dummy.o
  ibm_essl_dummy.o: CPPFLAGS  +=
  # The fft (and other) math libraries to link based on whether DUAL_FFTW is
  # turned on or off
  ifeq ($(DUAL_FFTW), -DDUAL_FFTW_OFF)
    LDLIBS   += -lrfftw -lfftw
  else
    LDLIBS   += -ldrfftw -ldfftw
  endif

#===============================================================================
#------------- Variables to control how regression tests are set up --------------#
# The location of the water (32M,10Ry) dataset that will be used in the tests
w3210       = $(data)/water_32M_10Ry
