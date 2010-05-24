#===============================================================================
#------------------------- Path to Charm++ and FFTW ----------------------------
  # Location of the charm installation
  CHARMBASE	= $(HOME)/charm
  # Location of the FFTW library installation
  FFT_HOME	= /soft/apps/fftw-2.1.5-double

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
  OPT       = -O3
  # What flags do we use when compiling the fragile portions of piny
  OPT_CARE  = -O2 -qstrict
  CPPFLAGS += $(DUAL_FFTW) -DFORTRANUNDERSCORE_OFF -DCMK_OPTIMIZE=1 \
	      -I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib 
  FFLAGS   += $(OPT)
  CFLAGS   += $(OPT)
  CXXFLAGS += $(OPT) -qarch=450 -qtune=450
  # Just to remove dcmf.h (included elsewhere) from the list of dependencies
  INCDIRS  += /bgsys/drivers/ppcfloor/comm/include
  DEPSTRIPDIRS += /bgsys/drivers/ppcfloor/comm/include

#-------------------------------------------------------------------------------
#------------------------------ Flags for linking ------------------------------
  LDFLAGS  += -L$(FFT_HOME)/lib -L/soft/apps/LAPACK \
	      -L/soft/apps/ibmcmp/xlf/bg/11.1/lib \
	      -L/bgsys/ibm_essl/sles10/prod/opt/ibmmath/lib \
	      -L/bgsys/drivers/ppcfloor/gnu-linux/powerpc-bgp-linux/lib
  LDLIBS   += -module CkMulticast -module comlib -lconv-util \
	      -lesslbg -lmass -lmassv -lxlfmath \
	      -lxlf90_r -lxlomp_ser -lrt -lpthread \
		  xerbla.o -lm -lz

# @note: Empty target specific appends (+=) hide previous global values for
#        make version 3.80 on bgp. Hence uncomment any of the following 
#        target-specific variables settings only if you actually want to 
#        modify the content.

#-------------------------------------------------------------------------------
#----------------- Flags and settings just for the driver code -----------------
$(libdriver):  CPPFLAGS += -I. -I$(driver) -I$(base) -I$(base)/include \
			   -I$(STANDARD_INC)
#$(libdriver):  FFLAGS   +=
#$(libdriver):  CFLAGS   +=
#$(libdriver):  CXXFLAGS +=

#-------------------------------------------------------------------------------
#----------------- Flags and settings just for the physics code ----------------
$(libphysics): CPPFLAGS += -I$(STANDARD_INC)
#$(libphysics): FFLAGS   +=
#$(libphysics): CFLAGS   +=
#$(libphysics): CXXFLAGS +=

# Where should we look for standard_include.h
STANDARD_INC             = $(physics)/include/ibm_noessl_par

#-------------------------------------------------------------------------------
#------------------ Flags and settings just for the math libs ------------------
#$(libmath):    CPPFLAGS +=
#$(libmath):    FFLAGS   +=
$(libmath):    CFLAGS   += -seq
$(libmath):    CXXFLAGS += -seq

#===============================================================================
#-------------------------------------------------------------------------------
  # Should we use dual fft or not
  DUAL_FFTW	= -DDUAL_FFTW_OFF
  # Which math library sources (that OpenAtom lugs around) need to be compiled
  src_math	= $(src_xerbla) $(src_eispack)
  # Special optimization options to used for compiling fastadd.C
  fastadd.o:     CXXFLAGS += -O3 -qarch=450d -qhot-simd -qhot=vector
  # Should we pass -D_IBM_ESSL_ as a preprocessor flag when building
  # ibm_essl_dummy.o
  ibm_essl_dummy.o: CPPFLAGS  += -D_IBM_ESSL_
  # The fft (and other) math libraries to link based on whether DUAL_FFTW is
  # turned on or off
  ifeq ($(DUAL_FFTW), -DDUAL_FFTW_OFF)
    LDLIBS   += -lrfftw -lfftw
  else
    LDLIBS   += -ldrfftw -ldfftw
  endif

#===============================================================================
#------------- Variables to control how scaling runs are set up --------------#
# The min num of cores to run on
procStart      = 64
# The max num of cores to run on
procEnd        = 4096
# The wall time to request for each run
walltime       = 15
# The number of cores (used) per node
ppn            = 4

#-------------------------------------------------------------------------------
# The location of the molecule database
perfDatabase = $(molDbase)
# The location of the dataset to use for the scaling runs
perfDataset = $(data)/water_32M_70Ry
# The location of the input files within the dataset
perfInpDir = $(perfDataset)/scaling/bgp

#-------------------------------------------------------------------------------
# The following placeholders can be used in anything that is part of the job
# submission command
# @C - num cores
# @N - num nodes
# @T - wall time
# @X - eXecutable

#-------------------------------------------------------------------------------
# The parallel config file to use 
perfConfig  = cpaimd_config.topo.@C
# The simulation config file to use 
perfInput   = water.input
# Any other arguments required. Placeholders can be used here too
perfOtherArgs   = 
# The command line to submit to the BGP scheduler. All arguments will be appended 
submitLine     = qsub -A CharmRTS --mode vn -n @N -t @T -O proc@C @X

