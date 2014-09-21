# Location of the charm installation
CHARMBASE     = $(HOME)/charm/pamilrts-bluegeneq-async-smp-xlc-prodz
#CHARMBASE     = $(HOME)/charm/mpi-bluegeneq-prod
# Location of the FFTW library installation
FFT_HOME      = /soft/libraries/alcf/current/xl/FFTW2/


#---------------------------------------------------------------
# Flags, include paths, libraries etc. on a per-target basis

# CPPFLAGS - Flags used for all preprocessing, compilation and linking
# FFLAGS   - Flags used for compiling and linking fortran code
# CFLAGS   - Flags used for compiling and linking C code
# CXXFLAGS - Flags used for compiling and linking C++ code
# LDFLAGS  - Flags used only for the link stage
# LDLIBS   - Extra libraries to be linked in


#---------------------------------------------------------------
#--------- Flags for the whole code ---------#
               # Optimization level and debug (Dont add other flags to OPT)
               OPT       = -O3
               # What flags do we use when compiling the fragile portions of piny
               OPT_CARE  = -O2 -DFORTRANUNDERSCORE_OFF
               CPPFLAGS += $(DUAL_FFTW) -DFORTRANUNDERSCORE_OFF -I$(FFT_HOME)/include
               FFLAGS   += $(OPT)
               CFLAGS   += $(OPT)
               CXXFLAGS += $(OPT) 
               # Just to remove dcmf.h (which gets tacked on somewhere) from the list of dependencies
               INCDIRS  += /bgsys/drivers/ppcfloor/comm/sys-fast/include
               INCDIRS += /bgsys/drivers/ppcfloor/spi/include/kernel/cnk
               DEPSTRIPDIRS += /bgsys/drivers/ppcfloor/comm/sys-fast/include
               DEPSTRIPDIRS  += /bgsys/drivers/ppcfloor
               DEPSTRIPDIRS  += spi/include/kernel
               DEPSTRIPDIRS  += /bgsys/drivers/ppcfloor/spi
               DEPSTRIPDIRS += /bgsys/drivers/ppcfloor/spi/include/kernel/cnk
               DEPSTRIPILES += kernel_impl.h



#---------------------------------------------------------------
#--------- Flags for linking ---------#
              LDFLAGS  += -L$(FFT_HOME)/lib -L/soft/libraries/alcf/current/xl/ZLIB/lib \
                        -L/soft/compilers/ibmcmp-may2014/xlf/bg/14.1/bglib64 \
                         -L/bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux/lib \
                         -L/soft/libraries/essl/5.1.1-1/essl/5.1/lib64
              LDLIBS   += -module CkMulticast -module comlib -lconv-util -lz -lesslbg -lesslsmpbg -lm -lxlfmath -lxlf90_r -lxl -lpthread -lrt -ldl

# @note: Empty target specific appends (+=) hide previous global values for
#        make version 3.80 on bgp. Hence uncomment any of the following 
#        target-specific variables settings only if you actually want to 
#        modify the content.

#---------------------------------------------------------------
#--------- Flags and settings just for the driver code ---------#
$(libdriver):  CPPFLAGS += -I. -I$(driver) -I$(base) -I$(base)/include -I$(STANDARD_INC)
#$(libdriver):  FFLAGS   +=
#$(libdriver):  CFLAGS   +=
#$(libdriver):  CXXFLAGS +=


#---------------------------------------------------------------
#--------- Flags and settings just for the physics code ---------#
$(libphysics): CPPFLAGS += -I$(STANDARD_INC) -DFORTRANUNDERSCORE_OFF
#$(libphysics): FFLAGS   += -Mnosecond_underscore
$(libphysics): FFLAGS   += -fno-second-underscore -fno-underscoring
#$(libphysics): CFLAGS   +=
#$(libphysics): CXXFLAGS +=

# Where should we look for standard_include.h
STANDARD_INC             = $(physics)/include/ibm_par


#---------------------------------------------------------------
#--------- Flags  and settings just for the math libs ---------#
#$(libmath):    CPPFLAGS +=
$(libmath):    FFLAGS   += -fno-second-underscore -fno-underscoring
$(libmath):    CFLAGS   += -seq -DFORTRANUNDERSCORE_OFF=1
$(libmath):    CXXFLAGS += -seq -DFORTRANUNDERSCORE_OFF=1

# Should we use dual fft or not
DUAL_FFTW                = -DDUAL_FFTW
#DUAL_FFTW                = 
# Which math library sources (that OpenAtom lugs around) need to be compiled
#src_math                 = $(src_xerbla) $(src_eispack) $(src_lapack) $(src_blas)
src_math                 = $(src_xerbla) $(src_eispack)
# Special optimization options to used for compiling fastadd.C
fastadd.o:     CXXFLAGS += -O3
# Should we pass -D_IBM_ESSL_ as a preprocessor flag when building ibm_essl_dummy.o
ibm_essl_dummy.o: CPPFLAGS  += -D_IBM_ESSL_
# The fft (and other) math libraries to link based on whether DUAL_FFTW is turned on or off
          ifeq ($(DUAL_FFTW), -DDUAL_FFTW_OFF)
               LDLIBS   += -lrfftw -lfftw
          else
               LDLIBS   += -ldrfftw -ldfftw
          endif

