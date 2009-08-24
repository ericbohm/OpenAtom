# Flags, include paths, libraries etc. on a per-target basis

# CPPFLAGS - Flags used for all preprocessing, compilation and linking
# FFLAGS   - Flags used for compiling and linking fortran code
# CFLAGS   - Flags used for compiling and linking C code
# CXXFLAGS - Flags used for compiling and linking C++ code
# LDFLAGS  - Flags used only for the link stage
# LDLIBS   - Extra libraries to be linked in


#--------- Flags for linking ---------#
# Linker flags (Used for all linking)
LDFLAGS  += -L$(FFT_HOME)/lib
# Additional libraries
LDLIBS   += -module CkMulticast -module comlib
ifeq ($(DUAL_FFTW), -DDUAL_FFTW_OFF)
  LDLIBS += -lrfftw -lfftw
else
  LDFLAGS+= $(MATH_LIB_PATH)
  LDLIBS += -ldrfftw -ldfftw $(MATH_LIB)
endif
LDLIBS   += -lz -lconv-util


#--------- Flags for the whole code ---------#
               CPPFLAGS += $(DUAL_FFTW) -DFORTRANUNDERSCORE -DCMK_OPTIMIZE=1
               FFLAGS   += -O3
               CFLAGS   += -O3
               CXXFLAGS += -O3


#--------- Flags just for the driver code ---------#
$(libdriver):  CPPFLAGS +=
$(libdriver):  FFLAGS   +=
$(libdriver):  CFLAGS   +=
$(libdriver):  CXXFLAGS +=


#--------- Flags just for the physics code ---------#
$(libphysics): CPPFLAGS += -I$(STANDARD_INC)
$(libphysics): FFLAGS   +=
$(libphysics): CFLAGS   += 
$(libphysics): CXXFLAGS += 


#--------- Flags just for the math libs ---------#
$(libmath):    CPPFLAGS +=
$(libmath):    FFLAGS   +=
$(libmath):    CFLAGS   += -seq
$(libmath):    CXXFLAGS += -seq

