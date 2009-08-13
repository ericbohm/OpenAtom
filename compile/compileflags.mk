####### Target and machine specific flags for compilation and linking #######

# Specify the flags, include paths, libraries etc. here on a per-target basis
# Preprocessor flags are used for all preprocessing, compilation and linking
# Compiler flags are used while compiling and linking
# Linker flags are used only for the link stage


# List of modules 
libdriver  := libCharmDriver.a
libphysics := libPinyInterface.a
libmath    := libMyMathLib.a
#--------- Flags for linking ---------#
# Linker flags (Used for all linking)
LDFLAGS  += -L$(FFT_HOME)/lib -memory gnu 
# Additional libraries
LDLIBS   += -module CkMulticast -module comlib
ifeq ($(DUAL_FFTW), -DDUAL_FFTW_OFF)
  LDLIBS += -lrfftw -lfftw
else
  LDLIBS += -ldrfftw -ldfftw
endif
LDLIBS   += -lz -lm -lconv-util


#--------- Flags for the driver ---------#
$(libdriver):  CPPFLAGS +=
$(libdriver):  FFLAGS   += $(FOPTS)
$(libdriver):  CFLAGS   += $(OPTS)
$(libdriver):  CXXFLAGS += $(OPTS) 


#--------- Flags for the physics code ---------#
$(libphysics): CPPFLAGS +=
$(libphysics): FFLAGS   += $(FOPTS) 
$(libphysics): CFLAGS   += 
$(libphysics): CXXFLAGS += 


#--------- Flags for the math libs ---------#
$(libmath):    CPPFLAGS +=
$(libmath):    FFLAGS   += $(FOPTS) 
$(libmath):    CFLAGS   += -seq $(OPTS)
$(libmath):    CXXFLAGS += -seq $(OPTS)


##### Other target-specific modifications ####
fastadd.o:        CXXFLAGS += $(MATH_OPTS)
ibm_essl_dummy.o: CXXFLAGS += $(IBM_ESSL_DUMMY) 

