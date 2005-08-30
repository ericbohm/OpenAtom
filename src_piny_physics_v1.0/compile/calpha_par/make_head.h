#-----------------------------------
# general directories for full code
#-----------------------------------

include ../../../compile/Makefile.config.gen
include ../../../compile/Makefile.config.charm

#-------------------------------
# user specific directories piny
#-------------------------------

MYHOME     = ../../../src_piny_physics_v1.0
CODE       = $(MYHOME)
DCODE      = $(MYHOME)
ECODE      = $(MYHOME)/../binary

#------------------------------
# general directories for piny
#-----------------------------
INCLUDES = $(DCODE)/include/pentium_par
EXE      = $(ECODE)/piny_md_calpha_par.x
LIB      = $(ECODE)/libPinyInterface.a

#--------------------------
# compiler
#--------------------------
CHARMC   = $(CHARMBASE)/bin/charmc $(OPTS)

FC       = $(CHARMC) -f77
CC       = $(CHARMC) 
LC       = $(CHARMC)
OPT      = 
OPT_CARE = 
OPT_GRP  = 
TFLAGS   = -I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib 
#CFLAGS   = -DLINUX -fno-second-underscore 
#FFLAGS   = -DLINUX -fno-second-underscore 
CFLAGS   = 
FFLAGS   = 
LIBS     = -lm 


