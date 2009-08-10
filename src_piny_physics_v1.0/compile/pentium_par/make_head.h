#-----------------------------------
# general directories for full code
#-----------------------------------
include ../../../compile/Makefile.config

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
EXE      = $(ECODE)/piny_md_pentium_par.x
LIB      = libPinyInterface.a

#--------------------------
# compiler
#--------------------------
CHARMC   = $(CHARMBASE)/bin/charmc 

FC       = $(CHARMC) -f77 -fno-second-underscore
CC       = $(CHARMC) 
LC       = $(CHARMC)
OPT_FULL = $(OPTS) 
OPT_CARE = -O2 $(DBG_FLAG)
OPT_GRP  = -O2 $(DBG_FLAG)
TFLAGS   = -I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib 
CFLAGS   = 
FFLAGS   =  
LIBS     = -lm 


