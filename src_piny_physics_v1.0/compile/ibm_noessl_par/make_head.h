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
INCLUDES = $(DCODE)/include/ibm_noessl_par
EXE      = $(ECODE)/piny_md_ibm_par.x
LIB      = $(ECODE)/libPinyInterface.a

#--------------------------
# compiler
#--------------------------
CHARMC   = $(CHARMBASE)/bin/charmc 

FC       = $(CHARMC) -f77
CC       = $(CHARMC) 
LC       = $(CHARMC)
OPT_FULL = $(OPTS)
OPT_CARE = -O2 $(DBG_FLAG) -qstrict
OPT_GRP  = -O2 $(DBG_FLAG) -qstrict
TFLAGS   = -I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib 
CFLAGS   = 
FFLAGS   = 
LIBS     =  


