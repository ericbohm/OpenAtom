#--------------------------
# user specific directories
#--------------------------
CHARM      = ../../../../sp/charm
MYHOME     = ../../../src_piny_physics_v1.0
CODE       = $(MYHOME)
DCODE      = $(MYHOME)
ECODE      = $(MYHOME)/../binary
CHARM_HOME = $(CHARM)
FFT_HOME   = $(CHARM_HOME)/work_cpaimd/versions/fftw

#--------------------------
# general directories
#--------------------------
INCLUDES = $(DCODE)/include/ibm_par
EXE      = $(ECODE)/piny_md_ibm_par.x
LIB      = $(ECODE)/libPinyInterface.a

#--------------------------
# compiler
#--------------------------
CHARMC   = $(CHARM_HOME)/bin/charmc 

FC       = $(CHARMC) -f77
CC       = $(CHARMC) 
LC       = $(CHARMC)
OPT      = $(OPTS)
OPT_CARE = -O2 $(DBG_FLAG) -qstrict
OPT_GRP  = -O2 $(DBG_FLAG) -qstrict
CFLAGS   = 
FFLAGS   = 
LIBS     = -lm


