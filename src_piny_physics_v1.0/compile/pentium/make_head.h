#--------------------------
# user specific directories
#--------------------------
MYHOME     = $(HOME)/leanMD/PHYSICS_PINY
CHARM_HOME = /expand8/home/glenn
CODE       = $(MYHOME)
DCODE      = $(MYHOME)
ECODE      = $(DCODE)/runable
FFT_HOME   = $(CHARM_HOME)/work_cpaimd/versions/fftw

#--------------------------
# general directories
#--------------------------
INCLUDES = $(DCODE)/include/pentium
EXE      = $(ECODE)/piny_md_pentium.x

#--------------------------
# compilers
#--------------------------
CHARM_OFF = 1
CHARMC    = $(CHARM_HOME)/charm/bin/charmc $(OPTS)

FC       = g77  
CC       = g++ 
LC       = g77
OPT      = -O2 
OPT_CARE = -O2 
OPT_GRP  = -O2
CFLAGS   = -DLINUX -fno-second-underscore 
FFLAGS   = -DLINUX -fno-second-underscore 
LIBS     = -lm  -lstdc++


