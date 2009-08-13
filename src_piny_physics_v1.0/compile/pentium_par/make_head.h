VPATH        += $(physics)/include/pentium_par 
CPPFLAGS     += -I$(physics)/include/pentium_par

OPT_FULL      = $(OPTS) 
OPT_CARE      = -O2 $(DBG_FLAG)
FFLAGS       += -fno-second-underscore

LDLIBS       += -lm 

