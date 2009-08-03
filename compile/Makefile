#===============================================================================
# include definitions
#-------------------------------------------------------------------------------
include ./Makefile.config

#===============================================================================
# general definitions
#-------------------------------------------------------------------------------
BINDIR = ../binary
TARGET = $(BINDIR)/OpenAtom
ifeq ($(DUAL_FFTW), -DDUAL_FFTW_OFF)
  FFT_INC = -lrfftw -lfftw
else
  FFT_INC = -ldrfftw -ldfftw
endif
    
#===============================================================================
all: full_code

#===============================================================================
docs:
	(cd doc; make;)

#===============================================================================
setup: 
	mkdir -p $(BINDIR)

cpaimd_driver: setup
	(cd $(MAKE_CPAIMD_DRV_INC); make driverlib OPTS="$(OPTS)";)

#===============================================================================
physics: setup
	test -n "$(MAKE_PHYS_INC)" && cd $(MAKE_PHYS_INC) &&  $(MAKE) PINYLIB 

#===============================================================================
doxygen:
	doxygen ../doc/Doxyfile

#===============================================================================
mymathlib: setup
	(cd $(MYMATH_HOME_INC); make OPTS="$(OPTS)";)

#===============================================================================
full_code: setup physics cpaimd_driver mymathlib
	cd $(CPAIMD)/compile; \
	$(CHARMC) -o $(TARGET) \
		  -language charm++ $(CHARMLIB) \
                  -L$(BINDIR) -lCharmDriver -lPinyInterface -lMyMathLib \
                  -L$(FFT_HOME)/lib $(FFT_INC) -lconv-util \
		  $(MATH_LIB) $(LNK_OPTS)
	cp charmrun $(BINDIR)

#===============================================================================
proj_code: setup physics cpaimd_driver mymathlib
	cd $(CPAIMD)/compile; \
	$(CHARMC) -o $(TARGET).prj \
		  -language charm++ $(CHARMLIB) \
                  -L$(BINDIR) -lCharmDriver -lPinyInterface -lMyMathLib \
	          -tracemode projections -lconv-util \
		  -L$(FFT_HOME)/lib $(FFT_INC) $(MATH_LIB) $(LNK_OPTS)

#===============================================================================
sum_code: setup physics cpaimd_driver mymathlib
	cd $(CPAIMD)/compile; \
	$(CHARMC) -o $(TARGET).sum \
		  -language charm++ $(CHARMLIB) \
                  -L$(BINDIR) -lCharmDriver -lPinyInterface -lMyMathLib \
	          -tracemode summary -lconv-util \
		  -L$(FFT_HOME)/lib -lrfftw -lfftw $(MATH_LIB) $(LNK_OPTS)

#===============================================================================
clean: 	clean_cpaimd_driver clean_physics clean_mymathlib
	/bin/rm -rf $(TARGET) $(TARGET).prj
#===============================================================================
#	/bin/rm -rf $(BINDIR)
#===============================================================================

clean_cpaimd_driver:
	cd $(MAKE_CPAIMD_DRV_INC); make clean; 

clean_physics:
	cd $(MAKE_PHYS_INC); make clean; 

clean_mymathlib:
	cd $(MYMATH_HOME_INC); make clean; 

valid_water: full_code
	(cd $(VALID_DIR); ./tidy water; $(RUNCOMMAND) $(AUTOTEST)/cpaimd_config.water $(AUTOTEST)/water.input >valid.result)
	perl $(VALID_EXEC)/valid.pl $(VALID_DIR)/ANSWERS_TRUE_PHYSICS $(VALID_DIR)/valid.result 

#===============================================================================

