#======================================================================
# include definitions
include ./Makefile.config.gen
include ./Makefile.config.piny
include ./Makefile.config.charm

#======================================================================
# general definitions

TARGET = binary/cpaimd

#======================================================================

all: full_code

#======================================================================
docs:
	(cd doc; make;)

#======================================================================
cpaimd_driver:	
	(cd $(MAKE_CPAIMD_DRV_INC); make driverlib OPTS="$(OPTS)";)

#======================================================================
physics:
	test -n "$(MAKE_PHYS_INC)" && cd $(MAKE_PHYS_INC) &&  $(MAKE) PINYLIB 

#======================================================================
doxygen:
	doxygen ../doc/Doxyfile

#======================================================================

mymathlib:
	(cd $(MYMATH_HOME); make OPTS="$(OPTS)";)

#======================================================================
full_code: physics cpaimd_driver mymathlib
	cd $(CPAIMD)/compile; \
	$(CHARMC) -o ../binary/cpaimd.x \
		  -language charm++ $(CHARMLIB) \
                  -L../binary -lCharmDriver -lPinyInterface -lMyMathLib\
                   -lconv-util -lrfftw -lfftw $(MATH_LIB)  


#======================================================================
#======================================================================
proj_code: physics cpaimd_driver mymathlib
	cd $(CPAIMD)/compile; \
	$(CHARMC) -o ../binary/cpaimd.x \
		  -language charm++ $(CHARMLIB) \
                  -L../binary -lCharmDriver -lPinyInterface -lMyMathLib\
		-tracemode summary \
                   -lconv-util -lrfftw -lfftw $(MATH_LIB)


#======================================================================
fulldp_code: physics cpaimddp_driver mymathlib
	cd $(CPAIMD)/compile; \
	$(CHARMC) -o ../binary/cpaimd_dp.x \
                  -language charm++ $(CHARMLIB) \
                  -L../binary -lCharmDPDriver -lPinyInterface -lMyMathLib\
                   -lconv-util -lrfftw -lfftw $(MATH_LIB)


#======================================================================
clean: 	clean_cpaimd_driver clean_physics clean_mymathlib
	/bin/rm -rf ../binary/cpaimd.x

cleandp: clean_cpaimddp_driver clean_physics clean_mymathlib
	/bin/rm -rf ../binary/cpaimd.x

clean_cpaimd_driver:
	cd $(MAKE_CPAIMD_DRV_INC); make clean; 

clean_cpaimddp_driver:
	cd $(MAKE_CPAIMDDP_DRV_INC); make clean; 

clean_physics:
	cd $(MAKE_PHYS_INC); make clean; 

clean_mymathlib:
	cd $(MYMATH_HOME); make clean; 

valid_water: full_code
	(cd $(VALID_DIR); ./tidy water; $(RUNCOMMAND) $(AUTOTEST)/cpaimd_config.water $(AUTOTEST)/water.input >valid.result)
	perl $(VALID_EXEC)/valid.pl $(VALID_DIR)/ANSWERS_TRUE_PHYSICS $(VALID_DIR)/valid.result 

#======================================================================


