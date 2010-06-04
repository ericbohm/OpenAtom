# Makefile for pairCalcTests

CHARMDIR = $(HOME)/charm
CHARMC = $(CHARMDIR)/bin/charmc $(OPTS)
FCC = gfortran
FFT_HOME = $(HOME)/fftw

base := ../..
driver := ../../src_charm_driver
physics = $(base)/src_piny_physics_v1.0
math = $(base)/src_mathlib
STANDARD_INC = $(physics)/include/pentium_par

LDFLAGS  += -L$(FFT_HOME)/lib -memory gnu
LDLIBS   += -module CkMulticast -module comlib -lz -lconv-util -lm -lfftw

CPPFLAGS += -I. -I$(driver) -I$(base) -I$(base)/include -I$(physics)/include -I$(physics)/interface -I$(STANDARD_INC)
CPPFLAGS += -DFORTRANUNDERSCORE #-DCMK_OPTIMIZE=1

BINARIES = pairCalcTest

all: ${BINARIES}

pairCalcTest: configure.o interface_hand.o CLA_Matrix.o piny_malloc.o friend_lib.o InstanceController.o ckPairCalculator.o CP_State_GSpacePlane.o ortho.o PeList.o orthoBuilder.o pcCreationManager.o pcCommManager.o dgemm.o lsame.o xerbla.o fastadd.o pcBuilder.o pcSectionManager.o MapTable.o MapFile.o pairCalcTest.o
	$(CHARMC) -language charm++ $(LDFLAGS) -o $@ *.o $(LDLIBS)
	
pairCalcTest.o: pairCalcTest.C pairCalcTest.h CPcharmParaInfo.h CPcharmParaInfoGrp.h configure.h pairCalcTest.decl.h pairCalcTest.def.h gStatePlane.decl.h gStatePlane.def.h
	${CHARMC} -o $@ pairCalcTest.C $(CPPFLAGS)
	
CP_State_GSpacePlane.o: CP_State_GSpacePlane.C CP_State_GSpacePlane.h gStatePlane.decl.h gStatePlane.def.h pairCalcTest.decl.h pairCalcTest.def.h configure.h gStatePlane.decl.h #MyGStateSlab.h
	${CHARMC} -o $@ CP_State_GSpacePlane.C $(CPPFLAGS)
	
InstanceController.o: InstanceController.C InstanceController.h instanceController.decl.h instanceController.def.h startupMessages.decl.h startupMessages.def.h gStatePlane.decl.h
	${CHARMC} -o $@ InstanceController.C $(CPPFLAGS)

fastadd.o: $(math)/fastadd.C
	${CHARMC} -o $@ $(math)/fastadd.C $(CPPFLAGS)

friend_lib.o: $(physics)/friend_lib/friend_lib.C
	${CHARMC} -o $@ $(physics)/friend_lib/friend_lib.C $(CPPFLAGS)

piny_malloc.o: $(physics)/friend_lib/piny_malloc.C
	${CHARMC} -o $@ $(physics)/friend_lib/piny_malloc.C $(CPPFLAGS)

interface_hand.o: $(physics)/interface/handle/interface_hand.C
	${CHARMC} -o $@ $(physics)/interface/handle/interface_hand.C $(CPPFLAGS)

dgemm.o: $(math)/dgemm.f
	${FCC} -c $(math)/dgemm.f
	
lsame.o: $(math)/lsame.f
	${FCC} -c $(math)/lsame.f
	
xerbla.o: $(math)/xerbla.f
	${FCC} -c $(math)/xerbla.f
	
configure.o: configure.h configure.C
	${CHARMC} -o $@ configure.C $(CPPFLAGS)
	
pcSectionManager.o: $(driver)/orthog_ctrl/pcSectionManager.C $(driver)/orthog_ctrl/pcSectionManager.h
	${CHARMC} -o $@ $(driver)/orthog_ctrl/pcSectionManager.C $(CPPFLAGS)

pcBuilder.o: $(driver)/paircalc/pcBuilder.C $(driver)/paircalc/pcBuilder.h
	${CHARMC} -o $@ $(driver)/paircalc/pcBuilder.C $(CPPFLAGS)

ckPairCalculator.o: $(driver)/paircalc/ckPairCalculator.C $(driver)/paircalc/ckPairCalculator.h ckPairCalculator.decl.h ckPairCalculator.def.h
	${CHARMC} -o $@ $(driver)/paircalc/ckPairCalculator.C $(CPPFLAGS)

pcCommManager.o: $(driver)/cp_state_ctrl/pcCommManager.C $(driver)/cp_state_ctrl/pcCommManager.h
	${CHARMC} -o $@ $(driver)/cp_state_ctrl/pcCommManager.C $(CPPFLAGS)

PeList.o: $(driver)/load_balance/PeList.C $(driver)/load_balance/PeList.h
	${CHARMC} -o $@ $(driver)/load_balance/PeList.C $(CPPFLAGS)

MapTable.o: $(driver)/load_balance/MapTable.C $(driver)/load_balance/MapTable.h
	${CHARMC} -o $@ $(driver)/load_balance/MapTable.C $(CPPFLAGS)

MapFile.o: $(driver)/utility/MapFile.C $(driver)/utility/MapFile.h
	${CHARMC} -o $@ $(driver)/utility/MapFile.C $(CPPFLAGS)
		
orthoBuilder.o: $(driver)/orthog_ctrl/orthoBuilder.C $(driver)/orthog_ctrl/orthoBuilder.h
	${CHARMC} -o $@ $(driver)/orthog_ctrl/orthoBuilder.C $(CPPFLAGS)	

pcCreationManager.o: $(driver)/main/pcCreationManager.C $(driver)/main/pcCreationManager.h startupMessages.decl.h startupMessages.def.h
	${CHARMC} -o $@ $(driver)/main/pcCreationManager.C $(CPPFLAGS)	

ortho.o: $(driver)/orthog_ctrl/ortho.C $(driver)/orthog_ctrl/CLA_Matrix.h ortho.decl.h ortho.def.h
	${CHARMC} -o $@ $(driver)/orthog_ctrl/ortho.C $(CPPFLAGS)

CLA_Matrix.o: $(driver)/orthog_ctrl/CLA_Matrix.C $(driver)/orthog_ctrl/CLA_Matrix.h CLA_Matrix.decl.h CLA_Matrix.def.h ortho.decl.h ortho.def.h ckPairCalculator.decl.h ckPairCalculator.def.h
	${CHARMC} -o $@ $(driver)/orthog_ctrl/CLA_Matrix.C $(CPPFLAGS)

CLA_Matrix.decl.h CLA_Matrix.def.h:  $(driver)/orthog_ctrl/CLA_Matrix.ci
	$(CHARMC) $(driver)/orthog_ctrl/CLA_Matrix.ci

ortho.decl.h ortho.def.h:  $(driver)/orthog_ctrl/ortho.ci
	$(CHARMC) $(driver)/orthog_ctrl/ortho.ci

ckPairCalculator.decl.h ckPairCalculator.def.h:  $(driver)/paircalc/ckPairCalculator.ci
	$(CHARMC) $(driver)/paircalc/ckPairCalculator.ci

startupMessages.decl.h startupMessages.def.h:  $(driver)/main/startupMessages.ci
	$(CHARMC) $(driver)/main/startupMessages.ci

pairCalcTest.decl.h pairCalcTest.def.h instanceController.decl.h instanceController.def.h: pairCalcTest.ci
	$(CHARMC) pairCalcTest.ci
	
gSpaceDriver.def.h gSpaceDriver.decl.h gStatePlane.decl.h gStatePlane.def.h: gspace.ci
	$(CHARMC) gspace.ci

clean:
	/bin/rm -f ${BINARIES} a.out core *.o
	/bin/rm -f *.decl.h *.def.h charmrun

run: all
	/bin/rm -f *.out
	./pairCalcTest cpaimd_config.p1
