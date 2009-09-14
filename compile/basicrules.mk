# Setp the basic compiler variables, flags, implicit rules etc
########### This stuff should be able take care of itself ############

# Required directory structure within the charm installation
CHARMBIN  = $(CHARMBASE)/bin
CHARMINC  = $(CHARMBASE)/include

# Compilers
CC        = $(CHARMBIN)/charmc
CXX       = $(CHARMBIN)/charmc
FC        = $(CHARMBIN)/charmc
DOXYGEN   = doxygen

# Basic compiler/linker flags
CFLAGS   += -language charm++
CXXFLAGS += -language charm++
FFLAGS   += -f77

####### Pattern rules
# Rule to generate dependency information for C++ source files
%.d: %.C
	$(info Generating dependencies for $<)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMBIN)/dep.pl $(CHARMINC) $(DEPSTRIPDIRS) \
	| sed 's|$*.o[ :]*|$*.o $@ : |g' > $@

%.d: %.cpp
	$(info Generating dependencies for $<)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMBIN)/dep.pl $(CHARMINC) $(DEPSTRIPDIRS) \
	| sed 's|$*.o[ :]*|$*.o $@ : |g' > $@

%.d: %.cxx
	$(info Generating dependencies for $<)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMBIN)/dep.pl $(CHARMINC) $(DEPSTRIPDIRS) \
	| sed 's|$*.o[ :]*|$*.o $@ : |g' > $@

# Rule to generate dependency information for C source files
%.d: %.c
	$(info Generating dependencies for $<)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMBIN)/dep.pl $(CHARMINC) $(DEPSTRIPDIRS) \
	| sed 's|$*.o[ :]*|$*.o $@ : |g' > $@

# Rule to generate dependency info for charm++ interface (ci) definition files
# @note: Need to handle pathological cases like multi-line module declarations 
%.di: %.ci
	$(info Generating dependencies for $<)
	@grep -oE "(extern[ ]+)?module[ ]+\w+" $< | \
	awk ' function printExternDeps(nExt,module,externs) { if (nExt>0) { printf "%s.decl.h: ",module; for (i=0;i<nExt;i++) printf "%s.decl.h ",externs[i]; printf "\n\n" } }   { if ($$1 ~ /extern/) { externs[nExt++] = $$3 } else { printExternDeps(nExt,modules[cnt-1],externs); nExt=0; modules[cnt++] = $$2 } }    END { printExternDeps(nExt,modules[cnt-1],externs); for (i=0;i<cnt;i++) printf "%s.decl.h %s.def.h ",modules[i],modules[i]; printf ": $<\n\t$(CXX) -c $$<\n\n" }' > $@


# Pattern rules copied from the built-in make rules
%.o: %.C
	$(COMPILE.C) $(OUTPUT_OPTION) $<

%.o: %.cpp
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<

%.o: %.cc
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

%.o: %.c
	$(COMPILE.c) $(OUTPUT_OPTION) $<

%.o: %.f
	$(COMPILE.f) $(OUTPUT_OPTION) $<

%.o: %.F
	$(COMPILE.F) $(OUTPUT_OPTION) $<

%.f: %.F
	$(PREPROCESS.F) $(OUTPUT_OPTION) $<

