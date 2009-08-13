# Define the base directory of the openatom tree
base := ../..
makedir := $(base)/compile

# Include basic make rules
include $(makedir)/srcdirs.mk
include $(makedir)/basicrules.mk
include $(makedir)/Makefile.config

# Directory structure and vpath related stuff
# Overriding values set in srcdirs.mk for a shorter, cleaner paths
CPAIMD        = $(base)
DEPSTRIPDIRS += 
VPATH        += $(mathlib)

# Specify the compilation target
TARGET        = libMyMathLib.a

# The relevant source files for this project
SRC           = \
                fastadd.C altRandom.C lst_sort_clean.C ibm_essl_dummy.C \
				random_num_gen.f genmatmul.f $(src_math)
HDR           = 
OBJ           = $(addsuffix .o, $(basename $(SRC)) )
INTF          = 
#src_math     = $(src_blas) $(src_lapack) $(src_eispack) $(src_xerbla)
src_blas      = dble_blas_stuff.f lsame.f xerbla.f zdotu.f zgemm.f zgemv.f dgemm.f dgemv.f
src_xerbla    = xerbla.f
src_eispack   = math_rs.f gefa_gesl.f
src_lapack    = dble_lapack_stuff.f dopgtr.f dorg2l.f dorg2r.f dspev.f \
                dspmv.f dspr2.f dsptrd.f dsteqr.f dsterf.f

########### Specify the flags, include paths, libraries etc.###########
# Preprocessor flags (Used for all preprocessing, compilation and linking)
CPPFLAGS += 
# Compiler flags (Used for compilation and linking)
CFLAGS    = -seq $(OPTS)
CXXFLAGS  = -seq $(OPTS)
FFLAGS    = -f77 $(FOPTS)
# Add some compile flags just for these targets
fastadd.o:        CXXFLAGS += $(MATH_OPTS)
ibm_essl_dummy.o: CXXFLAGS += $(IBM_ESSL_DUMMY) 
# Linker flags (Used for all linking)
LDFLAGS  += 
# Additional libraries
LDLIBS   += 

.PHONY: all clean again 

all: $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o *.a) $(TARGET) charmrun 

again: 
	$(MAKE) clean; $(MAKE)


# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(SRC)) ) )
ifneq ($(MAKECMDGOALS),clean)
-include $(depFiles)
-include $(INTF:.ci=.di)
endif

