# The relevant source files for this module
libmath_src   = \
                fastadd.C altRandom.C lst_sort_clean.C ibm_essl_dummy.C \
				random_num_gen.f genmatmul.f $(src_math)
libmath_obj   = $(addsuffix .o, $(basename $(libmath_src)) )
libmath_intf  = 
# The actual sources for this module, that are decided on a per-machine basis
#src_math     = $(src_blas) $(src_lapack) $(src_eispack) $(src_xerbla)
src_blas      = dble_blas_stuff.f lsame.f xerbla.f zdotu.f zgemm.f zgemv.f dgemm.f dgemv.f
src_xerbla    = xerbla.f
src_eispack   = math_rs.f gefa_gesl.f
src_lapack    = dble_lapack_stuff.f dopgtr.f dorg2l.f dorg2r.f dspev.f \
                dspmv.f dspr2.f dsptrd.f dsteqr.f dsterf.f

# Specify the list of directories whose contents should be stripped from prerequisite lists 
# during dependency generation
DEPSTRIPDIRS += 
# Add appropriate directory(ies) to the vpaths for the source file types present in this module
# so that they can be located from the build directory. This is a small effort to avoid swamping
# VPATH with a long list of directories hurting the build times that we hope to improve
fileTypes     = $(sort $(suffix $(libmath_src) $(libmath_intf) ) )
$(foreach suf, $(fileTypes), $(eval vpath %$(suf) $(mathlib)) )

# The primary target for this module
$(libmath): $(libmath_obj) 
	$(info-ar)
	$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

#Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(libmath_src)) ) )
ifneq ($(MAKECMDGOALS),clean)
-include $(depFiles)
-include $(libmath_intf:.ci=.di)
endif

