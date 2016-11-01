# Directory hierarchy for the openatom sources

# variable 'base' will need to be defined by any makefile that includes this, and then the 
# complete source tree structure will be available thro the following variables

# Utilities for path manipulation: Requires perl
PERL          = /usr/bin/perl
rel2abs       = $(shell $(PERL) -e "use File::Spec; print File::Spec->rel2abs('$(1)','$(2)');")
abs2rel       = $(shell $(PERL) -e "use File::Spec; print File::Spec->abs2rel('$(1)','$(2)');")
realpath      = $(shell $(PERL) -e "use Cwd; print Cwd::realpath('$(1)');")

################## Modify only if directory structure changes ####################
# Top-level directories
makedir       = $(base)/makefiles
driver        = $(base)/src_charm_driver
physics       = $(base)/src_piny_physics_v1.0
mathlib       = $(base)/src_mathlib
interoplib    = $(base)/src_interop
topinclude    = $(base)/include
docs          = $(base)/doc
data          = $(base)/data
tests         = $(base)/tests

# Charm driver directory structure
density       = $(driver)/cp_density_ctrl
state         = $(driver)/cp_state_ctrl
largesp       = $(driver)/cp_largesp_ctrl
fftslab       = $(driver)/fft_slab_ctrl
loadbal       = $(driver)/load_balance
main          = $(driver)/main
ortho         = $(driver)/orthog_ctrl
paircalc      = $(driver)/paircalc
strfact       = $(driver)/structure_factor
uber          = $(driver)/uber
util          = $(driver)/utility
alldriverdirs = $(main) $(density) $(state) $(fftslab) $(loadbal) $(ortho) $(paircalc) $(strfact) $(uber) $(util) $(largesp)

# PINY physics directory structure
pinymake      = $(physics)/compile
pinyinc       = $(physics)/include
abinitio      = $(physics)/abinito_physics
classical     = $(physics)/classical_physics
interface     = $(physics)/interface

abinitio_dirs = cp_nonlocal cp_local cp_xc_fnctl cp_integrate cp_ions cp_orthog
classical_dirs= integrate output vel_smpl_atms
interface_dirs= coords_cp cp_ewald search_base vel_sampl_cp vps_params parse \
                handle sim_params coords surf_params inter_params vel_sampl_class \
                lists path_integral intra_params intra_params/set_params mol_params \
				mol_params/set_params
pinysrcdirs   = $(abinitio_dirs:%=$(abinitio)/%) \
                $(classical_dirs:%=$(classical)/%) \
                $(interface_dirs:%=$(interface)/%) \
                $(physics)/piny_to_driver $(physics)/friend_lib $(physics)/mathlib

# Tests directory
testregr      = $(tests)/regression
testunit      = $(tests)/unit

# Directories usefule for regression tests
molDbase      = $(data)/DATABASE
regrInpDir    = $(w3210)/regression

# List of modules and libraries to be built (yuck,,,)
executable   := OpenAtom
moddriver    := CharmDriver
modphysics   := PinyInterface
modmath      := MyMathLib
modinterop   := Interop
libdriver    := lib$(moddriver).a
libphysics   := lib$(modphysics).a
libmath      := lib$(modmath).a
libinterop   := lib$(modinterop).a

