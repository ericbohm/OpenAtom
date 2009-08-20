# Directory hierarchy for the openatom sources

# variable 'base' will need to be defined by any makefile that includes this, and then the 
# complete source tree structure will be available thro the following variables

################## Modify only if directory structure changes ####################
# Top-level directories
makedir   = $(base)/compile
driver    = $(base)/src_charm_driver
physics   = $(base)/src_piny_physics_v1.0
mathlib   = $(base)/src_mathlib
topinclude= $(base)/include
docs      = $(base)/doc

# Build directories
build     = $(base)/binary
build_driver = $(build)
build_physics = $(build)
build_mathlib= $(build)

# Charm driver directory structure
density   = $(driver)/cp_density_ctrl
state     = $(driver)/cp_state_ctrl
largesp   = $(driver)/cp_largesp_ctrl
fftslab   = $(driver)/fft_slab_ctrl
loadbal   = $(driver)/load_balance
main      = $(driver)/main
ortho     = $(driver)/orthog_ctrl
paircalc  = $(driver)/paircalc
strfact   = $(driver)/structure_factor
uber      = $(driver)/uber
util      = $(driver)/utility
alldriverdirs = $(main) $(density) $(state) $(fftslab) $(loadbal) $(ortho) $(paircalc) $(strfact) $(uber) $(util) $(largesp)

# PINY physics directory structure
pinymake  = $(physics)/compile
pinyinc   = $(physics)/include

# List of modules and libraries to be built (yuck,,,)
moddriver  := CharmDriver
modphysics := PinyInterface
modmath    := MyMathLib
libdriver  := lib$(moddriver).a
libphysics := lib$(modphysics).a
libmath    := lib$(modmath).a

