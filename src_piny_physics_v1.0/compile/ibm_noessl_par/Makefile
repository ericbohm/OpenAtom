#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================
# MAKEFILE FOR PI_MD
#
# This makefile can be used from anywhere NFSd to the directory.
# The objs defined in ../make_defs/obj_files.h, are placed in the same 
# directory as the Makefile. You may build the executable by typing ``make''.
# You may delete all of the object files by typing "make clean"
# You may delete some of the object files by typing "make clean_care"
# or "make clean_group". Compilation instructions are in 
# ../make_defs/make_*.h files and make_spec.h. Machine specific routines 
# are in spec_files.h.  Machine specific directory tree and compiler options 
# are given in make_head.h.
#
#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================


#==========================
# General unix instructions
#==========================
ECHO    = /bin/echo making
RM      = /bin/rm -f
CP      = /bin/cp
MV      = /bin/mv -f
TOUCH   = touch

#==========================
# Machine specific info
#==========================
include make_head.h

#==========================
# Compiler commands
#==========================
FOBJ            = $(FC) $(FFLAGS) $(OPT_FULL) -I$(INCLUDES) -c -o $@
COBJ            = $(CC) $(CFLAGS) $(OPT_FULL) -I$(INCLUDES) $(TFLAGS) -c -o $@
COBJ_TEST       = $(CC) $(CFLAGS) $(OPT_FULL) -I$(INCLUDES) $(TFLAGS) -c -o $@
COBJ_FUSS       = $(CC) $(CFLAGS) $(DBG_FLAG) -I$(INCLUDES) $(TFLAGS) -c -o $@
COBJ_NOOPT      = $(CC) $(CFLAGS) $(DBG_FLAG) -I$(INCLUDES) $(TFLAGS) -c -o $@
COBJ_CARE       = $(CC) $(CFLAGS) $(OPT_CARE) -I$(INCLUDES) $(TFLAGS) -c -o $@
COBJ_GRP        = $(CC) $(CFLAGS) $(OPT_GRP)  -I$(INCLUDES) $(TFLAGS) -c -o $@
COBJ_DECL       = $(CC) $(CFLAGS) 
COBJ_FUSS_DECL  = $(CC) $(CFLAGS) 
COBJ_NOOPT_DECL = $(CC) $(CFLAGS) 
COBJ_CARE_DECL  = $(CC) $(CFLAGS) 
COBJ_GRP_DECL   = $(CC) $(CFLAGS) 
LD              = $(LC) $(OPT_FULL) $(CFLAGS) 

#==========================
# Make the executable 
#==========================
include spec_files.h
include ../make_defs/obj_files.h

$(EXE) : $(OBJS_COMP)
	$(ECHO) $(EXE)
	mkdir -p $(ECODE)
	$(RM) $(EXE)
	$(LD)  $(OBJS) -o  $(EXE) $(LIBS)

PINYLIB: $(LIB)

$(LIB) : $(LIB_OBJS_COMP)
	$(ECHO) $(LIB)
	mkdir -p $(ECODE)
	$(RM) $(LIB)
#	ar -qc $(LIB) $(LIB_OBJS) 
	$(LD)  $(LIB_OBJS) -o  $(LIB) $(LIBS)

#==========================
# Compilation instructions
#==========================
include ../make_defs/proto_files.h
include ../make_defs/make_friend.h
include ../make_defs/make_interface.h
include ../make_defs/make_interface_cp.h
include ../make_defs/make_interface_intra.h
include ../make_defs/make_interface_mol.h
include ../make_defs/make_main.h
include ../make_defs/make_math.h
include ../make_defs/make_abinitio_physics.h
include ../make_defs/make_classical_physics.h
include make_spec.h

#==========================
#remove all objs and libs
#==========================
clean       :
	/bin/rm -f $(LIB) $(LIB_DECLS) $(LIB_OBJS) $(DECLS) $(OBJS) $(ULIBS) *.o charmrun 

#==========================
#remove care objs and libs
#==========================
clean_care  :
	/bin/rm -f  $(DECLS) $(OBJS_CARE) $(ULIBS) *.o  charmrun

#==========================
#remove group objs and libs
#==========================
clean_group :
	/bin/rm -f $(DECLS) $(OBJS_GRP) $(ULIBS) *.o charmrun

#==========================
#remove group objs and libs
#==========================
clean_piny_to_driver :
	/bin/rm -f $(DECLS) $(LIB_FILES)

#===========================================================================
