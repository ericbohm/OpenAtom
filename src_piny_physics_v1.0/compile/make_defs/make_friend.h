#=================================================================
#               FRIEND_FILES
#=================================================================


#=================================================================
piny_malloc.o     :      $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
                         $(CODE)/friend_lib/piny_malloc.C
	$(ECHO) $@
	$(COBJ) $(DUAL_FFTW) -I$(FFT_HOME)/include $(CODE)/friend_lib/piny_malloc.C

#-----------------------------------------------------------------
piny_pup.o     :          $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
                         $(CODE)/friend_lib/piny_pup.C
	$(ECHO) $@
	$(COBJ) $(CODE)/friend_lib/piny_pup.C

#-----------------------------------------------------------------
friend_lib.o   :         $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
                         $(CODE)/friend_lib/friend_lib.C
	$(ECHO) $@
	$(COBJ) $(CODE)/friend_lib/friend_lib.C

#=================================================================
