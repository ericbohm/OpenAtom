CHARMDIR = /Users/minjung/packages/Charm

FLAGS = #-D_CP_DEBUG_UTIL_VERBOSE_  -DDEBUG_FFT

CHARMC = $(CHARMDIR)/bin/charmc $(OPTS) -g $(FLAGS)

FFTPATH = /sw/lib
FFTLIB = -L$(FFTPATH) -lfftw3 -lm
FFTINC = -I/sw/include


default: all
all: states


states : my_fftw.o fft_routines.o piny_pup.o main.o states.o friend_lib.o interface_hand.o piny_malloc.o configure_gwbse.o
	$(CHARMC) -language charm++ -o states main.o states.o my_fftw.o fft_routines.o  friend_lib.o interface_hand.o piny_malloc.o piny_pup.o configure_gwbse.o -lz $(FFTINC) $(FFTLIB)

main.o : main.C piny_pup.C main.h main.decl.h main.def.h states.decl.h 
	$(CHARMC) -o main.o main.C

friend_lib.o : friend_lib.C
	$(CHARMC) -o friend_lib.o friend_lib.C

interface_hand.o : interface_hand.C
	$(CHARMC) -o interface_hand.o interface_hand.C

piny_malloc.o : piny_malloc.C
	$(CHARMC) -o piny_malloc.o piny_malloc.C

piny_pup.o : piny_pup.C
	$(CHARMC) -o piny_pup.o piny_pup.C

configure_gwbse.o :configure_gwbse.C 
	$(CHARMC) -o configure_gwbse.o configure_gwbse.C

main.decl.h main.def.h : main.ci
	$(CHARMC) main.ci

my_fftw.o: my_fftw.C my_fftw.h
	$(CHARMC) -o my_fftw.o my_fftw.C $(FFTINC) $(FFTLIB)

fft_routines.o: fft_routines.C fft_routines.h
	$(CHARMC) -o fft_routines.o fft_routines.C


states.o : states.C states.h my_fftw.h fft_routines.h states.decl.h states.def.h main.decl.h
	$(CHARMC) -o states.o states.C


states.decl.h states.def.h : states.ci
	$(CHARMC) states.ci


clean:
	rm -f main.decl.h main.def.h main.o
	rm -f states.decl.h states.def.h states.o
	rm -f configure_gwbse.o friend_lib.o interface_hand.o piny_malloc.o piny_pup.o
	rm -f fft_routines.o my_fftw.o
	rm -f states charmrun

