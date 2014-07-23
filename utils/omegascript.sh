#!/bin/sh
cd $PBS_O_WORKDIR
~/gennodelist.pl $PBS_NODEFILE $PBS_JOBID
NODEFILE=.nodelist.$PBS_JOBID
module load Libraries/FFTW/2.1.5
module load Compilers/Intel/icsxe/2013.1.046
./tidy water; ../../build-O3/charmrun ++nodelist $NODEFILE ../../build-O3/OpenAtom regression/cpaimd_config.p1 regression/water.input.min.ees-nl1l1 +p1
rm $NODEFILE
