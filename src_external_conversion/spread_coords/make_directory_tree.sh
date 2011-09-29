#==========================================
# $1 is head directory $1 is number of beads-1
#==========================================
#!/bin/bash
#
mkdir -p $1
for i in `seq 0 $2`;
 do
   mkdir $1/Bead."$i"_Temper.0
done
#
#==========================================
