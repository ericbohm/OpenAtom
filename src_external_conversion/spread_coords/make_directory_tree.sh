#!/bin/bash
#==========================================
# $1 is head directory , $2 is number of beads-1
#==========================================
#
mkdir -p $1
for i in `seq 0 $2`;
 do
   mkdir $1/Bead."$i"_Temper.0
done
#
#==========================================
