#!/bin/bash

filename=$1
iterfreq=$2
max=$3
let loopmax="$max - $iterfreq"
for iter in $( seq "$iterfreq" "$loopmax")
do
    grep "Iter \[$iter\]" $filename | grep -v 'Rot' | grep -v 'Exit' | sed "s/Iter \[$iter\]/Iter c/" >co
    let iternext="$iter + $iterfreq"
    echo Comparing $iter to $iternext
    grep "Iter \[$iternext\]" $filename | grep -v 'Rot' | grep -v 'Exit' | sed "s/Iter \[$iternext\]/Iter c/" >cn
    diff co cn > d
    if [ $? -ne 0 ]; then
	echo "mismatch iteration $iter vs $iternext"
	cat d
	exit;
    fi
done
