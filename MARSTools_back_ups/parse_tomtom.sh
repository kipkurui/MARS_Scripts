#!/usr/bin/env bash
mkdir /tmp/test
tf=$1
tom_out=$2
mode=$3
for thr in 5 1 0.5 0.05 0.005 0.0005 0.00005 0.000005
	do
	tomtom -oc /tmp/test-tomtom -min-overlap 5 -dist $mode -evalue -thresh $thr -no-ssc $tf $tf
	
	#echo $thr >/tmp/test/f$thr
	for mot in $(grep "MOTIF" $tf  |cut -f2 -d" ");
	do 
		cut -f1 /tmp/test-tomtom/tomtom.txt | grep -c $mot  >>/tmp/test/f$thr
done
done

grep "MOTIF" $tf | cut -f2 -d" " >/tmp/test/000
paste /tmp/test/* >$tom_out
rm /tmp/test/*
rmdir /tmp/test/
