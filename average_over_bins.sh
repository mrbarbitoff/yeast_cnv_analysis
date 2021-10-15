#!/bin/bash

BINSIZE=$1

for i in *per_base.tsv
do
	./average_over_bins.py $i $BINSIZE > ${i%%.per_base.tsv}.bin.tsv &
done
wait

rm all.bin.tsv 2> /dev/null
cat *.bin.tsv > all.bin.tsv
