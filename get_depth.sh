#!/bin/bash

for i in bams/*sorted.bam
do
	TAG=$( basename $i )
	samtools depth -a -Q 10 $i > coverage/${TAG%%.sorted.bam}.per_base.tsv &
done
