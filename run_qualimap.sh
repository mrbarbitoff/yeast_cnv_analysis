#!/bin/bash


for i in bams/*.IV.sorted.bam
do
	TAG=$( basename $i )
	qualimap bamqc -bam $i -c -outdir qmap_IV/${TAG%%.sorted.bam}
done
