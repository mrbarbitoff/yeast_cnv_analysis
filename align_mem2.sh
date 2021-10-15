#!/bin/bash

for i in fastqs/*R1*fastq.gz
do
	TAG=${i##fastqs/}
	bwa mem -t 40 -R "@RG\tID:${TAG%%.R1*}\tSM:${TAG%%.R1.*}\tLB:1\tPL:ILLUMINA" /media/array/yeast_proj/chromosomal_mutatnts/74D_merged_sorted.fasta $i ${i/R1/R2} | samtools view -bS - > bams/${TAG%%.R1*}.bam
done
