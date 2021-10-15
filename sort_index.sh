#!/bin/bash

#sleep 14400

for i in bams/*.bam
do
	samtools sort -@ 40 -o ${i%%.bam}.sorted.bam $i
done

for i in bams/*.sorted.bam
do
	samtools index $i
done
