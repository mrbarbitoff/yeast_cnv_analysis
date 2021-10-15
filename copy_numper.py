#!/usr/bin/env python3
import argparse
import numpy as np
import subprocess
import re
import itertools
import statistics

parser = argparse.ArgumentParser(description='Files for alignment')
parser.add_argument('-ref', type=str, help='Reference alignment')
parser.add_argument('-fwr',  type=str, help='Forward reads')
parser.add_argument('-rwr', type=str, help='Reverse reads')
parser.add_argument('-tr', type=str, help='Transgen scaffolds')
parser.add_argument('-d', type=float, help='Depth')
parser.add_argument('-threads', type=int, help='Number of threads')
args = parser.parse_args()    
transgen_head = subprocess.call(f"minimap2 --secondary=no -ax sr -t {args.threads} {args.tr} {args.fwr} {args.rwr} | awk '$6 !~ /\*/' > aln.sam", shell=True)
transgen_dict={}
with open("aln.sam") as file:
    for line in file:
        if line.startswith('@'):
            continue
        tr_list = line.split("\t")	
        iD = tr_list[0]
        Flag = tr_list[1]
        if bin(int(Flag))[-7]=='1':
            iD = (f"{iD}_1")
        else:
            iD = (f"{iD}_2")    
        Chromosome = tr_list[2]
        Coordinate = tr_list[3]
        CIGAR = tr_list[5]
        transgen_dict[iD]=[Chromosome, Flag, Coordinate, CIGAR]
with open("QNAMES.txt", "w") as f:
    for iD in transgen_dict:
         f.write(f"{iD[:-2]}\n")
subset_1 = subprocess.call(f"~/bbmap/filterbyname.sh -Xmx60g in={args.fwr} out=filter_fwr.fastq.gz names=QNAMES.txt include=t", shell=True)
subset_2 = subprocess.call(f"~/bbmap/filterbyname.sh -Xmx60g in={args.rwr} out=filter_rwr.fastq.gz names=QNAMES.txt include=t", shell=True)
reference_head = subprocess.call(f"minimap2 --secondary=no -ax sr -t {args.threads} {args.ref} filter_fwr.fastq.gz filter_rwr.fastq.gz | awk '$6 !~ /\*/' > ref.sam", shell=True)
reference_dict={}
with open("ref.sam") as file:
    for line in file:
        if line.startswith('@'):
            continue
        ref_list = line.split("\t")
        ID = ref_list[0]
        flag = ref_list[1]
        if bin(int(flag))[-7]=='1':
            ID = (f"{ID}_1")
        else:
            ID = (f"{ID}_2")
        chromosome = ref_list[2]
        coordinate = ref_list[3]
        cigar = ref_list[5]
        reference_dict[ID]=[chromosome, flag, coordinate, cigar]
with open("aln.sam") as tr_file, open("copy_ref_to_transgen.sam", "w") as copy_file:   
    for tline in tr_file:
        if tline.startswith('@'):
            copy_file.write(f"{tline}")
            continue
        tr_list = tline.split("\t")
        tr_ID = tr_list[0]
        flag = tr_list[1]
        if bin(int(flag))[-7]=='1':
            tr_ID = (f"{tr_ID}_1")
        else:
            tr_ID = (f"{tr_ID}_2")
        if tr_ID in reference_dict:
            copy_file.write(f"{tline}")
ref_bam = subprocess.call("samtools view -S -b copy_ref_to_transgen.sam > copy_ref_to_transgen.bam", shell=True)
coordinates_rt = subprocess.call("samtools view copy_ref_to_transgen.bam | bedtools bamtobed -i copy_ref_to_transgen.bam > Co_ref_to_transgen.bed", shell=True)
sort_rt_coordinate = subprocess.call("sort -k1,1 -k2,2n Co_ref_to_transgen.bed > sCo_ref_to_transgen.bed", shell=True)
merge_rt_coordinates = subprocess.call("bedtools merge -i sCo_ref_to_transgen.bed > Comerge_ref.bed", shell=True)
transgen_bam = subprocess.call("samtools view -S -b aln.sam > Coaln.bam", shell=True)
transgen_bam_sort = subprocess.call("samtools sort Coaln.bam > sCoaln.bam" , shell=True)
transgen_depth = subprocess.call("samtools depth sCoaln.bam > transgen_depth.txt", shell=True)
depth_to_bed = subprocess.call("awk 'BEGIN { FS=\"\t\"; OFS=\"\t\"; }{print $1, $2-1, $2, $3}' transgen_depth.txt | sort -k1,1 -k2,2n > coordinates_tr_depth.bed", shell=True)
filter_intersection = subprocess.call("bedtools intersect -a coordinates_tr_depth.bed -b Comerge_ref.bed -v > only_transgen_reads.txt", shell=True)
transgen_copy_number={}
with open("only_transgen_reads.txt") as copy_number_file:
    for line in copy_number_file:
        copy_number_list = line.strip().split("\t")
        transgen = copy_number_list[0]
        if transgen not in transgen_copy_number:
            list_depth = []
            list_depth.append(int(copy_number_list[3]))
            transgen_copy_number[transgen] = list_depth
        if transgen in transgen_copy_number:
            transgen_copy_number[transgen].append(int(copy_number_list[3]))
with open("report_transgen_copy_number.txt", "w") as tr_depth:
    tr_depth.write('Trangene\tCopy_number\tStandart_deviation\n')
    for node in transgen_copy_number:
        depth = (transgen_copy_number[node])
        mean_depth = np.mean(depth)/args.d 
        dev_depth = np.array(depth)/args.d     
        st_dev = np.std(dev_depth)
        tr_depth.write(f"{node}\t{mean_depth}\t{st_dev}\n")

rename = subprocess.call(f"mv report_transgen_copy_number.txt {args.rwr}_copy_number_{args.tr}.txt", shell=True)
remove = subprocess.call('rm aln.sam ref.sam filter_fwr.fastq.gz filter_rwr.fastq.gz QNAMES.txt sCo_ref_to_transgen.bed Comerge_ref.bed Co_ref_to_transgen.bed copy_ref_to_transgen.bam copy_ref_to_transgen.sam Coaln.bam sCoaln.bam transgen_depth.txt coordinates_tr_depth.bed only_transgen_reads.txt', shell=True)
     	
     	
     	
     	
     	

