#!/usr/bin/env python

import sys
import numpy as np

tsv, bin_size = sys.argv[1:]
bin_size = int(bin_size)

sname = tsv.split('_')[0]

with open(tsv, 'r') as tsv_file:
    base_counter = 0
    bin_counter = 0
    coverages = []
    previous_tig = ''
    for line in tsv_file:
        content = line.strip().split('\t')
        if previous_tig == '':
            previous_tig = content[0]
        if bin_counter == bin_size or content[0] != previous_tig:
            mean_cov = np.mean(coverages)
            sd_cov = np.std(coverages)
            bin_counter = 0
            print(f'{sname}\t{previous_tig}\t{base_counter}\t{mean_cov}\t{sd_cov}')
            previous_tig = content[0]
            coverages = []
        coverages.append(int(content[2]))
        base_counter += 1
        bin_counter += 1

mean_cov = np.mean(coverages)
sd_cov = np.std(coverages)
bin_counter = 0
previous_tig = content[0]
print(f'{sname}\t{previous_tig}\t{base_counter}\t{mean_cov}\t{sd_cov}')



