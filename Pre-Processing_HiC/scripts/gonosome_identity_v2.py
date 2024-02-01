#!/usr/local/bin/python3

import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
import csv
import Bio
import numpy as np
import pandas as pd
import os #importing operating system commands
import glob #importing glob which helps read files

# read in files
chrlen_file = sys.argv[1] # chromosome length file
clen_df = pd.read_csv(chrlen_file, delimiter='\t', names=['chr', 'len'])
mappedfile = sys.argv[2] # mapped reads per chromosome file
chrom_l = []
mapped_l = []
with open(mappedfile) as f:
        for line in f:
            line = line.strip()
            if "chr" in line:
                chrom_l.append(line)
            else:
                mapped_l.append(line)
mapped_df = pd.DataFrame({'chr':chrom_l,
                              'mapped':mapped_l})

#claculate coverage per chromosome
print("chromosome mapped_coverage")
for c in chrom_l:
    N = mapped_df.loc[mapped_df['chr'] == c, 'mapped']
    G = clen_df.loc[clen_df['chr'] == c, 'len']
    cov = int(N) / int(G)
    print(c,cov)
