#!/usr/bin/env python

import sys
import tqdm
import requests
import pandas as pd

ena_md5sums = []
with open(sys.argv[1]) as fh:
    for line in fh:
        line = line.strip().split('\t')
        r1 = line[0] + "_1.fastq.gz"
        r1_md5 = line[1].split(';')[0]
        r2 = line[0] + "_2.fastq.gz"
        r2_md5 = line[1].split(';')[1]
        ena_md5sums.append((r1, r1_md5))
        ena_md5sums.append((r2, r2_md5))
ena_md5sums = set(ena_md5sums)

down_md5sums = []
with open(sys.argv[2]) as fh:
    for line in fh:
        line = line.strip().split()
        read = line[1].split('/')[-1]
        md5 = line[0]
        down_md5sums.append((read, md5))

down_md5sums = set(down_md5sums)

print(f"Downloaded: {len(down_md5sums)}")
print(f"All : {len(ena_md5sums)}")
print(down_md5sums - ena_md5sums)

print(ena_md5sums - down_md5sums)
#failed = ena_md5sums - down_md5sums
#for seq in failed:
#    print(seq[0].split('_')[0])
#

print(ena_md5sums == down_md5sums)


