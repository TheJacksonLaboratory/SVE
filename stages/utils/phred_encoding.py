#!/usr/bin/env python

#guess the base quality encoding from a 1E6 randomly selected alignment records
#used as: samtools view -Sh old.bam | SVE/stages/utils/phred_encoding.py 1E6 ./old.valid
import sys
import numpy as np

args = sys.argv
n     = int(round(float(args[1])))
valid = args[2] 
MIN,MAX = 255,0
RANGES = {
    'Sanger':       (33, 73),
    'Solexa':       (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (67, 104)
}
for line in sys.stdin:
    if line.startswith('@') and np.random.choice([True,False],p=[1.0/n,1.0-1.0/n]): 
        raw  = line.split('\t')
        for i in raw[10]:
            j = ord(i)
            if j < MIN: MIN = j
            if j > MAX: MAX = j
if MIN <= 33: #sanger adds no new inputs for our validation pipeline
    with open(valid,'a') as f:     f.write('\nbase quality encoding = Sanger:33-73\n')
    sys.exit()
elif MIN > 33 and MAX > 73:
    if MIN < 64:
        with open(valid,'a') as f: f.write('\nbase quality encoding = Solexa:59-104\n')
        sys.exit()
    elif MIN < 67:
        with open(valid,'a') as f: f.write('\nbase quality encoding = Illumina-1.3:64-104\n')
        sys.exit()
    else:
        with open(valid,'a') as f: f.write('\nbase quality encoding = Illumina-1.5:67-104\n')
        sys.exit()