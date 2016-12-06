#!/usr/bin/env python

#guess the base quality encoding from a 1E6 randomly selected alignment records
#used as: samtools view -Sh old.bam | SVE/stages/utils/phred_encoding.py 1E6 > phred.encoding
import sys
import numpy as np

args = sys.argv
n = int(round(float(args[1]))) 
MIN,MAX = 255,0
RANGES = {
    'Sanger': (33, 73),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (67, 104)
}
for line in sys.stdin:
    if line.startswith('@') and np.random.choice([True,False],p=[n,1.0-n]): 
        raw  = line.split('\t')
        for i in raw[10]:
            j = ord(i)
            if j < MIN: MIN = j
            if j > MAX: MAX = j
    if MIN <= 33:
        print 'Sanger:33-73'
        sys.exit()
    elif MIN > 33 and MAX > 73:
        if MIN < 64:
            print 'Solexa:59-104'
            sys.exit()
        elif MIN < 67:
            print 'Illumina-1.3:64-104'
            sys.exit()
        else:
            print 'Illumina-1.5:67-104'
            sys.exit()