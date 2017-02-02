#!/usr/bin/env python

import sys

args = sys.argv
n,m     = int(round(float(args[1]))),0
d = int(1E3)
valid = args[2] 
MIN,MAX = 255,0
RANGES = {
    'Sanger':       (33, 73),
    'Solexa':       (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (67, 104)
}
for line in sys.stdin:
    if not line.startswith('@'):
        m += 1
        if m%d==0: print('MAX=%s\tMIN=%s'%(MAX,MIN))
        raw  = [ord(j) for j in line.split('\t')[10]]
        for r in raw:
            if MIN > r: MIN = r
            if MAX < r: MAX = r
    if MIN <= 33: #sanger adds no new inputs for our validation pipeline
        with open(valid,'a') as f:     f.write('[base quality encoding] = Sanger:33-73\n')
        sys.exit()
    elif MIN > 33 and MAX > 73:
        if MIN < 64:
            with open(valid,'a') as f: f.write('[base quality encoding] = Solexa:59-104\n')
            sys.exit()
        elif MIN < 67:
            with open(valid,'a') as f: f.write('[base quality encoding] = Illumina-1.3:34-104\n')
            sys.exit()
        else:
            with open(valid,'a') as f: f.write('[base quality encoding] = Illumina-1.5:67-104\n')
            sys.exit()
    if m>n: break
with open(valid,'a') as f:             f.write('[base quality encoding] = Unkown Illumina:34-104\n')
sys.exit()