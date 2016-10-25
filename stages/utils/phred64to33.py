#!/usr/bin/env python

#base quality conversion Phred64 to Phred33 scale
import sys

for line in sys.stdin:
    if line.startswith('@'): #header
        print line
    else:                    #alignments
        raw = line.split('\t')
        print '\t'.join(raw[:10]+[''.join([chr(ord(i)-33) for i in raw[10]])]+raw[11:])