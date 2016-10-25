#!/usr/bin/env python

#HTSeq quality conversion Phred64 to Phred33 scale
import sys

for line in sys.stdin:
    if line.startswith('@'): #header
        print line
    else:                    #alignments
        raw = line.split('\t')
        print raw[:10]+[chr(ord(raw[10])-33)]+raw[11:]
