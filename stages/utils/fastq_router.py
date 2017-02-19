#!/usr/bin/env python
import argparse
import sys
import os
import itertools as it
import pysam
des = """
SVE:FASTQ PE routing, splitting and interleveing tool\n
fastq_router.py -f /data/fastq/test_1.fq.gz,/data/fastq/test_2.fq.gz -s 4 | bwa mem -M -p...
"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-f', '--fastq',type=str, help='single fastq or comma separated\t[None]')
parser.add_argument('-s', '--split_index',type=int, help='splitting index\t[4]')
args = parser.parse_args()
fastqs,split = args.fastq.split(','),args.split_index
if not all([os.path.exists(f) for f in fastqs]): raise IOError
#main iterator reads in two file handles and outputs interleved
with pysam.FastxFile(fastqs[0],'rb') as f1:
    with pysam.FastxFile(fastqs[1],'rb') as f2:
        x,values = 0, []
        for reads in it.izip(*[f1,f2]):
            if (x+split)%4==0: #default 4 lines of fastq format
                for i in range(2):
                    values += [reads[i].name,reads[i].sequence,
                               reads[i].comment,reads[i].quality]
                values = ['+' if v is None else v for v in values]
                print('\n'.join(values)+'\n') #could do name check assertion here
                values = []
            x += 1
sys.exit()