#HTSeq quality conversion Phred64 to Phred33 scale
import argparse
import sys
import time
import os
import HTSeq as ht

des = """phred64 to phred 33 base quality conversion"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-i', '--in_bam_path',type=str, help='input bam file\t[None]')
parser.add_argument('-o', '--out_bam_path',type=str, help='output bam file\t[None]')
args = parser.parse_args()

if args.in_bam_path is None or args.out_bam_path is None:
    print('bam files not specified!')
    raise IOError
else:
    in_bam_path  = args.in_bam_path
    out_bam_path = args.out_bam_path

start = time.time()
i = 0
bam_reader = ht.BAM_Reader(in_bam_path)
bam_writer = ht.BAM_Writer.from_BAM_Reader(out_bam_path,bam_reader)
for aln in bam_reader:
    if i%1E6==0: print('-')
    aln.read.qual -= 31
    bam_writer.write(aln)
stop = time.time()
print('%s record written in %s sec'%(i,int(round(stop-start))))
