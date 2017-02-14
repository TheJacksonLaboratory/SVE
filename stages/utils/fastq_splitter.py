#!/usr/bin/env python
#fastq splitting tool using the pysam API
import argparse
import pysam
import os
parser = argparse.ArgumentParser(description='FASTQ splitting tool')
parser.add_argument('-f', '--fastq',type=str, help='input fastq file')
parser.add_argument('-o', '--out_dir',type=str, help='output directory')
parser.add_argument('-s', '--split_factor',type=int, help='splitting factor')
args = parser.parse_args()
fastq_file = args.fastq
base_name = fastq_file.rsplit('/')[-1].rsplit('.')[0]
out_dir    = args.out_dir
split       = args.split_factor
x,f = 0,{} #base split counter
if not os.path.exists(out_dir): os.makedirs(out_dir)
for i in range(split):f[i] = open(out_dir+'/'+base_name+'_'+str(i)+'.fq','a')
with pysam.FastxFile(fastq_file) as in_fq:
    for entry in in_fq: #[1] get the fastq records
        values = ['+' if v is None else v for v in [entry.name,entry.sequence,entry.comment,entry.quality]]
        f[x%split].write('\n'.join(values)+'\n')
        x += 1
for i in range(split): f[i].close()
 
   