#!/usr/bin/env python
#fastq splitting tool using the pysam API, || by fastq file
import argparse
import os
import io
import gzip
import time
import multiprocessing as mp
import pysam

parser = argparse.ArgumentParser(description='FASTQ splitting tool')
parser.add_argument('-f', '--fastq',type=str, help='single fastq or comma separated\t[None]')
parser.add_argument('-o', '--out_dir',type=str, help='output directory\t[same as input on fastq]')
parser.add_argument('-s', '--split_factor',type=int, help='splitting factor\t[4]')
parser.add_argument('-c','--comp',action='store_true',help='gzip compressed output\t[True]')
args = parser.parse_args()

fastq_files = args.fastq.split(',')
out_dir     = args.out_dir
split       = args.split_factor
comp        = args.comp
if not os.path.exists(out_dir): os.makedirs(out_dir)

result_list = []
def collect_results(result):
    result_list.append(result)
    
def splitter(fastq_file,split,out_dir,comp=True):
    start = time.time()
    x,y,f,base_name = 0,4,{},fastq_file.rsplit('/')[-1].rsplit('.fq')[0] #or .fq.gz
    if not comp: #do not use gzip compression on output streams----------------------------------------------------
        for i in range(split): f[i] = open(out_dir+'/'+base_name+'.'+str(i)+'.fq','a')
        with pysam.FastxFile(fastq_file) as in_fq:
            for entry in in_fq:
                values = ['+' if v is None else v for v in [entry.name,entry.sequence,entry.comment,entry.quality]]
                f[x%split].write('\n'.join(values)+'\n')
                x += 1
        for i in range(split): f[i].close()
    else: #use the gzip compression on the output sreams-----------------------------------------------------------
        for i in range(split): f[i] = io.BufferedWriter(gzip.open(out_dir+'/'+base_name+'.'+str(i)+'.fq.gz','ab'))
        with io.BufferedReader(gzip.open(fastq_file,'rb')) as in_fq:
            values = []
            for line in in_fq: #y start at 0
                y = (y+1)%4
                values += [line]
                if y == 3:
                    f[x%split].write('\n'.join(values)+'\n')
                    x,values = x+1,[]
        for i in range(split): f[i].close()
    stop = time.time()
    return [fastq_file,stop-start]

#|| by number of fastq files presented
if __name__ == '__main__':
    print('found %s fastq files'%len(fastq_files))
    print('using gzip compression on outputs = %s'%comp)     
    start = time.time()
    p = mp.Pool(processes=len(fastq_files))
    for fastq_file in fastq_files:
        p.apply_async(splitter,
                      args=(fastq_file,split,out_dir),
                      callback=collect_results)
    p.close()
    p.join()
    for result in result_list: print(result)
    stop = time.time()
    print('completed in %s sec'%round(stop-start,0))
#---------------------------------------------------------------------------------------------------------------    