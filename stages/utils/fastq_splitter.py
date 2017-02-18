#!/usr/bin/env python
import argparse
import os
import pysam #will read regular fastq or read fq.gz files

parser = argparse.ArgumentParser(description='SVE:FASTQ splitting pipe tool')
parser.add_argument('-i', '--in_f',type=str, help='single fastq\t[None]')
parser.add_argument('-o', '--out_p',type=str, help='output pipe name\t[None]')
parser.add_argument('-s', '--split_factor',type=int, help='splitting factor\t[4]')
parser.add_argument('-j', '--jump',type=str, help='rows to jump for split\t[4]')
args = parser.parse_args()
#parse agruments and do some file checks---------------
in_f,out_p,s = args.in_f,args.out_p,args.split_factor
#do one pipe at a time, but can open multiple files
#file reading on f, pipe writting on p s is the modulatoing offset
def pipe_writer(in_f,out_p,s):
    x = 0
    with pysam.FastxFile(in_f,'rb') as f: #can read compressed inputs
        print('opening: in_f=%s\tout_p=%s\ts=%s'%(in_f,out_p,s))
        with open(out_p,'w') as p:
            print('%s past blocking W call'%out_p)
            for entry in f: #[1] get the fastq records
                if (x+s)%4==0: #this is the jump index for number of fastq rows
                    line = '\n'.join(['+' if v is None else v for v in \
                             [entry.name,entry.sequence,
                              entry.comment,entry.quality]])+'\n'
                    p.write(line)
                x += 1
    return 'wrote all data into pipe %s'%out_p
    
#|| by number of fastq files presented
if __name__ == '__main__':
    pipe_writer(in_f,out_p,s)
#---------------------------------------------------------------------------------------------------------------    