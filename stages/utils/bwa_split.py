#!/usr/bin/env python
import argparse
import os
import time
import subprocess32 as subprocess
import multiprocessing as mp

parser = argparse.ArgumentParser(description='SVE:|| FASTQ splitting,alignment,mapping test pipeline 0.4.1')
parser.add_argument('-r', '--ref_path',type=str, help='reference fasta file to align to\t[None]')
parser.add_argument('-f', '--fastq',type=str, help='single fastq or comma separated\t[None]')
parser.add_argument('-o', '--out_dir',type=str, help='output directory\t[same as input on fastq]')
parser.add_argument('-s', '--split_factor',type=int, help='splitting factor\t[4]')
args = parser.parse_args()
#parse agruments and do some file checks---------------
ref_path    = args.ref_path
fastq_files = args.fastq.split(',')
out_dir     = args.out_dir
split       = args.split_factor
if os.path.exists(out_dir): subprocess.call(['rm','-rf',out_dir])
os.makedirs(out_dir)
software = os.path.dirname(os.path.abspath(__file__))+'/../../../' #/software/SVE/stages/utils/../../../
router        = software+'/SVE/stages/utils/fastq_router.py'
bwa_mem       = software+'/bwa-master/bwa mem' 
samtools_view = software+'/samtools-1.3/samtools view'
#--------------------------------------------------------------
#if you need to pass back some statuses ect..
result_list = []
def collect_results(result):
    result_list.append(result)
                
#one of the split values to assemble...
def fastq_bwa_mem_piped(fastqs,i,out_dir,ref_path):
    output = ''
    bam = out_dir+'/'+fastqs[0].rsplit('/')[-1].rsplit('.fq')[0]+'.%s.bam'%i
    piped_mem = [router,'-f',','.join(fastqs),'-s','|',
                 bwa_mem,'-M','-p',ref_path,'|',   #-t for threads
                 samtools_view,'-Sb','-','>',bam]  #-@ for threads here
    try:#bwa mem call here-------------------------------------------
        output += subprocess.check_output(' '.join(piped_mem),
                                          stderr=subprocess.STDOUT,
                                          executable='/bin/bash',
                                          shell=True)
    except Exception as E:
        output += str(E) #error on the call-------------------------
    return [output]
    
#|| by number of fastq files presented
if __name__ == '__main__':
    print('found %s fastq files'%len(fastq_files))
    p = mp.Pool(processes=split)
    for i in range(split):
        print('dispatching bwa mem process %s'%i)
        p.apply_async(fastq_bwa_mem_piped,
                       args=(fastq_files,i,out_dir,ref_path),
                       callback=collect_results)
    p.close()
    p.join()
    print('all processes have been joined!!!')
    for i in result_list: print(i)
#---------------------------------------------------------------------------------------------------------------
#/software/SVE/stages/utils/bwa_split.py -r /data/human_g1k_v37_decoy/human_g1k_v37_decoy.fa -f /data/fastq/test_1.fq.gz,/data/fastq/test_2.fq.gz -o /data/test/ -s 4    