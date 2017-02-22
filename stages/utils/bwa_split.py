#!/usr/bin/env python
import argparse
import os
import time
import subprocess32 as subprocess
import multiprocessing as mp

parser = argparse.ArgumentParser(description='SVE:|| FASTQ splitting,alignment,mapping test pipeline 0.4.1')
parser.add_argument('-r', '--ref_path',type=str, help='reference fasta file to align to\t[None]')
parser.add_argument('-f', '--fastq',type=str, help='comma separated fastq\t[None]')
parser.add_argument('-o', '--out_dir',type=str, help='output directory\t[same as input on fastq]')
parser.add_argument('-s', '--split_factor',type=int, help='splitting factor\t[4]')
parser.add_argument('-t', '--threads',type=int,help='threads per split\t[4]')
parser.add_argument('-R', '--read_group',type=int, help='bwa read group string: "@RG\tID:foo\tSM:bar"[None]')
args = parser.parse_args()
#parse agruments and do some file checks---------------
if args.ref_path is not None:
    ref_path    = args.ref_path
else:
    print('--ref_path not specified')
    raise IOError
if args.fastq is not None:
    fastq_files = args.fastq.split(',')
else:
    print('--fastq files not located')
    raise IOError
if args.out_dir is not None:
    out_dir     = args.out_dir
else:
    print('--out_dir not specified')
    out_dir = '/'.join(fastq_files[0].rsplit('/')[:-1]) #use input directory
if args.split_factor is not None:
    split       = args.split_factor
else:
    split = 4 #default split factor
if args.threads is not None:
    t = args.threads
else:
    t = 4
if args.read_group is not None:
    read_group  = args.read_group
else: #take the fq name after stripping off all the '.','_','-' seperators
    base = fastq_files[0].rsplit('/')[-1].rsplit('.')[0].rsplit('_')[0].rsplit('-')[0]
    SM,ID,LB,PL = base,base+'_RG',base+'_RG','illumina'
    read_group = '@RG\tID:%s\tSM:%s\tLB:%s\tPL:%s'%(ID,SM,LB,PL)
    print('using read group tag:"%s" for bwa mem alignment'%read_group)

if os.path.exists(out_dir): subprocess.call(['rm','-rf',out_dir])
os.makedirs(out_dir)
software = os.path.dirname(os.path.abspath(__file__))+'/../../../' #/software/SVE/stages/utils/../../../
route            = software+'/SVE/stages/utils/fastq_route.py'
bwa_mem          = software+'/bwa-master/bwa mem' 
samtools_view    = software+'/samtools-1.3/samtools view'
sambamba_merge    = software+'/sambamba/sambamba merge'
sambamba_markdup = software+'/sambamba/sambamba markdup'
sambamba_sort    = software+'/sambamba/sambamba sort'
sambamba_index   = software+'/sambamba/sambamba index'
#--------------------------------------------------------------
#if you need to pass back some statuses ect..
result_list = []
def collect_results(result):
    result_list.append(result)
                
#one of the split values to assemble...
def fastq_bwa_mem_piped(fastqs,i,j,t,rg,out_dir,ref_path):
    output = ''
    bam = out_dir+'/'+fastqs[0].rsplit('/')[-1].rsplit('.fq')[0].rsplit('_')[0]+'.%s.bam'%j
    piped_mem = [bwa_mem,'-M','-t %s'%t,'-R',r"'%s'"%rg,ref_path,
                 '<(%s -f %s -i %s -j %s)'%(route,fastqs[0],i,j),
                 '<(%s -f %s -i %s -j %s)'%(route,fastqs[1],i,j),
                 '|',samtools_view,'-Sb','-','-@ %s'%t,
                 '|',sambamba_sort,'-t %s'%t,'--tmpdir=%s/temp'%out_dir,'-o',bam,'/dev/stdin']  #-@ for threads here
    try:#bwa mem call here-------------------------------------------
        output += subprocess.check_output(' '.join(piped_mem),
                                          stderr=subprocess.STDOUT,
                                          executable='/bin/bash',
                                          shell=True)
        output += subprocess.check_output(['rm','-rf','%s/temp'%out_dir])
    except Exception as E:
        output += str(E) #error on the call-------------------------
    return bam #list of bam files to merge into next step

    
#|| by number of fastq files presented
if __name__ == '__main__':
    start = time.time()
    print('found %s fastq files'%len(fastq_files))
    p = mp.Pool(processes=split)
    for j in range(split):
        print('dispatching bwa mem process %s'%j)
        p.apply_async(fastq_bwa_mem_piped,
                       args=(fastq_files,split,j,t,read_group,out_dir,ref_path),
                       callback=collect_results)
    p.close()
    p.join()
    print('all processes have been joined!!!')
    for i in result_list: print(i)
    
    #merge now
    output = ''
    bam = out_dir+'/'+fastq_files[0].rsplit('/')[-1].rsplit('.fq')[0].rsplit('_')[0]+'.bam'
    print('merging into one bam file %s'%bam)
    if len(result_list)>1:
        merge_bams = [sambamba_merge,'-t %s'%(split*t),'-l 9',bam]+result_list
        try:#bwa mem call here-------------------------------------------
            output += subprocess.check_output(' '.join(merge_bams),
                                              stderr=subprocess.STDOUT,
                                              executable='/bin/bash',
                                              shell=True)
            for i in result_list:
                output += subprocess.check_output(' '.join(['rm',i+'*']),
                                                  stderr=subprocess.STDOUT,
                                                  executable='/bin/bash',
                                                  shell=True)
        except Exception as E:
            output += str(E)
    else:
        try:
            output += subprocess.check_output(['mv',result_list[0],bam])
            output += subprocess.check_output(['mv',result_list[0]+'.bai',bam+'.bai'])
        except Exception as E:
            output += str(E)
    mark_dups = [sambamba_markdup,'-t %s'%(split*t),'-l 9','--tmpdir=%s/temp'%out_dir,
                 bam,bam.rsplit('.bam')[0]+'.dup.bam']
    #mark dups now
    try:#bwa mem call here-------------------------------------------
        output += subprocess.check_output(' '.join(mark_dups),
                                          stderr=subprocess.STDOUT,
                                          executable='/bin/bash',
                                          shell=True)
        output += subprocess.check_output(['rm','-rf','%s/temp'%out_dir])
        output += subprocess.check_output(' '.join(['rm',bam+'*']),
                                          stderr=subprocess.STDOUT,
                                          executable='/bin/bash',
                                          shell=True)
    except Exception as E:
        output += str(E)
    #now you can merge, mark duplicates, sort and index...
    stop = time.time()
    print('total time was %s sec'%round(stop-start,0))
    print(subprocess.check_output(['ls','-lh',out_dir]))
#---------------------------------------------------------------------------------------------------------------
#/software/SVE/stages/utils/bwa_split.py -r /data/human_g1k_v37_decoy/human_g1k_v37_decoy.fa -f /data/fastq/test_1.fq.gz,/data/fastq/test_2.fq.gz -o /data/test/ -s 4    