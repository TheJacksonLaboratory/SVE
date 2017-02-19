#!/usr/bin/env python
import argparse
import os
import time
import itertools as it
import subprocess32 as subprocess
import multiprocessing as mp
import pysam

parser = argparse.ArgumentParser(description='SVE:FASTQ splitting,alignment,mapping test pipeline V-4.0')
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
if not os.path.exists(out_dir): os.makedirs(out_dir)
software = os.path.dirname(os.path.abspath(__file__))+'/../../../' #/software/SVE/stages/utils/../../../
print(software)#-----------------------------------------------
bwa_mem       = software+'/bwa-master/bwa mem' 
samtools_view = software+'/samtools-1.3/samtools view'
#--------------------------------------------------------------
#if you need to pass back some statuses ect..
result_list = []
def collect_results(result):
    result_list.append(result)
    
#we are assuming that there is a pair of fastq files here for bwa mem input
#make 2xsplit number of output pipes avaible to the system
def make_pipes(fq,split,out_dir):
    P,base_name = {},fq.rsplit('/')[-1].rsplit('.')[0].rsplit('_')[0]
    for i in range(split):
        for j in range(2):
            if P.has_key(i):
                P[i][j] = out_dir+'/'+base_name+'.pipe.%s.'%j+str(i)
            else:
                P[i] = {j:out_dir+'/'+base_name+'.pipe.%s.'%j+str(i)}
    for i in range(split):
        for j in range(2):
            if os.path.exists(P[i][j]): os.remove(P[i][j])
            os.mkfifo(P[i][j])
    return P

#do one pipe at a time, but can open multiple files
#file reading on f, pipe writting on p s is the modulatoing offset
def pipe_writer(in_f,out_p,s):
    x = 0
    with pysam.FastxFile(in_f,'rb') as f: #can read compressed inputs
        print('opening: in_f=%s\tout_p=%s\ts=%s'%(in_f,out_p,s))
        with open(out_p,'w',0) as p:
            print('%s past blocking W call'%out_p)
            for entry in f: #[1] get the fastq records
                if (x+s)%4==0: #this is the jump index for number of fastq rows
                    line = '\n'.join(['+' if v is None else v for v in \
                             [entry.name,entry.sequence,
                              entry.comment,entry.quality]])+'\n'
                    p.write(line)
                x += 1
    return 'wrote all data into pipe %s'%out_p

def pipe_reader(P,i,ref_path):
    out,output = P[i][0].rsplit('.pipe.')[0]+'.%s.fq'%i,None
    with open(P[i][0],'w') as p1:
        print('%s past blocking R call'%P[i][0])
        with open(P[i][1],'w') as p2:
            print('%s past blocking R call'%P[i][1])
            with open(out,'w') as f:
                x,values = 0,{0:[],1:[]}
                for lines in it.izip(*[p1,p2]):
                    print lines[0],
                    print lines[1],
                    if x%4==3:
                        values[0] += [lines[0]]
                        values[1] += [lines[1]]
                        f.write(''.join(values[0]+values[1]))
                        values[0],values[1] = [],[]
    return 'finished writing and %s have been released'%P[i]
                
#No try to get bwa mem working on these pipes...
def paired_fastq_pipe_bwa_reader(P,i,ref_path):
    out,output = P[i][0].rsplit('.pipe.')[0]+'.%s.bam'%i,None #paired coding
    piped_mem = [bwa_mem,'-M',ref_path,'<(cat %s)'%P[i][0],'<(cat %s)'%P[i][1],'>',out]#,'|',
                #samtools_view,'-Sb','-','>',out]
    try:#bwa mem cll here-------------------------------------------
        output = subprocess.Popen(' '.join(piped_mem),
                                   shell=True,
                                   executable='/bin/bash')
        output.wait()
    except Exception as E:
        output = str(E) #error on the cal-------------------------
    return output
    
#|| by number of fastq files presented
if __name__ == '__main__':
    print('found %s fastq files'%len(fastq_files))
    F = make_pipes(fastq_files[0],split,out_dir)
    
    p = mp.Pool(processes=3*split)
    print('starting process forking and IPC')
    for i in F:
#        print('dispatching pipe reader %s'%F[i])
#        p.apply_async(paired_fastq_pipe_bwa_reader,
#                       args=(F,i,ref_path),
#                       callback=collect_results)
#        time.sleep(1)
        #now lets try one reader for every two write pipes:
        for j in F[i]: #fq1_pipe + fq2_pipe => paired_fastq_pipe_reader
            print('dispatching pipe writer %s'%F[i][j])
            p.apply_async(pipe_writer,
                          args=(fastq_files[j],F[i][j],i),
                          callback=collect_results)
            time.sleep(0.25)
    p.close()
    p.join()
    print('all processes have been joined!!!')
    for result in result_list: print(result)
    
    ready = []
    for i in F:
        for j in F[i]:
            os.remove(F[i][j])
            ready += [os.path.exists(F[i][j])]
    print('named pipes deleted = %s'%(not any(ready)))
#---------------------------------------------------------------------------------------------------------------    