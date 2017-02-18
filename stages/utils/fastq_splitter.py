#!/usr/bin/env python
import argparse
import os
import time
import subprocess32 as subprocess
import multiprocessing as mp
import pysam #will read regular fastq or read fq.gz files

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
    P,base_name,output = {},fq.rsplit('/')[-1].rsplit('.')[0].rsplit('_')[0],''
    for i in range(split):
        for j in range(2):
            if P.has_key(i):
                P[i][j] = out_dir+'/'+base_name+'.pipe.%s.'%j+str(i)
            else:
                P[i] = {j:out_dir+'/'+base_name+'.pipe.%s.'%j+str(i)}
    return P
    
#take down the pipes
def remove_pipes(P):
    output = ''
    for i in P:
        for j in P[i]:
            try:
                output += subprocess.check_output(' '.join(['rm', P[i][j]]),
                                                  stderr=subprocess.STDOUT,shell=True)
            except Exception as E: output += str(E)
    return output
    
#P{i:{1:'/tmp/pipe.1.1',2:'/tmp/pipe.1.2'}} 
#is one pair of fifo pipels for use with one instance of bwa mem
#takes in a references fasta and a pair of named pipes that feed in reads
#output gets piped to an instance of samtools view to make a BAM file
def bwa_mem_piped_bam(P,i,ref_path):
    start = time.time()
    bam_name,output = P[i][0].rsplit('.pipe.')[-1]+'.'+str(i)+'.bam',''
    piped_mem = [bwa_mem,'-M','-t 4',ref_path,P[i][0],P[i][1],'|',
                 samtools_view,'-Sb','-','>',bam_name]
    try:#this should start and block with any luck----------------------------
        output += subprocess.check_output(' '.join([piped_mem]),
                                          stderr=subprocess.STDOUT,
                                          shell=True)
    except Exception as E:
        output += str(E)
    stop = time.time()
    return output

#do one pipe at a time, but can open multiple files
#file reading on f, pipe writting on p s is the modulatoing offset
def pipe_writer(in_f,out_p,s):
    x = 0
    if os.path.exists(out_p): os.remove(out_p)
    os.mkfifo(out_p)
    with pysam.FastxFile(in_f,'rb') as f: #can read compressed inputs
        print('opening: in_f=%s\tout_p=%s\ts=%s'%(in_f,out_p,s))
        with open(out_p,'w') as p:
            print('%s past blocking W call'%out_p)
            for entry in f: #[1] get the fastq records
                if (x+s)%4==0:
                    #print('found interleved skip == %s\n'%s)
                    line = '\n'.join(['+' if v is None else v for v in \
                             [entry.name,entry.sequence,
                              entry.comment,entry.quality]])+'\n'
                    #print(line)
                    p.write(line)
                x += 1
        print('wrote all data into pipe %s'%out_p)
    return True

#read in from a single pipe, TO DO read from two pipes...    
def pipe_reader(P,i,j):
    out = P[i][j].rsplit('.pipe.')[0]+'.%s.%s.fq'%(j,i)
    print('reading from pipe %s'%P[i][j])
    with open(P[i][j],'r') as p:
        print('past the blocking read for pipe %s'%P[i][j])
        with open(out,'a') as f:
            x,values = 0,[]
            print('opened file %s'%out)
            for line in p:
                print line,
                values += [line]
                if x%4==3:
                    f.write(''.join(values))
                    values = []
                x += 1
    return True

#this is the old working splitter code that made files.....    
#def splitter(fastq_file,out_dir):
#    start = time.time()
#    x,f,base_name = 0,{},fastq_file.rsplit('/')[-1].rsplit('.fq')[0] #or .fq.gz
#    for i in range(split): f[i] = open(out_dir+'/'+base_name+'.'+str(i)+'.fq','a')
#    with pysam.FastxFile(fastq_file) as in_fq:
#        for entry in in_fq: #[1] get the fastq records
#            values = ['+' if v is None else v for v in [entry.name,entry.sequence,entry.comment,entry.quality]]
#            f[x%split].write('\n'.join(values)+'\n')
#            x += 1
#    for i in range(split): f[i].close()
#    stop = time.time()
#    return [fastq_file,stop-start]

#|| by number of fastq files presented
if __name__ == '__main__':
    print('found %s fastq files'%len(fastq_files))
    F = make_pipes(fastq_files[0],split,out_dir)
    
    print('starting process forking and IPC')
    p = mp.Pool(processes=2)
    for i in F:
        for j in F[i]:
            print('dispatching writer %s'%F[i][j])
            p.apply_async(pipe_writer,
                          args=(fastq_files[j],F[i][j],i),
                          callback=collect_results)
            time.sleep(0.25)
            print('dispatching reader %s'%F[i][j])
            p.apply_async(pipe_reader,
                          args=(F,i,j),
                          callback=collect_results)
            time.sleep(0.25)
    p.close()
    p.join() #should now all be blocking until a reader flushes thme out?
    print('all processes have been joined!!!')
    print(result_list)
    
    ready = []
    for i in F:
        for j in F[i]:
            os.remove(F[i][j])
            ready += [os.path.exists(F[i][j])]
    print('named pipes deleted = %s'%(not any(ready)))
#---------------------------------------------------------------------------------------------------------------    