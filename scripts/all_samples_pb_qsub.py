#!/usr/bin/env python
import argparse
import os
import glob
import subprocess32 as subprocess #to call qsub a bunch of times

des = """
script that auto generates PBS scripts for runing a bwa mem|aln sampe on a folder of .fq via qsub
saving the resulting bams,bais to a single root directory"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-r', '--ref_path',type=str, help='reference fasta')

parser.add_argument('-i', '--fqs_input_dir',type=str, help='directory for fastq files')
sample_to_lane_help = """
sample to multiple lane mapping string
samples are associated to lanes with the ':' symbol
lanes are seperated by the ',' symbol
samples are sperated by the ';' symbol
[EX SE/PE] --sample_to_lane NA12878:ERR1205,ERR1206NA19238:ERR1104
"""
parser.add_argument('-s', '--sample_to_lane',type=str, help=sample_to_lane_help)
parser.add_argument('-f', '--fqs_pattern',type=str, help='fqs search pattern [_1.fq.gz,_2.fq.gz]')

parser.add_argument('-g', '--replace_rg',action='store_true', help='replace reads groups')
parser.add_argument('-d', '--mark_duplicates',action='store_true', help='mark duplicate reads')
parser.add_argument('-a', '--algorithm', type=str, help='aln, mem, [piped_mem]')
parser.add_argument('-o', '--output_dir',type=str, help='outputdirectory to save .bam,.bai/ into')
parser.add_argument('-w', '--wall_time',type=str,help='wall time requested from cluster')
parser.add_argument('-m', '--memory',type=str,help='radom access memory needed')
parser.add_argument('-p', '--processors',type=int,help='processors needed')
parser.add_argument('-e', '--email_address',type=str,help='cluster email results to this email address')
args = parser.parse_args()

if args.ref_path is not None:
    ref_path = args.ref_path
else:
    print('reference fasta not found')
    raise IOError
    
if args.fqs_pattern is not None:
    fqs_pattern = {p:set([]) for p in args.fqs_pattern.split(',')}
else:
    fqs_pattern = {'_1.fq.gz':set([]),'_2.fq.gz':set([])} #default apptern here

if args.sample_to_lane is not None and args.fqs_input_dir is not None:
    s = args.sample_to_lane
    if s.find(':')!=-1 and s.find(',')!=-1: #associate se or pe fastq lanes to a sample name
        S = {s.split(':')[0]:{k:[] for k in s.split(':')[1].split(',')} for s in args.sample_to_lane.split(';')}
        F = {p.rsplit('/')[-1]:p for p in glob.glob(args.fqs_input_dir+'*')}
        for p in fqs_pattern:
            for f in F:
                if f.endswith(p):
                    fqs_pattern[p].add(f.split(p)[0])
        for p in fqs_pattern:
            for q in fqs_pattern:
                fqs_pattern[p].intersection_update(fqs_pattern[q])       
        for k in S:                     #for each sample
            for r in S[k]:              #the lane that should be attached to each sample
                for p in fqs_pattern:
                    if r in fqs_pattern[p]:
                        S[k][r] += [F[r+p]]
            for r in S[k]:
                S[k][r] = sorted(S[k][r])
        print('sample to lane mapping:\n')
        print(S)
    else:
        print(' --sample_to_lane: valid seperators not detected')
        raise IOError
elif args.fqs_input_dir is not None:
    print('--sample_to_lane arguments did not lead to a valid association')
    print('looking at fqs pattern inputs, assuming each pair is a single sample')
    F,S = {p.rsplit('/')[-1]:p for p in glob.glob(args.fqs_input_dir+'*')},{}
    for p in fqs_pattern:
        for f in F:
            if f.endswith(p):
                fqs_pattern[p].add(f.split(p)[0])
    for p in fqs_pattern:
        for q in fqs_pattern:
            fqs_pattern[p].intersection_update(fqs_pattern[q])
    for p in sorted(fqs_pattern):
        s = list(fqs_pattern[p])
        for i in s:
            if S.has_key(i):
                S[i][i] += [F[i+p]]
            else:
                S[i] = {i:[F[i+p]]}
else:
    print('--fqs_input_dir not specified')
    raise IOError

#use bwa mem as the default
if args.algorithm is not None:
    algo = args.algorithm 
else:
    algo = 'piped_mem' #default algoritm is the easiest/best

if args.output_dir is not None:
    out_dir = args.output_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
else:
    print('output directory not specified')
    raise IOError

if args.wall_time is not None:
    walltime = args.wall_time
else:
    walltime = '24:00:00'

if args.memory is not None:
    ram = args.memory
else:
    ram = '32gb' #default

if args.processors is not None:
    cpus = args.processors
else:
    cpus = 8     #default

if args.email_address is not None:
    email = args.email_address
else:
    print('email not valid')
    raise AttributeError
    
#make a temp directory for pbs scripts
pbs_dir = out_dir+'PBS'
if not os.path.exists(pbs_dir): os.makedirs(pbs_dir)
#write out the .pbs scripts for each set of reads
software_path = os.path.dirname(os.path.abspath(__file__))+'/../../'
PBS = []
python = software_path+'/anaconda/bin/python'
pb     = software_path+'/SVE/scripts/prepare_bam.py'
ml = 'module load '
modules = [ml+'perl/5.16.3',ml+'gcc/4.9.2',ml+'Root/v5.34.18',ml+'samtools/1.2'] #everything that you need here...
for sample in S:
    for rg in S[sample]:
        group_pbs = pbs_dir+'/prepare_bam_'+sample+'_'+rg+'.pbs'
        PBS += [group_pbs]
        with open(group_pbs,'w') as pbs:
            call = [python,pb,'-r',ref_path,'-f',','.join(S[sample][rg]),
                    '-s',sample,'-a',algo,'-o',out_dir]
            pbs.write('#!/bin/bash\n'+'\n'.join(modules)+'\n'+' '.join(call))
#execute qsub with the scripts, getting the jids back (can display these or attach to further monitor progress)
output,err = '',{}
for pbs in PBS: #test with one of these and a fast caller on a small file...
    print('processing %s'%pbs)
    try:
        command = ['qsub','-l','procs=%s'%cpus+',walltime=%s'%walltime+',mem=%s'%ram,
                   '-m','e','-M',email,'-o',pbs[0:-4]+'.log','-j oe',pbs]
        #print(' '.join(command)) #don't run it yet!
        output += subprocess.check_output(' '.join(command),stderr=subprocess.STDOUT,shell=True)
    #catch all errors that arise under normal call behavior
    except subprocess.CalledProcessError as E:
        print('call error: '+E.output)        #what you would see in the term
        err['output'] = E.output
        #the python exception issues (shouldn't have any...
        print('message: '+E.message)          #?? empty
        err['message'] = E.message
        #return codes used for failure....
        print('code: '+str(E.returncode))     #return 1 for a fail in art?
        err['code'] = E.returncode
    except OSError as E:
        print('os error: '+E.strerror)        #what you would see in the term
        err['output'] = E.strerror
        #the python exception issues (shouldn't have any...
        print('message: '+E.message)          #?? empty
        err['message'] = E.message
        #the error num
        print('code: '+str(E.errno))
        err['code'] = E.errno
print('output:\n'+output)

#remove/delete the intermediate pbs scripts TO DO...
