#!/usr/bin/env python
import argparse
import os
import glob
import subprocess32 as subprocess #to call qsub a bunch of times

des = """
script that auto generates PBS scripts for runing a single caller on a folder of bams via qsub
saving the results to a single root directory with sample-name subfolders for configuration/intermediate files"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-b', '--bam_dir',type=str, help='bam file or directory')
parser.add_argument('-o', '--output_dir',type=str, help='outputdirectory to save ...bam.bai/ into')
parser.add_argument('-t', '--wall_time',type=str,help='wall time requested from cluster')
parser.add_argument('-m', '--email_address',type=str,help='cluster email results to this email address')
args = parser.parse_args()

if args.wall_time is not None:
    walltime = args.wall_time
else:
    walltime = '04:00:00'

if args.bam_dir is not None:
    if args.bam_dir.endswith('.bam'):
        bam_list = [args.bam_dir]
    else:
        bam_list = glob.glob(args.bam_dir+'/*.bam')
    if len(bam_list)<1:
        print('bam files were not found')
        raise IOError
    samples = {}
    samples = {i.rsplit('/')[-1]:i for i in bam_list} #file name from path
else:
    print('bam file or directory not found')
    raise IOError

if args.output_dir is not None:
    out_dir = args.output_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
else:
    print('out put directory not specified')
    out_dir = args.bam_dir+'/'
    pass

if args.email_address is not None:
    email = args.email_address
else:
    print('email not specified')
    email = 'timothy.becker@jax.org' #default here
    
    
#make a temp directory for pbs scripts
pbs_dir = out_dir+'PBS'
if not os.path.exists(pbs_dir): os.makedirs(pbs_dir)
    
#write out the .pbs scripts
PBS = []
for k in samples:
    samtools = '/home/tbecker/software/samtools-1.2/samtools'
    sample_pbs = pbs_dir+'/'+k+'.pbs' #name of pbs script for sample k
    PBS += [sample_pbs]
#    if not os.path.exists(out_dir+k): os.makedirs(out_dir+k)
    with open(sample_pbs,'w') as pbs:
        index = [samtools, 'index', samples[k]]
        pbs.write('#!/bin/bash\n'+' '.join(index))
#execute qsub with the scripts, getting the jids back (can display these)
output,err = '',{}
for pbs in PBS:
    print('processing %s'%pbs)
    try:
        command = ['qsub','-l','walltime=%s'%walltime+',mem=16gb','-m','e','-M',email,'-o','index.log','-j oe',pbs]
        print(' '.join(command)) #don't run it yet!
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