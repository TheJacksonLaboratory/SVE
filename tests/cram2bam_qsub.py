import argparse
import os
import sys
import csv
import datetime
import subprocess32 as subprocess #to call qsub a bunch of times

des = """
script that auto generates PBS scripts for simple cram to bam conversion and submits them via qsub"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-r', '--ref_path',type=str, help='reference fasta')
parser.add_argument('-c', '--cram_list',type=str, help='csv cram file list')
parser.add_argument('-o', '--output_dir',type=str, help='outputdirectory to save ...bam/ into')
parser.add_argument('-m', '--email_address',type=str,help='cluster email results to this email address')
args = parser.parse_args()

#check ftp file list and make it a dict sample_name:ftp-URL 
if args.ref_path is not None:
    ref_path = args.ref_path
else:
    print('reference fasta not found')
    raise IOError

if args.cram_list is not None:
    cram_list = args.cram_list.split(',')
    samples = {}
    samples = {i.rsplit('/')[-1]:i for i in cram_list}
else:
    print('bam csv list not formed correctly')
    raise IOError

if args.output_dir is not None:
    out_dir = args.output_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
else:
    print('out put directory not specified')
    raise IOError

if args.email_address is not None:
    email = args.email_address
else:
    print('email not valid')
    raise AttributeError
    
#make a temp directory for pbs scripts
pbs_dir = out_dir+'PBS'
if not os.path.exists(pbs_dir): os.makedirs(pbs_dir)
    
#write out the .pbs scripts
PBS = []
for k in samples:
    python = '/home/tbecker/software/anaconda/bin/python'
    varp   = '/home/tbecker/software/SVE/tests/variant_processor.py'
    sample_pbs = pbs_dir+'/'+k+'.pbs' #name of pbs script for sample k
    PBS += [sample_pbs]
    if not os.path.exists(out_dir+k): os.makedirs(out_dir+k)
    with open(sample_pbs,'w') as pbs:
        bam = [python, varp, '-d','jax','-r',ref_path,'-b',samples[k],'-s','cram2bam','-o',out_dir]
        pbs.write('#!/bin/bash\n'+' '.join(bam))
#execute qsub with the scripts, getting the jids back (can display these)
output,err = '',{}
for pbs in PBS:
    print('processing %s'%pbs)
    try:
        command = ['qsub','-l','walltime=24:00:00,mem=12gb','-m','e','-M',email,'-o','cram.log','-j oe',pbs]
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