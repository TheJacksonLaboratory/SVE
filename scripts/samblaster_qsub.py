#!/usr/bin/env python
import argparse
import os
import glob
import subprocess32 as subprocess #to call qsub a bunch of times

#parse commandline arguments for usage
des = """
Automated PBS job resequencing (speedseq) generator for low coverage g1k sample data"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-b', '--bam_dir',type=str, help='input data directory bams')
parser.add_argument('-o', '--out_dir',type=str, help='outputdirectory to use with samblaster')
parser.add_argument('-w', '--wall_time',type=str,help='wall time requested from cluster')
parser.add_argument('-m', '--memory',type=str,help='radom access memory needed')
parser.add_argument('-p', '--cpus',type=str,help='|| processors needed')
parser.add_argument('-e', '--email_address',type=str,help='cluster email results to this email address')
args = parser.parse_args()

if args.bam_dir is not None:
    bam_dir = args.bam_dir
else:
    print('no input')
    raise IOError
if args.out_dir is not None:
    out_dir = args.out_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
else:
    print('no output')
    raise IOError

if args.wall_time is not None:
    walltime = args.wall_time
else:
    walltime = '12:00:00'
if args.memory is not None:
    ram = args.memory
else:
    ram = '64gb'
if args.cpus is not None:
    cpus = args.cpus
else:
    cpus = '16'
if args.email_address is not None:
    email = args.email_address
else:
    print('email not valid')
    raise AttributeError

samples = {}
mapped_bams   = glob.glob(bam_dir+'/*.bam')
print(mapped_bams)
samples = {i.rsplit('/')[-1].split('.')[0]:[i] for i in mapped_bams}

#make a temp directory for pbs scripts
PBS,pbs_dir = [],out_dir+'PBS'
if not os.path.exists(pbs_dir): os.makedirs(pbs_dir)
software_path = os.path.dirname(os.path.abspath(__file__)) + '/../../'
samtools   = software_path+'/samtools-1.2/samtools'
samblaster = software_path+'/speedseq/bin/samblaster'
for k in samples:
    job_pbs = pbs_dir+'/'+k+'.realign.pbs'
    PBS += [job_pbs]
    ml = 'module load'
    modules =    [ml+' samtools/1.2\n'+ml+' gcc/4.9.2\n']
    #samtools view -h samp.bam | samblaster [-a] [-e] [-d samp.disc.sam] [-s samp.split.sam] [-u samp.umc.fasta] -o /dev/null    
    blast = [samtools,'view','-Sh',samples[k][0],'|',samblaster,'-e',
             '-d',out_dir+'/'+k+'.disc.sam',
             '-s',out_dir+'/'+k+'.split.sam'
             '-u',out_dir+'/'+k+'.um.fa',
             '-o','/dev/null\n',
             samtools,'view','-Sb',out_dir+'/'+k+'.disc.sam','>',out_dir+'/'+k+'.disc.bam\n',
             samtools,'index',out_dir+'/'+k+'.disc.bam\n',
             samtools,'view','-Sb',out_dir+'/'+k+'.split.sam','>',out_dir+'/'+k+'.split.bam\n',
             samtools,'index',out_dir+'/'+k+'.disc.bam\n',
             'rm',out_dir+'/'+k+'.disc.sam',out_dir+'/'+k+'.split.sam\n']
    with open(job_pbs,'w') as pbs:
        pbs.write('#!/bin/bash\n'+\
                  ' '.join(modules)+'\n'+\
                  ' '.join(blast)+'\n')
#execute qsub with the scripts, getting the jids back (can display these or attach to further monitor progress)
output,err = '',{}
for pbs in PBS: #test with one of these and a fast caller on a small file...
    print('processing %s'%pbs) #mark something to save
    try:
        command = ['qsub','-l','walltime=%s,mem=%s,procs=%s'%(walltime,ram,cpus),'-m','e','-M',email,
                   '-o',pbs[0:-4]+'.log','-j oe',pbs]
        print(' '.join(command))
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
