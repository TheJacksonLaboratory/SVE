import argparse
import os
import math
import socket
import random
import hashlib
import subprocess32 as subprocess #to call qsub a bunch of times

#parse commandline arguments for usage
des = """
Automated PBS job varsim generator for multi sample human genome illumina PE simulation"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-r', '--ref_path',type=str, help='reference fasta needed to write vcf g1k output files')
parser.add_argument('-i', '--varsim_data_path',type=str, help='varsim data directory with variation files')
parser.add_argument('-o', '--out_dir',type=str, help='outputdirectory to save into')
parser.add_argument('-c', '--total_coverage',type=int, help='target coverage for WGS simulation')
parser.add_argument('-n', '--nlanes',type=int, help='number of sequencing lanes used to achieve coverage depth')
parser.add_argument('-R', '--read_length',type=int, help='length of the reads')
parser.add_argument('-M', '--mean_fragment_size',type=int, help='mean fragment(insert) size')
parser.add_argument('-S', '--sd_fragment_size',type=int, help='standard deviation of fragement(insert) size')

parser.add_argument('-j', '--j',type=int,help='number of job samples to simulate\t[1]')
parser.add_argument('-w', '--wall_time',type=str,help='wall time requested from cluster')
parser.add_argument('-m', '--memory',type=str,help='radom access memory needed')
parser.add_argument('-e', '--email_address',type=str,help='cluster email results to this email address')
args = parser.parse_args()

if args.out_dir is not None:
    out_dir = args.out_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
else:
    print('no output')
    raise IOError
if args.ref_path is not None:
    ref_path = args.ref_path
else:
    ref_path = ''
    print('ref_path not provided, will not be able to write a vcf and g1k file')
if args.varsim_data_path is not None:
    data_path = args.varsim_data_path
else:
    ref_path = ''
    print('varsim_data_path not provided, will not be able to simulate')
    
if args.total_coverage is not None: total_coverage = args.total_coverage
else:                               total_coverage = 1

if args.nlanes is not None: nlanes = args.nlanes
else:                       nlanes = 1

if args.read_length is not None: read_length = args.read_length
else:                            read_length = 100

if args.mean_fragment_size is not None: mean_fragment_size = args.mean_fragment_size
else:                                   mean_fragment_size = 350

if args.sd_fragment_size is not None: sd_fragment_size = args.sd_fragment_size
else:                                 sd_fragment_size = 50

if args.j is not None: jobs = args.j
else:                  jobs = 1

if args.wall_time is not None:
    walltime = args.wall_time
else:
    walltime = '12:00:00'
if args.memory is not None:
    ram = args.memory
else:
    ram = '32gb'
if args.email_address is not None:
    email = args.email_address
else:
    print('email not valid')
    raise AttributeError

def get_identifier(length=1000):
    l = int(round(math.log(length),0))+1
    return ''.join(random.sample(hashlib.md5(socket.gethostname()).hexdigest(),min(l,hashlib.md5().digestsize)))
    
#make a temp directory for pbs scripts
pbs_dir = out_dir+'PBS'
pe_reads_dir = out_dir+'pe_reads/'
vcf_dir = out_dir+'vcf/'
if not os.path.exists(pbs_dir): os.makedirs(pbs_dir)
if not os.path.exists(pe_reads_dir): os.makedirs(pe_reads_dir)
if not os.path.exists(vcf_dir): os.makedirs(vcf_dir)
#write out the .pbs scripts
PBS = []
python = '/home/tbecker/software/anaconda/bin/python'
varsim = '/home/tbecker/software/varsim/varsim.py'
art    = '/home/tbecker/software/art_bin_VanillaIceCream/art_illumina'
params = ['--vc_in_vcf',data_path+'/All.vcf.gz','--sv_insert_seq',data_path+'/insert_seq.txt',
          '--sv_dgv', data_path+'/GRCh37_hg19_supportingvariants_2013-07-23.txt','--reference',ref_path,
          '--read_length',str(read_length),'--mean_fragment_size',str(mean_fragment_size),
          '--sd_fragment_size',str(sd_fragment_size),'--nlanes',str(nlanes),'--total_coverage',str(total_coverage),
          '--vc_num_snp 0','--vc_num_ins 0','--vc_num_del 0',
          '--vc_num_mnp 0','--vc_num_complex 0','--vc_percent_novel 0',
          '--vc_min_length_lim 0','--vc_max_length_lim 49',
          '--sv_num_ins 500','--sv_num_del 3000','--sv_num_dup 500','--sv_num_inv 1000','--sv_percent_novel 0.9',
          '--sv_min_length_lim 50', '--sv_max_length_lim 500000']
for j in range(jobs): #add --id sample
    s_id = get_identifier(jobs*10).upper()
    job_pbs  = pbs_dir+'/job_'+s_id+'.pbs'
    sample_dir = out_dir+'/'+s_id+'/'
    if not os.path.exists(sample_dir): os.makedirs(sample_dir)
    log_dir  = sample_dir+'/log/' 
    if not os.path.exists(log_dir): os.makedirs(log_dir)
    work_dir = sample_dir+'/work/'
    if not os.path.exists(work_dir): os.makedirs(work_dir)
    PBS += [job_pbs]
    with open(job_pbs,'w') as pbs:
        seed = str(int(s_id,base=16))
        run = [python,varsim]+params+['--seed',seed,'--simulator_executable',art,'--simulator','art','--id',s_id,
                                      '--out_dir',sample_dir,'--log_dir',log_dir,'--work_dir',work_dir]
        clean     = ['rm','-rf',log_dir,work_dir] #clean out the temp files
        rename1   = ['mv',sample_dir+'lane0.read1.fq.gz',sample_dir+s_id+'_1.fq.gz']
        rename2   = ['mv',sample_dir+'lane0.read2.fq.gz',sample_dir+s_id+'_2.fq.gz']
        movereads = ['mv',sample_dir+s_id+'*fq.gz',pe_reads_dir]
        movevcf   = ['mv',sample_dir+'*.truth.*',vcf_dir]
        rm        = ['rm','-rf',sample_dir] #:::TO DO::: test this one first
        #:::TO DO mv and clean the output files to only those that are needed
        pbs.write('#!/bin/bash\n'+'\n'+ \
                  ' '.join(run)+'\n'+ \
                  ' '.join(clean)+'\n'+ \
                  ' '.join(rename1) + '\n' + \
                  ' '.join(rename2) + '\n' + \
                  ' '.join(movereads) + '\n' + \
                  ' '.join(movevcf) + '\n')
#execute qsub with the scripts, getting the jids back (can display these or attach to further monitor progress)
output,err = '',{}
for pbs in PBS: #test with one of these and a fast caller on a small file...
    print('processing %s'%pbs) #mark something to save
    try:
        command = ['qsub','-l','walltime=%s'%walltime+',mem=%s'%ram,'-m','e','-M',email,'-o',pbs[0:-4]+'.log','-j oe',pbs]
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

#remove/delete the intermediate pbs scripts TO DO...