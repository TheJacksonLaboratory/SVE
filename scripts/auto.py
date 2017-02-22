#!/usr/bin/env python

#auto.py --ref_path|--ref_file --fqs --out_dir --sample --database
#prep and variant calling steps___________________________________________________________________________________
#(A) --pr_cpus --ref_path (optional if indexes are not found)
#(B) --pb_cpus --pb_threads --pb_mem | uses the bwa_split alogorithm
#(C) --vp_cpus --vp_threads --vp_mem | default ordering with a few || steps
#(D) [1A: breakdancer, 1B: breakseq, 1C: cnmops, 2A:cnvnator, 2B:hydra, 4 lumpy, 5 delly, 6 GS, 7 GATK (optional)]
#merging steps----------------------------------------------------------------------------------------------------
#(E) --fsv_cpus --in_dir --out_dir --model
#(F) --vp_cpus --vp_threads --vp_mem [8 tigra-ext ith fsv input for --target]
#(G) --fsv_cpus --in_dir --contig_dir --model --cross_map_chain --out_dir
#----------------------------------------------------------------------------------------------------------------- 
import argparse
import os
import sys
import glob
import socket
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import stage_utils as su
import svedb
import stage

def path(path):
    return os.path.abspath(path)[:-1]

#[1]parse command arguments
des = """
Autonomous Structural Variation Engine:
Given a .fa reference file and a pair: NA12878_1.fq.gz,NA12878_2.fq.gz, 
produce a FusorSV VCF file all_samples_merge.vcf with comprehensive high quality SV calls.
#[USAGE] auto.py --ref_path|--ref_file --fqs --out_dir --sample --database
#prep and variant calling steps___________________________________________________________________________________
#(A) --pr_cpus --ref_path (optional if indexes are not found)
#(B) --pb_cpus --pb_threads --pb_mem | uses the bwa_split alogorithm
#(C) --vp_cpus --vp_threads --vp_mem | default ordering with a few || steps
#(D) [1A: breakdancer, 1B: breakseq, 1C: cnmops, 2A:cnvnator, 2B:hydra, 4 lumpy, 5 delly, 6 GS, 7 GATK (optional)]
#merging steps----------------------------------------------------------------------------------------------------
#(E) --fsv_cpus --in_dir --out_dir --model
#(F) --vp_cpus --vp_threads --vp_mem [8 tigra-ext ith fsv input for --target]
#(G) --fsv_cpus --in_dir --contig_dir --model --cross_map_chain --out_dir
#----------------------------------------------------------------------------------------------------------------- """
parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-o', '--out_dir', type=str, help='output directory to store resulting files\t[None]')
parser.add_argument('-r', '--ref', type=str, help='fasta reference path (if indexes are not found, run (A))\t[None]')
parser.add_argument('-d','--database',type=str, help='database configuration file\t[SVE/data]')
fqs_help = """
fq comma-sep file path list\t[None]
[EX PE] --fqs ~/data/sample1_FWD.fq,~/data/sample1_REV.fq"""
parser.add_argument('-f', '--fqs',type=str, help=fqs_help)
parser.add_argument('-m', '--model',type=str, help='data fusion model\t[FusorSV/data/models/g1k_v37_decoy.P3.pickle]')
parser.add_argument('-P','--cpus',type=str, help='number of cpus for alignment and sorting, ect\t[1]')
parser.add_argument('-T','--threads',type=str, help='number of threads per CPU\t[4]')
parser.add_argument('-M','--mem',type=str, help='ram in GB units to use for processing per cpu/thread unit\t[4]')
args = parser.parse_args()

#read the database configuration file
dbc = {'srv':'','db':'','uid':'','pwd':''}
if args.database is not None:
    with open(args.database, 'r') as f:
        params = f.read().split('\n') #newline seperated configuration file
    try:
        dbc['srv']  = params[0].split('srv=')[-1]
        dbc['db'] = params[0].split('srv=')[-1]
        dbc['uid'] = params[0].split('srv=')[-1]
        dbc['pwd'] = params[0].split('srv=')[-1]
    except Exception:
        print('invalid database configuration')
        print('running the SVE without the SVEDB')
        pass
else:
    with open(os.path.dirname(os.path.abspath(__file__))+'/../data/svedb.config', 'r') as f:
        params = f.read().split('\n') #newline seperated configuration file
    try:
        dbc['srv'] = params[0].split('srv=')[-1]
        dbc['db']  = params[1].split('db=')[-1]
        dbc['uid'] = params[2].split('uid=')[-1]
        dbc['pwd'] = params[3].split('pwd=')[-1]
        schema = {}
        with svedb.SVEDB(dbc['srv'], dbc['db'], dbc['uid'], dbc['pwd']) as dbo:
            dbo.embed_schema()   #check the schema for a valid db
            schema = dbo.schema
        if len(schema)<1:
            print('dbc:%s' % [c + '=' + dbc[c] for c in dbc])
            print('invalid database configuration')
            print('running the SVE without the SVEDB')
        else:
            print('dbc:%s'%[c+'='+dbc[c] for c in dbc])
            print('valid database configuration found')
            print('running the SVE with the SVEDB')
    except Exception:
        print('invalid database configuration')
        print('running the SVE without the SVEDB')
        pass

host = socket.gethostname()
directory = path('~/'+host+'/') #users base home folder as default plus hostname
if args.out_dir is not None:    #optional reroute
    directory = args.out_dir
if not os.path.exists(directory): os.makedirs(directory)
if args.ref is not None:
    ref_fa_path = args.ref
else:
    print('no ref pattern found')
    ref_fa_path = ''

refbase = ref_fa_path.rsplit('/')[-1].split('.fa')[0]
if args.fqs is not None:
    reads = args.fqs.split(',') #CSL returns a list of 1+
else:
    print('no fqs pattern found')
    reads = []
if not all([os.path.exists(r) for r in reads]):
    print('fastq files not found!')
    raise IOError
if args.bam is not None:
    bam = args.bam
else:
    print('no bam pattern found')
    bam = ''

if args.sample is not None:
    SM = args.sample
else:
    SM = None

if args.cpus is not None:
    cpus = int(args.cpus)
else:
    cpus = 1
if args.threads is not None:
    threads = args.threads
else:
    threads = 4
if args.mem is not None:
    mem = int(args.mem)
else:
    mem = 4
if args.algorithm is not None:
    algorithm = args.algorithm
else:
    algorithm = 'piped_split'
    
#(A) check for the reference preparations in that directory

#(B) prepare bam files using the fastest method

#(C) start as many variant callers as you can (in series ...)

#
    
     
