#!/usr/bin/env python
#BAM input metacalling pipeline
#given a reference directory (same name as the refrence) that already contains all of the
#preprossesing steps such as indexing, masking, chrom spliting, etc...

#simple_sve.py ref_path bam_path out_dir
import argparse
import os
import sys
import glob
import socket
import time
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import stage_utils as su
import svedb
import stage
import subprocess32 as subprocess

def path(path):
    return os.path.abspath(path)[:-1]

def Dedup_Sort(in_bam, mem, threads):
    # Mark duplications
    st = stage.Stage('picard_mark_duplicates',dbc)
    params = st.get_params()
    params['-t']['value'] = threads
    params['-m']['value'] = mem
    st.set_params(params)
    dedup_bam = st.run(run_id,{'.bam':[sorted_bam]})
    if (dedup_bam == False):
        print "ERROR: picard_mark_duplicates fails"
    else:
        subprocess.call(["mv",dedup_bam,sorted_bam])
        # Sort bam
        st = stage.Stage('sambamba_index',dbc)
        st.set_params(params)
        st.run(run_id,{'.bam':[sorted_bam]})


#[1]parse command arguments
des = """
Processes Illumina PE .fq reads into a .bam file given a reference .fa file as input.
[Note] mem, sub-pipelines assume reads > 75bp in length"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-a','--algorithm', type=str, help='aln|mem|speed_seq')
parser.add_argument('-g', '--replace_rg',action='store_true', help='replace reads groups\t[False]')
parser.add_argument('-m', '--mark_duplicates',action='store_true', help='mark duplicate reads\t[False]')
parser.add_argument('-s', '--sample',type=str, help='sample name\t[input]')
parser.add_argument('-o', '--out_dir', type=str, help='output directory to store resulting files\t[None]')
parser.add_argument('-r', '--ref', type=str, help='fasta reference file path\t[None]')
parser.add_argument('-d','--database',type=str, help='database configuration file\t[SVE/data]')
parser.add_argument('-A','--realign',action='store_true', help='Realign')
fqs_help = """
fq comma-sep file path list\t[None]
[EX PE] --fqs ~/data/sample1_FWD.fq,~/data/sample1_REV.fq"""
parser.add_argument('-f', '--fqs',type=str, help=fqs_help)
parser.add_argument('-b', '--bam',type=str, help='bam file path\t[None]')
parser.add_argument('-t','--threads',type=int, help='number of threads per CPU\t[4]')
parser.add_argument('-M','--mem',type=int, help='ram in GB units to use for processing per cpu/thread unit\t[4]')
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
    algorithm = 'speed_seq'
    
#take in bam file(s) run
with svedb.SVEDB(dbc['srv'], dbc['db'], dbc['uid'], dbc['pwd']) as dbo:
    dbo.embed_schema()
    print('\n<<<<<<<<<<<<<USING HOST %s>>>>>>>>>>>>>>>=\n')%host
    print('using reference name = %s'%refbase)
    ref_id = -1
    try:
        ref_id = dbo.get_ref_id(refbase)
    except IndexError:
        print('unkown reference: run starting with -1')
    print('using ref_id=%s'%str(ref_id))
    dbo.new_run('illumina',host,ref_id)
    run_id = dbo.get_max_key('runs')
    print('starting run_id = %s'%run_id)
    
    stage_meta = su.get_stage_meta()
    ids = su.get_stage_name_id(stage_meta)
   
    if args.bam is not None:
        if args.realign:
            a_start = time.time()
            aligner_params = {'.fa':[ref_fa_path],'.bam':bam,'out_dir':[directory]}
            st = stage.Stage('speedseq_realign',dbc)
            aligner_stage_params = st.get_params()
            aligner_stage_params['-t']['value'] = threads
            aligner_stage_params['-m']['value'] = mem
            st.set_params(aligner_stage_params)
            outs = st.run(run_id,aligner_params)
            a_stop = time.time()
            print('SVE:picard_mark_duplicates time was % hours'%round((a_stop-a_start)/(60**2),1))
        if args.mark_duplicates:
            d_start = time.time()
            st = stage.Stage('picard_mark_duplicates',dbc)
            outs = st.run(run_id,{'.bam':[bam]})
            d_stop = time.time()
            print('SVE:picard_mark_duplicates time was % hours'%round((d_stop-d_start)/(60**2),1))
        if args.replace_rg: #set to do one at a time only for now...
            r_start = time.time()
            base = bam.rsplit('/')[-1].rsplit('.')[0].rsplit('_')[0].rsplit('-')
            if SM is None: SM = base
            st = stage.Stage('picard_replace_rg',dbc)
            outs = st.run(run_id,{'.bam':[bam],
                                  'platform_id':['illumina'],
                                  'SM':[SM]})
            r_stop = time.time()
            print('SVE:picard_replace_rg time was %s sec'%round((r_stop-r_start)/(60**2),1))
    else:
        a_start = time.time()
        base = su.get_common_string_left(reads).rsplit('/')[-1].rsplit('.')[0]
        if SM is None: SM = base
        aligner_params = {'.fa':[ref_fa_path],'.fq':reads,'platform_id':['illumina'],'SM':[SM],'out_dir':[directory]}
        # Set the stage's parameters
        st = stage.Stage('fq_to_bam_piped',dbc)
        aligner_stage_params = st.get_params()
        aligner_stage_params['-t']['value'] = threads
        aligner_stage_params['-m']['value'] = mem
        if algorithm == 'mem':
            st = stage.Stage('fq_to_bam_piped',dbc)
            st.set_params(aligner_stage_params)
            # outs will receive ".sorted.bam"
            sorted_bam = st.run(run_id,aligner_params)
            Dedup_Sort(sorted_bam, mem, threads)
        elif algorithm == 'speed_seq':
            st = stage.Stage('speedseq_align',dbc)
            st.set_params(aligner_stage_params)
            outs = st.run(run_id,aligner_params)
        elif algorithm == 'aln':
            st = stage.Stage('bwa_aln',dbc)
            st.set_params(aligner_stage_params)
            # outs will receive ".sorted.bam"
            sorted_bam = st.run(run_id,aligner_params)
            Dedup_Sort(sorted_bam, mem, threads)
            
        a_stop = time.time()
        print('SVE:BAM:%s was completed in %s hours'%(algorithm,round((a_stop-a_start)/(60.0**2),4)))
