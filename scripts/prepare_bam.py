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
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import stage_utils as su
import svedb
import stage

def path(path):
    return os.path.abspath(path)[:-1]

#[1]parse command arguments
des = """
Processes .fq reads into a .bam file given a reference sequence as input"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-a','--algorithm', type=str, help='aln|mem|piped_mem')
parser.add_argument('-g', '--replace_rg',action='store_true', help='replace reads groups')
parser.add_argument('-m', '--mark_duplicates',action='store_true', help='mark duplicate reads')
parser.add_argument('-s', '--sample',type=str, help='sample name')
parser.add_argument('-o', '--out_dir', type=str, help='output directory to store resulting files')
parser.add_argument('-r', '--ref', type=str, help='fasta reference file path')
parser.add_argument('-d','--database',type=str, help='database configuration file')
fqs_help = """
fq comma-sep file path list
[EX PE] --fqs ~/data/sample1_FWD.fq,~/data/sample1_REV.fq"""
parser.add_argument('-f', '--fqs',type=str, help=fqs_help)
parser.add_argument('-b', '--bam',type=str, help='bam file path')
parser.add_argument('-P','--cpus',type=str, help='number of cpus for alignment and sorting, ect')
parser.add_argument('-M','--mem',type=str, help='ram in GB units to use for processing per cpu')
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
if args.mem is not None:
    mem = int(args.mem)
else:
    mem = 4
    
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
    
    #:::TO DO::: make args the pick which stage -s bwa_mem,samtools_sort,picard_replace_rg,samtools_index
    #:::TO DO::: defailt is all
    if args.mark_duplicates:
       st = stage.Stage('picard_mark_duplicates',dbc)
       outs = st.run(run_id,{'.bam':[bam]})
    if args.replace_rg: #set to do one at a time only for now...
        base = bam.rsplit('/')[-1].rsplit('.')[0]
        if SM is None: SM = base
        st = stage.Stage('picard_replace_rg',dbc)
        outs = st.run(run_id,{'.bam':[bam],
                              'platform_id':['illumina'],
                              'SM':[SM]})
    else:
        if args.algorithm == 'mem': #standard 75+bp illuminam PE read
            base = su.get_common_string_left(reads).rsplit('/')[-1].rsplit('.')[0]
            if not os.path.exists(directory+base+'.sam'):
                if SM is None: SM = base
                #alignment, mapping and add read groups
                st = stage.Stage('bwa_mem',dbc)
                bwa_mem_params = st.get_params()
                bwa_mem_params['-t']['value'] = cpus
                st.set_params(bwa_mem_params)
                outs = st.run(run_id,{'.fa':[ref_fa_path],'.fq':reads,
                                      'platform_id':['illumina'],
                                      'SM':[SM],
                                      'out_dir':[directory]}) 
            if not os.path.exists(directory+base+'.bam'): #look for the finished sorted bam file
                #picard sam to bam, sort, mark duplicates and index
                st = stage.Stage('picard_sam_convert',dbc)
                outs = st.run(run_id,{'.sam':[directory+base+'.sam']})
        elif args.algorithm == 'piped_mem':
            base = su.get_common_string_left(reads).rsplit('/')[-1].rsplit('.')[0]
            if SM is None: SM = base
            st = stage.Stage('fq_to_bam_piped',dbc)
            bwa_mem_params = st.get_params()
            bwa_mem_params['-t']['value'] = cpus
            bwa_mem_params['-m']['value'] = mem
            st.set_params(bwa_mem_params)
            outs = st.run(run_id,{'.fa':[ref_fa_path],'.fq':reads,
                                  'platform_id':['illumina'],
                                  'SM':[SM],
                                  'out_dir':[directory]})
        elif args.algorithm == 'aln':
            #index gapped/ungapped alnment files
            st = stage.Stage('bwa_aln',dbc)
            bwa_aln_params = st.get_params()
            bwa_aln_params['-t']['value'] = cpus
            st.set_params(bwa_aln_params)
            #this is more of the cascading code we want to use:::::::::::::::::::::::::::::::::::::::
            #fix the JSON so that outs has more or less the cascaded file you want to pass into the next stage
            for read in reads:
                outs = st.run(run_id,{'.fa':[ref_fa_path],
                                      '.fq':[read]})
                #print('stage d params: '+str(st.params))
            #paired end alignment :::TO DO::: need to update bwa_sampe to add the read groups                      
            st = stage.Stage('bwa_sampe',dbc)
            l,r = reads[0].strip('.fq'),reads[1].strip('.fq')
            outs = st.run(run_id,{'.fa':[ref_fa_path],
                                  '.fq':[l+'.fq',r+'.fq'],
                                  '.sai':[l+'.sai',r+'.sai']}) 
            
            #picard conversion and RG additions needed for gatk VC...
            st = stage.Stage('picard_sam_convert',dbc)
            outs = st.run(run_id,{'.sam':[l+ids['bwa_sampe']+'.sam']})
    
            """
            #divet conversion
            st = stage.Stage('mrfast_divet')
            outs = st.run(run_id,{'.fa': [base+'node_data/'+ref_base+'.fa'],
                                  'L.fq':[l+'.fq'],'R.fq':[r+'.fq']}) #same as r
            #print(outs)""" 
