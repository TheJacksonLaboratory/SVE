#!/usr/bin/env python
#BAM input metacalling SV pipeline
#given a reference directory (same name as the refrence) that already contains all of the
#preprossesing steps such as indexing, masking, chrom spliting, etc...

import argparse
import os
import sys
import glob
import socket
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import stage_utils as su
import read_utils as ru
import svedb
import stage

def path(path):
    return os.path.abspath(path)[:-1]

#[1]parse command arguments
des = """
Shortened Version of the SVE that uses pre-made .bam files
Allthough .bam files are not compatible with all callers such as Variation Hunter"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-o', '--out_dir', type=str, help='output directory to store resulting files')
parser.add_argument('-r', '--ref', type=str, help='fasta reference file path')
bam_help = """
BAM comma-sep bam file path list, assuming matching bam.bai is in place...
[EX] --bam ~/data/sample1.bam,~/data/sample2.bam,~/data/sample3.bam"""
parser.add_argument('-b', '--bams',type=str, help=bam_help)
parser.add_argument('-D', '--read_depth',type=int, help='Average Depth of Coverage in ROI')
parser.add_argument('-L', '--read_length',type=int, help='Read Length')
parser.add_argument('-t', '--targets',type=str, help='target region files for input: vcf,bd,1gk formats')
parser.add_argument('-s','--stages',type=str, help='stage name list')
parser.add_argument('-c','--chroms',type=str, help='chrom name list')
parser.add_argument('-p','--partitions',type=int, help='number of partitions or ||L')
parser.add_argument('--debug',action='store_true',help='save result/error data to db')
parser.add_argument('-d','--database',type=str, help='database configuration file')
parser.add_argument('-e','--erase_db',action='store_true',help='reset and clear the schema bound db')
#
parser.add_argument('-v','--verbose',action='store_true',help='be verbose with caller stdout/stderr')
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

#db clearing option
with svedb.SVEDB(dbc['srv'], dbc['db'], dbc['uid'], dbc['pwd']) as dbo:
    if args.erase_db: dbo.new()  # reset db if needed

host = socket.gethostname()
directory = path('~/'+host+'/') #users base home folder as default plus hostname

if args.stages is not None:
    sids = []
    with svedb.SVEDB(dbc['srv'], dbc['db'], dbc['uid'], dbc['pwd']) as dbo:
        stage_meta = su.get_stage_meta()
        sids = su.get_stage_name_id(stage_meta)
#if args.meta_call.upper()=='ALL':
#need to check the stage_id and that type like 'variant%'
    try: #get just the callers name:stage_id
        staging = {c:sids[c] for c in args.stages.split(',')}
        print('processing with : %s'%(sorted(staging.keys()),))
    except Exception:
        print('unknown processor name used as input argument: %s'%args.stages)
        print('availble stages are:\n------------------------------\n%s\n' % '\n'.join(sorted(sids.keys())))
        raise KeyError
else:
    stage_meta = su.get_stage_meta()
    sids = su.get_stage_name_id(stage_meta)
    print('missing value for caller stage_id list')
    print('availble stages are:\n------------------------------\n%s\n' % '\n'.join(sorted(sids.keys())))
    raise AttributeError

if args.out_dir is not None:    #optional reroute
    directory = args.out_dir
if not os.path.exists(directory): os.makedirs(directory)

if args.bams is not None:
    bams = args.bams.split(',') #CSL returns a list of 1+
else:
    print('bam input file error')
    raise IOError 
    
#collect read depth and length estimates for RD analysis
RD,RL = None,None
if args.read_depth is not None:
    RD = args.read_depth #CSL returns a list of 1+
else:
    print('read_depth not given')
if args.read_length is not None:
    RL = args.read_length #CSL returns a list of 1+
else:
    print('read_length not given')
#collect read depth and length estimates for RD analysis
if not(RD is None) and not(RL is None):
    auto_RD_RL = False
else:
    auto_RD_RL = True

if args.targets is not None:
    targets = args.targets.split(',') #CSL returns a list of 1+
else:
    targets = None
    print('no targets')

if args.ref is not None and os.path.exists(args.ref):
    ref_fa_path = args.ref
    #names = ru.get_fasta_seq_names(ref_fa_path) #this will throw an exception
else:
    #escape this if bam_split is the only stage
    noref = staging.has_key('bam_split_simple') or staging.has_key('bam_split_all') or \
            staging.has_key('bam_clean') or staging.has_key('bam_stats') or\
            staging.has_key('samtools_merge') or staging.has_key('picard_merge')
    if len(staging)<=1 and noref:
        ref_fa_path = '/missing_ref.fa'
    else:
        print('reference fasta error')
        raise IOError
refbase = ref_fa_path.rsplit('/')[-1].split('.fa')[0]

if args.chroms is not None:
    chroms = args.chroms.split(',') #CSL returns a list of 1+
    #allow some kind of ALL tag for future use
else:
    print('chr selection not used scanning reference')
    if args.ref is not None and ref_fa_path != '/missing_ref.fa':
        chroms = ru.get_fasta_seq_names(ref_fa_path)
    elif len(staging)<=1 and noref:
        print('using only the bam file')
        #have to have a bam file and read the header here
    else:
        print('missing chrom and reference')
        raise AttributeError

if args.partitions is not None:
    partitions = args.partitions
else:
    partitions = 1

if args.verbose is not None:
    verbose = True
#[3] start DB and execute the Pipeline

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
    dbo.new_run('illumina',host,ref_id,debug=args.debug)
    run_id = dbo.get_max_key('runs')
    print('starting run_id = %s'%run_id)
    
    print('processing bam list')
    print(bams)

    if staging.has_key('cram2bam'):
        st = stage.Stage('cram2bam',dbc)
        outs = st.run(run_id,{'.fa':[ref_fa_path],'.cram':bams,'out_dir':[directory]})
        if verbose: print(outs)
    
    if staging.has_key('cram2bam_split_all'):
        st = stage.Stage('cram2bam_split_all',dbc)
        split_params = st.get_params()
        split_params['p']['value'] = 1
        st.set_params(split_params)
        outs = st.run(run_id,{'.fa':[ref_fa_path],'.cram':bams,
                              'chroms':chroms,'out_dir':[directory]})
        if verbose: print(outs)
    
    #good for QC or for auto estimating RD,RL
    #leaving RD callers blank with imply bam_stats first
    #bam_stats will imply bam_clean to ensure full SV compatibilty, IE GATK, GenomeSTRiP will fail...
    if staging.has_key('bam_clean') or staging.has_key('bam_stats') or\
       staging.has_key('cnmops') and auto_RD_RL or \
       staging.has_key('cnvnator') and auto_RD_RL or \
       staging.has_key('genome_strip') and auto_RD_RL:
           
        #check for the *_S3 files first
        in_stats = glob.glob(directory+'*_S'+sids['bam_stats'])
        h,v = False,False
        for i in range(len(in_stats)):
            if in_stats[i].endswith('.header'): h = True
            if in_stats[i].endswith('.valid'):  v = True 
        if h and v: #run bam_clean
            st = stage.Stage('bam_clean',dbc)
            bam_clean_params = st.get_params()
            bam_clean_params['-t'] = 4 #default threads
            st.set_params(bam_clean_params)
            outs = st.run(run_id,{'.header':[in_stats[0]],'.valid':[in_stats[1]],
                                  '.bam':bams,'out_dir':[directory]})
        else:
            print('---------------running stats and conditional cleaning-------------------')
            st = stage.Stage('bam_stats',dbc)
            outs = st.run(run_id,{'.bam':bams,'out_dir':[directory]})
            if not type(outs) is None and len(outs)>0:
                outs = outs[0].split('read statistics\n')[-1].split('\n')
            try: #pull out the positions here
                RD = int(round(float(outs[1].split(' = ')[-1]),0))  #average depth
                RL = int(round(float(outs[24].split(' = ')[-1]),0)) #average length
            except Exception:
                print('--------------RD,RL not determined, setting default values-----------------')
                RD,RL = 30,100
            st = stage.Stage('bam_clean',dbc)
            in_stats = glob.glob(directory+'*'+sids['bam_stats']+'.header') +\
                       glob.glob(directory+'*'+sids['bam_stats']+'.valid')
            header,valid = [],[]
            for i in range(len(in_stats)):
                if in_stats[i].endswith('.header'): header = in_stats[i]
                if in_stats[i].endswith('.valid'):  valid  = in_stats[i]
            if len(header)>0 and len(valid)>0: #run bam_clean
                bam_clean_params = st.get_params()
                bam_clean_params['-t'] = 4 #default threads
                bam_clean_params['-m'] = 8 #default memory
                st.set_params(bam_clean_params)
                #print the .valid file for debugging-------
                print('-------------valid file output-----------------')
                with open(valid) as f: print(f.readlines())
                outs += st.run(run_id,{'.header':[header],'.valid':[valid],
                                       '.bam':bams,'out_dir':[directory]})
                #:::TO DO::: check on last time for validation...
        if verbose: print(outs)

    if staging.has_key('bam2cram'):
        st = stage.Stage('bam2cram',dbc)
        outs = st.run(run_id,{'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)
    
    if staging.has_key('bam_split_all'):
        st = stage.Stage('bam_split_all',dbc)
        outs = st.run(run_id,{'.bam':bams,'chroms':chroms,'out_dir':[directory]})
        if verbose: print(outs)
            
    if staging.has_key('bam_split_simple'):
        st = stage.Stage('bam_split_simple',dbc)
        outs = st.run(run_id,{'.bam':bams,'chroms':chroms,'out_dir':[directory]})
        if verbose: print(outs)
    
    if staging.has_key('samtools_index'):
        st = stage.Stage('samtools_index',dbc)
        outs = st.run(run_id,{'.bam':bams})
        if verbose: print(outs)
            
    if staging.has_key('samtools_merge'):
        st = stage.Stage('samtools_merge',dbc)
        outs = st.run(run_id,{'.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)
            
    if staging.has_key('picard_merge'):
        st = stage.Stage('picard_merge',dbc)
        outs = st.run(run_id,{'.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)
        
    if staging.has_key('lumpy'):
        #lumpy    
        st = stage.Stage('lumpy',dbc)
        outs = st.run(run_id, {'.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)

    if staging.has_key('cnmops'):
        st = stage.Stage('cnmops',dbc)
        cnmops_params = st.get_params()
        cnmops_params['window']['value'] = ru.expected_window(depth=RD,length=RL,target=100)
        if len(bams)<=1:
            cnmops_params['mode']['value']   = 3  
        elif len(bams)==2:
            cnmops_params['mode']['value']   = 1
        else:
            cnmops_params['mode']['value']   = 0
        cnmops_params['normal']['value'] = 3      #poisson normalization
        cnmops_params['cir_seg']['value'] = True
        cnmops_params['cores']['value']  = 1      #bibc threads
        st.set_params(cnmops_params)
        outs = st.run(run_id, {'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})    
        if verbose: print(outs)
    
    if staging.has_key('delly'):     
        #delly
        st = stage.Stage('delly',dbc)
        outs = st.run(run_id, {'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)
    
    if staging.has_key('breakdancer'):
        #breakdancer
        st = stage.Stage('breakdancer',dbc)
#        bd_params = st.get_params()
#        bd_params['-l']['value'] = True
#        st.set_params(bd_params)
        outs = st.run(run_id, {'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)    

    if staging.has_key('breakseq'):
        #breakseq
        brkptlib_path = '/'.join(ref_fa_path.rsplit('/')[0:-1])+'/'+refbase+sids['breakseq']
#        print(brkptlib_path)
        #check that it exists and swap it out if need be....
        st = stage.Stage('breakseq',dbc)
        bs_params = st.get_params()
        bs_params['window']['value']   = 2*RL
        bs_params['junction']['value'] = 4*RL
        st.set_params(bs_params)
        outs = st.run(run_id, {'.fa':[ref_fa_path],'.gff':[brkptlib_path+'.brkptlib.gff'],
                               '.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)
    
    if staging.has_key('hydra'):
        #hydra multi
        st = stage.Stage('hydra',dbc)
        outs = st.run(run_id, {'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
        if args.verbose: print(outs)
  
    if staging.has_key('cnvnator'):
        #file dump issue related to paramiko ENV variable and the root system in cnvnator
        #cnvnator VC
        st = stage.Stage('cnvnator',dbc)
        cnvnator_params = st.get_params()    #automatically get the depth and length
        cnvnator_params['window']['value'] = ru.expected_window(depth=RD,length=RL,target=100)
        st.set_params(cnvnator_params)
        outs = st.run(run_id, {'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)
    
    if staging.has_key('genome_strip'):
        #genomestrip
        gs_ref_path = '/'.join(ref_fa_path.rsplit('/')[0:-1])+'/'+refbase+sids['genome_strip_prepare_ref']
        st = stage.Stage('genome_strip',dbc)
        outs = st.run(run_id, {'.fa':[gs_ref_path+'.fa'],
                               '.fa.svmask.fasta':[gs_ref_path+'.fa.svmask.fasta'],
                               '.ploidymap.txt':  [gs_ref_path+'.ploidymap.txt'],
                               '.rdmask.bed': [gs_ref_path+'.rdmask.bed'],
                               '.gcmask.fasta': [gs_ref_path+'.gcmask.fasta'],
                               '.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)

    if staging.has_key('gatk_haplo'):
        #gatk Haplotyper VC
        st = stage.Stage('gatk_haplo',dbc)
        outs = st.run(run_id,{'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
        if verbose: print(outs)
        #gatk VairantReCalibration using dbsnp, known:Hapmap, 1Kgenomes, etc...
    
    if staging.has_key('tigra') and targets is not None:
        st = stage.Stage('tigra',dbc)
        tigra_params = st.get_params()     #automatically get the depth and length
        tigra_params['p']['value'] = 1     #cpus for sorting, leave at 1  
        tigra_params['F']['value'] = 0     #bedIntersect flanking bp
        tigra_params['L']['value'] = 1000  #bp away from the breakpoint
        tigra_params['A']['value'] = 1000  #bp into the breakpoint
        tigra_params['Q']['value'] = 2     #min quality
        tigra_params['P']['value'] = 1000  #max read depth
        tigra_params['H']['value'] = 200   #max nodes
        st.set_params(tigra_params)
        outs = st.run(run_id,{'.fa':[ref_fa_path],'.bam':bams,
                              '.calls':targets,'.vcf':targets,'out_dir':[directory]})
        if verbose: print(outs)
#    need to clear out the directory when completed here...
#    
#    if callers.has_key('variationhunter'):
    #variation hunter
#    st = stage.Stage('variationhunter',dbc)
#    outs = st.run(1,{'.DIVET.vh': [base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S26.DIVET.vh',
#                                   base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S26.DIVET.vh']})
    #pindel
    
    #svseq2
    
    #gindel
    
    #svelter

#now you can run FusorSV for one sample if needed                       