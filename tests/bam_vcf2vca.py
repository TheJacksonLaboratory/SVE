#BAM input metacalling pipeline
#given a reference directory (same name as the refrence) that already contains all of the
#preprossesing steps such as indexing, masking, chrom spliting, etc...

#simple_sve.py ref_path bam_path out_dir
import argparse
import os
import sys
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
Shortened Version of the SVE that uses pre-made .bam files
Allthough .bam files are not compatible with all callers such as Variation Hunter"""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-o', '--out_dir', type=str, help='output directory to store resulting files')
parser.add_argument('-r', '--ref', type=str, help='fasta reference file path')
bam_help = """
BAM comma-sep bam file path list, assuming matching bam.bai is in place...
[EX] --bam ~/data/sample1.bam,~/data/sample2.bam,~/data/sample3.bam"""
parser.add_argument('-b', '--bams',type=str, help=bam_help)
args = parser.parse_args()

#[2] apply and bind args to variables
host = socket.gethostname()
directory = path('~/'+host+'/') #users base home folder as default plus hostname
if args.out_dir is not None:    #optional reroute
    directory = args.out_dir
if not os.path.exists(directory): os.makedirs(directory)
if args.ref is not None:
    ref_fa_path = args.ref
else:
    raise IOError
refbase = ref_fa_path.rsplit('/')[-1].split('.fa')[0]
bams = args.bams.split(',') #CSL returns a list of 1+
if args.bams is not None:
    bams = args.bams.split(',') #CSL returns a list of 1+
else:
    raise IOError
    
#[3] start DB and execute the Pipeline
#try a new schema here?....???
srv = 'arc-gis.ad.engr.uconn.edu'
db  = 'sve' #make a new DB for human tests?
uid = 'sv_calibrator'
pwd = 'sv_calibrator'
#take in bam file(s) run
with svedb.SVEDB(srv, db, uid, pwd) as dbo:
    dbo.embed_schema()
    print('\n<<<<<<<<<<<<<USING HOST>>>>>>>>>>>>>>>=%s\n')%host
    print('using reference name = %s'%refbase)
    ref_id =  dbo.get_ref_id(refbase)  
    print('using ref_id=%s'%str(ref_id))
    dbo.new_run('illumina',host,ref_id)
    run_id = dbo.get_max_key('runs')
    print('starting run_id = %s'%run_id)
    
    stage_meta = su.get_stage_meta()
    ids = su.get_stage_name_id(stage_meta)
    
    print('processing bam list')
    print(bams)
    
    #samtools_bam_stats
    st = stage.Stage('bam_stats')
    outs = st.run(run_id,{'.bam':bams,'out_dir':[directory]})
    print(outs)

    #gatk Haplotyper VC
    st = stage.Stage('gatk_haplo')
    outs = st.run(run_id,{'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
    print(outs)
    #gatk VairantReCalibration using dbsnp, known:Hapmap, 1Kgenomes, etc...
   
    #delly
    st = stage.Stage('delly')
    outs = st.run(run_id, {'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
    print(outs)
    
    #breakdancer
    st = stage.Stage('breakdancer')
    outs = st.run(run_id, {'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
    print(outs)    
    
    #lumpy    
    st = stage.Stage('lumpy')
    outs = st.run(run_id, {'.bam':bams,'out_dir':[directory]})
    print(outs)
    
#    #hydra multi
#    st = stage.Stage('hydra')
#    outs = st.run(run_id, {'.bam':bams,'out_dir':[directory]})
#    print(outs)
    
    #file dump issue related to paramiko ENV variable and the root system in cnvnator
    #cnvnator VC
    st = stage.Stage('cnvnator')
    cnvnator_params = st.get_params()
    cnvnator_params['window']['value'] = 100
    st.set_params(cnvnator_params)
    outs = st.run(run_id, {'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})
    print(outs)
    
    #cnmops VC still not handling no data output issues fully
    #cnmops.R needs some more robust error handling:
    #Reference sequence:  1
    #Error in checkForRemoteErrors(val) :
    #  one node produced an error: missing value where TRUE/FALSE needed
    #Calls: singlecn.mops ... clusterApply -> staticClusterApply -> checkForRemoteErrors
    #Execution halted
    
    st = stage.Stage('cnmops')
    cnmops_params = st.get_params()
    cnmops_params['window']['value'] = 100
    if len(bams)<=1:
        cnmops_params['mode']['value']   = 3  
    elif len(bams)==2:
        cnmops_params['mode']['value']   = 1
    else:
        cnmops_params['mode']['value']   = 0
    cnmops_params['normal']['value'] = 3      #poisson normalization
    cnmops_params['cir_seg']['value'] = False
    cnmops_params['cores']['value']  = 1      #bibc threads
    st.set_params(cnmops_params)
    outs = st.run(run_id, {'.fa':[ref_fa_path],'.bam':bams,'out_dir':[directory]})    
    print(outs)
    
    #genomestrip
    gs_ref_path = '/'.join(ref_fa_path.rsplit('/')[0:-1])+'/'+refbase+ids['genome_strip_prepare_ref']
    st = stage.Stage('genome_strip')
    outs = st.run(run_id, {'.fa':[gs_ref_path+'.fa'],
                           '.fa.svmask.fasta':[gs_ref_path+'.fa.svmask.fasta'],
                           '.ploidymap.txt':  [gs_ref_path+'.ploidymap.txt'],
                           '.rdmask.bed': [gs_ref_path+'.rdmask.bed'],
                           '.bam':bams,'out_dir':[directory]})
    print(outs)
    
#    need to clear out the directory when completed here...
#    
    #variation hunter
#    st = stage.Stage('variationhunter')
#    outs = st.run(1,{'.DIVET.vh': [base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S26.DIVET.vh',
#                                   base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S26.DIVET.vh']})
    #pindel have this one
    
    #svseq2 have this one
    
    #gindel
    
    #WES specific VC
    #ExomeDepth
    #contra
    #xmmer
    
    #somatic mutation callers...
    #snpsniper VC
    #mutec VC

    #combined tech callers
    #multi-break-SV

    #do a select all on table t...
    #data = dbo.select_all(t) #bool=>bit(1) are coming back as int, datetime.datetime are datetime...
                     
    #sharded for each node...NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN                         