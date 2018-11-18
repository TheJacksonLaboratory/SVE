#!/usr/bin/env python
import argparse
import glob
import gc
import os
import time
import multiprocessing as mp
#pip installed libs
import numpy as np
#local libs
import svu_utils as su
import read_utils as ru
import fusor_utils as fusor

des = """
FusorSV - SV call (VCF4.0+) fusion framework
Timothy James Becker, PhD candidate, UCONN 05/25/2016-03/18/2017, Version 0.1.2"""
parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-r', '--ref_path',type=str, help='reference fasta needed to write vcf g1k output files\t[None]')
des = """
input directory with sample named folders and caller id tagged vcf files
[EX] /path/sample/sample_S4.vcf implies that there are calls of id 4 for sample\t[None]"""
parser.add_argument('-i', '--in_dir',type=str, help=des)
parser.add_argument('-a', '--ctg_dir',type=str,help='assembly contig directory\t[None]')
parser.add_argument('-c', '--chroms',type=str,help='comma seperated chrom listing\t[1,2,...22,X,Y,MT]')
parser.add_argument('-o', '--out_dir',type=str, help='outputdirectory to save ...bam.bai/ into\t[None]')
parser.add_argument('-f', '--apply_fusion_model_path',type=str, help='apply a model for statistics')
parser.add_argument('-k', '--k_fold',type=int,help='k_fold cross validation, k=0 implies no validation\t[0]')
parser.add_argument('-n', '--n_k',type=int,help='number of times to do k_fold\t[1]')
parser.add_argument('-p', '--cpus',type=int,help='number of cores to use in ||\t[1]')
parser.add_argument('-b', '--bins',type=int,help='number of requested discriminating features\t[14]')
parser.add_argument('-e', '--obs',type=int,help='number of observations needed per bin\t[1000]')
parser.add_argument('-g', '--min_g',type=float,help='minimum group expectation contribution\t[0.0]')
parser.add_argument('-l', '--over_m',type=float,help='overlap allowed before removal in merge step\t[0.0]')
parser.add_argument('-m', '--pre_cluster',action='store_true', help='cluster the calls for all samples first\t[False]')
parser.add_argument('-s', '--brkpt_smoothing',action='store_true', help='brkpt_smoothing algo\t[False]')
parser.add_argument('-u', '--ucsc_filter',action='store_true', help='apply stringent UCSC filter tracks\t[False]')
stage_mapping = """
1:1 mapping of caller ids to stage names (and back):
stage_map_json_file -> {0:'G1K-P3',-1:'fusorSV',1:'MetaSV',4:'BreakDancer',9:'cnMOPS',10:'CNVnator',
                        11:'Delly',13:'GATK',14:'GenomeSTRiP',17:'Hydra',18:'Lumpy',35:'BreakSeq',
                        36:'Pindel',38:'Tigra'}\t[../data/stage_map.json]
"""
parser.add_argument('-S', '--stage_map_json_file',type=str,help=stage_mapping)
parser.add_argument('-E', '--stage_exclude_list',type=str,help='comma seperated id list to exclude from test/training\t[1,13,36]')
parser.add_argument('-M', '--cluster_overlap',type=float,help='reciprocal overlap needed for clustering\t[0.5]')
parser.add_argument('-L', '--lift_over',type=str,help='liftover chain file path or default\t[./data/hg19ToHg38.over.chain.gz]')
parser.add_argument('-C', '--clean', action='store_true',help='keep all kfold run data and print extra details\t[False]')
args = parser.parse_args()

if args.in_dir is not None:
    in_dir = args.in_dir
else:
    print('no VCF input')
    raise IOError
if args.out_dir is not None:
    out_dir = args.out_dir
else:
    print('no output directory specified')
    raise IOError
if args.ctg_dir is not None:
    ctg_dir = args.ctg_dir
else:
    print('no contig directory specified')
    ctg_dir = None
if args.ref_path is not None:
    ref_path = args.ref_path
else:
    ref_path = ''
    print('ref_path not provided, will not be able to write a vcf and g1k file')
write_stats   = True   #new default
write_vcf_g1k = True   #new default
write_model   = True   #new default
if args.k_fold is not None: cross_fold = args.k_fold
else:                       cross_fold = 1
if args.n_k is not None and cross_fold > 1: n_cross = args.n_k
else:                                       n_cross = 1
if args.cpus is not None: cpus = args.cpus
else:                     cpus = 1
if args.bins is not None: beta = args.bins
else:                     beta = 16
if args.obs is not None: obs = args.obs
else:                    obs = 1000
if args.min_g is not None: min_g = args.min_g
else:                      min_g = 0.0
if args.over_m is not None: over_m = args.over_m
else:                       over_m = 1E-10        #1bp overlap in a human genome = 3E9
if args.chroms is not None:
    chroms = args.chrom.split(',')
else:
    chroms = [str(i) for i in range(1,23)]+['X','Y','MT'] #limit to chroms of interest
if args.apply_fusion_model_path is not None: apply_fusion_model_path = args.apply_fusion_model_path
else:                                        apply_fusion_model_path = None
if args.cluster_overlap is not None: cluster_overlap = args.cluster_overlap
else:                                cluster_overlap = 0.0
if args.lift_over is not None:
    if args.lift_over.upper()=='DEFAULT':        #need to do -L default to get default
        lift_over = ru.get_local_path('liftover/hg19ToHg38.over.chain.gz')    #default
    else: lift_over = args.lift_over                 #lift over file for downstream analysis
else: lift_over = None
#load a user defined stage mapping id or use the one provided
if args.stage_map_json_file is not None:
    callers= ru.get_stage_map(args.stage_map_json_file)
    if callers == {}:
        callers = ru.get_stage_map(ru.get_local_path('stage_map.json'))
else:
    callers = ru.get_stage_map(ru.get_local_path('stage_map.json'))
    
stage_exclude_list = [1,36] #default exclude list
if args.stage_exclude_list is not None:
    try:
        stage_exclude_list = [int(i) for i in args.stage_exclude_list.rsplit(',')]
    except Exception as E:
        print('error parsing stage id exclude list')
        print('using default stage id exclude list: %s'%stage_exclude_list)

    # make sure all excluded ids are valid callers
    for c in stage_exclude_list:
        if not c in callers:
            stage_exclude_list = [1, 36]
            print('invalid caller id in stage id exclude list: %d'%c)
            print('using default stage id exclude list: %s'%stage_exclude_list)
            break
else:
    print('using default stage id exclude list: %s'%stage_exclude_list)

result_list = [] #async queue to put results for || stages
def collect_results(result):
    result_list.append(result)

#[1] File Partitioning------------------------------------------------------------------  
#read and convert to SVULTB                                     all SV caller VCF inputs
#saving each partition to a seperate pickle file               || by sample << partition
def partition_call_sets(sample,k,O,R,B,chroms,flt,flt_exclude,caller_exclude):
    sname = sample[sample.rfind('/')+1:]                      #extract sample identifier
    print('reading sample %s'%sname)
    sname_partition_path = out_dir+'/svul/'+sname                            #build path
    S,V = su.vcf_glob_to_svultd(sample+'/*vcf',chroms,O,flt=flt,flt_exclude=flt_exclude)
    S = su.filter_call_sets2(S,R,exclude=flt_exclude)                    #filter svmasks
    Q = fusor.slice_samples([[sname,S]])                                         #legacy
    P = fusor.partition_sliced_samples(Q,B,exclude=caller_exclude)            #partition
    success = fusor.write_partitions_by_sample(sname_partition_path,P)    #write to disk
    return [sname,success]                                   #report back to async queue
#[1] File Partitioning------------------------------------------------------------------

#[2] Pool Partitions--------------------------------------------------------------------
#this is || by partition: t, b
def merge_partitions(callers,t,b):
    partition_path = out_dir+'/svul/'
    P = fusor.read_partitions_by_caller(partition_path,callers,t,b,False)
    return [sname,P]
#[2] Pool Partitions--------------------------------------------------------------------

#[3a] Fit the Model Partition----------------------------------------------------------------------------
def prior_model_partition(snames,t,b,k,callers,caller_exclude,min_g,brkpt_smoothing):
    print('starting model partition:\tt=%s\tb=%s'%(t,b))
    start = time.time()
    P = fusor.read_partitions_by_caller(out_dir+'/svul/',callers,caller_exclude,t,b,False) 
    #P,snames = fusor.pre_cluster_samples(P,r=0.9),['CLUSTER']                       #optional preclustering 
    T = fusor.target_by_sample(P,k)                                          #extract the partitioned target
    J = fusor.all_samples_all_pairs_magnitudes(P,snames)                        #pool all feature magnitudes
    K = fusor.all_samples_all_callers_bkpts(P,snames)                                  #pool all breakpoints
    D,NN = fusor.pooled_distance(J)                       #this is done inside the all_group_weights as well
    W = fusor.all_group_weights(J,k,mode='j')                             #pooled D,NN and all group weights
    E = fusor.select_groups(W,min_g)                                      #gamma is a group selection cutoff
    A = fusor.pileup_group_by_sample(P,E,(k,))                              #now its just one partition here
    if brkpt_smoothing:
        print('optimizing model partition with brkpt smoothing:\tt=%s\tb=%s'%(t,b))  
        alpha = fusor.target_filter_cutoff_exhaustive_brkpt(A,P,T,E,K)
        stop = time.time()
        print('fusion model partition:\tt=%s\tb=%s\t%s sec\tcappa=%s'%(t,b,round(stop-start,2),alpha[t][b]))
    else:
        print('optimizing model partition:\tt=%s\tb=%s'%(t,b))   
        alpha = fusor.target_filter_cutoff_exhaustive(A,E,T)                        #optimal cutoff location 
        stop = time.time()
        print('fusion model partition:\tt=%s\tb=%s\t%s sec\talpha=%s'%(t,b,round(stop-start,2),alpha[t][b]))
    return [(t,b),J,D,E,alpha,K]
#[3a] Fit the Model Partition----------------------------------------------------------------------------

#[3b] Posterior Estimation Partition---------------------------------------------------------------------
def post_model_partition(apply_fusion_model_path,snames,t,b,k,callers,caller_exclude,min_g):
    print('starting posterior estimate on partition:\tt=%s\tb=%s'%(t,b))
    start = time.time()
    #[1] load prior model values----------------------------------------
    B,J,D,E,alpha,n,K = fusor.import_fusion_model(apply_fusion_model_path)      #import the existing model
    #[2] load new input data partitions
    P = fusor.read_partitions_by_caller(out_dir+'/svul/',callers,caller_exclude,t,b,False)   #all samples
    J_new = fusor.all_samples_all_pairs_magnitudes(P,snames)                 #pool all feature magnitudes
    #[3] construct the posterior estimator using:        the prior data, new data and imputed true row==k
    J_post = fusor.additive_magnitude_smoothing(J,J_new,k)        #k is used to swap a row J_prime into J
    D_post,NN_post = fusor.pooled_distance(J_post)                      #get the new data distance matrix
    W_post = fusor.all_group_weights(J_post,k,mode='j')  #calculate the pooled D,NN and all group weights
    E_post = fusor.select_groups(W_post,min_g)                         #gamma is a group selection cutoff 
    alpha_post = fusor.post_filter_cutoff(E,E_post,alpha)                       #updated filter estimates
    stop = time.time()
    print('posterior estimate on partition:\tt=%s\tb=%s\t%s sec\talpha=%s'%(t,b,round(stop-start,2),
                                                                            alpha_post[t][b]))
    return [(t,b),J_post,D_post,E_post,alpha_post,K]
#[3b] Posterior Estimation Partition--------------------------------------------------------------------

#[4] Apply Model To Samples--------------------------------------------------------------------------------
#supports both single model training and posterior estimate via additive smoothing and diagnostics:  mantel
def apply_model_to_samples(sample,ref_path,chroms,types,bins,callers,O,
                           model_path,apply_fusion_model_path,k,f_id,
                           over_m,r=0.5,smoothing=False,detail=False,verbose=False,IDX=6):
    sname = sample[sample.rfind('/')+1:]                                          #extract sample identifier
    print('starting fusorSV discovery on sample %s'%sname)
    ref_seq = {'.'.join(ref_path.rsplit('/')[-1].rsplit('.')[0:-1]):ru.read_fasta(ref_path,True)}      #~3GB
    hist,C = {},{}
    cross_fold_stats,detailed_stats = {c:{} for c in callers},{}                         #each single caller
    for t in types:
        for c in cross_fold_stats:
            if not t in cross_fold_stats[c]: cross_fold_stats[c][t] = []
    B,J,J_post,D,D_post,alpha,alpha_post,K,n,n_post = {},{},{},{},{},{},{},{},0,0
    #[1] apply the model here------------------------------------------------------------------------------
    if apply_fusion_model_path is None:#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        if verbose: print('loading base fusion model partitions for %s'%sname)
        B,J,D,E,alpha,n,K = fusor.import_fusion_model(model_path)                    #import the base model
        P = fusor.read_partitions_by_sample(partition_path,sname)            #read all partitions for sname
        Q = fusor.unpartition_sliced_samples(P)                              #unpartition for merging later
        if verbose: print('projection of all call sets on %s'%sname)
        A = fusor.pileup_group_by_sample(P,E,(k,))                    #projection of all calls for a sample
        F = fusor.filter_pileup_by_sample(A,alpha,leave_in=False)   #filter using optimal cutof in the mode
        if smoothing:
            #now do breakpoint smoothing algorithm---------------------------------------------------------
            F = fusor.best_smooth_brkpt_samples(F,K,P)
            #now do breakpoint smoothing algorithm---------------------------------------------------------
        fusor.merge_filtered_samples(Q,F,f_id,snames,[],over_m)
    else:#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                                                           
        if verbose: print('loading base and posterior estimate partitions for %s'%sname)        
        B,J,D,E,alpha,n,K = fusor.import_fusion_model(apply_fusion_model_path)          #import prior model
        B,J_post,D_post,E_post,alpha_post,n_post,K = fusor.import_fusion_model(model_path) #import new data
        P = fusor.read_partitions_by_sample(partition_path,sname)            #read all partitions for sname
        Q = fusor.unpartition_sliced_samples(P)                              #unpartition for merging later
        A = fusor.pileup_group_by_sample(P,E,(k,))                    #projection of all calls for a sample                    
        F = fusor.filter_pileup_by_sample(A,alpha,E_post,leave_in=False)         #filter cutoff in the mode
        if smoothing:
            #now do breakpoint smoothing algorithm---------------------------------------------------------
            F = fusor.best_smooth_brkpt_samples(F,K,P)
            #now do breakpoint smoothing algorithm---------------------------------------------------------
        fusor.merge_filtered_samples(Q,F,f_id,snames,[],over_m)      #expectation priorty merge back result
    #[2] do some scoring and write out the sample results, returning for global performance----------------
    start = time.time()
    for t in Q:                                  #give the metrics for each sample and write out the result
        for c in set(cross_fold_stats).difference(set([f_id,k])):
            if verbose: print('%s%s'%(callers[c],''.join(['-' for i in range(80)])))   #other callers first
            cross_fold_stats[c][t] += [fusor.pretty_stats(Q,types,t,k,c,sname,r,verbose)]
        if verbose: print('fusorSV%s'%(''.join(['-' for i in range(80)])))                    #then fusorSV
        cross_fold_stats[f_id][t] += [fusor.pretty_stats(Q,types,t,k,f_id,sname,r,verbose)]
        C[t] = []
        if Q[t].has_key(f_id) and Q[t][f_id].has_key(sname):
            C[t] = Q[t][f_id][sname]
            for i in cross_fold_stats[f_id][t][-1][6]:
                C[t][i][IDX][k] = [-1]                             #flag the idx with the target key and -1
    if verbose: print('writing VCF for %s'%sname)
    G = su.svult_to_genome(C,O)                                               #start conversion back to VCF
    hist[sname] = su.genome_to_vcf(G,ref_seq,types,chroms,callers,
                                   out_dir+'/vcf/'+sname+'_S'+str(f_id)+'.vcf',sname,target_key=k)     #VCF
    if detail:
        for c in callers:
            detailed_stats[c] = {}
            for t in types:
                detailed_stats[c][t] = {}
                for i in range(len(B[t])-1):
                    detailed_stats[c][t][bins[t][i]] = []      
        for t in Q:
            for c in set(detailed_stats.keys()).difference(set([f_id,k])):            #do the other callers
                if verbose: print('%s%s'%(callers[c],''.join(['-' for x in range(80)])))
                T = fusor.pretty_bin_stats(Q,types,t,B,bins,k,c,sname,r,verbose)
                for i in range(len(B[t])-1):
                    detailed_stats[c][t][bins[t][i]] = [T[bins[t][i]]]
            if verbose: print('fusorSV%s'%(''.join(['-' for i in range(80)])))
            T = fusor.pretty_bin_stats(Q,types,t,B,bins,k,f_id,sname,r,verbose)                 #no fusorSV
            for i in range(len(B[t])-1):
                detailed_stats[f_id][t][bins[t][i]] = [T[bins[t][i]]]
    
    stop = time.time()
    if verbose: print('scoring completed for %s in %s sec'%(sname,round(stop-start,2)))
    return [sname,cross_fold_stats,hist,detailed_stats]
#[4] Apply Model To Samples--------------------------------------------------------------------------------                 

if __name__ == '__main__':        
    full_start = time.time()

    if not os.path.exists(out_dir):                   os.makedirs(out_dir)
    if not os.path.exists(out_dir+'/tigra_ctg_map/'): os.makedirs(out_dir+'/tigra_ctg_map/')
    if not os.path.exists(out_dir+'/g1k/'):           os.makedirs(out_dir+'/g1k/')
    if not os.path.exists(out_dir+'/vcf/'):           os.makedirs(out_dir+'/vcf/')
    if not os.path.exists(out_dir+'/svul/'):          os.makedirs(out_dir+'/svul/')
    if not os.path.exists(out_dir+'/models/'):        os.makedirs(out_dir+'/models/')
    if not os.path.exists(out_dir+'/visual/'):        os.makedirs(out_dir+'/visual/')
    if not os.path.exists(out_dir+'/visual/bed/'):    os.makedirs(out_dir+'/visual/bed/')
    #will need file preparation documentation or use the full SVE to do this...
    files = glob.glob(in_dir+'*') #get all sample directories
    file_exclude = ['PBS','NA19240'] #27 samples exclude these folders
    samples,snames = [],[]
    for i in files:
        sname = i.rsplit('/')[-1]
        if not sname in file_exclude: #skips ids here
            samples += [i]
            snames  += [i.rsplit('/')[-1]]      
    #snames,samples = snames[0:2],samples[0:2] #testing line
    
    print('processing samples %s\n for chroms %s'%(samples,chroms))  
    #automate the coordinate construction for use with alternate assemblies
    O = ru.get_coordinate_offsets('human_g1k_v37_decoy_coordinates.json') #load the reference offset map
    R = []
    #mask these regions------------------------------------------------------------------------------------
    R += ru.get_mask_regions('human_g1k_v37_decoy_svmask.json',O)               #svmask from ref complexity
#    R += ru.get_mask_regions('human_g1k_v37_decoy_ucsc_known_genes.json',O)          #known gene locations
#    R += ru.get_mask_regions('human_g1k_v37_decoy_ucsc_gap.json',O)                      #gaps in assembly
#    R += ru.get_mask_regions('human_g1k_v37_decoy_ucsc_microsatellite.json',O)            #microsatellites
#    R += ru.get_mask_regions('human_g1k_v37_decoy_ucsc_segmental_dups.json',O)             #segmental dups
#    R += ru.get_mask_regions('human_g1k_v37_decoy_ucsc_simple_repeats.json',O)             #simple repeats
#    R += ru.get_mask_regions('human_g1k_v37_decoy_ucsc_repeat_masker.json',O)        #repeat masker tracks
    print('merging the svmask regions')
    start = time.time()
    R = ru.flatten_mask_regions(R,O,complement=False)                                       #single IRanges
    stop = time.time()
    print('svmask regions merged in %s sec'%(round(stop-start,2)))
    #mask these regions------------------------------------------------------------------------------------
    k,flt,f_id,m_id = 0,0,-1,1             #k=true_id,flt=filter 0 is .,PASS,f_id=fusorSV_id,m_id=metaSV_id
    exclude_callers = stage_exclude_list            #exclude caller options put any id here to exclude
    B = {t:[1,100,250,500,1000,5000,10000,50000,100000,1000000] for t in range(0,8)}
#   B = fusor.distribute_bins(Q,k,n_b=beta,m_b=obs,lower=None,upper=None,event=False) #equal power distribution
    B[1] = [1,50,100,1000,1000000]        
    B[2] = [1,50,100, 400, 600, 950, 1250, 1550, 1950, 2250, 2950, 3650, 4800, 6150, 9000, 18500, 100000, 10000000]
    B[3] = [1,50,500,1000,5000,10000,50000,250000,10000000]
    B[5] = [1,50,100,250,500,1000,2500, 3500, 45000, 80000, 115000, 180000, 260000, 300000, 500000, 1000000]
    types = {0:'SUB',1:'INS',2:'DEL',3:'DUP',4:'CNV',5:'INV',6:'TRA',7:'BND'}
    bins  = {t:su.pretty_ranges(B[t],'') for t in B}
    partition_path = out_dir+'/svul/'
    total_partitions = len(glob.glob(partition_path+'*.pickle'))

#    sname = samples[0][samples[0].rfind('/')+1:]                      #extract sample identifier
#    print('reading sample %s'%sname)
#    sname_partition_path = out_dir+'/svul/'+sname                            #build path
#    S,V = su.vcf_glob_to_svultd(samples[0]+'/*vcf',chroms,O,flt=flt,flt_exclude=[])
#    S = su.filter_call_sets2(S,R,exclude=[])                    #filter svmasks
#    Q = fusor.slice_samples([[sname,S]])                                         #legacy
#    P = fusor.partition_sliced_samples(Q,B,exclude=[])            #partition
#    success = fusor.write_partitions_by_sample(sname_partition_path,P)
    
    
    #||||||||||||||||||||||||||||||||||||||BY SAMPLE|||||||||||||||||||||||||||||||||||||||||||||

    #[1] read, parse, structure, select, partition and write out data for each sample if not done
    if total_partitions<1: #skip this step if you have already read and partitioned the data set      
        print('reading, parsing, partitioning and writing sample VCFs')
        start = time.time()
        p1 = mp.Pool(processes=cpus)
        for sample in samples:
            p1.apply_async(partition_call_sets,
                           args=(sample,k,O,R,B,chroms,flt,[],exclude_callers),
                           callback=collect_results)
            time.sleep(0.25)
        p1.close()
        p1.join()
        L = []
        for i in result_list:
            if i is not None: L+=[i]
        result_list = []
        gc.collect()
        snames = [i[0] for i in L] #passing list of sample names
        #only have to read in the samples once
        stop = time.time()
        total_partitions = len(glob.glob(partition_path+'*.pickle'))
        print('finished reading %s out of %s samples generating %s partitions in %s sec'%\
              (len(snames),len(samples),total_partitions,round(stop-start,2)))
    #||||||||||||||||||||||||||||||||||||||BY SAMPLE|||||||||||||||||||||||||||||||||||||||||||||
    if n_cross>1 and cross_fold>1: #for each run will permute using the k_fold divisor to partition
            print('employing crossfold validation measures...runs=%s\tkfold=%s'%(n_cross,cross_fold))
    for n_k in range(n_cross): #default is 1 and cross_fold = 0
        n_start = time.time()
        tst_ids = []
        if cross_fold>1:
            tst_ids = sorted(list(np.random.choice(range(len(samples)),len(samples)/cross_fold,replace=False)))
        trn_ids = sorted(list(set(range(len(samples))).difference(set(tst_ids))))
        tst_str = ''.join([hex(i)[2:].upper() for i in tst_ids]) #get a id sorted hex string of the ids used
        trn_str = ''.join([hex(i)[2:].upper() for i in trn_ids]) #get a id sorted hex string of the ids used
        
        #||||||||||||||||||||||||||||||||||||||BY PARTITION||||||||||||||||||||||||||||||||||||||||||
        #[2]train or apply the model
        #load the data and build a model if one isn't already available in ||
        if apply_fusion_model_path is None: #train a new model assuming k_fold=0 here
            model_path = out_dir+'/models/'+'.'.join(ref_path.rsplit('/')[-1].rsplit('.')[0:-1])+\
                         '.'+in_dir[0:-1].rsplit('/')[-1]+trn_str+'.pickle'
            if not os.path.exists(model_path): #write a model if it hasn't been done yet
                start = time.time()            #now in || for faster performance and less RAM
                p1 = mp.Pool(processes=cpus)        
                for t in types:
                    for b in range(len(B[t])-1):
                        p1.apply_async(prior_model_partition,
                                       args=([snames[i] for i in trn_ids],t,b,k,
                                             callers,exclude_callers,min_g,False),
                                       callback=collect_results)
                        time.sleep(0.25)
                p1.close()
                p1.join()
                L = []
                for i in result_list:
                    if i is not None: L+=[i]
                result_list = []
                gc.collect()
                stop = time.time()
                print('finished modeling in %s sec'%round(stop-start,2))    
                J,D,E,alpha,n,K = fusor.assemble_model(L)
                fusor.export_fusion_model(B,J,D,E,alpha,len(snames),K,model_path)
                L = []
                gc.collect()
            else: #can clip this one the || sample application is completed
                B,J,D,E,alpha,n,K = fusor.import_fusion_model(model_path)
            
        else: #apply an existing model and then use additive smoothing with the all new input data
            B,J,D,E,alpha,n,K = fusor.import_fusion_model(apply_fusion_model_path)
            model_path = model_path = out_dir+'/models/'+'.'.join(ref_path.rsplit('/')[-1].rsplit('.')[0:-1])+\
                         '.'+in_dir[0:-1].rsplit('/')[-1]+trn_str+'.post.pickle'
            #now look at the new data and make a model for it, minus the true (estimate it)
            if not os.path.exists(model_path): #write a model if it hasn't been done yet
                start = time.time()            #now in || for faster performance and less RAM
                p1 = mp.Pool(processes=cpus)        
                for t in types:
                    for b in range(len(B[t])-1):
                        p1.apply_async(post_model_partition,
                                       args=(apply_fusion_model_path,snames,t,b,k,
                                             callers,exclude_callers,min_g),
                                       callback=collect_results)
                        time.sleep(0.25)
                p1.close()
                p1.join()
                L = []
                for i in result_list:
                    if i is not None: L+=[i]
                result_list = []
                gc.collect()
                stop = time.time()
                print('finished estimation in %s sec'%round(stop-start,2))    
                J_post,D_post,E_post,alpha_post,n_post,K = fusor.assemble_model(L)
    
                fusor.export_fusion_model(B,J_post,D_post,E_post,alpha_post,len(snames),K,model_path)
                L = []
                gc.collect()
            else: #can clip this one the || sample application is completed
                B,J_post,D_post,E_post,alpha_post,n_post,K = fusor.import_fusion_model(model_path)
        #||||||||||||||||||||||||||||||||||||||BY PARTITION||||||||||||||||||||||||||||||||||||||||||
              
        #||||||||||||||||||||||||||||||||||||||BY SAMPLE|||||||||||||||||||||||||||||||||||||||||||||
        print('apply fusion model to sample inputs and generating fusorSV ouput')
        if n_cross>1 and cross_fold>1:
            print('scoring the crossfold run %s out of...runs=%s\tkfold=%s'%(n_k,n_cross,cross_fold))
        else:
            tst_ids,tst_str = trn_ids,trn_str
        start = time.time()
        p1 = mp.Pool(processes=max(1,cpus/2))
        for sample in [samples[i] for i in tst_ids]:
            p1.apply_async(apply_model_to_samples,
                           args=(sample,ref_path,chroms,types,bins,callers,O,
                                 model_path,apply_fusion_model_path,k,f_id,
                                 over_m,0.5,args.brkpt_smoothing,True,True,6),
                           callback=collect_results)
            time.sleep(0.25)
        p1.close()
        p1.join()
        L = []
        for i in result_list:
            if i is not None: L+=[i]
        result_list = []
        gc.collect()
        #only have to read in the samples once
        stop = time.time()
        print('finished reading samples in %s sec'%round(stop-start,2))
        #||||||||||||||||||||||||||||||||||||||BY SAMPLE|||||||||||||||||||||||||||||||||||||||||||||         
        cross_fold_stats,hist,detailed_stats = fusor.assemble_stats(L)
        ref_seq = {'.'.join(ref_path.rsplit('/')[-1].rsplit('.')[0:-1]):ru.read_fasta(ref_path,True)}
        #compute cross_fold averages
        #if not apply_fusion_model_path is None: 
        for c in cross_fold_stats:
            print('%s--------------------------------------------------------------'%callers[c])
            for t in cross_fold_stats[c]:
                if len(cross_fold_stats[c][t])>0:#have all samples here to plot/look at!
                    #idea[1] look at scatterplot with a line
                    prec = round(np.mean([i[2] for i in cross_fold_stats[c][t]]),2)
                    rec  = round(np.mean([i[3] for i in cross_fold_stats[c][t]]),2)
                    f1   = round(np.mean([i[4] for i in cross_fold_stats[c][t]]),2)
                    j    = round(np.mean([i[5] for i in cross_fold_stats[c][t]]),2)
                    n    = round(np.mean([i[7] for i in cross_fold_stats[c][t]]),2)
                    m    = round(np.mean([i[8] for i in cross_fold_stats[c][t]]),2)
                    l_mu = round(np.mean([i[9] for i in cross_fold_stats[c][t]]),2)
                    r_mu = round(np.mean([i[10] for i in cross_fold_stats[c][t]]),2)
                    print('average for t=%s\tprec=%s\trec=%s\tf1=%s\tj=%s\tn=%s\tm=%s\tl_mu=%s\tr_mu=%s'%\
                          (types[t],prec,rec,f1,j,n,m,l_mu,r_mu))
        #CHECK-SCORES-----------------------------------------------------------------------------------------------
        fusor.export_caller_performance(cross_fold_stats,callers,
                                        out_dir+'/visual/cross_fold_stats.'+tst_str+'.tsv')
        fusor.export_detailed_performance(detailed_stats,callers,
                                          out_dir+'/visual/detailed_stats.'+tst_str+'.tsv')
        fusor.export_caller_by_type_and_bin(E,alpha,callers,types,bins,
                                            out_dir+'/visual/callers_tbj.'+tst_str+'.tsv')
        fusor.export_distance_matrix(D,callers,types,bins,
                                     out_dir+'/visual/sim_matrix.'+tst_str+'.tsv',sim=True)
        #(a) check for an Rscript engine
        #(b) use the command parser to fire up the Rscript and check
        #(c) for the drawing libraries to ggplot, ect to use
        #-------------------------------------------------------------------

        vcf_glob = out_dir+'/vcf/*_S-1.vcf' #fusorSV VCFS only
        # now cluster the calls
        if cluster_overlap > 0.0:
            print('cluster merging tool processing samples')
            out_vcf  = out_dir+'/vcf/all_samples_genotypes.'+tst_str+'.vcf'
            print('completed = %s'%su.fusorSV_vcf_multi_sample_merge(vcf_glob,out_vcf,
                                                                     ref_seq[ref_seq.keys()[0]],
                                                                     overlap=cluster_overlap))
        if lift_over is not None:
            #now do a liftover to partition the calls with a possible new reference space
            su.fusorSV_vcf_liftover_samples(out_dir+'/vcf/all_samples_genotypes*.vcf*',ref_path,lift_over) #default is on
        if n_cross>1 and cross_fold>1 and not args.clean: #can clean up the VCF and model files...
            print('cleaning interim data for run %s out of...runs=%s\tkfold=%s'%(n_k, n_cross,cross_fold))
            os.remove(model_path)
            for vcf in glob.glob(vcf_glob): os.remove(vcf)
        n_stop = time.time()
        print('run %s in %s sec'%(n_k,round(n_stop-n_start,2)))
        print(''.join([':::' for i in range(40)]))
