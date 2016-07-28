#svedb_test.py
import os
import sys
import socket
import subprocess32 as subprocess
import time
import HTSeq as ht
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import generate_variants as gv
import read_utils as sr
import svedb
import stage

def path(path):
    return os.path.abspath(path)[:-1]

#srv = 'arc-gis.ad.engr.uconn.edu'
srv = 'arc-gis.ad.engr.uconn.edu'
db  = 'sve'
uid = 'sv_calibrator'
pwd = 'sv_calibrator'

with svedb.SVEDB(srv, db, uid, pwd) as dbo:
    #setup shared and sharded file paths and directories
    host = socket.gethostname()
    directory = os.path.dirname(os.path.abspath('~'))+'/node_data/'#+socket.gethostname()+'/'
    
    rnd_test = True
    #ref_base = 'hg19_chr13'
    ref_base = 'rg1'    
    ref_fa_path = directory+ref_base+'.fa'
    mut_fa_path = directory+host+'/'+ref_base+'_R1_CASE.fa'
    ctrl_fa_path = directory+host+'/'+ref_base+'_R1_CTRL.fa'
    #vca_path = directory+host+'/rg1__R1.pkl'
    if not os.path.exists(directory): os.makedirs(directory)
    if not os.path.exists(directory+host): os.makedirs(directory+host)
    #move over real refs when needed at this step
    if not rnd_test:
        mv_path = os.path.dirname(os.path.abspath('~'))+'/refs/hg19/hg19_chr13.fa*'
        subprocess.check_output(' '.join(['cp',mv_path,directory]),shell=True)
    #generate test data0000000000000000000000000000000000000000000000000000000000000000000000000
    
    #make a new random reference fasta------------------------------------------------------
    if rnd_test:
        #can have multiple chroms
        start = time.time()
        chrom_names = ['chr'+str(i) for i in range(1,5)]
        #chrom_names += ['chrX','chrY']             #can simulate sex/ploidy here...
        L = int(1E5)  #L is sum(l1,l2,...ln)   
        chroms = gv.gen_chrom_lens(L,chrom_names) #{'chr1':200,'chr2':100, ... ,'chrn':10}
        print(chroms)
        seqs = gv.gen_random_seqs(chroms,method='fast')
        seqs = sorted(seqs,key=lambda x: len(x.seq),reverse=True)
        print('random reference written: %s'%sr.write_fasta(seqs,ref_fa_path))
        print('random reference written: %s'%sr.write_fasta(seqs,ctrl_fa_path))    
        stop = time.time()
        print('real cpu time: %s seconds'%(stop-start))
        #make a new random reference fasta------------------------------------------------------    
    else:
        #use an existing fasta and make edits on it.............................................
        seqs = sr.read_fasta(directory+ref_base+'.fa')
        #use an existing fasta and make edits on it.............................................
    print('control fasta written: %s'%sr.write_fasta(seqs,ctrl_fa_path))     
    
    #read a ref and generate mutations,limited to one mutation type/len per run==========================    
    mut_len,mut_rate,mut_type,muts,vcas = int(1E2),0.1,'DEL',[],[]#{'DUP':8}  #mut len, rate, and type
    chroms = sr.read_fasta(ref_fa_path,dictionary=True) #read all at once into memory->dict...
    chr_lens = [len(chroms[k]) for k in chroms]
    print('checking mut pos requests...')
    print({k:len(chroms[k]) for k in chroms})
    for k in chroms:  #for each chrom: 'chr1','chr2',...'chrN'        
        mut_pos = gv.gen_mut_pos(len(chroms[k]),mut_len,mut_rate)
        vca = gv.gen_var_calls(ref_fa_path,k,mut_len,mut_type,mut_pos)
        mut_pos = [vc.pos for vc in vca] #correct the pos_list filtering out overlapps, etc...
        vcs = ''
        try:
            vcs = gv.apply_var_calls(chroms[k].seq,vca)
        except KeyError as e:
            print('variant encoding mismatch, application aborted...')        
        mut = [gv.str2seq(vcs, k)]
        muts += mut
        vcas += vca
        print('applied mut_lens match target: %s\n'%gv.check_mut_lens(chroms[k],vca,mut_type,mut))
        print('chr_len=%s, mut_mag=%s, mut_len=%s, mut_rate=%s'%(len(chroms[k]),len(mut_pos),mut_len,mut_rate))
    sr.write_fasta(muts,mut_fa_path)    
    if type(mut_type) is dict: mut_type = mut_type.keys()[0] #swap to is key for db entry FIX ME LATER
    #generate test data0000000000000000000000000000000000000000000000000000000000000000000000000
    print('\n:::::: %s VARIATIONS SYNTHESIZED ::::::\n'%len(vcas))
    
    #start db entriesDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
    #drop and create all
    dbo.new() #reset everything here...
    dbo.embed_schema() #embed the db schema into the dbo object
    #get table dict
    tables = dbo.select_tables()
    for t in tables['names']: print t #good you can get table names
    #get field info on a table
    t = 'runs'#tables['names'][2]
    fields = dbo.select_fields(t) #you can get field info from a table name
    print fields
    
    s = dbo.obj_to_blob(vcas,True)
    sil = [1]
    silb = dbo.obj_to_blob([1],True)
    
    t = 'refs'
    rsb = dbo.obj_to_blob(chroms)
    v = {'ref_id':1,'name':ref_base,'ref_len':sum(chr_lens),
         'seq_names':','.join(chroms.keys()),
         'seq_lens':','.join([str(l) for l in chr_lens]),'seqs':rsb,'url':''} 
    dbo.insert(t,v)
    
    t = 'runs'
    start = dbo.time() #
    v = {'run_id':1,'platform_id':'illumina','node_id':socket.gethostname(),'calibrating':True,
         'ref_id':1,'mut_mag':len(vcas),'mut_len':mut_len,'mut_type':mut_type,
         'mut_true_vc':s,'stage_id_list':silb,'stage_depth':len(sil),'start':start}
    dbo.insert(t,v)
    #time.sleep(5)
    #start db entriesDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


    #do this once for each ref ... RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR    
    #test out stages here...
    base = path('~')
    #split by chrom here and index for random access...
    ss = sr.read_fasta(base+'node_data/'+ref_base+'.fa')
    chrom_names = sr.write_fasta_by_chrom(ss,base+'node_data','')
    ss_len = len(ss)
    ss = [s.name for s in ss] #swap data for names
    for s in ss:
        #index for alignment
        st = stage.Stage('bwa_index')
        outs = st.run(1,{'.fa':[base+'node_data/'+s+'.fa']})
        print(outs)
        st = stage.Stage('picard_dict')
        outs = st.run(1,{'.fa':[base+'node_data/'+s+'.fa']})
        print(outs)
    
    #index for alignment
    st = stage.Stage('bwa_index')
    outs = st.run(1,{'.fa':[base+'node_data/'+ref_base+'.fa']})
    st = stage.Stage('picard_dict')
    outs = st.run(1,{'.fa':[base+'node_data/'+ref_base+'.fa']})
    """
    #genomestrip index and genome masking
    st = stage.Stage('genome_strip_prepare_ref')
    outs = st.run(1, {'.fa':[base+'node_data/'+ref_base+'.fa']})
    print(outs)
    """
    
    #mrfast/variation hunter fasta ref indexing
    st = stage.Stage('mrfast_index')
    outs = st.run(1,{'.fa':[base+'node_data/'+ref_base+'.fa']})
    print(outs)
    
    #do this once for each ref ... RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
    

    #sharded for each node...NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN     
    #PREPARE BAM FILES BAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAM    
    #CASE DATA SET------------------------------------------------------------------------    
    #simulation 
    st = stage.Stage('art_illumina')
    art_params = st.get_params()
    art_params['-f']['value'] = 5
    art_params['-l']['value'] = 75
    art_params['-m']['value'] = 300
    art_params['-s']['value'] = 20
    st.set_params(art_params)
    outs = st.run(1,{'.fa':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE.fa']})
    reads = []
    print("\nouts= %s\n"%(outs,))
    for o in outs:
        if o[-3:]=='.fq': reads.append(o)
    print(reads)
    print('stage a params: '+str(st.params))
    
    #index gapped/ungapped alnment files
    st = stage.Stage('bwa_aln')
    bwa_aln_params = st.get_params()
    bwa_aln_params['-t']['value'] = 12
    st.set_params(bwa_aln_params)
    for read in reads:
        outs = st.run(1,{'.fa':[ref_fa_path],
                         '.fq':[read]})
        #print('stage d params: '+str(st.params))
    #paired end alignment
    st = stage.Stage('bwa_sampe')
    l,r = base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1L',base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1R'
    outs = st.run(1,{'.fa':[base+'node_data/'+ref_base+'.fa'],
                     '.fq':[l+'.fq',r+'.fq'],'.sai':[l+'.sai',r+'.sai']}) 
    
    #initial compression
    st = stage.Stage('samtools_view')
    outs = st.run(1,{'.sam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1.sam']})
    #print(outs)    
    outs = st.run(1,{'.sam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4.sam']})
    #print(outs)
    #sorting of compressed entries
    st = stage.Stage('samtools_sort')
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S5.bam']})
    #print(outs)
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S5.bam']})
    #print(outs)
    #indexing of compressed and indexed .bam
    st = stage.Stage('samtools_index')
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S5.bam']})
    #print(outs)
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S5.bam']})
    #print(outs)
    
    #picard conversion and RG additions needed for gatk VC...
    st = stage.Stage('picard_sam_convert')
    outs = st.run(1,{'.sam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4.sam']})
    #print(outs)    
    st = stage.Stage('picard_replace_rg')
    outs = st.run(1,{'noRG.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8noRG.bam'],
                     'platform_id':['illumina']})
    st = stage.Stage('picard_index')
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8.bam']})
    
    #divet conversion
    st = stage.Stage('mrfast_divet')
    outs = st.run(1,{'.fa': [base+'node_data/'+ref_base+'.fa'],
                     'L.fq':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1L.fq'],
                     'R.fq':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1R.fq']})
    #print(outs)
    #CASE DATA SET------------------------------------------------------------------------ 
    
    #CTRL DATA SET------------------------------------------------------------------------
    #simulation 
    st = stage.Stage('art_illumina')
    st.set_params(art_params)
    outs = st.run(1,{'.fa':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL.fa']})
    reads = []
    print("\nouts= %s\n"%(outs,))
    for o in outs:
        if o[-3:]=='.fq': reads.append(o)
    print(reads)
    print('stage a params: '+str(st.params))
    
    #index gapped/ungapped alnment files
    st = stage.Stage('bwa_aln')
    st.set_params(bwa_aln_params)
    for read in reads:
        outs = st.run(1,{'.fa':[ref_fa_path],
                         '.fq':[read]})
        #print('stage d params: '+str(st.params))
    #paired end alignment
    st = stage.Stage('bwa_sampe')
    l,r = base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1L',base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1R'
    outs = st.run(1,{'.fa':[base+'node_data/'+ref_base+'.fa'],
                     '.fq':[l+'.fq',r+'.fq'],'.sai':[l+'.sai',r+'.sai']}) 
    
    #initial compression
    st = stage.Stage('samtools_view')
    outs = st.run(1,{'.sam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1.sam']})
    #print(outs)    
    outs = st.run(1,{'.sam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4.sam']})
    #print(outs)
    #sorting of compressed entries
    st = stage.Stage('samtools_sort')
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S5.bam']})
    #print(outs)
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S5.bam']})
    #print(outs)
    #indexing of compressed and indexed .bam
    st = stage.Stage('samtools_index')
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S5.bam']})
    #print(outs)
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S5.bam']})
    #print(outs)
    
    
    #picard conversion and RG additions needed for gatk VC...
    st = stage.Stage('picard_sam_convert')
    outs = st.run(1,{'.sam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4.sam']})
    #print(outs)    
    st = stage.Stage('picard_replace_rg')
    outs = st.run(1,{'noRG.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S8noRG.bam'],
                     'platform_id':['illumina']})
    #print(outs)
    st = stage.Stage('picard_index')
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S8.bam']})
    
    #divet conversion
    st = stage.Stage('mrfast_divet')
    outs = st.run(1,{'.fa': [base+'node_data/'+ref_base+'.fa'],
                     'L.fq':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1L.fq'],
                     'R.fq':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1R.fq']})
    #print(outs)
    #CTRL DATA SET------------------------------------------------------------------------ 
    #PREPARE BAM FILES BAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAMBAM 
    
    #samtools_bam_stats
    st = stage.Stage('bam_stats')
    outs = st.run(1,{'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S5.bam']})
    print('CHECK STATS ON SIMULATED VERSUS KNOWN DELLY TEST:')
    print('---------------SIMULATED------------------------')
    print(outs)
    """
    st = stage.Stage('bam_stats')
    outs = st.run(1,{'.bam':[base+'/software/delly-master/test/DEL.bam']})
    print('---------------DELLY TEST-----------------------')
    print(outs)
    """
    
    """
    #samtools VC
    st = stage.Stage('samtools_snp')
    outs = st.run(1,{'.fa':[base+'node_data/'+ref_base+'.fa'],
                     '.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S5.bam']})
    print(outs)
    st = stage.Stage('vcftools_filter')
    outs = st.run(1,{'.vcf':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S5_S11.vcf']})
    print(outs)
        
    #gatk Haplotyper VC
    st = stage.Stage('gatk_haplo')
    outs = st.run(1,{'.fa':[base+'node_data/'+ref_base+'.fa'],
                     '.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8.bam']})
    print(outs)
    #gatk VairantReCalibration using dbsnp, known:Hapmap, 1Kgenomes, etc...
          
    #delly
    st = stage.Stage('delly')
    outs = st.run(1, {'.fa':[base+'node_data/'+ref_base+'.fa'],
                      '.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S8.bam',
                              base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8.bam'],
                      'out_dir':[base+'node_data/'+host+'/']})
    print(outs)
    
    #breakdancer
    st = stage.Stage('breakdancer')
    outs = st.run(1, {'.fa':[base+'node_data/'+ref_base+'.fa'],
                      '.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S8.bam',
                              base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8.bam'],
                      'out_dir':[base+'node_data/'+host+'/']})
    print(outs)    
    
    #lumpy    
    st = stage.Stage('lumpy')
    outs = st.run(1, {'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S8.bam',
                              base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8.bam'],
                      'out_dir':[base+'node_data/'+host+'/']})
    print(outs)
    
    #hydra multi
    st = stage.Stage('hydra')
    outs = st.run(1, {'.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S8.bam',
                              base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8.bam'],
                      'out_dir':[base+'node_data/'+host+'/']})
    print(outs)
    
    #cnvnator VC
    st = stage.Stage('cnvnator')
    cnvnator_params = st.get_params()
    cnvnator_params['window']['value'] = 100
    st.set_params(cnvnator_params)
    outs = st.run(1, {'.fa':[base+'node_data/'+ref_base+'.fa'],
                      '.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S8.bam',
                              base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8.bam']})
    print(outs)

    #cnmops VC
    st = stage.Stage('cnmops')
    cnmops_params = st.get_params()
    cnmops_params['window']['value'] = 100
    cnmops_params['mode']['value']   = 0      #CASE/CTRL WGS analysis
    cnmops_params['normal']['value'] = 3      #poisson normalization
    cnmops_params['cores']['value']  = 2      #bibc threads
    st.set_params(cnmops_params)
    outs = st.run(1, {'.fa':[base+'node_data/'+ref_base+'.fa'],
                      '.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S8.bam',
                              base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8.bam']})    
    print(outs)
    
    #genomestrip
    st = stage.Stage('genome_strip')
    outs = st.run(1, {'.fa':[base+'node_data/'+ref_base+'_S18'+'.fa'],
                      '.fa.svmask.fasta':[base+'node_data/'+ref_base+'_S18'+'.fa.svmask.fasta'],
                      '.ploidymap.txt':  [base+'node_data/'+ref_base+'_S18'+'.ploidymap.txt'],
                      '.rdmask.bed': [base+'node_data/'+ref_base+'_S18'+'.rdmask.bed'],
                      '.bam':[base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S1_S4_S8.bam',
                              base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S1_S4_S8.bam'],
                      'out_dir':[base+'node_data/'+host+'/']})
    print(outs)
    """
    
    #variation hunter
#    st = stage.Stage('variationhunter')
#    outs = st.run(1,{'.DIVET.vh': [base+'node_data/'+host+'/'+ref_base+'_R1_CTRL_S26.DIVET.vh',
#                                   base+'node_data/'+host+'/'+ref_base+'_R1_CASE_S26.DIVET.vh']})
    #pindel
    
    #svseq2
    
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