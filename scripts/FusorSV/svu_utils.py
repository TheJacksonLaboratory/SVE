#interfacing with SVU class and forming core data structures
import glob
import os
import csv
import copy
import re
import datetime
#pip installed libs
import numpy as np
import HTSeq as ht  # :::TO DO::: refactor out this one
#local libs
import fusion_utils as fu
import read_utils as ru
from structural_variant_unit import SVU
import crossmap as cs #do conditional or function import
import mygene         #do conditional or functional import

def get_sv_types():
    return SVU().get_sv_types()

def pretty_ranges(B,units):
    s,size_map = [],{0:'',1:'K',2:'M',3:'G',4:'T',5:'P',6:'E'}
    for b in B:
        x,m = [b,0],0
        if x[0]/1000>0:
            while x[0]/1000>0:
                x[1] = x[0]%1000 
                x[0],m = x[0]/1000,m+1
            d = ''
            if x[1]>0: d = '.'+str(x[1])[:1]
            s += [str(x[0])+d+size_map[m]+units]
        else:
            s += [str(b)+units]
    t = []
    for i in range(len(s)-1):
        t += [s[i]+'-'+s[i+1]]
    return t
    
#first beta vcf attempt, uses the sorting order of chroms
#CHROM, POS, ID, REF, ALT, QUAL, FILTER INFO
#QUAL is Phred-scale quality score for the assertion made in ALT 10log(prob(call in ALT is wrong))
def svul_to_genome(S,O):
    L = {O[k]:k for k in O} #offsets as keys
    #ctg = max([len(L[k]) for k in L])
    C,B = [],sorted(L.keys())
    for i in range(len(S)):
        #find each chrom bin and subtract the offset
        x1,x2,t,ys,wx,end = S[i][0],S[i][1],S[i][2],S[i][3],S[i][4],S[i][6:] #can have flt on index 7
        xo,chrx = 0,''
        for j in range(0,len(B)-1):
                if B[j] <= x1 < B[j+1]: xo,chrx = B[j],L[B[j]]
        for y in ys:
            wy,yo,chry = S[i][5],0,''
            for j in range(0,len(B)-1):
                if B[j] <= y[0] < B[j+1]: yo,chry = B[j],L[B[j]]
            #expand each svu into a source destination with chrom tags
            C += [[chrx,int(x1)-xo,int(x2)-xo,wx,chry,int(y[0])-yo,int(y[1])-yo,wy,t]+end]
    #C = sorted(C,key=lambda x: (x[0].rjust(ctg),x[1])) #SORTING ISSUES::::::::::::::::::::::::::::::
    return C

#G = [chrx, posx1, posx2, wx, chry, posy1, posy2, wy, type, {idx}]
def svult_to_genome(S,O):
    L = {O[k]:k for k in O} #offsets as keys
    ctg = max([len(L[k]) for k in L])
    G = []
    for t in S:
        G += svul_to_genome(S[t],O)
    G = sorted(G,key=lambda x: (x[0].rjust(ctg),x[1])) #SORTING ISSUES HERE::::::::::::::::::::::::::::::
    return G

#G = [chrx, posx1, posx2, wx, chry, posy1, posy2, wy, type, {idx}]
def svult_to_bed(S,O,bed_path,cid,rgb=[255,255,255]):
    data,i = [],1
    #fields = ['chrom','chromStart','chromEnd','name','score','strand','itemRgb'] #8bit RGB?
    G = svult_to_genome(S,O)
    for v in G:
        data += [[v[0],v[1],v[2],str(i)+'-'+cid,1.0,'.',','.join([str(c) for c in rgb])]]
        i += 1
    with open(bed_path, 'wb') as bed:
        csv_w = csv.writer(bed,delimiter='\t')
        #csv_w.writerow(fields)
        for row in data:
            csv_w.writerow(row)     

#bed file generator for the base call set C
def s2bed(S,O,bed_base_path):
    for s in S:
        types = list(set(list(S[s][:,2])))
        for t in types:
            i = np.where(S[s][:,2]==t)[0]
            cid = str(t)+'C'+str(s)
            svult_to_bed(S[s][i],O,bed_base_path+cid+'.bed',cid)

#:::TO DO::: take in all calls and use the alpha to add the PASS or lowqual status
#G = [chrx,x1,x2,wx,chry,y1,y2,wy,t,{idx}]
#what is needed to generate the VCF header and how can it be gathered
#before this steps happens, IE how to collect the metadata
#ref should be a dict with {'ref_name':{'chr1':len(chr1)}} taken from the info file
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO (SAMPLE)
def genome_to_vcf(D,ref_seq,types,chroms,callers,out_path,sname,
                  target_key=None,header_path=None,fltr=1): #:::TO DO::: flt needs work in the mergeing step
    refname = ref_seq.keys()[0] #get the ref name
    C = {k:len(ref_seq[refname][k]) for k in ref_seq[refname]} #get the chrom names
    ctg = max([len(k) for k in C])
    cs = sorted(C.keys(), key =lambda x: x.rjust(ctg)) #SORTING ISSUES::::::::::::::::::::::::::::::
    if header_path is None: #default vcf_header is in the data directory of fusionSVU
        path = os.path.dirname(os.path.abspath(__file__))+'/data/header_template.vcf'
    else:
        path = header_path
    raw = []
    with open(path,'r') as f: raw = f.readlines() #load the header template
    #assert a vcf header fiel was loaded....needs headers. etc could just use the vcf validator tool!
    header,data = [],[] #the vcf header and the data atbel portion of the final file
    hist = {t:{0:[],1:[]} for t in types}
    if len(raw)>10 and raw[0].startswith('##fileformat=VCFv'): #header has at least 10 lines
        #clear the header template comments done in //
        for row in raw:
            if not row.startswith('//'): #clear newlines
                if row.startswith('##'): header += [row.replace('\n','')]
        for i in range(len(header)):
            if header[i].startswith('##fileDate='):
                header[i] += str(datetime.date.today()).replace('-','')
            if header[i].startswith('##reference='):
                header[i] += refname #name
        #construct a name length pair for each contig,,, '##contig=<ID=,len=,>'
        for k in cs:
            header += [''.join(['##contig=<ID=',k,',len=',str(C[k]),'>'])]
        s = ','.join([str(k)+':'+callers[k] for k in sorted(callers.keys())])
        header += ['##INFO=<ID=SVMETHOD,Number=.,Type=String,Description="'+s+'">']
        header += ['\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sname])] #GT 1/1
        #header complete now parse and generate the data rows before writting data
        #G = [chrx, posx1, posx2, wx, chry, posy1, posy2, wy, type, {idx}]
        #:::TO DO::: take in all calls and use the alpha to add the PASS or lowqual status
        for i in range(len(D)):
            if D[i][0] in chroms:
                idx = D[i][9] #{11:set([2])} //can we dig out the origunal ID field to attach here?
                for k in idx: idx[k] = list(idx[k])
                
                target = 1
                if target_key in idx:
                    target = 0
                    idx.pop(target_key)
                
                hist[D[i][8]][target] += [D[i][3]]

                svmethod = str(D[i][9]).replace(' ','')
                svmethod = svmethod.replace('{','')
                svmethod = svmethod.replace('}','')
                svmethod = svmethod.replace('[','')
                svmethod = svmethod.replace(']','')
                row = [D[i][0],                                #CHROM
                       str(D[i][1]),                           #POS
                       'fusorSV_'+str(i+1),                    #ID
                       ref_seq[refname][D[i][0]].seq[D[i][1]], #REF
                       '<'+types[D[i][8]]+'>',                 #ALT
                       '.',                                    #QUAL need to calculate
                       'PASS',                                 #devise a pre and post filter
                       'SVTYPE='+types[D[i][8]]+';SVLEN='+str(D[i][2]-D[i][1])+';'+\
                       'END='+str(D[i][2])+';CHR2='+str(D[i][4])+';IMPRECISE;SVEX='+str(D[i][3])+\
                       ';SVMETHOD='+svmethod+';TARGET='+str(target),
                       'GT','0/1']
                data += [row]
        #some string conversion now
        vcf = '\n'.join(header)+'\n'
        for i in data:
            vcf += '\t'.join(i)+'\n'
        with open(out_path,'w') as f: #looks good can just write now and run with the vcf validator tool
            f.write(vcf)
        return hist
        
#from tigra readme DOCS, tab-seperated
#CHR     
#START_OUTER     
#START_INNER     
#END_INNER       
#END_OUTER       
#TYPE_OF_EVENT   
#SIZE_PREDICTION 
#MAPPING_ALGORITHM       
#SEQUENCING_TECHNOLOGY   
#SAMPLEs 
#TYPE_OF_COMPUTATIONAL_APPROACH  
#GROUP   
#OPTIONAL_ID
#1       829757  829757  829865  829865  DEL       116     MAQ     SLX     NA19238,NA19240    RP      WashU
def genome_to_g1k(D,types,chroms,sname,out_path,target_key=None,ex_flt=0.0,
                  map_alg='BWA',sq_tech='SLX',comp='FUSION',grp='JAXGM'):
    data,g1k = [],''
    if target_key is None: #this is a default for when you want to check all the calls
        for i in range(len(D)):
            if D[i][0] in chroms and D[i][8] in types and D[i][7]>=ex_flt:
                data += [[D[i][0],str(D[i][1]),str(D[i][1]),str(D[i][2]),str(D[i][2]),
                          types[D[i][8]],str(D[i][2]-D[i][1]),
                          map_alg,sq_tech,sname,comp,grp,'fusorSV_'+str(i+1)]]
    else: #this means that the call did not overlap with a target call
        for i in range(len(D)):
            idx = D[i][9] #{0:[-1]} flag in idx field for overlap with target
            for k in idx: idx[k] = list(idx[k])
            if D[i][0] in chroms and D[i][8] in types and D[i][7]>=ex_flt and not target_key in idx:
                data += [[D[i][0],str(D[i][1]),str(D[i][1]),str(D[i][2]),str(D[i][2]),
                          types[D[i][8]],str(D[i][2]-D[i][1]),
                          map_alg,sq_tech,sname,comp,grp,'fusorSV_'+str(i+1)]]
    for i in data:
        g1k += '\t'.join(i)+'\n'
    with open(out_path,'w') as f:
        f.write(g1k)
    return data

#t:[[x1,x2,t,[y],wx,wy,{idx}]] => t:[[chrx,x1,x2,t,[[chry,y1,y2]],wx,wy,{idx}]]
def svult_to_glt(S,O):
    L = {O[k]:k for k in O}    #offsets as keys
    B,G = sorted(L.keys()),{}  #sort the offsets
    for t in S:
        C = []             
        for i in range(len(S[t])):
            #find each chrom bin and subtract offset
            x1,x2,y   = S[t][i][0],S[t][i][1],S[t][i][3]
            wx,wy,idx = S[t][i][4],S[t][i][5],S[t][i][6]
            chrx,ys   = '',[[] for k in y]
            for j in range(0,len(B)-1):
                if B[j] <= x1 < B[j+1]:
                    chrx,x1,x2 = L[B[j]],int(x1)-B[j],int(x2)-B[j]
                for k in range(len(y)):
                    if B[j] <= y[k][0] < B[j+1]:
                        ys[k]  = [L[B[j]],int(y[k][0])-B[j],int(y[k][1])-B[j]]
            C += [[chrx,x1,x2,t,ys,wx,wy,idx]]
        G[t] = C
    return G

#svua is [x1,x2,t,y1,y2,wx,wy]
#svul is [x1,x2,t,y1,y2,wx,wy,idx={}]
def svua_to_svul(S):
    C = []
    for i in range(len(S)):
        C += [list(S[i])+[{}]]
    return C
    
def svul_to_svua(C):
    if len(C) > 0:
        S = np.zeros((len(C),4),dtype='u4')
        for i in range(len(C)):
            S[i] = C[i]
            S[i][3] = 0 #hard coded y ------------------------------------------------
    else:
        S = np.array([],dtype='u4')
    return S
#mark for updates----------------------------------------------

#given: Q[type][f_id][sname][i] return tigra indecies to lookup ctgs
def tigra_ids(Q,sname,idx=6,f_id=-1,t_id=38):
    I = {}
    for t in Q:
        if Q[t].has_key(f_id) and Q[t][f_id].has_key(sname) and len(Q[t][f_id][sname])>0:
            for i in range(len(Q[t][f_id][sname])):
                if Q[t][f_id][sname][i][idx].has_key(t_id):
                    I[i] = {k:{} for k in Q[t][f_id][sname][i][idx][t_id]}
    return I

def info_to_idx(info):
    svmethod = info.split('SVMETHOD=')[-1].split(';')[0]
    raw = svmethod.split(',')
    idx,prev = {},0    
    for i in range(len(raw)):
        if raw[i].find(':')>0:
            split = raw[i].split(':')
            prev = int(split[0])
            idx[prev] = set([int(split[1])])
        else:
            idx[prev].add(int(raw[i]))
    return idx

def idx_to_str(idx):
    s = ''
    for k in sorted(idx):
        s += str(k)+':'+','.join(sorted(list(idx[k])))+'|'
    return s[:-1]

def get_info_v(info,k,d=';'):
    #match the key in start position, or ; delimited or with a whitespace in front
    p = '\A'+k+'=|['+d+']'+k+'=|[\s]'+k+'='      
    m = re.search(p,info)
    if m is None: v = ''
    else:         v = info[m.end():].split(d)[0]
    return v

def info_to_svtype(info):
    return get_info_v(info,'SVTYPE')

def info_to_end(info):
    return int(get_info_v(info,'END'))

def info_to_len(info):
    return int(info.split('SVLEN=')[-1].split(';')[0])

def info_to_consensus(info):
    return get_info_v(info,'CONSENSUS')

def info_to_svex(info):
    return float(info.split('SVEX=')[-1].split(';')[0])

def info_to_target(info):
    return [int(i.split('_')[0]) for i in info.split('TARGET=')[-1].split(';')[0].split(',')]

def info_to_svmethod(info):
    return info.split('SVMETHOD=')[-1].split(';')[0]

def g1kP3_info_to_fusorSV_info(info,k,i):
    return info+';SVEX=1.0;SVMETHOD=%s:%s;TARGET=1'%(k,i)

def lift_tuple_same_strand(lift_tuple):
    return True

def g1kP3_liftover(vcf_path,ref_path,chain_path,
                     CHR=0,POS=1,ID=2,REF=3,ALT=4,QUAL=5,FILT=6,INFO=7,FORMAT=8,SAMPLE=9,
                     add_chr=True):
    return True

#given a multi sample delly file, query
def delly_vcf_reader(vcf_glob,out_vcf,reference,samples=['HG00513','HG00733','NA19240'],
                     VCF_CHR=0,VCF_POS=1,VCF_ID=2,VCF_REF=3,VCF_ALT=4,
                     VCF_QUAL=5,VCF_FILT=6,VCF_INFO=7,VCF_FORMAT=8,VCF_SAMPLE=9):
    header,data,snames = [],[],[]
    for vcf in glob.glob(vcf_glob):
        with open(vcf,'r') as f:
            sname,i,k = vcf.rsplit('/')[-1].split('.vcf')[0],0,0
            for line in f:
                if not line.startswith('#'):
                    r = line.split('\n')[0].split('\t')
                    r[VCF_ID] += '_'+sname
                    r[VCF_POS] = int(r[VCF_POS])
                    ln = info_to_end(r[VCF_INFO])-r[VCF_POS] #
                    r[VCF_INFO] = g1kP3_info_to_fusorSV_info(r[VCF_INFO],k,i)
                    data += [[r[VCF_CHR],r[VCF_POS],info_to_end(r[VCF_INFO]),ln]+r[VCF_ID:]]
                    i += 1
    vcfs = glob.glob(vcf_glob) #grab first header only
    with open(vcfs[0],'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.split('\n')[0].split('\t')]
    snames = header[-1][VCF_FORMAT+1:]                      
    #now sort and cluster calls for genotyping
    data = coordinate_sort(data)
    return True

#
def slice_merged(flat_data):
    mag = 0
    B = [(1,100),(100,250),(250,500),(500,1000),(1000,2500),
         (2500,5000),(5000,10000),(10000,50000),(50000,100000),
         (100000,1000000),(1000000,100000000)]
    sliced_data = {'INS':{b:0 for b in B},'DEL':{b:0 for b in B},
                   'DUP':{b:0 for b in B},'INV':{b:0 for b in B}}
    for row in flat_data:
        t = row[4].replace('<','').replace('>','')
        l = info_to_end(row[7])-int(row[1])
        mag += l
        for b in B:
            if l>b[0] and l<=b[1]:
                sliced_data[t][b] += 1
    return sliced_data
        
    
def g1kP3_vcf_multi_sample_merge(vcf_glob,out_vcf,reference,overlap=0.5,
                                 VCF_CHR=0,VCF_POS=1,VCF_ID=2,VCF_REF=3,VCF_ALT=4,
                                 VCF_QUAL=5,VCF_FILT=6,VCF_INFO=7,VCF_FORMAT=8,VCF_SAMPLE=9):
    header,data,snames = [],[],[]
    for vcf in glob.glob(vcf_glob):
        with open(vcf,'r') as f:
            sname,i,k = vcf.rsplit('/')[-1].split('_S0.vcf')[0],0,0
            snames += [sname]
            for line in f:
                if not line.startswith('#'):
                    r = line.split('\n')[0].split('\t')
                    r[VCF_ID] += '_'+sname
                    r[VCF_POS] = int(r[VCF_POS])
                    ln = info_to_end(r[VCF_INFO])-r[VCF_POS] #
                    r[VCF_INFO] = g1kP3_info_to_fusorSV_info(r[VCF_INFO],k,i)
                    data += [[r[VCF_CHR],r[VCF_POS],info_to_end(r[VCF_INFO]),ln]+r[VCF_ID:]]
                    i += 1
    vcfs = glob.glob(vcf_glob) #grab first header only
    with open(vcfs[0],'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.split('\n')[0].split('\t')]
    header[-1]  = '\t'.join(header[-1][0:9]+['FORMAT']+snames) #fix up the sample colums                        
    raw = []
    for i in range(0,len(data),2):
        raw += [data[i]]
    data = raw
    #now sort and cluster calls for genotyping
    data = coordinate_sort(data)
    cluster = coordinate_cluster(data,overlap)
    flat_data = clusters_to_flattened_str(cluster,snames,reference)
    vcf = '\n'.join([''.join(h) for h in header])+'\n'
    for i in flat_data:
        vcf += '\t'.join(i)+'\n'
    with open(out_vcf,'w') as f: #looks good can just write now and run with the vcf validator tool
        f.write(vcf)
        return True
    return False                                   

#downstream analysis check of lifted coordinates
def fusorSV_liftover(vcf_path,ref_path,chain_path,
                     CHR=0,POS=1,ID=2,REF=3,ALT=4,QUAL=5,FILT=6,INFO=7,FORMAT=8,SAMPLE=9,
                     add_chr=True):
    if add_chr: chrom = 'chr'
    else:       chrom = ''        
    I,header,data,err = {},[],[],[]
    with open(vcf_path,'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.split('\n')[0].split('\t')]
            else:
                data += [line.split('\n')[0].split('\t')]
    #now do the coordinate to coordinate liftover search
    mapTree,targetChromSizes, sourceChromSizes = cs.read_chain_file(chain_path)
    for row in data:
        try:
            q_chr = chrom + row[CHR]
            q_start = int(row[POS])
            q_stop  = info_to_end(row[INFO])
            fusorSV_id = int(row[ID].replace('fusorSV_',''))
            I[fusorSV_id] = cs.map_coordinates(mapTree,q_chr,q_start,q_stop)
        except Exception:
            err += [row]
    return True

#read the fasta files and parse out the correct locations to liftover
def breakseq2_fasta_parser(fasta_path,CHROM=0,POS=1,SET=2,TYPE=3,ID=4):
    raw = []
    with open(fasta_path,'r') as f: raw += f.readlines()
    raw = ''.join(raw).split('\n')
    while raw[-1]=='': raw.pop() #trim any empty ends here
    data = []
    for i in range(0,len(raw),2):
        row = raw[i].split(':')
        chrom,pos,st,typ,ids = row[0:ID+1]
        start,stop = [int(j) for j in pos.split('-')]
        chrom = chrom.replace('>','')
        end = []
        if len(row)>ID+1: end += row[ID+1:]
        data += [[chrom,start,stop,st,typ,ids]+end+[raw[i+1]]]
    return data

#lift coordinates for a breakseq2 library
def breakseq2_liftover(brkpt_path, source_ref_path, destination_ref_path, chain_path, add_chr=True):
    source_ref_path = '/Users/tbecker/Desktop/TEMP/human_g1k_v37_decoy/human_g1k_v37_decoy.fa'
    source_ref = ru.read_fasta(source_ref_path,True)
    destination_ref_path = '/Users/tbecker/Desktop/TEMP/human_g1k_v38_decoy_hla/human_g1k_v38_decoy_hla.fa'
    destination_ref = ru.read_fasta(destination_ref_path,True)
    fasta_path = '/Users/tbecker/Desktop/TEMP/human_g1k_v37_decoy/human_g1k_v37_decoy_S35.brkptlib.fna'
    chain_path = '/Users/tbecker/Desktop/TEMP/fusionSVU/data/liftover/hg19ToHg38.over.chain.gz'
    return True

def add_chrom(data,CHR=0):
    for row in data:
        if row[CHR] in [str(x) for x in range(23)]+['Y','X','MT']:
            row[CHR] = 'chr'+row[CHR]
    return data
    
#cluster a set of fusorSV per sample VCF files
def fusorSV_vcf_multi_sample_merge(vcf_glob,out_vcf,reference,overlap=0.5,add_chr=True,
                                   VCF_CHR=0,VCF_POS=1,VCF_ID=2,VCF_REF=3,VCF_ALT=4,
                                   VCF_QUAL=5,VCF_FILT=6,VCF_INFO=7,VCF_FORMAT=8,VCF_SAMPLE=9):
    header,data,snames = [],[],[]
    for vcf in glob.glob(vcf_glob):
        with open(vcf,'r') as f:
            sname = vcf.rsplit('/')[-1].split('_S-1.vcf')[0]
            snames += [sname]
            for line in f:
                if not line.startswith('#'):
                    r = line.split('\n')[0].split('\t')
                    r[VCF_ID] += '_'+sname
                    r[VCF_POS] = int(r[VCF_POS])
                    data += [[r[VCF_CHR],r[VCF_POS],info_to_end(r[VCF_INFO]),info_to_len(r[VCF_INFO])]+\
                              r[VCF_ID:VCF_SAMPLE+1]]
    vcfs = glob.glob(vcf_glob) #grab first header only
    with open(vcfs[0],'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.split('\n')[0].split('\t')]
    header[-1]  = '\t'.join(header[-1][0:9]+snames) #fix up the sample colums                        

    #cluster calls for genotyping
    cluster = coordinate_cluster(data,overlap)
    flat_data = clusters_to_flattened_str(cluster,snames,reference)
    if add_chr: flat_data = add_chrom(flat_data)
    #remove_C = ['C_12596', 'C_20039', 'C_17709', 'C_5949', 'C_6014', 'C_21055', 'C_17312', 'C_9188', 'C_2567', 'C_18819', 'C_18814']
    #remove_S = list(set(snames).difference(set(['NA19017','NA12878','HG00419','NA19238','NA19239','NA19625','NA18525']))) 
    #flat_data = filter_fusorSV_vcf_multi_sample_merge(flat_data,remove_C,remove_S)
    #write the header and flat data to a new vcf file
    vcf = '\n'.join([''.join(h) for h in header])+'\n'
    for i in flat_data:
        vcf += '\t'.join(i)+'\n'
    with open(out_vcf,'w') as f: #looks good can just write now and run with the vcf validator tool
        f.write(vcf)
        return True
    return False

def filter_fusorSV_vcf_multi_sample_merge(flat_data,remove_C,VCF_ID=2):
    filtered_data = []
    for row in flat_data:
        if row[VCF_ID] not in remove_C:
            filtered_data += [row]
    return filtered_data

#read in individual fusorSV sample VCF files, multi_sample VCF sample/id/validation table
def fusorSV_multi_sample_merge_query(fusorSV_vcf_dir,sample_id_validation):
    samples = glob.glob(fusorSV_vcf_dir+'*_S-1.vcf') #get fusorSV VCF files out
    snames  = [s.rsplit('/')[-1].split('_S-1.vcf')[0] for s in samples]
    merged_samples  = fusorSV_vcf_dir+'/all_samples_genotypes_liftover.mapped.vcf'
    raw,header,merged_data = [],[],[]
    with open(merged_samples,'r') as f:
        raw = f.readlines()
    for line in raw:
        if line.startswith('#'):
            header += [line]
            header[-1] = header[-1].replace('\n','')
        else:                    
            merged_data += [line.split('\t')]
            merged_data[-1][-1] = merged_data[-1][-1].replace('\n','')
    header_key,i = header[-1].split('\t'),0
    header_key[0] = header_key[0].replace('#','')
    while i < len(header_key) and header_key[i]!='FORMAT':i += 1
    if i < len(header_key): i += 1
    n = len(header_key[i:]) #number of samples to look at
    contig_len = max([len(i[0]) for i in merged_data])
    merged_data = sorted(merged_data,key=lambda x: (x[0].zfill(contig_len),int(x[1])))
    return merged_data

#need the bedtools intersect output for refseq
def fusorSV_bed_gene_converter(refseq_bed,gene_counts,REFSEQ_ID=3):
    data,genes = [],{}
    with open(refseq_bed,'r') as f:
        for line in f:
            data += [line.split('\t')]
    mg = mygene.MyGeneInfo()
    for row in data:
        q = str(mg.query(row[REFSEQ_ID])['hits'][0]['symbol'])
        if not genes.has_key(q): genes[q]  = 1
        else:                    genes[q] += 1
    s = 'total genes=%s\taverage exons per gene=%s\n'%(len(genes),sum([genes[g] for g in genes])*1.0/len(genes))
    for g in sorted(genes):
        s += '%s,'%g
    with open(gene_counts,'w') as f:
        f.write(s[:-2])
        return True
    return False
    
def fusorSV_multi_sample_merge_query_write(query_data,header,out_vcf):
    s = '\n'.join(header)+'\n'
    s += '\n'.join(['\t'.join(i) for i in query_data])+'\n'
    with open(out_vcf,'w') as f:
        f.write(s)
        return True
    return False
    
#using string INFO field: IE 'DEL', 'DUP', ect
def query_svtype(data,svtype,INFO=7):
    result = []
    for row in data:
        if info_to_svtype(row[INFO])==svtype:
            result += [row]
    return result
    
#get the frequency out of genotyped rows
def query_frequency(data,c,freq,GT=9):
    result = []
    for row in data:
        if c == '<' and sum([1 if i=='0/1' else 0 for i in row]) < freq:
            result += [row]
        elif c == '<=' and sum([1 if i=='0/1' else 0 for i in row]) <= freq:
            result += [row]
        elif c == '>' and sum([1 if i=='0/1' else 0 for i in row]) > freq:
            result += [row]
        elif c == '>=' and sum([1 if i=='0/1' else 0 for i in row]) >= freq:
            result += [row]
        elif c == '==' and sum([1 if i=='0/1' else 0 for i in row]) == freq:
            result += [row]
        elif c == '!=' and sum([1 if i=='0/1' else 0 for i in row]) != freq:
            result += [row]
    return result

def query_svex(data,c,svex,INFO=7):
    result = []
    for row in data:
        if c == '<' and info_to_svex(row[INFO]) < svex:
            result += [row]
        elif c == '<=' and info_to_svex(row[INFO]) <= svex:
            result += [row]
        elif c == '>' and info_to_svex(row[INFO]) > svex:
            result += [row]
        elif c == '>=' and info_to_svex(row[INFO]) >= svex:
            result += [row]
        elif c == '==' and info_to_svex(row[INFO]) == svex:
            result += [row]
        elif c == '!=' and info_to_svex(row[INFO]) != svex:
            result += [row]
    return result

#returns all target values that all agree,
#set agree=False and you get out the conflicting calls
def query_target(data,t=0,agree=True,INFO=7):
    result = []
    for row in data:
        if agree and all([t==j for j in info_to_target(row[INFO])]):
            result += [row]
        elif not agree and (any([t!=j for j in info_to_target(row[INFO])]) and \
                            any([t==j for j in info_to_target(row[INFO])])):
                   result += [row]
    return result

#returns where a caller is present for all the samples
def query_caller_presence(data,t_id=38,p=0.5,INFO=7):
    result = []
    for row in data:
        M = info_to_svmethod(row[INFO]).split('|')
        S = svmethod_to_sample(M)
        x = [1 if t_id in S[k] else 0 for k in S]
        if 1.0*sum(x)/len(x)>=p:
            result += [row]
    return result

#returns calls that have a certain number of supporting algorithms
#agree means that all the samples have to have more than c callers
def query_caller_number(data,c='<',x=3,agree=True,INFO=7):
    results = []
    for row in data:
        M = info_to_svmethod(row[INFO]).split('|')
        S = svmethod_to_sample(M)
        if agree:
            if c=='<' and all([len(S[k])<x for k in S]):
                results += [row]
            elif c=='<=' and all([len(S[k])<=x for k in S]):
                results += [row]
            elif c=='>' and all([len(S[k])>x for k in S]):
                results += [row]
            elif c=='>=' and all([len(S[k])>=x for k in S]):
                results += [row]
        else:
            if c=='<' and any([len(S[k])<x for k in S]) and any([len(S[k])>=x for k in S]):
                results += [row]
            elif c=='<=' and any([len(S[k])<=x for k in S]) and any([len(S[k])>x for k in S]):
                results += [row]
            elif c=='>' and any([len(S[k])>x for k in S]) and any([len(S[k])<=x for k in S]):
                results += [row]
            elif c=='>=' and any([len(S[k])>=x for k in S]) and any([len(S[k])<x for k in S]):
                results += [row]
    return results
    
#for a list of samples return the calls that they are genotyped possitive
def query_sample_presence(data,s_ids,INFO=7):
    result = []
    for row in data:
        M = info_to_svmethod(row[INFO]).split('|')
        S = svmethod_to_sample(M)
        for s_id in s_ids:
            if s_id in S.keys():
                result += [row]
                break
    return result

#returns each samples caller group
def svmethod_to_sample(svmethod):
    S = {}
    for i in svmethod:
        c,s = i.split(':')
        s = set([x.split('_')[-1] for x in s.split(',')])
        for j in s:
            if S.has_key(j): S[j] += [int(c)]
            else:            S[j]  = [int(c)]
    for j in S:
        S[j] = tuple(sorted(S[j]))
    return S
    
#given a a list of FusorSV rows, do a downstream analysis
#on the files where all group frequencies are returned for each sample
def get_svmethod_gfreq(data,INFO=7):
    G = {}
    for row in data:
        M = info_to_svmethod(row[INFO]).split('|')
        S = svmethod_to_sample(M)
        for i in [S[k] for k in S]:
            for j in i:
                if G.has_key(j): G[j] += 1
                else:            G[j]  = 1
    return G       
    
def get_row_freq(row,INFO=7):
    svmethods,S = info_to_svmethod(row[INFO]).split('|'),{}
    C = {int(i.split(':')[0]):{j.split('_')[-1]:int(j.split('_')[0]) for j in i.split(':')[-1].split(',')} for i in svmethods}
    for c in C:
        for s in C[c]:
            if not S.has_key(s):
                S[s]  = {c:[C[c][s]]}
            else: 
                if not S[s].has_key(c):
                    S[s][c]  = [C[c][s]]
                else:
                    S[s][c] += [C[c][s]]
    return S



#group frequency of a given set of data
#samples = {sname:{f_id:0},...}
def get_group_validation_frequency(data,samples,INFO=7):
    F = {}
    for row in data:
        S = get_row_freq(row,INFO)
        for s in S:
            if s in samples:
                g = tuple(sorted(S[s].keys()))
                if not F.has_key(g):
                    F[g] = {}
    return F
        
 
#custom search tool
#Cluster Type SVEX TARGET SVMETHOD fusroSV_id
def three_list_fusorSV_id_search(id_path,repaired_vcf,flat_data):
    #new calls we are looking for
    with open(id_path,'r') as f:
        raw = f.readlines()
    data1 = [r.strip('\n').split('\t') for r in raw]
    #old ones that were repaired
    header,data2,err = [],[],[]
    with open(repaired_vcf,'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.split('\n')[0].split('\t')]
            else:
                data2 += [line.split('\n')[0].split('\t')]
    filtered,removed = [],[]
    for e in data1:
        if e[5] not in [f[2] for f in data2]:
            filtered += [e]
        else:
            removed += [e]
            
    ids,N = [e[2] for e in data2],{}
    for i in ids:
        for row in flat_data:
            fids = row[7].split('FUSORSVID=')[-1].split(';')[0].split(',')
            if i in fids and len(fids) == 1:
                N[i] = row[2]
    1.0*sum([N[k] for k in N])/len(N) #average allele frequency
    return []

#sort rows by chrom, type, start, length
def coordinate_sort_type(data,MAX=100,CHR=0,TYPE=6,R0=1,LEN=3):
    return sorted(data,key=lambda x: (x[TYPE],x[CHR].zfill(MAX),x[R0],x[LEN]))            

def coordinate_sort_pos(data,MAX=100,CHR=0,TYPE=6,R0=1,LEN=3):
    return sorted(data,key=lambda x: (x[CHR].zfill(MAX),x[R0],x[LEN])) 
    
def coordinate_overlap(a,b,r0=1,r1=2,chr1=0,svtype=6): #add s1, s2 later
    if a[svtype]!=b[svtype] or a[chr1]!=b[chr1]:
        x = 0.0
    else:
        l = (a[r1]-a[r0]+1)+(b[r1]+b[r0]+1)
        u = float(min(l,max(a[r1],b[r1])-min(a[r0],b[r0])+1))
        i = 1.0*abs(a[r0]-b[r0])+abs(a[r1]-b[r1])
        x = max(0.0,u-i)/u
    return x

def coordinate_partition_type(data,TYPE=6):
    types = set([x[TYPE] for x in data])
    T = {t:[] for t in types}
    for row in data:
        T[row[TYPE]] += [row]
    return T 
    
#fast LR overlap clustering
#for each contig and for each svtype
#for every two rows in a cluster, they have at least >= overlap amount
def coordinate_cluster(data,overlap=0.5,MAX=100):
    d = coordinate_partition_type(coordinate_sort_pos(data,MAX=MAX))
    C = {t:{} for t in d}
    for t in d:
        i,j,n = 0,0,len(d[t])
        while i < n-1:
            j = i+1
            while coordinate_overlap(d[t][i],d[t][j])>=overlap and j < n-1:
                j += 1
            C[t][i] = d[t][i:j]
            i = j
        if not C[t].has_key(j) and j < n-1: C[t][j] = [d[j]]
    F,i = {},0
    for t in C:
        for j in sorted(C[t]):
            F[i] = C[t][j]
            i += 1
    return F
    
#recursive overlap clustering, following 1000 genomes phase 3 supplements
def coordinate_cluster_g1k_p3_style(data,overlap=0.5):
    C = {}
    #will use a helper function to pass in C and overlap only
    return C
    
def clusters_to_flattened_str(cluster,snames,reference,average=True,
                              CHR=0,POS=1,END=2,ID=4,REF=5,ALT=6,QUAL=7,
                              FILT=8,INFO=9,FORMAT=10,SAMPLE=11,
                              VCF_CHR=0,VCF_POS=1,VCF_ID=2,VCF_REF=3,VCF_ALT=4,
                              VCF_QUAL=5,VCF_FILT=6,VCF_INFO=7,VCF_FORMAT=8,VCF_SAMPLE=9):
    data,seq_name_len = [],max([len(k) for k in reference])
    for k in cluster:
        row = [cluster[k][0][CHR],cluster[k][0][POS],'C_'+str(k),cluster[k][0][REF],
               cluster[k][0][ALT],'.','PASS',cluster[k][0][INFO],'GT','0/1']
        if average:
            pos = int(round(1.0*sum([e[POS] for e in cluster[k]])/len(cluster[k]),0)) #average pos
            end = int(round(1.0*sum([e[END] for e in cluster[k]])/len(cluster[k]),0)) #avergae end
            svlen = end-pos+1 #average svlen
            svex = str(sum([info_to_svex(e[INFO]) for e in cluster[k]])/len(cluster[k])) #average svex
        idx,targets = {},[] #merged idx now has caller_id:set([row+samplename]) mapping
        f_ids = [cluster_fusorSV_id_to_sample(e[ID]) for e in cluster[k]]
        fusorSV_ids = [e[ID] for e in cluster[k]]
        for e in cluster[k]:
            new = info_to_idx(e[INFO])
            sname = cluster_fusorSV_id_to_sample(e[ID])
            targets += [str(info_to_target(e[INFO])[0])+'_'+sname] #now target flag and sample
            for i in new:
                if not idx.has_key(i):
                    idx[i] = set([str(j)+'_'+sname for j in new[i]])
                else:
                    idx[i] = idx[i].union(set([str(j)+'_'+sname for j in new[i]]))
        svmethod = idx_to_str(idx)
        row[VCF_POS] = str(pos)
        row[VCF_INFO] = cluster_info_update(row[VCF_INFO],svlen,end,svex,svmethod,targets,fusorSV_ids)
        row[VCF_REF]  = reference[row[VCF_CHR]].seq[pos] #new average pos
        row[VCF_SAMPLE] = cluster_to_samples(snames,f_ids)
        data += [row]
    data = sorted(data,key=lambda x: (x[VCF_CHR].zfill(seq_name_len),int(x[VCF_POS])))
    return data

def cluster_to_samples(snames,f_ids):
    samples = ''
    for s in snames:
        if s in f_ids:
            samples += '0/1'+'\t'
        else:
            samples += '0/0'+'\t'
    return samples[:-1]   

#verage up and update the row            
def cluster_info_update(info,svlen,end,svex,svmethod,targets,
                        fusorSV_ids=None,
                        remove_destination_coordinates=True):
    if remove_destination_coordinates:
        info = info.split('CHR2=')[0]+';'.join(info.split('CHR2=')[-1].split(';')[1:])
    l = info.split('SVLEN=')[0]
    r = ';'.join(info.split('SVLEN=')[-1].split(';')[1:])
    info = l+'SVLEN='+str(svlen)+';'+r
    l = info.split('END=')[0]
    r = ';'.join(info.split('END=')[-1].split(';')[1:])
    info = l+'END='+str(end)+';'+r
    l = info.split('SVEX=')[0]
    r = ';'.join(info.split('SVEX=')[-1].split(';')[1:])
    info = l+'SVEX='+str(svex)+';'+r
    l = info.split('SVMETHOD=')[0]
    r = ';'.join(info.split('SVMETHOD=')[-1].split(';')[1:])
    info = l+'SVMETHOD='+str(svmethod)+';'+r
    l = info.split('TARGET=')[0]
    if fusorSV_ids is not None: #fusorSV_ids here
        info = l+'FUSORSVID='+','.join(fusorSV_ids)+';'+'TARGET='+','.join(targets)
    else:
        info = l+'TARGET='+','.join(targets)    
    return info
            

def cluster_fusorSV_id_to_sample(fusorSV_id):
    return fusorSV_id.rsplit('_')[-1]

#downstream analysis check of lifted coordinates
def fusorSV_fix_merged_samples(vcf_in_path,vcf_out_path):       
    header,data = [],[]
    with open(vcf_in_path,'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.split('\n')[0].split('\t')]
            else:
                r = line.split('\n')[0].split('\t')
                r[2] += '_'+r[-1]
                data += [r[0:8]]
    header[-1] = ['\t'.join(header[-1][0:8])] #fix the colums specs
    s = '\n'.join([h[0] for h in header])+'\n'
    s += '\n'.join(['\t'.join(row) for row in data])+'\n'
    with open(vcf_out_path,'w') as f:
        f.write(s)
        return True
    return False

def fusorSV_vcf_liftover(vcf_in_path,ref_path,chain_path):
    mapTree,targetChromSizes, sourceChromSizes = cs.read_chain_file(chain_path)
    cs.crossmap_vcf_file(mapTree,vcf_in_path,chain_path,ref_path)
    return True

def fusorSV_vcf_liftover_samples(sample_dir_glob,ref_path,chain_path):
    vcfs = glob.glob(sample_dir_glob)
    for vcf in sorted(vcfs):
        fusorSV_vcf_liftover(vcf,ref_path,chain_path)

#reads and parses the fusorSV VCF file here to attach supporting ids
def fusorSV_support_ids(vcf_path,ID=2,INFO=7,s_id=[38]):
    I,header,data,err = {s:{} for s in s_id},[],[],[]
    with open(vcf_path,'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.split('\n')[0].split('\t')]
            else:
                data += [line.split('\n')[0].split('\t')]
    for row in data:
        try:
            fusorSV_id = int(row[ID].replace('fusorSV_',''))
            idx = info_to_idx(row[INFO])
            for s in s_id:
                if idx.has_key(s):
                    I[s][fusorSV_id] = {k:{} for k in idx[s]}      
        except Exception:
            err += [row]
            pass
    return I

#given a fusorSV_id mapps in the original VCF line of a supporting VCF row
def support_id_map(M,V,ID=2,s_id=[],callers=None):
    U,K = {s:{} for s in s_id},{s:{} for s in s_id}
    for s in U:
        if V.has_key(s):
            for j in range(len(V[s])):
                k = tuple(V[s][j].svu[0])
                if U[s].has_key(k): U[s][k] += [j]
                else:               U[s][k]  = [j]
            K[s] = {tuple(U[s][k]):list(k) for k in U[s]}    
    for s in M:
        for i in M[s]:
            for j in M[s][i]:
                for k in K[s]:
                    if j in k:
                        row = V[s][j].as_vcf_row()
                        if callers is None:
                            row = row[0:ID]+['S'+str(s)+'_'+str(j)]+row[ID+1:]
                        else:
                            row = row[0:ID]+[callers[s]+'_'+str(j)]+row[ID+1:]
                        for x in range(len(row)):
                            if not type(row[x]) is str:
                                row[x] = str(row[x])
                        M[s][i][j] = [str(x) for x in row]
    return M

#SM already has the s_ids that are attached...
def support_id_search(SM,fusorSV_vcf,ID=2,INFO=7):
    IS,header,data,err = {},[],[],[]
    with open(fusorSV_vcf,'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.split('\n')[0].split('\t')]
            else:
                data += [line.split('\n')[0].split('\t')]
    header[-1] = header[-1][0:-2] #just go up to the INFO feild for now
    for row in data:
        try:
            fusorSV_id = int(row[ID].replace('fusorSV_',''))
            idx = info_to_idx(row[INFO])
            for s in SM:
                if idx.has_key(s):
                    if IS.has_key(fusorSV_id):
                        IS[fusorSV_id] += [SM[s][fusorSV_id][k] for k in idx[s]]
                    else:
                        IS[fusorSV_id]  = [SM[s][fusorSV_id][k] for k in idx[s]] 
        except Exception:
            err += [row]
            pass
    return IS,header

#write the supporting VCF file for downstream analyis
def write_support_id_vcf(IS,header,support_vcf_path):
    S = '\n'.join(['\t'.join(i) for i in header])+'\n'
    for i in sorted(IS):
        S += '\n'.join(['\t'.join(j) for j in IS[i]])+'\n'
    with open(support_vcf_path,'w') as f:
        f.write(S)
        return True
    return False
    
#replace original vcf id field with the fusorSV supporting ids instead
#replace ALT with the ALT contig with the longest tigraSV CTG if one exists    
        
#given V and fusorSV idex into it, look back and merge
#INFO;CTG= fields to retrieve the CTG ids to look into the fasta file in tigra_Y
def tigra_id_to_ctg_map(M,V,t_id=38):
    #i will be the highest index, look back into i-x entries to merge x records
    U = {}
    for j in range(len(V[t_id])):
        k = tuple(V[t_id][j].svu[0])
        if U.has_key(k): U[k] += [j]
        else:            U[k]  = [j]
    K = {tuple(U[k]):list(k) for k in U}
    for s in M:
        for i in M[s]: #list of t_id keys
            for j in M[s][i]:
                for k in K:
                    if j in k: #pull out the tigra CTG from the info field
                        M[s][i][j] = {V[t_id][x].info.split('CTG=')[-1].split(';')[0]:'' for x in k}
    return M
    
#given a tigra ctg map M, a coordinate offset map O and the fasta file ctg_fasta
#dig out the contig sequences in relation to the chrom start/stop position
#and do something with them (DP sliding alignment with affine gap?)
def tigra_ctg_search(M,ctg_fasta,t_id=38):
    seq,raw,err = {},'',0
    with open(ctg_fasta,'r') as f:
        raw = f.readlines()
    for i in range(0,len(raw),2): #tigra fasta is even idex = >ctg, odd idex = sequence
        if raw[i].startswith('>') and i < len(raw)+1:
            ctg = raw[i].upper().replace('>','').split(' ')[0]
            seq[ctg] = raw[i+1].replace('\n','')
    for s in M:
        for i in M[s]:
            for j in M[s][i]:
                for ctg in M[s][i][j]:
                    if seq.has_key(ctg): M[s][i][j][ctg] = seq[ctg]
                    else:                err += 1
    print('%s parsed and searched with %s errors'%(ctg_fasta,err))
    return M

#given the result of the tigra_ids->togra_id_to_ctg_map->tigra_ctg_search
#make one easy to search tsv file with:
#f_id t_id t_ctg fasta_style_seq
def write_tigra_ctg_map(M,tigra_tsv_path,t_id=38):
    S = '\t'.join(['fusorSV_id','tigra_id','ctg_id','fasta_seq'])+'\n'
    for s in M:
        for i in sorted(M[s].keys()):         #fusorSV call id
            for j in sorted(M[s][i].keys()):  #tigra call id attached to fusorSV call id i
                for c in M[s][i][j]: #tigra ctg attached to the fusorSV call id location and type
                    S += '\t'.join([str(i),str(j),c,M[s][i][j][c]])+'\n'
    with open(tigra_tsv_path,'w') as f:
        f.write(S)
        return True
    return False

#lift over analysis

#TO DO this will need to be updated once the final data structure is set
#construct a filtered SVD of FLTR==PASS
def construct_svult(vcr,chroms,offset_map,s_id,flt=0,upper=int(500E3)):
    sx,vc_i,vx,k,j = {},{},[],[],0 #filtered container, and j is the originating row
    for vc in vcr: #iterate on the variant call record
        vx += [SVU(vc,offset_map)]
        if vx[-1].chrom in chroms and vx[-1].filter >= flt and vx[-1].svlen < upper:
            sx[tuple(vx[-1].svu[0])] = j
            #if len(vx[-1].svu)>1:
            #   sx[tuple(vx[-1].svu[1])] = j 
        j += 1
    svua = np.zeros((len(sx),7),dtype='u4')
    k = sorted(sx.keys()) #sorted svua key (x1,x2,t,y1,y2,wx,wy):j for row
    for i in range(len(k)):
        svua[i] = k[i]     #sorted, unique entries
        vc_i[i] = sx[k[i]] #get j back
    svult = {}
    types = sorted(list(set(list(svua[:,2]))))
    for t in types:
        l = np.where(svua[:,2]==t)[0]
        L = []
        for i in l:
            x = list(svua[i])
            L += [x[0:3]+[[[x[3],x[4]]],np.float64(x[5]),np.float64(x[6]),{s_id:{vc_i[i]}}]] #build it out
        svult[t] = L
    return svult,vx

def print_svult(C):
    print('x1\tx2\tt\t[y]\twx\twy\t{idx}')
    for t in C:
        for i in C[t]:
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s'%(str(i[0]),str(i[1]),str(i[2]),
                                                str(i[3]),str(i[4]),str(i[5]),str(i[6])))
            
#given a bash type wildcard glob path, reads the vcf as svult
#do a flt param map for each caller id =>{s_id:flt_val}
def vcf_glob_to_svultd(path_glob,chroms,offset_map,flt=0,flt_exclude=[]):        
    vcfs,S,V = glob.glob(path_glob),{},{}
    for vcf in vcfs:
        vcr = ht.VCF_Reader(vcf)
        s_id = id_trim(vcf)
        if s_id in flt_exclude:
            S[s_id],V[s_id] = construct_svult(vcr,chroms,offset_map,s_id,-1)
        else:
            S[s_id],V[s_id] = construct_svult(vcr,chroms,offset_map,s_id,flt)
    return S,V

#given a vcf with SVCP naming convention, trim to an int value
def id_trim(s):
    i = 0
    try:
        i = int(s.rsplit('/')[-1].rsplit('_S')[-1].split('.')[0])
    except Exception:
        print('not _S*.vcf named...')
    return i 

#ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
                    
#input is a set of call sets S and list of list ref regions R
#ouput is a new set of sets T that does not overlap with any of Rs elements
def filter_call_sets(S,R,exclude=[]):
    T = {}
    types = set([i for l in [S[k].keys() for k in S] for i in l])
    for k in S:
        N = {}
        for t in types:
            if S[k].has_key(t):
                if k in exclude: N[t] = S[k][t]
                else:            N[t] = fu.filter_regions(S[k][t],R)
        T[k] = N
    return T

def filter_call_sets2(S,R,exclude=[]):
    T = {}
    types = set([i for l in [S[k].keys() for k in S] for i in l])
    for k in S:
        N = {}
        for t in types:
            if S[k].has_key(t):
                if k in exclude: N[t] = S[k][t]
                else:            N[t] = filter_regions2(S[k][t],R)
        T[k] = N
    return T

def filter_regions2(C,R):
    D = [-1]
    i,j,upper,n,m = 0,0,0,len(C)+1,len(R)+1 #boundries here
    if   n>1 and m<=1: upper = C[-1][1]
    elif m>1 and n<=1: upper = R[-1][1]
    elif n>1 and m>1:  upper = max(C[-1][1],R[-1][1])
    C += [[upper+2,upper+2,0,[],0,0,{}],[upper+4,upper+4,0,[],0,0,{}]] #pad out the end of C
    R += [[upper+2,upper+2,0,[],0,0,{}],[upper+4,upper+4,0,[],0,0,{}]] #pad out the end of R
    
    while i+j < n+m:  #pivioting dual ordinal indecies scan left to right on C1, C2        
        a = long(C[i][0])-long(R[j][0])
        b = long(C[i][0])-long(R[j][1])
        c = long(C[i][1])-long(R[j][0])
        d = long(C[i][1])-long(R[j][1])
        if  a==0 and d==0:     #[7] C[i] and R[j] are equal on x
            #print('equal to\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]
            i += 1
            j += 1 
        elif  c<0:               #[1] C[i] disjoint of left of R[j]
            #print('disjoint left\ti=%s\tj=%s'%(i,j))
            i += 1    
        elif  b>0:               #[6] C[i] disjoint right of R[j]
            #print('disjoint right\ti=%s\tj=%s'%(i,j))        
            j += 1 
        elif  a<0 and d<0:       #[2] C[i] right overlap to R[j] left no envelopment
            #print('right overlap\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]      
            i += 1 
        elif  a>0 and d>0:       #[4] C[i] left overlap of R[j] right no envelopment
            #print('left overlap\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]
            j += 1
        elif  a<=0 and d>=0:     #[3] C[i] envelopment of R[j]
            #print('envelopment\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]
            j += 1
        elif  a>=0 and d<=0:     #[5] C[i] enveloped by R[j]
            #print('being eveloped\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]
            i += 1 
        if i>=n: i,j = n,j+1 #sticky indecies wait for eachother
        if j>=m: j,i = m,i+1 #sticky indecies wait for eachother
    while len(C) > 0 and C[-1][0] > upper:  C.pop()    
    while len(R) > 0 and R[-1][0] > upper:  R.pop()
    return [C[i] for i in sorted(set(range(len(C))).difference(set(D[1:])))]    