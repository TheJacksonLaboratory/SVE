#this is genomestrip DEL call splitter
import argparse
import os

def read_genomestrip_del_genotypes(gs_path,cutoff=1.0):
    header,data = [],[] #init variables
    with open(gs_path,'r') as f:
        for line in f:
            if not line.startswith('##GATK'):
                if line.startswith('#'):
                    header += [line.replace('\n','').split('\t')]
                else:
                    row = line.replace('\n','').split('\t')
                    row[2] = row[2].replace('DEL','GS')
                    data += [row]
            if line.startswith('#CHROM'):
                i,h = 0,line.replace('\n','').split('\t')
                samples = {}
                for j in range(len(h)): #should be 9 and up...
                    if h[j]=='FORMAT': i = j+1
                j = len(h)
                h = h[i:]
                for x in range(len(h)):
                    samples[h[x].rsplit('_')[0]] = x+9
    H = []
    for h in header: #take out the format
        if not h[0].startswith('##FORMAT='): H += [h]
    H[-1] = header[-1][0:8] #clip off the sample names
    #scan through and look for the last INFO TAG and append a new row for ICN:8 info
    S = {k:[] for k in samples}
    for i in range(len(data)):
        for k in samples:
            cn = 2.0
            try:
                fltr = data[i][6]
                cn = float(data[i][samples[k]].split(':')[1])
                if cn <= cutoff:# and fltr=='PASS':
                    svtype = '<DEL>' #pick the correct SV TYPE
                    info = data[i][7].split(';')
                    info = ';'.join(info[2:3]+info[-3:])
                    info = info.replace('SVTYPE=DEL','SVTYPE=DEL')
                    info = info.replace('SVTYPE=CNV','SVTYPE=DEL') #and set the info tag too               
                    S[k] += [data[i][0:4]+[svtype]+data[i][5:6]+[fltr]+[info+';ICN='+str(int(round(cn,0)))]]
#                if cn > 2 and cn < 20:# and fltr=='PASS':
#                    svtype = '<DUP>' #pick the correct SV TYPE
#                    info = data[i][7].split(';')
#                    info = ';'.join(info[2:3]+info[-3:])
#                    info = info.replace('SVTYPE=DEL','SVTYPE=DUP')
#                    info = info.replace('SVTYPE=CNV','SVTYPE=DUP') #and set the info tag too               
#                    S[k] += [data[i][0:4]+[svtype]+data[i][5:6]+[fltr]+[info+';ICN='+str(cn)]]
            except Exception:
                pass
    return S,H
 
def read_genomestrip_cnv_genotypes(path,cutoff=3):
    header,data = [],[] #init variables
    with open(path,'r') as f:
        for line in f:
            if not line.startswith('##GATK'):
                if line.startswith('#'):
                    header += [line.replace('\n','').split('\t')]
                else:
                    data += [line.replace('\n','').split('\t')]
            if line.startswith('#CHROM'):
                i,h = 0,line.replace('\n','').split('\t')
                samples = {}
                for j in range(len(h)): #should be 9 and up...
                    if h[j]=='FORMAT': i = j+1
                j = len(h)
                h = h[i:]
                for x in range(len(h)):
                    samples[h[x]] = x+9
    H = []
    for h in header: #take out the format
        if not h[0].startswith('##FORMAT='): H += [h]
    H[-1] = header[-1][0:8] #clip off the FORMAT
    S = {k:[] for k in samples}
    for i in range(len(data)):
        for k in samples:
            cn = 2
            try: cn = int(data[i][samples[k]].split(':')[1])
            except Exception: pass
#            if cn < 2:
#                data[i][4] = '<DEL>' #pick the correct SV TYPE
#                data[i][7] = data[i][7].replace('SVTYPE=CNV','SVTYPE=DEL') #and set the info tag too
#                S[k] += [data[i][0:8]]
            if cn >= cutoff:
                data[i][4] = '<DUP>' #pick the correct SV TYPE
                data[i][7] = data[i][7].replace('SVTYPE=CNV','SVTYPE=DUP') #and set the info tag too
                S[k] += [data[i][0:8]]
    return S,H

def merge_genomestrip_calls(S1,S2,H1,H2):
    S,H,j = {},[],0
    for h in H1: H += [h]
    for i in range(len(H2)):
        if not H2[i] in H and not H2[i][0].startswith('##fileDate='):
            H = H[0:i+j]+[H2[i]]+H[i+j:]
            j += 1
    K = {0:[],1:[],2:[]} #find all the matching and single keys
    for k in S1:
        if S2.has_key(k): K[0] += [k]
        else:             K[1] += [k]
    for k in S2:
        if S1.has_key(k): K[0] += [k]
        else:             K[2] += [k]
    for k in K[1]: S[k] = S1[k]
    for k in K[2]: S[k] = S2[k]
    for k in K[0]: #sort and merge these ones
        i,j,n,m = 0,0,len(S1[k]),len(S2[k])
        while i+j < n+m:
            if i < n and j < m:
                if S1[k][i][0:2]<S2[k][j][0:2]:
                    if S.has_key(k):
                        S[k] += [S1[k][i]]
                    else:
                        S[k]  = [S1[k][i]]
                    i += 1
                else:
                    if S.has_key(k):
                        S[k] += [S2[k][j]]
                    else:
                        S[k]  = [S2[k][j]]
                    j += 1
            else:
                if i < n:
                    if S.has_key(k):
                        S[k] += [S1[k][i]]
                    else:
                        S[k]  = [S1[k][i]]
                    i += 1
                else:
                    if S.has_key(k):
                        S[k] += [S2[k][j]]
                    else:
                        S[k]  = [S2[k][j]]
                    j += 1
    return S,H

def write_vcfs(S,header,path):
    for k in S:
        with open(path+'/'+k+'_S14.vcf','w') as f:
            f.write('\n'.join(['\t'.join(row) for row in header])+'\n')
            f.write('\n'.join(['\t'.join(row) for row in S[k]])+'\n')
    print('done')

if __name__ == '__main__':
    des = """
    GenomeSTRiP2.0 SVDiscovery to DEL, CNVDiscovery to DUP genotype sample splitter.
    Reads genotyped VCF files and splits by sample merging the DEL calls and CNV calls
    with a user definable integer copy number threashold, del_cutoff, dup_cutoff"""
    parser = argparse.ArgumentParser(description=des)
    parser.add_argument('-d', '--del_vcf',type=str, help='SVDiscovery genotype VCF')
    parser.add_argument('-c', '--cnv_vcf',type=str, help='CNVDiscovery genotype VCF')
    parser.add_argument('-o', '--out_dir',type=str, help='split+infered+merged per sample VCF folder')
    parser.add_argument('-a', '--del_cutoff',type=float, help='float values less than or equal to a are DEL calls [1.0]')
    parser.add_argument('-b', '--cnv_cutoff',type=float, help='integer rounded values greater than or equal to b are DUP calls [3]')
    args = parser.parse_args()
    
    if args.del_vcf is not None:
        del_vcf = args.del_vcf
    else:
        print('SVDiscovery Step genotyped VCF file not found')
        raise IOError
    
    if args.cnv_vcf is not None:
        cnv_vcf = args.cnv_vcf
    else:
        print('CNVDiscovery Step genotyped VCF file not found')
        raise IOError
        
    if args.out_dir is not None:
        out_dir = args.out_dir
    else:
        print('output directory not specified')
        raise IOError
    
    if args.del_cutoff is not None:
        del_cutoff = args.del_cutoff
    else:
        del_cutoff = 1.0
    
    if args.cnv_cutoff is not None:
        cnv_cutoff = int(round(args.cnv_cutoff,0))
    else:
        cnv_cutoff = 3
    gs_del_path = del_vcf #'/Users/tbecker/Desktop/TEMP/vcf_experiments/population_g1k_50X.genotypes.vcf'
    gs_cnv_path = cnv_vcf #'/Users/tbecker/Desktop/TEMP/vcf_experiments/gs_cnv.genotypes.vcf'
    S1,H1 = read_genomestrip_del_genotypes(gs_del_path)
    S2,H2 = read_genomestrip_cnv_genotypes(gs_cnv_path)
    S,H   = merge_genomestrip_calls(S1,S2,H1,H2)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    vcf_path = out_dir #'/Users/tbecker/Desktop/TEMP/P3_gs_del_cnv_vcf/'
    write_vcfs(S,H,vcf_path)
