import os
import sys
import datetime
#relink = os.path.dirname(os.path.abspath(__file__))+'/../../'
#sys.path.append(relink) #go up one in the modules
#import stage_utils as su
#import read_utils as ru

#TIGRA-ext.pl pipeline intersect.bed to tigra.vcf
def tigra_ext_bed_to_vcf(bed_path,sname,ref_seq,vcf_path,header_path=None):
    raw,header,data,clean,final,s = [],[],[],{},[],''
    with open(bed_path,'r') as f:
        raw = [i.split('\t') for i in f.read().split('\n')]
    for i in raw:
        if len(i)==10:
            row = i[8].split('@') #switch ordering
            data += [[row[0],row[1],row[3],row[2]]+row[4:]]
    #do some metadata collection
    refname = ref_seq.keys()[0] #get the ref name
    C = {k:len(ref_seq[refname][k]) for k in ref_seq[refname]} #get the chrom names
    ctg = max([len(k) for k in C])
    cs = sorted(C.keys(), key =lambda x: x.rjust(ctg)) #SORTING ISSUES::::::::::::::::::::::::::::::
    #build a header
    if header_path is None: #default vcf_header is in the data directory of SVCP
        h_path = os.path.dirname(os.path.abspath(__file__))+'../../data/header_template.vcf'
    else:
        h_path = header_path
    with open(h_path,'r') as f: raw = f.readlines()
    if len(raw)>10 and raw[0].startswith('##fileformat=VCFv'): #header has at least 10 VCF lines
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
        header += ['\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sname])] #GT 1/1
        #header complete now parse and generate the data rows before writting data
        for row in data: #could clean duplicates if we want to or pass into fusion framework
            if clean.has_key(row[0]): clean[row[0]] += [row]
            else:                     clean[row[0]]  = [row]
        for c in cs:
            if clean.has_key(c):
                final += clean[c]     
        s = '\n'.join(header)+'\n'
        s += '\n'.join(['\t'.join(i) for i in final])
        with open(vcf_path,'w') as f:
            f.write(s)
            return True
    return False


#ref_seq = {'human_g1k_v37_decoy':ru.read_fasta('/Users/tbecker/Desktop/TEMP/human_g1k_v37_decoy/human_g1k_v37_decoy.fa',True)}
#header_path = '/Users/tbecker/Desktop/TEMP/SVCP/data/header_template.vcf'
#snames = ['HG00096', 'HG00268', 'HG00419', 'HG00759', 'HG01051', 'HG01112', 
#          'HG01500', 'HG01565', 'HG01583', 'HG01595', 'HG01879', 'HG02568', 
#          'HG02922', 'HG03006', 'HG03052', 'HG03642', 'HG03742', 'NA12878', 
#          'NA18525', 'NA18939', 'NA19017', 'NA19238', 'NA19239', 'NA19625', 
#          'NA19648', 'NA20502', 'NA20845']
#for sname in snames: #pull it from a train.py analysis
#    bed_path = '/Users/tbecker/Desktop/TEMP/tigra_Y/'+sname+'/intersect.bed'
#    vcf_path = '/Users/tbecker/Desktop/TEMP/tigra_Y/'+sname+'/'+sname+'_S38.vcf'
#    tigra_ext_bed_to_vcf(bed_path,sname,ref_seq,vcf_path,header_path)
