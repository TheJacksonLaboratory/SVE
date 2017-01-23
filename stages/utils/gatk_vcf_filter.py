#GATK 1G vcf filtering tool
import argparse
import glob
import os

def read_gatk_vcf(path):
    header,data = [],[] #init variables
    with open(path,'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.replace('\n','').split('\t')]
            else:
                data += [line.replace('\n','').split('\t')]
    return header,data

def filter_by_sv_len(raw,lower,upper):
    #i = {0:'CHROM',1:'POS',2:'ID',3:'REF',4:'ALT',5:'QUAL',6:'FILTER',7:'INFO',8:'FORMAT',9:'GENOTYPE'}
    data = []
    x = [1 if abs(len(r[3])-len(r[4]))>=lower and abs(len(r[3])-len(r[4]))<=upper else 0 for r in raw]
    for i in range(len(x)):
        if x[i]>0: data += [raw[i]]
    return data

def write_filtered_vcf(header,data,path):
    with open(path+'_S13.vcf','w') as f:
        f.write('\n'.join(['\t'.join(row) for row in header])+'\n')
        f.write('\n'.join(['\t'.join(row) for row in data]))
    print('done')

des ='GATK vcf filtering tool'
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-i', '--sample_dir',type=str, help='sample directory')
parser.add_argument('-o', '--out_dir',type=str, help='output directory')
parser.add_argument('-l', '--lower',type=int, help='lower filter cutoff in bp')
parser.add_argument('-u', '--upper',type=int, help='upper filter cutoff in bp')
args = parser.parse_args()

#now glob this on a folder
lower,upper = args.lower,args.upper
gatk_paths = glob.glob(args.sample_dir+'/*/*.vcf') #weak pattern to glob onto...
out_dir = args.out_dir
if not os.path.exists(out_dir): os.makedirs(out_dir)
for gatk_path in gatk_paths:
    sname = gatk_path.rsplit('/')[-2]
    file_name = gatk_path.rsplit('/')[-1]
    dir_path = gatk_path.replace(file_name,'')
    print('GATK filtering sample %s'%sname)
    header,raw = read_gatk_vcf(gatk_path)
    data = filter_by_sv_len(raw,lower,upper)
    write_filtered_vcf(header,data,out_dir+'/'+sname+'_'+str(lower)+'bp_'+str(upper)+'bp')
