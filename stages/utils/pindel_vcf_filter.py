#pindel 2G vcf filtering tool
import argparse
import glob
import os

def read_pindel_vcf(path):
    header,data = [],[] #init variables
    with open(path,'r') as f:
        for line in f:
            if line.startswith('#'):
                header += [line.replace('\n','').split('\t')]
            else:
                data += [line.replace('\n','').split('\t')]
    return header,data

#filter with a lower and upperbound
#can clean off the assembly to make a smaller file
def filter_by_sv_len(raw,lower,upper,clean=True):
    #i = {0:'CHROM',1:'POS',2:'ID',3:'REF',4:'ALT',5:'QUAL',6:'FILTER',7:'INFO',8:'FORMAT',9:'GENOTYPE'}
    data = []
    for i in range(len(raw)):
        x = 0
        try: x = int(raw[i][7].split('SVLEN=')[-1].split(';')[0])
        except Exception: pass
        if x>=lower and x <= upper:
            if clean: data += [raw[i][0:3]+['.','.']+raw[i][5:]]
            else:     data += [raw[i]]
    return data    

def write_filtered_vcf(header,data,path):
    with open(path+'_S36.vcf','w') as f:
        f.write('\n'.join(['\t'.join(row) for row in header]))
        f.write('\n'.join(['\t'.join(row) for row in data]))
    print('done')

des ='Pindel vcf filtering tool'
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-i', '--sample_dir',type=str, help='sample directory')
parser.add_argument('-o', '--out_dir',type=str, help='output directory')
parser.add_argument('-l', '--lower',type=int, help='lower filter cutoff in bp')
parser.add_argument('-u', '--upper',type=int, help='upper filter cutoff in bp')
args = parser.parse_args()

#now glob this on a folder
lower,upper = args.lower,args.upper
pindel_paths = glob.glob(args.sample_dir+'/*/*.vcf') #weak pattern to glob onto...
out_dir = args.out_dir
if not os.path.exists(out_dir): os.makedirs(out_dir)
for pindel_path in pindel_paths:
    sname = pindel_path.rsplit('/')[-2]
    file_name = pindel_path.rsplit('/')[-1]
    dir_path = pindel_path.replace(file_name,'')
    print('GATK filtering sample %s'%sname)
    header,raw = read_pindel_vcf(pindel_path)
    data = filter_by_sv_len(raw,lower,upper)
    write_filtered_vcf(header,data,out_dir+'/'+sname+'_'+str(lower)+'bp_'+str(upper)+'bp')
                   
                   
