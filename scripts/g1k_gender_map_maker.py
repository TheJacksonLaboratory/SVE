#takes the bam directory and pulls out the sample names
#and then makes a gender map for GenomeSTRiP2.0
import os
import sys
import glob

def build_gender_map(bam_dir,g1k_gender_file,gender_map_file):
    bams = glob.glob(bam_dir+'*.bam')
    samples = [b.split('/')[-1].split('.')[0] for b in bams]
    with open(g1k_gender_file,'r') as f:
        raw = f.readlines()[0] #specific to this file
    data = [row.split('\t') for row in raw.split('\r')]
    header,i,j,gender_map = data[0],0,0,''
    while header[i]!='Gender': i+=1
    while header[j]!='Individual ID': j += 1    
    data = data[1:]
    for row in data:
        if row[j] in samples:
            gender_map += row[j]+' '+row[i]+'\n'
    with open(gender_map_file,'w') as f:
        f.write(gender_map)
    return True
    
    
bam_dir = '/data/ch-lee-lab/g1k_low_cov_group1_bams/'
g1k_gender_file = '/Users/tbecker/Desktop/TEMP/SVE/data/g1k_gender_trio_20130606.tsv'
gender_map_file = '/Users/tbecker/Desktop/TEMP/SVE/data/g1k_low_cov_group1.map'
build_gender_map(bam_dir,g1k_gender_file,gender_map_file)