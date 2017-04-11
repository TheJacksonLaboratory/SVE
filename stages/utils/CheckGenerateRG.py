import os
import sys
import subprocess32 as subprocess
import string
import random

def CheckRG(samtools,bam,out_name,result):
    file = out_name+'.RG'
    command = [samtools,'view','-SH',bam,'>',file]
    print (' '.join(command))
    subprocess.check_output(' '.join(command),stderr=subprocess.STDOUT,shell=True)

    result=[]

    with open(file,'r') as f:
        header = f.readlines()
    for l in range(len(header)):
        i = 0
        RG = {x.split(':')[0]:x.split(':')[-1].replace('\n','') for x in header[l].split('@RG')[-1].split('\t')[1:]}
        result.append(RG)
    return result

def GenerateRG(sample):
    id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    pl = 'Illumina'
    return '@RG\tID:'+id+'\tSM:'+sample+'\tPL:'+pl
