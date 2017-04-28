import os
import sys
import subprocess32 as subprocess
import string
import random

def CheckRG(samtools,bam,out_name,result):
    file = out_name+'.bam.header'
    command = [samtools,'view','-SH',bam,'-o',file]
    print (' '.join(command))
    subprocess.check_output(' '.join(command),stderr=subprocess.STDOUT,shell=True)

    result=[]

    with open(file,'r') as f:
        header = f.readlines()
    for l in range(len(header)):
        i = 0
        if header[l][0:3] == "@RG":
          RG = {x.split(':')[0]:x.split(':')[-1].replace('\n','') for x in header[l].split('@RG')[-1].split('\t')[1:]}
          result.append(RG)

    clean = ['rm','-f',file]
    subprocess.check_output(' '.join(clean),stderr=subprocess.STDOUT,shell=True)
    return result

def GenerateRG(sample):
    if sample == '': sample = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    id = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    lb = 'Solexa_' + sample
    pl = 'Illumina'
    return r'\t'.join(['@RG', 'ID:'+id, 'LB:'+lb, 'PL:'+pl,
                         'PU:'+sample,'SM:'+sample])
