#path resolution tester for automated path construction
import os
import subprocess32 as subprocess

software = os.path.dirname(os.path.abspath(__file__))+'/../../../'
command = [software+'/samtools-1.0/bin/samtools','view','-SH','/data/ch-lee-lab/g1k_P4_illumina_80X/HG00512.high_cov.bam']
output = []
try:
    output = subprocess.check_output(command,stderr=subprocess.STDOUT)
except Exception as E:
    print('error occured')
    pass
print(output)