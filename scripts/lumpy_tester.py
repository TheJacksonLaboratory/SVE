import os
import sys
import glob

file_list = []
for i in range(len(file_list)):
    file_list[i] = file_list[i].replace('\n','')

#treat each */sname.* as a sample with three files: -B, -D, -S for lumpy express
samples = {}
for i in range(len(file_list)):
    stub  = file_list[i].rsplit('/')[-1]
    sname = stub.split('.')[0]
    d = stub.find('.disc.')
    s = stub.find('.split.')
    if samples.has_key(sname):
        if   d!=-1: samples[sname]['D'] = file_list[i]
        elif s!=-1: samples[sname]['S'] = file_list[i]
        else:       samples[sname]['B'] = file_list[i]
    else:
        if   d!=-1: samples[sname] = {'D':file_list[i]}
        elif s!=-1: samples[sname] = {'S':file_list[i]}
        else:       samples[sname] = {'B':file_list[i]}
#check there are files for each
B,D,S = [],[],[]
for sname in sorted(samples):
    if not (samples[sname].has_key('D') and samples[sname].has_key('S') and samples[sname].has_key('B')):
        print('files list malformed for sname:\t%s'%sname)
        samples.pop(sname)
    else:
        B += [samples[sname]['B']]
        D += [samples[sname]['D']]
        S += [samples[sname]['S']]
        
#samples is goo to go now
command = ['/data/ch-lee-lab/software/anaconda/bin/python',
           '/data/ch-lee-lab/software/lumpy-sv/bin/lumpyexpress',
           '-T /data/tbecker/cloudSV_R1/S18',
           '-B',','.join(B),'-D',','.join(D),'-S',','.join(S),
           '-m 2','-o /data/tbecker/cloudSV_R1/P4_all_S18.vcf']