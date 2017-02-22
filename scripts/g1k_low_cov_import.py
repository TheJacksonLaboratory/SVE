#!/usr/bin/env python
import time
import sys
import os
import multiprocessing as mp
import subprocess32 as subprocess

#wget from the ftp a specified sample, unless it is already in the download.check file
def wget(base_url,log_path,sample):
        print('starting sample %s'%sample)
        output,err = '',''
        #[1]unmapped index
        url = base_url+'/%s/alignment/%s.unmapped*.bam.bai'%(sample,sample)
        command = ['cd','/'.join(log_path.rsplit('/')[0:-1])+'/','&&','wget','-c',url]
        try:
            output += subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT, shell=True)
        except Exception:
            err += '\t'.join(command)+'\n'
            pass
        
        #[2]unmapped bam
        url = base_url+'/%s/alignment/%s.unmapped*.bam'%(sample,sample)
        command = ['cd','/'.join(log_path.rsplit('/')[0:-1])+'/','&&','wget','-c',url]
        try:
            output += subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT, shell=True)
        except Exception:
            err += '\t'.join(command)+'\n'
            pass
        
        url = base_url+'/%s/alignment/%s.mapped*.bam.bai'%(sample,sample)
        command = ['cd','/'.join(log_path.rsplit('/')[0:-1])+'/','&&','wget','-c',url]
        #[3]mapped index
        try:
            output += subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT, shell=True)
        except Exception:
            err += '\t'.join(command)+'\n'
            pass
        
        url = base_url+'/%s/alignment/%s.mapped*.bam'%(sample,sample)
        command = ['cd','/'.join(log_path.rsplit('/')[0:-1])+'/','&&','wget','-c',url]
        #[4]mapped bam
        try:
            output += subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT, shell=True)
        except Exception:
            err += '\t'.join(command)+'\n'
            pass
        
            
        print('output:\n'+output)
        #[3a]execute the command here----------------------------------------------------
        #[3b]do a os directory/data check or a ls type command to check the size
        #[3b]of the produced data to ensure that everything is fine...        
        if err == '':
            with open(log_path+'_'+sample,'w') as f:
                f.write(output)
            return output
        else:
            with open(log_path+'_'+sample,'w') as f:
                f.write(err)
            return 'error on sample %s'%sample

results = []
def collect_results(result):
    results.append(result)  

if __name__ == '__main__':
    #g1k low covergae downloader tool:
    cpus = 4 #number of simultaneous downloads
    base_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/'
    sample_list_path = os.path.dirname(os.path.abspath(__file__))+'/g1k_sample_list.txt'
    high_list_path = os.path.dirname(os.path.abspath(__file__))+'/high_cov_sample_list.txt'
    log_path = os.path.dirname(os.path.abspath(__file__))+'/g1k_download_log'
    
    print('using sample_list: %s'%sample_list_path)
    print('using high_list: %s'%high_list_path)
    print('writing logs to %s'%log_path)
    print('using %s cpus and wget commands'%cpus)
    
    raw_list = []
    high_list = []
    with open(sample_list_path,'r') as f:
        sample_list = f.readlines()
    with open(high_list_path,'r') as f:
        high_list = f.readlines()
    raw_list = [s.replace('\n','') for s in sample_list]
    high_list = [s.replace('\n','') for s in high_list]
    
    #take out the 27 samples and pick another 27*4 = 27*5 = 135 samples saving to a new file list
    #deterministic selection based on alphabetical ordering
    sample_list = []
    for sample in raw_list:
        if sample not in high_list:
            sample_list += [sample]
    pick_list = []
    for i in range(0,len(sample_list),len(sample_list)/108):
        pick_list += [sample_list[i]]
    pick_list = high_list+pick_list
    
    #pick_list = pick_list[0:4] #debug
    
    #start || wget calls
    p1 = mp.Pool(processes=cpus)
    for sample in pick_list: #for each sample calculate all posisble combinations and score them
        p1.apply_async(wget, args=(base_url,log_path,sample), callback=collect_results)
        time.sleep(1)
    p1.close()
    p1.join()
    
    L = []
    for i in results:
        if not i.startswith('error on sample'): L += [i]
        else: print(i)
    print('%s samples were successfully downloaded'%len(L))

