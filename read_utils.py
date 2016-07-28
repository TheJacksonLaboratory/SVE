#modifying this library to have the ability to skip the NNNNN regions of the regnome
#could also get into the alignibilty of a given read length onto a seq like genome strip...

import os
import json
import numpy as np
import HTSeq as ht


#takes in the multi-chrom fasta file and rads it by seq/chrom
#reads each sequence and then writes that portion to the chrom_fasta_dir
#using the chrom_base='' seq_name.fa
#write_fasta_by_chrom(path,path[0:path.rfind('/')],'')

def read_fasta_substring(fasta_path,chrom,pos,end):
    ss = ''
    for s in ht.FastaReader(fasta_path):
        if s.name==chrom:
            ss = s
            return ss[pos:end] #short circuit
    return ss
    
def read_fasta_chrom(fasta_path,chrom):
    ss = ''
    for s in ht.FastaReader(fasta_path):
        if s.name==chrom:
            ss = s
            return ss
    return ss

def read_fasta(fasta_path,dictionary=False,trimN=False):
    ref = None
    if dictionary:
        ref = dict( (s.name, s) for s in ht.FastaReader(fasta_path))
    else:
        ss = []
        for s in ht.FastaReader(fasta_path): ss+=[s]
        ref = ss    
    if trimN:
        if dictionary:
            for k in ref:
                ref[k].seq = ref[k].seq.replace('N','')
        else:
            for i in range(0,len(ref)):
                ref[i].seq = ref[k].seq.replace('N','')
    return ref

#read a genome mask file and convert to a region list of list object
#take out each sequence name and then the sequence and project to bits
#use the bit to integer ranges to calculate the integer ranges on the mask
#add the optional offsets to each of these using the offset map
#def read_fasta_mask(fasta_path,mask={'1':1,'0':0},offset_map=None):
#    M = {}
#    for s in ht.FastaReader(fasta_path):
#        M[s.name] = fu.str2ir(s.seq)
#    return M

def get_fasta_seq_names(fasta_path):
    ss = []
    for s in ht.FastaReader(fasta_path): ss+=[s.name]
    return ss

def get_fasta_seq_lens(fasta_path):
    ss = []
    for s in ht.FastaReader(fasta_path): ss+=[len(s)]
    return ss

def get_fasta_seq_names_lens(fasta_path):
    ss = {}
    for s in ht.FastaReader(fasta_path): ss[s.name]=len(s)
    return ss

def write_fasta(seqs, fasta_path):
    with open(fasta_path, 'w') as fasta:
        for seq in seqs: seq.write_to_fasta_file(fasta)
        return True

#ss is a HTSeq Sequence list?   
def write_fasta_by_chrom(ss, chrom_fasta_dir, chrom_base=''):
    names = []
    for s in ss:
        name = chrom_fasta_dir+'/'+chrom_base+s.name+'.fa'
        names += [name]
        with open(name, 'w') as fasta: s.write_to_fasta_file(fasta)
    return names

def write_fasta_mask(M,json_path):
    with open(json_path,'w') as f:
        json.dump(M,f)
    return True

#compute an expectation given randomly distributed short reads for the RD windows (hist bins)
def expected_window(depth=20,length=100,target=100):
    return int(round((1.0*target/depth)*2.0*length,0))           

def get_coordinate_offsets(json_name):
    path = os.path.dirname(os.path.abspath(__file__))+'/data/'
    info,O = {},{}
    with open(path+json_name,'r') as f:
        info = json.load(f)
    for i in info:
        O[str(i)] = info[i]
    return O

def write_coordinate_offsets(fasta_path,json_path):
    #read in the fasta lengths and then sort them by length into a json offset map
    L = get_fasta_seq_names_lens(fasta_path)
    l_i = list(np.argsort([L[k] for k in L]))[::-1] #sort by max length
    O,offset = {},0 #starts at zero
    for i in l_i:
        O[L.keys()[i]] = offset
        offset = offset+L[L.keys()[i]] #old + new + 1
    with open(json_path,'w') as f:
        json.dump(O,f)
    return True

def get_chrom_dict(json_name):
    path = os.path.dirname(os.path.abspath(__file__))+'/data/'
    info,I = {},{}
    with open(path+json_name,'r') as f:
        info = json.load(f)
    for i in info:
        I[str(i)] = int(info[i])
    return I
    
#input is json data store for svmask regions and  offset map O
#output is a sorted by ref start pos list of list to filter on
def get_mask_regions(json_name,O):
    path = os.path.dirname(os.path.abspath(__file__))+'/data/'
    M = {}
    with open(path+json_name,'r') as f:
        M = json.load(f) #load the svmask per sequence
    N = []    
    for k in M:
        for i in range(len(M[k])):
            N += [[np.uint32(O[k])+np.uint32(M[k][i][0]),np.uint32(O[k])+np.uint32(M[k][i][1])]] #apply offsets
    N = sorted(N,key=lambda x: x[0])
    return N

def write_mask_regions(json_name):
    return True

def bed_mask_to_json_mask(bed_path,json_path):
    bed_data = []
    with open(bed_path, 'r') as f:
        bed_data = [i.split('\t') for i in f.read().split('\n')]
    data = {}
    for row in bed_data: #chrom,start,stop
        if len(row)==3:
            if data.has_key(row[0]): data[row[0]] += [[int(row[1]),int(row[2])]]
            else:                    data[row[0]]  = [[int(row[1]),int(row[2])]]
    for k in data:
        data[k] = sorted(data[k], key = lambda x: x[0])
    #no coordinate sort per contig
    with open(json_path,'w') as f:
        json.dump(data,f)
    return True
    
def get_offsets(chroms,order):
    x,m = 0,{k:0 for k in chroms}
    for k in order:
        m[k] = x
        x += chroms[k]
    return m
   
        
    
    
    
    
    
    




    