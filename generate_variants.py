#refactor as variant_utils, include in cvm functions

import numpy as np
import HTSeq as ht
import time


def gen_random_seqs(seqs={'chr1':1E5},alphabet=['A','C','G','T'],
                   prob=[0.25,0.25,0.25,0.25],method='slow'):
    ss = []
    sigma = np.asarray(alphabet,dtype='S1') #do set difference for SNPs
    for k in seqs:
        if method=='slow': #slow method
            s = ''
            for i in range(0,seqs[k]):
                s += np.random.choice(sigma,size=1,p=prob).tostring()
            seq = ht.Sequence(s,k)
        elif method=='fast':                                 #~100x faster
            s = np.random.choice(sigma,size=seqs[k],p=prob)  #size param
            seq = ht.Sequence(s.tostring(),k)
        else: raise KeyError
        ss += [seq]
    return ss[::-1]
        
#[0] given a type and variant size interval VS --> {INS:1E0}

def gen_chrom_lens(L,chroms):
    c = len(chroms)
    ls = np.random.geometric(0.4,L) #partitioned lens of each chrom that will sum to L
    LS = np.histogram(ls,bins=list(set(ls)))  #bin the geometric dist
    chr_lens = list(LS[0][0:c])               #chop off largest c
    r = L-sum(chr_lens)                       #calculate leftovers
    chr_lens = [i+int(r/c) for i in chr_lens] #distribute leftovers
    if r%c!=0: chr_lens[0]+=L-sum(chr_lens)   #correct to sum to L
    return {chroms[i]:chr_lens[i] for i in range(0,c)}

#[1] generate uniform positions and attach up to 50% edit rate
#    this is done by taking the target length 1E8 and dividing by 2*VS
def gen_mut_pos(ref_len, mut_len, mut_rate):
    start = time.time()
    L = np.asarray(range(0,int(ref_len)))
    s = get_num_mut(ref_len, mut_len, mut_rate)
    if s > 1E6: mut_pos = np.sort(np.random.choice(L,size=s,replace=True,p=None))
    else:       mut_pos = np.sort(np.random.choice(L,size=s,replace=False,p=None))
    stop = time.time()
    print('Position list with %s variants generated'%s)
    print('Position list elapsed time was %s seconds\n' %(stop-start))
    return mut_pos

def get_num_mut(full_len, mut_len, mut_rate):
    return int(full_len*mut_rate/mut_len)

#take in a list of named sequences
#and copy a random section of length l
#to every other chrom for shard testing of bam files...
def translocate(ss,k,l):
    #get the ss for k
    for s in ss:
        if s.name == k:
            t = s
            break
    l = int(l)
    L = np.asarray(range(0,len(t.seq)))
    a = np.random.choice(L-l,size=1)
    copy = t.seq[a:max(a+l,len(t.seq))] #this will be length l
    for s in ss:
        if s.name != k:
            L = np.asarray(range(0,len(t.seq)))
            a = np.random.choice(L-l,size=1)
            s.seq = s.seq[0:a]+copy+s.seq[a:] #insert the copied sequence
    return ss
            
#[2] for each position in the uniform distribution apply the edit
#    and update the current position (DEL = -L, INS = +L, SUB = 0, INV = 0, DUP = L*CN)
#    the final result should throw and error or not tabulate the edit
#    if the simple merging step consumes preselected positions
#    fa_path=fasta file path, vca_path = variant call array pickle path
#    l = full ref length, m = var length, var_type = INS,DEL,SUB,INV,DUP, pos_list = int pos list
def gen_var_calls(ref_fa_path,chrom,mut_len,mut_type,mut_pos):
    refdict = dict((s.name,s) for s in ht.FastaReader(ref_fa_path))
    full_ref = refdict[chrom].seq                 #ref string in memory,use .tostring()        <-----
    l = len(full_ref)                             #original length
    #function map of the form ref,alt
    vt = {'INS':lambda s: ['',gen_alt_seq(s)],
          'DEL':lambda s: [s,''],
          'SUB':lambda s: [s,gen_alt_seq(s)],
          'INV':lambda s: [s,s[::-1]],
          'DUP':lambda s: [s['s'],dup(s['s'],s['i'])]} #set to do copy number : 1=>5
    
    start = time.time()
    i,vca,cur,info = 0,[],0,''     #init position list pointer
    #starting left and going right on the position list
    #apply each variant consuming intermediate positions along the way
    for pos in mut_pos:
        if cur<=pos and pos+mut_len <=l:
            cur = pos+mut_len
            if type(mut_type) is dict:#variants with a params
                if mut_type.has_key('DUP'):
                    [ref,alt] = vt['DUP']({'s':full_ref[pos:cur],'i':mut_type['DUP']})
                    info = 'TYPE= DUP:'+str(mut_type['DUP'])+'; END= '+str(cur)
            else:                     #vairants without params
                [ref,alt] = vt[mut_type](full_ref[pos:cur])
                info = 'TYPE= '+mut_type+'; END= '+str(cur)
            vca+=[gen_var_call(chrom,pos,ref,alt,info)]  
            i+=1
    stop = time.time()
    #now save the vca pickle but will return it for now...
    print('%s variants requested, %s generated in response'%(len(mut_pos),len(vca)))
    print('Variant generation elapsed time was %s seconds\n' %(stop-start))
    if len(vca)<1 and l>mut_len:
        mut_pos = gen_mut_pos(l,mut_len,(1.1*mut_len/l))
        return gen_var_calls(ref_fa_path, chrom, mut_len,mut_type,mut_pos)
    else:
        return vca         

#given a possible empty string = ref and a length = default is 0
#either produce a new string that has a different random nuclideotide
#for use as a substitution generator or just a full random insertion
#generator for insertion operations
def gen_alt_seq(ref,alphabet=['A','C','G','T'],length=0):
    alt = ''
    l = len(ref)
    sigma = set(alphabet) #do set difference for INS/SUB
    for i in range(0,l):         #have to do this for each position
        cand = np.sort(list(sigma-set([ref[i]]))) #candidates to make a sigle mut
        alt += np.random.choice(cand)
    return alt

#one sigma look-back
def gen_alt_seq2(ref,length=0):
    alt = ''
    l = len(ref)
    dna = set(['A','C','G','T']) #do set difference for INS/SUB
    if l>0: alt += np.random.choice(list(dna-set([ref[0]])))[0] #get one new one
    for i in range(1,l):         #have to do this for each position
        cand = np.sort(list(dna-set([ref[i]]))) #candidates to make a snp
        #some correction code maybe needed?
        prob = [0.475,0.475,0.475]     #one will be set to 0.2 to penalize similiar
        for j in range(0,3):     #slight random selection
            if alt[i-1]==cand[j]: prob[j]=0.05  #adjustment via single lookback chain
        alt += np.random.choice(cand) #get one new one
    if l<=0: alt = np.random.choice(list(dna),size=length).tostring()
    return alt
    
    
#[3] Keep track pf the applyed edits by making a HTSeq VariantCall
#    alt, chrom, filter, format, id, info, pos, qual, ref, samples
def gen_var_call(chrom,pos,ref,alt,info):
    vc = ht.VariantCall(chrom,pos,'.',ref,alt,None,None,info)
    return vc

def dup(s,i):
    return ''.join([s for j in range(0,i)])

def get_dup(ref,alt):
    a,r = len(alt),len(ref)
    if r>0 and a>1 and a%r==0: #a is non-negative integer multiple of r
        m = alt.split(ref)     #split the alt by the ref
        if all([s == '' for s in m]): return len(m)-1
    return 0

def check_mut_lens(full_ref,vca,mut_type,muts):
    tests = []
    for mut in muts:
        l,m = len(full_ref),len(mut)
        if type(mut_type) is dict:#variants with a params
            if mut_type.has_key('DUP'):
                delta = sum([len(vc.alt)-len(vc.ref) for vc in vca])
        else: #do all the other ones here...
            if mut_type == 'INS':
                delta = sum([len(vc.alt)-len(vc.ref) for vc in vca])
            if mut_type == 'DEL':
                delta = sum([len(vc.ref)-len(vc.alt) for vc in vca])
            if mut_type == 'SUB' or mut_type == 'INV': delta = 0
        test = [l+delta==m]
        if not test: print('issue found:%s'%([l,delta,m],))
        tests += [test]
    return all(tests)

#given the full refrence string = full_ref and the 
#variant call array = vca produce the mutated string = mut
#as quickly as possible...
def apply_var_calls(full_ref,vca):
    mut,pos = '',0 #mut is the full mutation sequence, pos is old until the end of the if...
    for i in range(0,len(vca)): #where to stat = pos, what the ref looks like, what the alternate is
        #copy out the inside sections up to the variant
        ref,alt = vca[i].ref,vca[i].alt
        if ref=='':                               #INS
            mut += full_ref[pos:vca[i].pos]+alt
            pos = vca[i].pos          #no skiping for ins in ref 
        elif alt=='':
            mut += full_ref[pos:vca[i].pos]+''    #DEL
            pos = vca[i].pos+len(ref) #skip the del section in ref
        elif alt==ref[::-1]:
            mut += full_ref[pos:vca[i].pos]+alt   #INV
            pos = vca[i].pos+len(alt) #move past the INV sectio
        elif len(alt)==len(ref) and alt!=ref:     #SUB
            mut += full_ref[pos:vca[i].pos]+alt
            pos = vca[i].pos+len(alt) #move past the SUB section
        elif get_dup(ref,alt)>0:
            mut += full_ref[pos:vca[i].pos]+alt   #DUP
            pos = vca[i].pos+len(ref) #move past the DUP section
        else: raise KeyError
    mut+=full_ref[pos:] # left over sections on RHS
    return mut

def str2seq(string,name):
    return ht.Sequence(string,name)