import re
import numpy as np
import fusion_utils as fu

#raw VCF reader that provides similar VCF reading function as HTSeq
#VCF_CHR=0,VCF_POS=1,VCF_ID=2,VCF_REF=3,VCF_ALT=4,VCF_QUAL=5,VCF_FILT=6,VCF_INFO=7,VCF_FORMAT=8,VCF_SAMPLE=9
# def VCF_Reader(vcf_path):
# 
#     # DEBUG
#     # print "Reading VCF: "+str(vcf_path)
# 
#     data,header,raw,i = [],[],[],0
#     with open(vcf_path,'r') as f:
#         for line in f:
#             if line.startswith('#'):
#                 header += [line.split('\n')[0].split('\t')]
#             else:
#                 raw += [line.split('\n')[0].split('\t')]
#                 data += [VariantCall(raw[i])]
#                 
#                 # DEBUG
#                 # print str(data)
#                 
#                 i += 1
#     return data
# 
# class VariantCall:
#     def __init__(self,row):
#         r = row + ['' for i in range(10-len(row))]
#         self.chrom      = r[0]
#         self.pos        = int(r[1])
#         self.id         = r[2]
#         self.ref        = r[3]
#         self.alt        = r[4]
#         self.qual       = r[5]
#         self.filter     = r[6]
#         self.info       = r[7]
#         self.format     = r[8]
#         self.gt         = r[9:]
#     def __enter__(self):
#         return self
# 
#     def __exit__(self,type,value,traceback):
#         return 0

class SVU:
    def __init__(self,vc=None,offset_map=None,ref_path=None,bam_path=None,conf_split=None):
        #parse the less complex fields first
        self.valid_svtypes = {'SUB':0,'RPL':0,
                              'INS':1,'INS:MEI':1,'INS:ME:ALU':1,'INS:ME:L1':1,
                              'DEL':2,'DEL:ME''DEL:ME:ALU':2,'DEL:ME:L1':2,
                              'DUP':3,'DUP:TANDEM':3,'ITX':3,
                              'CNV':4,'INV':5,'CTX':6,'TRA':6,'BND':7}
        # DEBUG
        # print "Breakpoint 1"
        
        self.ref_path = ref_path #can repair the ref consensus string
        self.bam_path = bam_path #can repair the alt consensus string
        self.parse_chrom(vc.chrom)      #string value, trim chrom or chr to chr1->1
        self.pos    = int(vc.pos.pos)            #uint here; changed from vc.pos.pos
        self.id     = vc.id             #string
        self.repair_id()
        self.ref    = vc.ref            #string
        self.repair_ref()        #can get from fasta
        self.alt    = vc.alt            #string
        self.qual   = vc.qual           #string initially
        self.filter = vc.filter.upper() #enforce uppercase
        self.info   = vc.info.upper()   #enforce uppercase
#         self.format = vc.format
#         self.gt     = vc.gt
        #more complex repairs and parsing to populate the object
        self.repair_info()              #clean out any wierd delimiters...
        self.parse_end()         #unint self.end 
        self.parse_svtype()      #string for type but could be int latter self.svtype
        self.parse_svlen()       #uint here from self.svlen
        self.parse_alt_chrom()   #look for alternate chrom tags in the info field
#        self.repair_alt()       #can get from bam or local alignments
        self.repair_qual()       #becomes float value
        self.repair_filter()     #becomes -1,0,1
        self.parse_conf()        #look for confidence values in the info field        
        self.parse_svu(offset_map,conf_split)  #this is the SVU version
        
    def __enter__(self):
        return self

    def __exit__(self,type,value,traceback):
        return 0
    
    def get_sv_types(self):
        return self.valid_svtypes
    
    #find info position invariant key k in s and return string value v
    def get_info_v(self,k,d=';'):
        #match the key in start position, or ; delimited or with a whitespace in front
        p = '\A'+k+'=|['+d+']'+k+'=|[\s]'+k+'='      
        m = re.search(p,self.info)
        if m is None: v = ''
        else:         v = self.info[m.end():].split(d)[0]
        return v
    
    #TO DO parsing Cx and Cy into this form
    def parse_chrom(self,chrom):
        chrom_tag = chrom.split('CHROM')[-1]
        chrom_tag = chrom_tag.split('CHR')[-1]
        chrom_tag = chrom_tag.split('chrom')[-1]
        # chrom_tag = chrom_tag.split('chr')[-1]
        self.chrom = chrom_tag
        #need to fix this alt chrom later x chrom, y chrom
        self.alt_chrom = chrom_tag #did this out of the TRA,BND or INFO tags
    
    #parse alternate chrom embeded in the info field
    def parse_alt_chrom(self):
        self.alt_chrom = self.chrom #TO DO di into the info field some more
    
    #some import parsing methods here to calculate the correct VCF values
    #get the reference ending position from either an info tage or by reading the ref string + pos    
    def parse_end(self):
        end = 1
        try: #need exact matching here
            end = int(self.get_info_v('END'))
        except Exception: #look at ref and alt string lengths
            end = self.pos+len(self.ref)
        self.end = end
        
    def svtype_map(self):
        return {self.valid_svtypes[k]:k for k in self.valid_svtypes}

    #get the svtype by either parsing the info or digging back into the ref->alt strings
    def parse_svtype(self):
        svtype = -1
        try:
            svtype = int(self.valid_svtypes[self.get_info_v('SVTYPE')])
        except Exception: #look at ref and alt
            try:
                svtype = int(self.valid_svtypes[self.get_info_v('MERGE_TYPE')])
            except Exception:
                n,m = len(self.ref),len(self.alt) #dig out from GATK small callers...
                if   n > m:  svtype = self.valid_svtypes['DEL']
                elif n < m:  svtype = self.valid_svtypes['INS']
                elif n == m: svtype = self.valid_svtypes['SUB']
        self.svtype = svtype
        
    #get the svlen from the info field or digging back into the ref->alt strings
    #DEL,INV,SUB should each have a 0 y value
    #DUP should be set to the same start pos for a tandum duplication
    #INS should be the insertion length
    def parse_svlen(self):
        svlen = 0
        v = self.get_info_v('SVLEN')
        #if there is no SVLEN tag, look for an END tag or ALT REF differences
        if v == '': #no SVLEN key in info and we are using the rectangular form for a svu
            if self.svtype != 1:
                svlen = abs(self.end-self.pos)
            elif self.svtype == 0:
                svlen = len(self.alt)
            else: #svtype is SUb or INS
                svlen = abs(len(self.ref)-len(self.alt))
        else:
            try:
                svlen = int(v)
                if svlen < 0: svlen *= -1
                #need to do BND types a bit more to make this work fully
            except Exception:
                svlen = abs(len(self.ref)-len(self.alt))#need sign here?
            if self.end < self.pos+svlen:
                self.end = self.pos+svlen
        self.svlen = svlen
    
    #try to parse out any CIPOS,CIEND,POSrange,ENDrange values from info
    def parse_conf(self):
        conf = [self.pos,self.pos,self.end,self.end]
        s = self.get_info_v('CIPOS')
        if len(s)>0:
            try: conf[0:2] = [int(x) for x in s.split(',')]
            except Exception: pass
            if conf[0]<=0: #-+style setected
                conf[0] = self.pos-conf[0]
                conf[1] = self.pos+conf[1]
        s = self.get_info_v('POSRANGE')
        if len(s)>0:
            try: conf[0:2] = [int(x) for x in s.split(',')]
            except Exception: pass
            if conf[0]<=0: #-+style setected
                conf[0] = self.pos-conf[0]
                conf[1] = self.pos+conf[1]
        s = self.get_info_v('CIEND')
        if len(s)>0:
            try: conf[2:4] = [int(x) for x in s.split(',')]
            except Exception: pass
            if conf[2]<=0: #-+style setected
                conf[2] = self.pos-conf[2]
                conf[3] = self.pos+conf[3]
        s = self.get_info_v('ENDRANGE')
        if len(s)>0:
            try: conf[2:4] = [int(x) for x in s.split(',')]
            except Exception: pass
            if conf[2]<=0: #-+style setected
                conf[2] = self.pos-conf[2]
        self.conf = conf
    
    def repair_id(self):
        self.id = self.id.replace(' ','')
        if self.id=='': self.id='.'
        
    #check and clean if not using the ; delimiter
    def repair_info(self):
        delims = ['\n','\t','\r','; '] #clean up some weird delims are here
        delim = delims[np.argmax([len(self.info.split(i)) for i in delims])]
        self.info = self.info.replace(delim,';')
        self.info = self.info.replace('<','')
        self.info = self.info.replace('>','')
        self.info = self.info.replace(';;',';') #weird end tag in BD
    
    #using the pos and end, repairs the VCF to have the reference string
    def repair_ref(self):
#        ref = ''
#        if self.ref =='' or self.ref =='.': ref = 'N'
#        elif self.svtype == self.valid_svtypes['DEL'] or self.svtype == self.valid_svtypes['SUB'] or \
#             self.svtype == self.valid_svtypes['DUP'] or self.svtype == self.valid_svtypes['INV'] or \
#             self.svtype == self.valid_svtypes['CNV'] or self.svtype == self.valid_svtypes['BND']:
#            try:
#                ref = ru.read_fasta_substring(self.ref_path,self.chrom,self.pos,self.end)
#            except Exception:
#                text = 'ref string repair failed with path=%s, chrom =%s, pos=%s, end=%s'
#                print(text%(self.ref_path,self.chrom,self.pos,self.end))
        self.ref = self.ref.replace(' ','')
        self.ref = self.ref.replace('<','')
        self.ref = self.ref.replace('>','')
        if self.ref =='' or self.ref =='.': self.ref = 'N'
                
    #using the .bam file, repairs the VCF to have the alternate consensus string for reads mapped
    def repair_alt(self):
        #this will have to wait a bit using pysam/pysamstats
        if type(self.alt) is list: self.alt = ''.join(self.alt)
    
    def repair_qual(self):
        qual = 0
        try:
            qual = float(self.qual)
        except ValueError: #got a . or ''
            qual = 0.0
        self.qual = qual
    
    #encode the filter as 1 = PASS, 0 = LowQua, etc
    def repair_filter(self):
        if self.filter=='PASS': fltr = 1
        elif self.filter=='LOWQUAL': fltr = -1
        else: fltr = 0
        self.filter = fltr
        
    #once you have repaired ref and alt consensus do edit distance
    def ref_alt_dist(self):
        self.dist = fu.edit_dist(self.ref,self.alt)
    
    def svtypes(self,t):
        svtypes = {0:'SUB',1:'INS',2:'DEL',3:'DUP',4:'CNV',5:'INV',6:'TRA',7:'BND'}
        return svtypes[t]
    
    def filters(self,f):
        filters = {-1:'LowQual',0:'.',1:'PASS'}
        return filters[f]
        
    #return a searchable dict of preprocessed and cleaned values
    def as_dict(self):
        alt = self.alt
        if len(alt) > 0: alt = alt[0]
        return {'chrom':self.chrom,'pos':self.pos,'end':self.end,
                'ref':self.ref,'alt':alt,'qual':self.qual,'id':self.id,
                'filter':self.filters(self.filter),'svtype':self.svtypes(self.svtype),
                'svlen':self.svlen,'info':self.info}

    def as_vcf_row(self):
        R = ['chrom','pos','id','ref','alt','qual','filter','info']
        D = self.as_dict()
        return [D[r] for r in R]
    
    #TO Do fix the Cx Cy desingnation for TRA, ITX,CTX,BND type SVs
    def parse_vcfu(self):
        self.vcfu = [self.chrom,self.pos,self.end,self.svtype,self.alt_chrom,self.svlen]#use chrom twice here
        
    #return a discrete geometric representation
    #svu = [x1=refpos,x2=refend,t=svtype,y=[[y1,y2]]:contextual with t, wx, wy,{idx}]
    #can use d(ref,alt) metrics in the future...
    #TO DO store more than just the svlength inside the [y]
    def parse_svu(self,O,conf_split=None):
        if conf_split is not None and conf_split:
            if self.conf == [self.pos,self.pos,self.end,self.end]:
                self.svu = [[O[self.chrom]+self.pos,O[self.chrom]+self.end,
                             self.svtype,0,0,1,1]]
            else:
                #outer first, then inner split structure
                self.svu = [[O[self.chrom]+self.conf[0],O[self.chrom]+self.conf[3],
                             self.svtype,O[self.alt_chrom]+self.svlen,O[self.alt_chrom]+self.svlen,1,1],
                            [O[self.chrom]+self.conf[1],O[self.chrom]+self.conf[2],
                             self.svtype,0,0,1,1]]
        else:
            self.svu = [[O[self.chrom]+self.pos,O[self.chrom]+self.end,
                        self.svtype,0,0,1,1]]
            
    def get_svu(self):
        return self.svu
    
    def array_pos_len(self):
        return [sv for sv in self.svu]
    
    def display(self):
        print('pos=\t%s\tend=\t%s'%(self.pos,self.end))
        print('ref=\t%s'%self.ref)
        print('alt=\t%s'%self.alt)
        print('svtype=\t%s\tsvlen=\t%s'%(self.svtype,self.svlen))
        