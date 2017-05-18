import copy
import csv
import gc
import json
import gzip
import glob
import itertools as it
import numpy as np
cimport numpy as np
cimport cython

# read all rle encoded gzip compressed freq files from a glob directory root
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def read_gzip_freq(str freq_rle_directory, list glob_inc = []):
    cdef str s,line
    cdef list G    #glob list of string paths
    cdef dict R,X  #result and intermediate result
    R,X,G = {},{},glob.glob(freq_rle_directory+'*.freq.gz') #all files ending with .freq.rle.gz
    for s in glob_inc: G.remove(s)
    for g in G:
        with gzip.GzipFile(g, 'r') as f:
            for line in f: #not sure what line is here... str?
               X  = json.loads(line)
               R[X.keys()[0]] = X[X.keys()[0]]
    return R

# write rle encoded gzip compressed freq files
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def write_gzip_freq(dict R, str refbase_directory):
    cdef str s,json_path
    for s in R:
        json_path = refbase_directory+'.'+s+'.freq.gz'
        with gzip.GzipFile(json_path, 'w') as f:
            f.write(json.dumps({s:R[s]}))

#each chrom is denoted with a >
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def read_fasta(str path, list chroms):
    cdef int i
    cdef str x,raw
    cdef list splits,lines,keys
    cdef dict ss
    raw,splits,lines,ss = '',[],[],{}
    with open(path,'r') as f: raw = f.read()
    splits = raw.split('>')[1:]
    raw = ''
    gc.collect()
    for i in range(len(splits)):
        lines = splits[i].split('\n')
        lines[0] = lines[0].split(' ')[0] #whitespace
        ss[lines[0]] = bytearray(''.join(lines[1:]))
    keys = ss.keys()
    for x in keys:
        if x not in chroms: ss.pop(x)
    return ss

#KMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMER
#[0] use the sqlite3 DBO to do this
#[1] generate the tabel with inserts of kmer as pkey and count as count
#[2] scan the table and return the count at each position (in ||)
#KMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMERKMER
 
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def edit_dist(str s1, str s2, list w):
    cdef unsigned int i,j,k,u,v
    u,v,k = len(s1),len(s2),2
    if u<v: u,v,s1,s2 = v,u,s2,s1 #flip to longest of the two
    cdef np.ndarray[unsigned int, ndim=2] D = np.zeros([u+1,k], dtype=np.uint32)
    for i in range(u+1):   D[i][0] = i
    for j in range(v%k+1): D[0][j] = j
    for j in range(1,v+1):
        for i in range(1,u+1):
            if s1[i-1] == s2[j-1]:    
                D[i][j%k] = D[i-1][(j-1)%k] #matching
            else:                           #mismatch, del, ins, sub 
                D[i][j%k] = min(D[i-1][j%k]+w[0],D[i][(j-1)%k]+w[1],D[i-1][(j-1)%k]+w[2])
    return D[u][v%k]
 
#BIT VECTOR01010101010101010101010101010101010101010101010101010101010101
#read from P and writes integer ranges to R and returns end point j   
@cython.boundscheck(False)
def bit2ir(np.ndarray[np.int64_t, ndim=1] P not None,
           np.ndarray[np.uint32_t,ndim=2] R not None):
    cdef np.uint32_t i,j
    cdef np.int64_t a
    j = 0    
    if P.shape[0] > 1:
        a = P[0]
        for i in range(0,P.shape[0]-1):
            if P[i+1]-P[i]>1:
                R[j][0] = <np.uint32_t>a
                if P[i] == a:
                    R[j][1] = <np.uint32_t>a
                else:
                    R[j][1] = <np.uint32_t>P[i]+1
                a = P[i+1]
                j += 1
        R[j][0] = <np.uint32_t>a
        if P[i+1] == a:
            R[j][1] = <np.uint32_t>a
        else:
            R[j][1] = <np.uint32_t>P[i+1]+1
        j += 1
    elif P.shape[0] > 0:
        R[j][0] = <np.uint32_t>P[0]
        R[j][1] = <np.uint32_t>P[0]
        j += 1
    return j

#hard coded '1'==1 is mask me please...
@cython.boundscheck(False)
def str2ir(str C):
    cdef np.uint32_t i,j,m
    m = len(C)
    cdef list N = []
    cdef np.ndarray B = np.zeros([m,1],dtype=np.uint8)
    cdef np.ndarray R = np.zeros([m/2+1,2], dtype=np.uint32)    
    for i in range(m):
        if C[i]=='1': B[i] = 1
    j = <np.uint32_t>bit2ir(np.where(B)[0],R)
    for i in range(0,j):
        N += [[<long>R[i][0],<long>R[i][1]]]
    return N

@cython.boundscheck(False)
@cython.nonecheck(False) 
def jaccard_score(list C1, list C2, bint self_merge=False):
    cdef long i
    cdef long double x,y,z
    if self_merge:    
        C1 = merge_regions(C1)
        C2 = merge_regions(C2)         
    I,U,D1,D2 = LR(C1,C2)
    x = np.sum([<long double>np.log(<long double>I[i][1]-<long double>I[i][0]+1) for i in range(len(I))])
    y = np.sum([<long double>np.log(<long double>U[i][1]-<long double>U[i][0]+1) for i in range(len(U))])
    if y <= 0.0: z = <long double>(0.0)
    else: z = x/y    
    return z

@cython.boundscheck(False)
@cython.nonecheck(False) 
def log1p_feature_magnitudes(list C1, list C2, bint self_merge=False):
    cdef long i
    cdef double f1,f2,f3,f4
    if self_merge:    
        C1 = merge_regions(C1)
        C2 = merge_regions(C2)           
    I,U,D1,D2 = LR(C1,C2)
    f1 = np.sum([np.log1p(np.abs(<double>I[i][1]-<double>I[i][0])) for i in range(len(I))])
    f2 = np.sum([np.log1p(np.abs(<double>U[i][1]-<double>U[i][0])) for i in range(len(U))])
    f3 = np.sum([np.log1p(np.abs(<double>D1[i][1]-<double>D1[i][0])) for i in range(len(D1))])
    f4 = np.sum([np.log1p(np.abs(<double>D2[i][1]-<double>D2[i][0])) for i in range(len(D2))])
    return [f1,f2,f3,f4]

@cython.boundscheck(False)
@cython.nonecheck(False) 
def feature_magnitudes(list C1, list C2, bint self_merge=False):
    cdef long i
    cdef double f1,f2,f3,f4
    if self_merge:    
        C1 = merge_regions(C1)
        C2 = merge_regions(C2)           
    I,U,D1,D2 = LR(C1,C2)
    f1 = np.sum([np.abs(<double>I[i][1]-<double>I[i][0])+1   for i in range(len(I))])
    f2 = np.sum([np.abs(<double>U[i][1]-<double>U[i][0])+1   for i in range(len(U))])
    f3 = np.sum([np.abs(<double>D1[i][1]-<double>D1[i][0])+1 for i in range(len(D1))])
    f4 = np.sum([np.abs(<double>D2[i][1]-<double>D2[i][0])+1 for i in range(len(D2))])
    return [f1,f2,f3,f4]

@cython.boundscheck(False)
@cython.nonecheck(False) 
def brkpt_score(list X):
    cdef long i
    cdef list Y,Z
    cdef double rm,rs,lm,ls
    Y,Z = [],[]
    for i in range(len(X)):
        if X[i] is not None:
            Y += [[X[i][0],X[i][1]]]
    if len(Y)>0:
        rm = <double>np.mean([<double>Y[i][0] for i in range(len(Y))])
        rs = <double>np.std([<double>Y[i][0] for i in range(len(Y))])
        lm = <double>np.mean([<double>Y[i][1] for i in range(len(Y))])
        ls = <double>np.std([<double>Y[i][1] for i in range(len(Y))])
        Z = [rm,rs,lm,ls]
    return Z
    
@cython.boundscheck(False)
@cython.nonecheck(False) 
cdef double overlap(long c_i0,long c_i1,long c_j0,long c_j1,bint check_x=True,bint check_y=False):
    cdef long l
    cdef double i,u,x
    if c_i1>=c_i0 and c_j1>=c_j0:
        l = (c_i1-c_i0+1)+(c_j1-c_j0+1)                     #total lengths
        u = <double>min(l,max(c_i1,c_j1)-min(c_i0,c_j0)+1)  #total union area
        i = 1.0*<double>abs(c_i0-c_j0)+abs(c_i1-c_j1)
        x = max(0.0,u-i)/u
    else:
        x = 0.0
    return x

@cython.boundscheck(False)
@cython.nonecheck(False) 
cdef double overlap2(long c_i0,long c_i1,long c_j0,long c_j1,bint check_x=True,bint check_y=False):
    cdef double a,b,c,d,i,u
    if c_i0<=c_i1: a,b = <double>c_i0,<double>c_i1
    else:          a,b = <double>c_i1,<double>c_i0    
    if c_j0<=c_j1: c,d = <double>c_j0,<double>c_j1
    else:          c,d = <double>c_j1,<double>c_j0
    i = abs(a-c)+abs(b-d)                           
    u = min((b-a+1)+(d-c+1),max(b,d)-min(a,b)+1) 
    return max(<double>0.0,(u-i)/u)

#number of n calls intersected by one of m calls with r% reciprocal overlap
#number of m calls intersected by one of n calls with r% reciprocal overlap
@cython.boundscheck(False)
@cython.nonecheck(False) 
def metric_score_zig(list C1,list C2,double r,bint self_merge=False):
    cdef long i,j,n,m
    cdef double x,y
    #cdef dict X
    if self_merge:
        C1 = merge_regions(C1)
        C2 = merge_regions(C2)
    X = {'l':[],'r':[]} #breakpoint diferential lists
    n,m = len(C1)-1,len(C2)-1
    i,j = 0,0
    x,y = 0.0, 0.0
    if n > 0 and m > 0:
        while i+j < n+m:
            if   C2[j][0]>C1[i][1]: i+=1
            elif C1[i][0]>C2[j][1]: j+=1
            elif overlap(<long>C1[i][0],<long>C1[i][1],<long>C2[j][0],<long>C2[j][1]) >= r: 
                x,y = x+1.0,y+1.0
                #X['l'] += [[<long>C1[i][0]-<long>C2[j][0]]]
                #X['r'] += [[<long>C1[i][1]-<long>C2[j][1]]]
                if   j<m and C1[i][1]>C2[j+1][0]: j+=1 
                elif i<n and C2[j][1]>C1[i+1][0]: i+=1           
                else: i,j = i+1,j+1
            else: #overlap but not close enough
                if   j<m and C1[i][1]>C2[j+1][0]: j+=1
                elif i<n and C2[j][1]>C1[i+1][0]: i+=1            
                else: i,j = i+1,j+1
            if i>=n: i,j = n,j+1
            if j>=m: j,i = m,i+1
        x = x/<double>n
        y = y/<double>m
    return {'n:m':x,'m:n':y}

@cython.boundscheck(False)
@cython.nonecheck(False) 
def metric_score(list C1,list C2,double r,bint self_merge=False,
                 bint d1_ids=False,bint d2_ids=False):
    cdef long i,j,n,m
    cdef double a,b,c,w
    cdef dict A,B
    cdef list X,W
    if self_merge:
        C1 = merge_regions(C1)
        C2 = merge_regions(C2)
    n,m = len(C1),len(C2)
    A,B,X,W,d1_i,d2_i = {},{},[None for i in range(n)],[0.0,0.0,0.0,0.0],[],[]
    a,b = 0.0,0.0
    if n > 0 and m > 0:
        for i in range(n):
            for j in range(m):
                c = overlap(<long>C1[i][0],<long>C1[i][1],<long>C2[j][0],<long>C2[j][1])
                if c >= r:
                    if not A.has_key(i) or A[i]<c:
                        A[i] = c
                        X[i] = [<long>C1[i][0]-<long>C2[j][0],<long>C1[i][1]-<long>C2[j][1]]
                    if not B.has_key(j) or B[j]<c:
                        B[j] = c 
        a = 1.0*len(A)/<double>n
        b = 1.0*len(B)/<double>m
        W = feature_magnitudes(C1,C2,True)       
        if d1_ids: d1_i = list(set(range(n)).difference(set(A.keys())))
        if d2_ids: d2_i = list(set(range(m)).difference(set(B.keys())))
    return {'n':n,'m':m,'n:m':a,'m:n':b,'j':W,'b':X,'d1_ids':d1_i,'d2_ids':d2_i}
    
@cython.boundscheck(False)
@cython.nonecheck(False)    
#remove entries that intersect with regions in R
#assume both are coordinate sorted?
def filter_regions(list C, list R, double r=0.0, check_x=True, check_y=False):
    cdef list A = []
    cdef int i,j,n,m
    n = len(C)
    m = len(R)
    for i in range(n):
        b = False
        for j in range(m):
            if overlap(<long>C[i][0],<long>C[i][0],<long>R[j][0],<long>R[j][1],check_x,check_y) > r:
                b = True
        if not b: A += [C[i]]
    return A

#fast if any element of C2 touches and element of C1
@cython.boundscheck(False)
@cython.nonecheck(False)
def filter_regions2(list C,list R):
    cdef long a,b,c,d,i,j,n,m,upper
    cdef list D = [-1]
    i,j,upper,n,m = 0,0,0,len(C)+1,len(R)+1 #boundries here
    if   n>1 and m<=1: upper = C[-1][1]
    elif m>1 and n<=1: upper = R[-1][1]
    elif n>1 and m>1:  upper = max(C[-1][1],R[-1][1])
    C += [[upper+2,upper+2,0,[],0,0,{}],[upper+4,upper+4,0,[],0,0,{}]] #pad out the end of C
    R += [[upper+2,upper+2,0,[],0,0,{}],[upper+4,upper+4,0,[],0,0,{}]] #pad out the end of R
    
    while i+j < n+m:  #pivioting dual ordinal indecies scan left to right on C1, C2        
        a = <long>C[i][0]-<long>R[j][0]
        b = <long>C[i][0]-<long>R[j][1]
        c = <long>C[i][1]-<long>R[j][0]
        d = <long>C[i][1]-<long>R[j][1]
        if  a==0 and d==0:     #[7] C[i] and R[j] are equal on x
            #print('equal to\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]
            i += 1
            j += 1 
        elif  c<0:               #[1] C[i] disjoint of left of R[j]
            #print('disjoint left\ti=%s\tj=%s'%(i,j))
            i += 1    
        elif  b>0:               #[6] C[i] disjoint right of R[j]
            #print('disjoint right\ti=%s\tj=%s'%(i,j))        
            j += 1 
        elif  a<0 and d<0:       #[2] C[i] right overlap to R[j] left no envelopment
            #print('right overlap\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]      
            i += 1 
        elif  a>0 and d>0:       #[4] C[i] left overlap of R[j] right no envelopment
            #print('left overlap\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]
            j += 1
        elif  a<=0 and d>=0:     #[3] C[i] envelopment of R[j]
            #print('envelopment\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]
            j += 1
        elif  a>=0 and d<=0:     #[5] C[i] enveloped by R[j]
            #print('being eveloped\ti=%s\tj=%s'%(i,j))
            if D[-1]!=i: D+=[i]
            i += 1 
        if i>=n: i,j = n,j+1 #sticky indecies wait for eachother
        if j>=m: j,i = m,i+1 #sticky indecies wait for eachother
    while len(C) > 0 and C[-1][0] > upper:  C.pop()    
    while len(R) > 0 and R[-1][0] > upper:  R.pop()
    return [C[i] for i in sorted(set(range(len(C))).difference(set(D[1:])))]

@cython.boundscheck(False)
@cython.nonecheck(False) 
def filter_by_mag(list C, long lower, long upper, bint self_merge=False):
    cdef unsigned int i,j,n
    cdef list D = []
    if self_merge: C = merge_regions(C)
    for i in range(len(C)):
        n = get_xy_mag(C[i])
        if n >= <unsigned int>lower and n < <unsigned int>upper:
            D += [C[i]]
    return D 

#given a svult C1 add calls from C2 that do not
#overlap the calls in C1
@cython.boundscheck(False)
@cython.nonecheck(False) 
def append_bin(list C1,list C2):
    cdef long i,j,n,m
    cdef double c
    cdef dict A 
    cdef list C3,S
    A,C3,S = {},[],[]
    n,m = len(C1),len(C2)
    if m > 0:
        for i in range(n):
            for j in range(m):
                c = overlap(<long>C1[i][0],<long>C1[i][1],<long>C2[j][0],<long>C2[j][1])
                if c > 0.0: A[j] = i
        S = sorted(list(set(range(m)).difference(set(A.keys()))))#take only non overlapping
        for i in range(n): C3 += [C1[i]]
        for j in S:        C3 += [C2[j]]
    else:
        C3 = C1
    return C3    

#start with highest accuracy bin first
#and then attempt to add each
#will redo this fuse for use with new mterics: n&m = (w*n:m+m:n/w)2, j = jaccard
#TO DO update for new metrics and keys: n&m, j
@cython.boundscheck(False)
@cython.nonecheck(False) 
def bin_fuse(dict N, dict A, str opt='n:m',bint priority=False, bint merge=False):
    cdef dict F,o
    cdef list P
    cdef long i,j,t,b
    cdef tuple g
    cdef str f
    F,o = {},{'n:m':2,'m:n':3}    
    for t in A:
        F[t] = []        #initialize the calls for type t
        P = A[t].keys()  #bins of A[t] in order
        if priority:     #prioritize bins by score expectation
            if opt=='n&m':
                P = list(np.argsort([A[t][b][opt][2]+A[t][b][opt][3] for b in A[t]])[::-1])
            else:
                P = list(np.argsort([A[t][b][opt][o[opt]] for b in A[t]])[::-1])
        for b in P:
            g,f = A[t][b][opt][0],A[t][b][opt][1]
            if g != (None,):
                if priority:
                    if merge: F[t] = append_bin(F[t],merge_regions(N[len(g)][g][f][t][b]))
                    else:     F[t] = append_bin(F[t],N[len(g)][g][f][t][b])
                else:
                    if merge: F[t] += merge_regions(N[len(g)][g][f][t][b])
                    else:     F[t] += N[len(g)][g][f][t][b]
    return F

@cython.boundscheck(False)
@cython.nonecheck(False)
cdef unsigned int get_y_mag(list y):
    cdef unsigned int i,t,x
    x = 0
    for i in range(len(y)):
        if type(y[i]) is list and len(y[i])>1:
            t = <unsigned int>(abs(<long>y[i][1]-<long>y[i][0]+1))
            if t > x: x = t
    return x

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef unsigned int get_xy_mag(list c):
    cdef unsigned int i,t,x,y
    t,y,x = 0,0,0
    x = <unsigned int>(abs(<long>c[1]-<long>c[0]+1))
#    try:
#        x = <unsigned int>(abs(<long>c[1]-<long>c[0]+1))
#    except Exception as E:
#        print('type is %s'%type(c))
#        if type(c) is list:
#            print('len is %s'%len(c))
#            if len(c)>2:
#                print('type c[0] is %s'%type(c[0]))
#                print('type c[1] is %s'%type(c[1]))
#        print(c)
    for i in range(len(c[3])):
        t = <unsigned int>(abs(<long>c[3][i][1]-<long>c[3][i][0]+1))
        if t > y: y = t
    return max(x,y)

@cython.boundscheck(False)
@cython.nonecheck(False)    
cdef unsigned int get_mag(list C):
    cdef int i,n
    cdef list L = []
    n = len(C)
    for i in range(n):
        #should be the average, sum ?
        L += [<unsigned int>max(abs(<long>C[i][1]-<long>C[i][0])+1,<long>get_y_mag(C[i][3]))]
    return L   

@cython.boundscheck(False)
@cython.nonecheck(False) 
def partition_by_mag(list C, list B, bint self_merge=False, bint interpolate=False):
    cdef unsigned int i,j,n,b
    cdef double v,w_a,w_b
    cdef list M
    b = len(B)
    M,P = [],{i:[] for i in range(b-1)}
    if self_merge: C = merge_regions(C)
    if interpolate:
        M = [<double>(B[i]+B[i+1])/2.0 for i in range(b-1)] #bin midpoints
        for i in range(len(C)):
            v = <double>get_xy_mag(C[i])
            if   v < M[0]: #check first and last
                P[0]   += [[C[i][0],C[i][1],C[i][2],C[i][3],1.0*C[i][4],1.0*C[i][5],C[i][6]]]
            elif v >= M[b-2]:
                P[b-1] += [[C[i][0],C[i][1],C[i][2],C[i][3],1.0*C[i][4],1.0*C[i][5],C[i][6]]]
            else:
                for j in range(b-2):
                    if v >= M[j] and v < M[j+1]: #bin center nearest neighbors
                        w_a = 1.0-(v-M[j])/(M[j+1]-M[j]) #linear interpolation
                        w_b = 1.0 - w_a                  #and complement
                        P[j]  += [[C[i][0],C[i][1],C[i][2],C[i][3],w_a*C[i][4],w_a*C[i][5],C[i][6]]]
                        P[j+1]+= [[C[i][0],C[i][1],C[i][2],C[i][3],w_b*C[i][4],w_b*C[i][5],C[i][6]]]
    else:
        for i in range(len(C)):
            n = get_xy_mag(C[i])
            for j in range(b-1):
                if n >= B[j] and n < B[j+1]:
                    P[j] += [C[i]]
    return P
    
@cython.boundscheck(False)
#alternate overlap calculation
def region2svult(list R, int t):
    cdef dict S
    cdef int i
    S = {t:[]}
    for i in range(len(R)):
        S[t] += [[R[i][0],R[i][1],t,[[R[i][0],R[i][1]]],1.0,1.0,{-1:set([i])}]]
    return S

#breakpoint code-----------------------------------------------------------------------
#pull out all original svus from V using the idx map in S[t][i]=>V[s_id][j]
@cython.boundscheck(False)
def idx_svus_from_svult(dict S, dict V, int t, int i):
    cdef dict sidx
    cdef int s_id,x
    sidx = {} #pull out all original svus
    for s_id in S[t][i][6]: #hardcoded 6=>{idx} = {s_id1:set([1,2,3]),sid2:set([1,2,3])
        sidx[s_id] = [V[s_id][x] for x in S[t][i][6][s_id]]
    return sidx

#have a summary chain: {s_id:svu}
@cython.boundscheck(False)
def conf_idx_svu(dict sidx):
    cdef dict L
    cdef int s_id,x
    L = {}
    for s_id in sidx:
        L[s_id] = [x.conf for x in sidx[s_id]]
    return L

#input is a s_d:idx_svu map
#output is a flattened left(x1) and right(x2) breakpoint list
@cython.boundscheck(False)
def xs_conf_idx_svu(dict sidx, bint unique=True):
    cdef dict L
    cdef int s_id,i
    L = {0:[],1:[]}
    for s_id in sidx:
        for i in range(len(sidx[s_id])):
            L[0] += sidx[s_id][i].conf[0:2]
            L[1] += sidx[s_id][i].conf[2:4]
    if unique:
        L[0] = list(set(L[0]))
        L[1] = list(set(L[1]))
    return L

#using the s_id mapped to the original entries in V
#get only the neirest left and right neighbors per s_id
#TO DO implement break point weighting by caller id which will drive expansion
#collapse operations downstream of the fusion segment
@cython.boundscheck(False)
def xs_knn_conf_idx_svu(dict S,int t,int i,dict sidx,int k=1,s_id_bp=None):
    cdef dict X
    cdef set knnl,knnr
    cdef list Knnl,L,Knnr,R
    cdef unsigned int l,r
    cdef int j,xs
    if s_id_bp is None: #this is a s_id priority break point weighting
        s_id_bp = {sid:1.0 for sid in sidx}
    X = {0:[],1:[]}
    for s_id in sidx:
        knnl,L,knnr,R = set([]),[],set([]),[]
        for xs in range(len(sidx[s_id])):
            for j in range(len(sidx[s_id][xs].svu)):
                knnl.add(sidx[s_id][xs].svu[j][0]) #grab unique
                knnr.add(sidx[s_id][xs].svu[j][1]) #entries of
        Knnl,Knnr = list(knnl),list(knnr)
        for l in knnl:
            if l <= S[t][i][0]: L += [l]
        for r in knnr:
            if r >= S[t][i][1]: R += [r]
        X[0] += sorted(L,key=lambda x: abs(S[t][i][0]-x))[0:min(len(L),k)]
        X[1] += sorted(R,key=lambda x: abs(S[t][i][1]-x))[0:min(len(R),k)]   
    return X

#input is a weighted left and right breakpoint list
#output is min,median,max,n,mean,std
@cython.boundscheck(False)
def stats_idx_svu(dict L, list A):
    cdef dict S
    cdef int i
    S = {0:[],1:[]} 
    for i in L:
        if len(L[i])>0:
            S[i] = {'min':np.min(L[i]),'med':int(round(np.median(L[i]),0)),'max':np.max(L[i]),
                    'outer':[np.min(L[i]) if i==0 else np.max(L[i])][0],
                    'num':len(L[i]),'mean':int(round(np.mean(L[i]),0)),'std':int(round(np.std(L[i]),0))}
        else:
            S[i] = {'min':A[i],'med':A[i],'max':A[i],'outer':A[i],
                    'num':0,'mean':A[i],'std':0}
    return S

#for each svu in the t:svul get the
#min,median,max,n,mean,std of l,r breakpoints
@cython.boundscheck(False)
def bp_stats(dict S, dict V,int t,bint unique=True):
    cdef dict sidx,lrbps,stats
    cdef list T
    cdef int i,n
    n,T = len(S[t]),[]
    for i in range(n):
        sidxs = idx_svus_from_svult(S,V,t,i)
        lrbps = xs_conf_idx_svu(sidxs,unique)
        stats = stats_idx_svu(lrbps)
        T += [stats]
    return T

#for each svu in the t:svul get the
#unique knn or knn left and right breakpoints
@cython.boundscheck(False)
def bp_lr_knn_stats(dict S,dict V,int t,int k=1):
    cdef dict sidx,lrbps,stats
    cdef list T
    cdef int i,n
    n,T = len(S[t]),[]
    for i in range(n):
        sidxs = idx_svus_from_svult(S,V,t,i)
        lrbps = xs_knn_conf_idx_svu(S,t,i,sidxs,k)
        stats = stats_idx_svu(lrbps,S[t][i]) #error here with original
        T += [stats]
    return T

#TO DO, look at clusters here to find alternate fusion segments
@cython.boundscheck(False)
def bp_kmeans(dict S,dict V,int t,int k,bint unique=True):
    cdef dict sidx,lrbps
    cdef list T
    cdef int i,n
    n,T = len(S[t]),[]
    for i in range(n):
        sidxs = idx_svus_from_svult(S,V,t,i)
        lrbps = xs_conf_idx_svu(sidxs,unique)
    return T

@cython.boundscheck(False)
def bp_analysis(dict S,dict V,int k=3,bint unique=True,stat='outer'):
    cdef dict T
    cdef list s
    cdef int i,w,t
    T = copy.deepcopy(S)
    for w in S:
        for t in S[w]:
            #do the stats
            #do the bp_lr_knn_stats
            s = bp_lr_knn_stats(S[w],V,t,k)
            for i in range(len(s)):
                T[w][t][i][0:2] = [s[i][0][stat],s[i][1][stat]]
            #do the bp_kmeans
    return T
    
#given a svul, join idxs from list j
@cython.boundscheck(False)
@cython.nonecheck(False)
################################################################################
#given a svul, join idxs from list j
cdef dict join_idx(list C,list j,long p):
    cdef long i,k,c
    cdef dict A,idx
    A = {}
    for idx in [C[i][p] for i in j]: #hard coded idx=6
        for k in idx:
            for c in idx[k]:
                if A.has_key(k): A[k].add(c)
                else:            A[k] =  {c}        
    return A

@cython.boundscheck(False)
@cython.nonecheck(False)
#given two idx merge them into one
cdef dict merge_idx(dict idx1, dict idx2):
    cdef long k,c
    cdef dict A
    A = {}
    for k in idx1:
        for c in idx1[k]:
            if A.has_key(k): A[k].add(c)
            else:            A[k] =  {c}
    for k in idx2:
        for c in idx2[k]:
            if A.has_key(k): A[k].add(c)
            else:            A[k] =  {c}
    return A

@cython.boundscheck(False)
@cython.nonecheck(False)
#given a svul, join the y ranges from list j
#if they overlap, otherwise leave them alone
cdef list join_y(list C,list x,long p):
    cdef long i,j,n,b
    cdef list A,B,F,R,k
    A,B,F,R = [],[],[],[]
    for i in x: #partition forward and reverse
        for k in C[i][p]: #hard coded y=3
            if k[0]>k[1]: R += [[k[1],k[0]]] #------>REV
            else:         F += [[k[0],k[1]]]
    #forward partitions
    F = sorted(F)
    n = len(F)
    if n > 0:
        i = 0
        while i < n-1:
            j = i+1
            b = F[i][1]
            while b+1 >= F[j][0] and j < n-1:
               if b < F[j][1]: b = F[j][1]
               j += 1
            A += [[F[i][0],b]]
            i = j
        if len(A) > 0 and A[-1][1]+1 >= F[i][0]:
            if A[-1][1] < F[i][1]:
                A[-1][1] = F[i][1]
        else:
            A += [[F[i][0],F[i][1]]]
    #reverse partitions
    R = sorted(R)
    n = len(R)
    if n > 0:
        i = 0
        while i < n-1:
            j = i+1
            b = R[i][1]
            while b+1 >= R[j][0] and j < n-1:
               if b < R[j][1]: b = R[j][1]
               j += 1
            B += [[R[i][0],b]]
            i = j
        if len(B) > 0 and B[-1][1]+1 >= R[i][0]:
            if B[-1][1] < R[i][1]:
                B[-1][1] = R[i][1]
        else:
            B += [[R[i][0],R[i][1]]]
    return A+[[B[i][1],B[i][0]] for i in range(len(B))]       #----->REV

@cython.boundscheck(False)
@cython.nonecheck(False)    
#if the y values fron Y1 and Y2 overlap, merge them together
#otherwise add them as unique entries and return list of lists
cdef list merge_y(list Y1,list Y2):
    cdef long i,j,n,b
    cdef list A,B,F,R
    A,B,F,R = [],[],[],[]
    for i in range(len(Y1)): #partition forward and reverse
        if Y1[i][0]>Y1[i][1]: R += [[Y1[i][1],Y1[i][0]]] #----->REV
        else:                 F += [[Y1[i][0],Y1[i][1]]] 
    for i in range(len(Y2)):
        if Y2[i][0]>Y2[i][1]: R += [[Y2[i][1],Y2[i][0]]] #----->REV
        else:                 F += [[Y2[i][0],Y2[i][1]]] 
    #forward partition    
    F = sorted(F)
    n = len(F)
    if n > 0:
        i = 0
        while i < n-1:
            j = i+1
            b = F[i][1]
            while b+1 >= F[j][0] and j < n-1:
               if b < F[j][1]: b = F[j][1]
               j +=1
            A += [[F[i][0],b]]
            i = j
        if len(A)>0 and A[-1][1]+1>=F[i][0]:
            if A[-1][1]<F[i][1]:
                A[-1][1] = F[i][1]
        else:
            A += [[F[i][0],F[i][1]]]
    #reverse partition    
    R = sorted(R)
    n = len(R)
    if n > 0:
        i = 0
        while i < n-1:
            j = i+1
            b = R[i][1]
            while b+1 >= R[j][0] and j < n-1:
               if b < R[j][1]: b = R[j][1]
               j +=1
            B += [[R[i][0],b]]
            i = j
        if len(B)>0 and B[-1][1]+1>=R[i][0]:
            if B[-1][1]<R[i][1]:
                B[-1][1] = R[i][1]
        else:
            B += [[R[i][0],R[i][1]]]
    return A+[[B[i][1],B[i][0]] for i in range(len(B))]       #----->REV

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
#join via method:{min,max,mu,f1}
def join_wx(list C,list x,long p):
    cdef long i    
    cdef double w
    cdef list M
    if len(x)>1:
        M = [np.abs(<double>C[i][1]-<double>C[i][0])+1.0 for i in x]
        w = np.sum([M[i]*C[i][p] for i in range(len(M))])/np.sum(M)
    else:
        w = C[x[0]][p]
    return w

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
#join via method:{min,max,mu,f1}
def join_wy(list C,list x,long p):
    cdef long i    
    cdef double w
    cdef list M
    if len(x)>1:
        M = [np.abs(<double>C[i][1]-<double>C[i][0])+1.0 for i in x]
        w = np.sum([M[i]*C[i][p] for i in range(len(M))])/np.sum(M)
    else:
        w = C[x[0]][p]
    return w
    
@cython.boundscheck(False)
@cython.nonecheck(False)    
#adpated for wieghted svul with idx
#need to update for y lists y1[3]-y2[4]
#C[i][4] is the wx which is averaged during a merge step
#C[i][5] is the wy which is averaged during a merge step
def merge_regions(list C,check_x=True,check_y=True):
    cdef long i,j,b,n
    cdef list M = []
    n = len(C)
    if n > 0:
        i,W = 0,[]
        while i < n-1:
            j = i+1
            b = C[i][1]
            while b+1 >= C[j][0] and j < n-1:
                if b < C[j][1]: b = C[j][1]
                j += 1
            M += [[C[i][0],b,C[i][2],
                   join_y(C,range(i,j),3),    #[3] will be y=[[y1,y2]]
                   join_wx(C,range(i,j),4),   #[4] will be wx
                   join_wy(C,range(i,j),5),   #[5] will be wy
                   join_idx(C,range(i,j),6)]] #[6] will be {idx}
            i = j                              #----------span of i to j here-------------
        if len(M)>0 and M[-1][1]+1>=C[i][0]:   #----------potential span of i-1 to i here-
            if M[-1][1]<C[i][1]:
                M[-1][1] = C[i][1]
                M[-1][3] = join_y(C,[i-1,i],3)
                M[-1][4] = join_wx(C,[i-1,i],4)
                M[-1][5] = join_wy(C,[i-1,i],5)
                M[-1][6] = join_idx(C,[i-1,i],6)
        else:                                  #------------only i is here----------------
            M += [[C[i][0],C[i][1],C[i][2],
                   join_y(C,[i],3),
                   join_wx(C,[i],4),
                   join_wy(C,[i],5),
                   join_idx(C,[i],6)]]
    return M

@cython.boundscheck(False)
@cython.nonecheck(False)
#[1] c1 disjoint of left of c2
cdef void orientation_1(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i] #no intersection
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]:
        U[-1][1] = C1[j1][1]
        U[-1][3] = merge_y(U[-1][3],C1[j1][3])
        U[-1][6] = merge_idx(U[-1][6],C1[j1][6])
    else:                         
        U += [[C1[j1][0],C1[j1][1],C1[j1][2],C1[j1][3],
               C1[j1][4],C1[j1][5],C1[j1][6]]]
    #[d1]--------------------[d1]
    if D1[-1][1]+1 >= C1[j1][0]:  #extend segment
        if D1[-1][1]+1!=C2[j2-1][0]:
            D1[-1][1] = C1[j1][1]
            D1[-1][3] = merge_y(D1[-1][3],C1[j1][3])
            D1[-1][6] = merge_idx(D1[-1][6],C1[j1][6])
    else:                         #new segment                        
        D1 += [[C1[j1][0],C1[j1][1],C1[j1][2],C1[j1][3],
               C1[j1][4],C1[j1][5],C1[j1][6]]]            
    #[d2]--------------------[d2] #no set two difference

@cython.boundscheck(False)
@cython.nonecheck(False)
#[2] c1 right overlap to c2 left no envelopment
cdef void orientation_2(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C2[j2][0]: 
        I[-1][1] = C1[j1][1] #was C2[j2][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3]))
        I[-1][6] = merge_idx(I[-1][6],merge_idx(C1[j1][6],C2[j2][6]))
    else:                       
        I += [[C2[j2][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],
               merge_idx(C1[j1][6],C2[j2][6])]]                 
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]: 
        U[-1][1] = C2[j2][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
        U[-1][6] =  merge_idx(U[-1][6],merge_idx(C1[j1][6],C2[j2][6]))
    else:                       
        U += [[C1[j1][0],C2[j2][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],
               merge_idx(C1[j1][6],C2[j2][6])]]
    #[d1]--------------------[d1]
    if D1[-1][1]+1 >= C1[j1][0]: 
        D1[-1][1] = C1[j1][1]
        D1[-1][3] = merge_y(D1[-1][3],C1[j1][3])
        D1[-1][6] = merge_idx(D1[-1][6],C1[j1][6])
        if D1[-1][1] > C2[j2][0]-1:
            D1[-1][1] = C2[j2][0]-1
            if D1[-1][1] < D1[-1][0]: D1.pop()
    else:
        D1 += [[C1[j1][0],C2[j2][0]-1,C1[j1][2],C1[j1][3],
                C1[j1][4],C1[j1][5],C1[j1][6]]]            
    #[d2]--------------------[d2] 
    D2 += [[C1[j1][1]+1,C2[j2][1],C2[j2][2],C2[j2][3],
            C2[j2][4],C2[j2][5],C2[j2][6]]]

@cython.boundscheck(False)
@cython.nonecheck(False)
#[3] c1 envelopment of c2
cdef void orientation_3(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C2[j2][0]: 
        I[-1][1] = C2[j2][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3]))
        I[-1][6] = merge_idx(I[-1][6],merge_idx(C1[j1][6],C2[j2][6]))            
    else:                       
        I += [[C2[j2][0],C2[j2][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],
               merge_idx(C1[j1][6],C2[j2][6])]]                        
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]: 
        U[-1][1] = C1[j1][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
        U[-1][6] = merge_idx(U[-1][6],merge_idx(C1[j1][6],C2[j2][6]))
    else:                       
        U += [[C1[j1][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],
               merge_idx(C1[j1][6],C2[j2][6])]]
    #[d1]--------------------[d1] 
    if D1[-1][1]+1 >= C1[j1][0]:
        D1[-1][1] = C2[j2][0]-1
        if C2[j2][1] < C1[j1][1]:
            D1 += [[C2[j2][1]+1,C1[j1][1],C1[j1][2],C1[j1][3],
                    C1[j1][4],C1[j1][5],C1[j1][6]]]
    elif D1[-1][1] >= C2[j2][0]:
        D1[-1][1] = C2[j2][0]-1
        if C2[j2][1] < C1[j1][1]:  #has a right side
            D1 += [[C2[j2][1]+1,C1[j1][1],C1[j1][2],C1[j1][3],
                    C1[j1][4],C1[j1][5],C1[j1][6]]]
    else:
        if C1[j1][0] < C2[j2][0]:  #has a left side
            D1 += [[C1[j1][0],C2[j2][0]-1,C1[j1][2],C1[j1][3],
                    C1[j1][4],C1[j1][5],C1[j1][6]]]
        if C2[j2][1] < C1[j1][1]:  #has a right side
            D1 += [[C2[j2][1]+1,C1[j1][1],C1[j1][2],C1[j1][3],
                    C1[j1][4],C1[j1][5],C1[j1][6]]]

@cython.boundscheck(False)
@cython.nonecheck(False)
#[4] c1 left overlap of c2 right no envelopment            
cdef void orientation_4(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C1[j1][0]: 
        I[-1][1] = C2[j2][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3]))
        I[-1][6] = merge_idx(I[-1][6],merge_idx(C1[j1][6],C2[j2][6]))          
    else:                       
        I += [[C1[j1][0],C2[j2][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],
               merge_idx(C1[j1][6],C2[j2][6])]]           
    #[u]----------------------[u]
    if U[-1][1]+1 >= C2[j2][0]: 
        U[-1][1] = C1[j1][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
        U[-1][6] = merge_idx(U[-1][6],merge_idx(C1[j1][6],C2[j2][6]))
    else:                       
        U += [[C2[j2][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],
               merge_idx(C1[j1][6],C2[j2][6])]]       
    #[d1]--------------------[d1]
    D1 += [[C2[j2][1]+1,C1[j1][1],C1[j1][2],C1[j1][3],
            C1[j1][4],C1[j1][5],C1[j1][6]]]            
    #[d2]--------------------[d2] 
    if D2[-1][1]+1 >= C2[j2][0]: 
        D2[-1][1] = C2[j2][1]
        D2[-1][3] = merge_y(D2[-1][3],C2[j2][3])
        D2[-1][6] = merge_idx(D2[-1][6],C2[j2][6])
        if D2[-1][1] > C1[j1][0]-1:
            D2[-1][1] = C1[j1][0]-1
            if D2[-1][1] < D2[-1][0]: D2.pop()
    else:                        
        D2 += [[C2[j2][0],C1[j1][0]-1,C2[j2][2],C2[j2][3],
                C2[j2][4],C2[j2][5],C2[j2][6]]]

@cython.boundscheck(False)
@cython.nonecheck(False)
#[5] c1 enveloped by c2
cdef void orientation_5(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C1[j1][0]: 
        I[-1][1] = C1[j1][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3])) 
        I[-1][6] = merge_idx(I[-1][6],merge_idx(C1[j1][6],C2[j2][6]))            
    else:                       
        I += [[C1[j1][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],
               merge_idx(C1[j1][6],C2[j2][6])]]
                                      
    #[u]----------------------[u]
    if U[-1][1]+1 >= C2[j2][0]: 
        U[-1][1] = C2[j2][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
        U[-1][6] = merge_idx(U[-1][6],merge_idx(C1[j1][6],C2[j2][6]))
    else:                       
        U += [[C2[j2][0],C2[j2][1],C2[j2][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C2[j2][4],C2[j2][5],
               merge_idx(C1[j1][6],C2[j2][6])]]
    #[d1]--------------------[d1] #no set one difference
    #[d2]--------------------[d2]
    if D2[-1][1]+1 >= C2[j2][0]: 
        D2[-1][1] = C1[j1][0]-1
        if C1[j1][1] < C2[j2][1]:
            D2 += [[C1[j1][1]+1,C2[j2][1],C2[j2][2],C2[j2][3],
                    C2[j2][4],C2[j2][5],C2[j2][6]]]
    elif D2[-1][1] >= C1[j1][0]:
        D2[-1][1] = C1[j1][0]-1
        if C1[j1][1] < C2[j2][1]:  #has a right side
            D2 += [[C1[j1][1]+1,C2[j2][1],C2[j2][2],C2[j2][3],
                    C2[j2][4],C2[j2][5],C2[j2][6]]]
    else:
        if C2[j2][0] < C1[j1][0]:  #has a left side
            D2 += [[C2[j2][0],C1[j1][0]-1,C2[j2][2],C2[j2][3],
                    C2[j2][4],C2[j2][5],C2[j2][6]]]
        if C1[j1][1] < C2[j2][1]:  #has a right side
            D2 += [[C1[j1][1]+1,C2[j2][1],C2[j2][2],C2[j2][3],
                    C2[j2][4],C2[j2][5],C2[j2][6]]]

@cython.boundscheck(False)
@cython.nonecheck(False)
#[6] c1 disjoint right of c2
cdef void orientation_6(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i] #no instersection
    if U[-1][1]+1 >= C2[j2][0]:   
        U[-1][1] = C2[j2][1]
        U[-1][3] = merge_y(U[-1][3],C2[j2][3])
        U[-1][6] = merge_idx(U[-1][6],C2[j2][6])
    else:                         
        U += [[C2[j2][0],C2[j2][1],C2[j2][2],C2[j2][3],
               C2[j2][4],C2[j2][5],C2[j2][6]]]
    #[d1]--------------------[d1] #no set one difference
    #[d2]--------------------[d2] 
    if D2[-1][1]+1 >= C2[j2][0]:
        if D2[-1][1]+1!=C1[j1-1][0]:
            D2[-1][1] = C2[j2][1]
            D2[-1][3] = merge_y(D2[-1][3],C2[j2][3])
            D2[-1][6] = merge_idx(D2[-1][6],C2[j2][6])
    else:                        
        D2 += [[C2[j2][0],C2[j2][1],C2[j2][2],C2[j2][3],
                C2[j2][4],C2[j2][5],C2[j2][6]]] 

@cython.boundscheck(False)
@cython.nonecheck(False)
#[7] c1 and c2 are equal on x
cdef void orientation_7(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C1[j1][0]:
        I[-1][1] = C1[j1][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3]))
        I[-1][6] = merge_idx(I[-1][6],merge_idx(C1[j1][6],C2[j2][6]))
    else:
        I += [[C1[j1][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],
               merge_idx(C1[j1][6],C2[j2][6])]]
               
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]:
        U[-1][1] = C1[j1][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
        U[-1][6] = merge_idx(U[-1][6],merge_idx(C1[j1][6],C2[j2][6]))
    else:
        U += [[C1[j1][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],
               merge_idx(C1[j1][6],C2[j2][6])]]
    #[d1]----------------------[d1]
    #[d2]----------------------[d2]

@cython.boundscheck(False)
@cython.nonecheck(False)
def LR(list C1,list C2):
    cdef long a,b,c,d,j1,j2,n1,n2
    cdef list I,U,D1,D2
    j1,j2,upper = 0,0,0         #initializations and padding
    n1,n2 = len(C1)+1,len(C2)+1 #boundries here
    I,U,  = [[-2,-2,0,[],0,0,{}]],[[-2,-2,0,[],0,0,{}]]
    D1,D2 = [[-2,-2,0,[],0,0,{}]],[[-2,-2,0,[],0,0,{}]]
    if n1 > 1 and n2 > 1:
        upper = max(C1[-1][1],C2[-1][1])
        C1 += [[upper+2,upper+2,0,[],0,0,{}],[upper+4,upper+4,0,[],0,0,{}]] #pad out the end of C1
        C2 += [[upper+2,upper+2,0,[],0,0,{}],[upper+4,upper+4,0,[],0,0,{}]] #pad out the end of C2
        while j1+j2 < n1+n2:  #pivioting dual ordinal indecies scan left to right on C1, C2
            a = <long>C1[j1][0]-<long>C2[j2][0]
            b = <long>C1[j1][0]-<long>C2[j2][1]
            c = <long>C1[j1][1]-<long>C2[j2][0]
            d = <long>C1[j1][1]-<long>C2[j2][1]
            if    C1[j1][0:2]==C2[j2][0:2]:    #[7] c1 and c2 are equal on x
                orientation_7(C1,j1,C2,j2,I,U,D1,D2)
                j1 += 1
                j2 += 1 
            elif  c<0:               #[1] c1 disjoint of left of c2
                orientation_1(C1,j1,C2,j2,I,U,D1,D2)                   
                j1 += 1    
            elif  b>0:               #[6] c1 disjoint right of c2
                orientation_6(C1,j1,C2,j2,I,U,D1,D2)             
                j2 += 1 
            elif  a<0 and d<0:       #[2] c1 right overlap to c2 left no envelopment
                orientation_2(C1,j1,C2,j2,I,U,D1,D2)           
                j1 += 1 
            elif  a>0 and d>0:       #[4] c1 left overlap of c2 right no envelopment
                orientation_4(C1,j1,C2,j2,I,U,D1,D2) 
                j2 += 1 
            elif  a<=0 and d>=0:     #[3] c1 envelopment of c2
                orientation_3(C1,j1,C2,j2,I,U,D1,D2)
                j2 += 1 
            elif  a>=0 and d<=0:     #[5] c1 enveloped by c2
                orientation_5(C1,j1,C2,j2,I,U,D1,D2)
                j1 += 1 
            if j1>=n1: j1,j2 = n1,j2+1 #sticky indecies wait for eachother
            if j2>=n2: j2,j1 = n2,j1+1 #sticky indecies wait for eachother
        #pop off extras for each features (at most two at the end)
        while len(C1) > 0 and C1[-1][0] > upper:  C1.pop()    
        while len(C2) > 0 and C2[-1][0] > upper: C2.pop()
        while len(I)  > 0 and I[-1][0]>upper:      I.pop()
        if len(I) > 0 and I[-1][1]>upper:         I[-1][1] = upper
        while len(U) > 0 and U[-1][0]>upper:      U.pop()
        if len(U) > 0 and U[-1][1]>upper:         U[-1][1] = upper
        while len(D1) > 0 and D1[-1][0]>upper:    D1.pop()
        if len(D1) > 0 and D1[-1][1]>upper:       D1[-1][1] = min(C2[-1][0]-1,C1[-1][1])
        while len(D2) > 0 and D2[-1][0]>upper:    D2.pop()
        if len(D2) > 0 and D2[-1][1]>upper:       D2[-1][1] = min(C1[-1][0]-1,C2[-1][1])
    else:
        if   n1==1:  
            if n2>1: U,D2 = U+C2,D2+C2
        elif n2==1:
            if n1>1: U,D1 = U+C1,D1+C1 
    return I[1:],U[1:],D1[1:],D2[1:] #cut the disjoint padding at the start of the features 
#L->R Merging Scan->->->->->->->->->->->->->->->->->->->->->->->->->->->->

@cython.boundscheck(False)
@cython.nonecheck(False)
#[1] c1 disjoint of left of c2
cdef void orientation_1_no_idx(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i] #no intersection
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]:
        U[-1][1] = C1[j1][1]
        U[-1][3] = merge_y(U[-1][3],C1[j1][3])
    else:                         
        U += [[C1[j1][0],C1[j1][1],C1[j1][2],C1[j1][3],
               C1[j1][4],C1[j1][5],C1[j1][6]]]
    #[d1]--------------------[d1]
    if D1[-1][1]+1 >= C1[j1][0]:  #extend segment
        if D1[-1][1]+1!=C2[j2-1][0]:
            D1[-1][1] = C1[j1][1]
            D1[-1][3] = merge_y(D1[-1][3],C1[j1][3])
    else:                         #new segment                        
        D1 += [[C1[j1][0],C1[j1][1],C1[j1][2],C1[j1][3],
               C1[j1][4],C1[j1][5],C1[j1][6]]]            
    #[d2]--------------------[d2] #no set two difference

@cython.boundscheck(False)
@cython.nonecheck(False)
#[2] c1 right overlap to c2 left no envelopment
cdef void orientation_2_no_idx(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C2[j2][0]: 
        I[-1][1] = C1[j1][1] #was C2[j2][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3]))
    else:                       
        I += [[C2[j2][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],C1[j1][6]]]                 
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]: 
        U[-1][1] = C2[j2][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
    else:                       
        U += [[C1[j1][0],C2[j2][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],C1[j1][6]]]
    #[d1]--------------------[d1]
    if D1[-1][1]+1 >= C1[j1][0]: 
        D1[-1][1] = C1[j1][1]
        D1[-1][3] = merge_y(D1[-1][3],C1[j1][3])
        D1[-1][6] = merge_idx(D1[-1][6],C1[j1][6])
        if D1[-1][1] > C2[j2][0]-1:
            D1[-1][1] = C2[j2][0]-1
            if D1[-1][1] < D1[-1][0]: D1.pop()
    else:
        D1 += [[C1[j1][0],C2[j2][0]-1,C1[j1][2],C1[j1][3],
                C1[j1][4],C1[j1][5],C1[j1][6]]]            
    #[d2]--------------------[d2] 
    D2 += [[C1[j1][1]+1,C2[j2][1],C2[j2][2],C2[j2][3],
            C2[j2][4],C2[j2][5],C2[j2][6]]]

@cython.boundscheck(False)
@cython.nonecheck(False)
#[3] c1 envelopment of c2
cdef void orientation_3_no_idx(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C2[j2][0]: 
        I[-1][1] = C2[j2][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3]))           
    else:                       
        I += [[C2[j2][0],C2[j2][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],C1[j1][6]]]                        
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]: 
        U[-1][1] = C1[j1][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
    else:                       
        U += [[C1[j1][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],C1[j1][6]]]
    #[d1]--------------------[d1] 
    if D1[-1][1]+1 >= C1[j1][0]:
        D1[-1][1] = C2[j2][0]-1
        if C2[j2][1] < C1[j1][1]:
            D1 += [[C2[j2][1]+1,C1[j1][1],C1[j1][2],C1[j1][3],
                    C1[j1][4],C1[j1][5],C1[j1][6]]]
    elif D1[-1][1] >= C2[j2][0]:
        D1[-1][1] = C2[j2][0]-1
        if C2[j2][1] < C1[j1][1]:  #has a right side
            D1 += [[C2[j2][1]+1,C1[j1][1],C1[j1][2],C1[j1][3],
                    C1[j1][4],C1[j1][5],C1[j1][6]]]
    else:
        if C1[j1][0] < C2[j2][0]:  #has a left side
            D1 += [[C1[j1][0],C2[j2][0]-1,C1[j1][2],C1[j1][3],
                    C1[j1][4],C1[j1][5],C1[j1][6]]]
        if C2[j2][1] < C1[j1][1]:  #has a right side
            D1 += [[C2[j2][1]+1,C1[j1][1],C1[j1][2],C1[j1][3],
                    C1[j1][4],C1[j1][5],C1[j1][6]]]

@cython.boundscheck(False)
@cython.nonecheck(False)
#[4] c1 left overlap of c2 right no envelopment            
cdef void orientation_4_no_idx(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C1[j1][0]: 
        I[-1][1] = C2[j2][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3]))         
    else:                       
        I += [[C1[j1][0],C2[j2][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],C1[j1][6]]]           
    #[u]----------------------[u]
    if U[-1][1]+1 >= C2[j2][0]: 
        U[-1][1] = C1[j1][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
    else:                       
        U += [[C2[j2][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],C1[j1][6]]]       
    #[d1]--------------------[d1]
    D1 += [[C2[j2][1]+1,C1[j1][1],C1[j1][2],C1[j1][3],
            C1[j1][4],C1[j1][5],C1[j1][6]]]            
    #[d2]--------------------[d2] 
    if D2[-1][1]+1 >= C2[j2][0]: 
        D2[-1][1] = C2[j2][1]
        D2[-1][3] = merge_y(D2[-1][3],C2[j2][3])
        if D2[-1][1] > C1[j1][0]-1:
            D2[-1][1] = C1[j1][0]-1
            if D2[-1][1] < D2[-1][0]: D2.pop()
    else:                        
        D2 += [[C2[j2][0],C1[j1][0]-1,C2[j2][2],C2[j2][3],
                C2[j2][4],C2[j2][5],C2[j2][6]]]

@cython.boundscheck(False)
@cython.nonecheck(False)
#[5] c1 enveloped by c2
cdef void orientation_5_no_idx(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C1[j1][0]: 
        I[-1][1] = C1[j1][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3]))           
    else:                       
        I += [[C1[j1][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],C1[j1][6]]]
                                      
    #[u]----------------------[u]
    if U[-1][1]+1 >= C2[j2][0]: 
        U[-1][1] = C2[j2][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
    else:                       
        U += [[C2[j2][0],C2[j2][1],C2[j2][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C2[j2][4],C2[j2][5],C1[j1][6]]]
    #[d1]--------------------[d1] #no set one difference
    #[d2]--------------------[d2]
    if D2[-1][1]+1 >= C2[j2][0]: 
        D2[-1][1] = C1[j1][0]-1
        if C1[j1][1] < C2[j2][1]:
            D2 += [[C1[j1][1]+1,C2[j2][1],C2[j2][2],C2[j2][3],
                    C2[j2][4],C2[j2][5],C2[j2][6]]]
    elif D2[-1][1] >= C1[j1][0]:
        D2[-1][1] = C1[j1][0]-1
        if C1[j1][1] < C2[j2][1]:  #has a right side
            D2 += [[C1[j1][1]+1,C2[j2][1],C2[j2][2],C2[j2][3],
                    C2[j2][4],C2[j2][5],C2[j2][6]]]
    else:
        if C2[j2][0] < C1[j1][0]:  #has a left side
            D2 += [[C2[j2][0],C1[j1][0]-1,C2[j2][2],C2[j2][3],
                    C2[j2][4],C2[j2][5],C2[j2][6]]]
        if C1[j1][1] < C2[j2][1]:  #has a right side
            D2 += [[C1[j1][1]+1,C2[j2][1],C2[j2][2],C2[j2][3],
                    C2[j2][4],C2[j2][5],C2[j2][6]]]

@cython.boundscheck(False)
@cython.nonecheck(False)
#[6] c1 disjoint right of c2
cdef void orientation_6_no_idx(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i] #no instersection
    if U[-1][1]+1 >= C2[j2][0]:   
        U[-1][1] = C2[j2][1]
        U[-1][3] = merge_y(U[-1][3],C2[j2][3])
    else:                         
        U += [[C2[j2][0],C2[j2][1],C2[j2][2],C2[j2][3],
               C2[j2][4],C2[j2][5],C2[j2][6]]]
    #[d1]--------------------[d1] #no set one difference
    #[d2]--------------------[d2] 
    if D2[-1][1]+1 >= C2[j2][0]:
        if D2[-1][1]+1!=C1[j1-1][0]:
            D2[-1][1] = C2[j2][1]
            D2[-1][3] = merge_y(D2[-1][3],C2[j2][3])
    else:                        
        D2 += [[C2[j2][0],C2[j2][1],C2[j2][2],C2[j2][3],
                C2[j2][4],C2[j2][5],C2[j2][6]]] 

@cython.boundscheck(False)
@cython.nonecheck(False)
#[7] c1 and c2 are equal on x
cdef void orientation_7_no_idx(list C1,long j1,list C2,long j2,
                        list I,list U,list D1,list D2):
    #[i]----------------------[i]
    if I[-1][1]+1 >= C1[j1][0]:
        I[-1][1] = C1[j1][1]
        I[-1][3] = merge_y(I[-1][3],merge_y(C1[j1][3],C2[j2][3]))
    else:
        I += [[C1[j1][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],C1[j1][6]]]
               
    #[u]----------------------[u]
    if U[-1][1]+1 >= C1[j1][0]:
        U[-1][1] = C1[j1][1]
        U[-1][3] = merge_y(U[-1][3],merge_y(C1[j1][3],C2[j2][3]))
    else:
        U += [[C1[j1][0],C1[j1][1],C1[j1][2],
               merge_y(C1[j1][3],C2[j2][3]),
               C1[j1][4],C1[j1][5],C1[j1][6]]]
    #[d1]----------------------[d1]
    #[d2]----------------------[d2]

@cython.boundscheck(False)
@cython.nonecheck(False)
def LR_no_idx(list C1,list C2):
    cdef long a,b,c,d,j1,j2,n1,n2
    cdef list I,U,D1,D2
    j1,j2,upper = 0,0,0         #initializations and padding
    n1,n2 = len(C1)+1,len(C2)+1 #boundries here
    I,U,  = [[-2,-2,0,[],0,0,{}]],[[-2,-2,0,[],0,0,{}]]
    D1,D2 = [[-2,-2,0,[],0,0,{}]],[[-2,-2,0,[],0,0,{}]]
    if n1 > 1 and n2 > 1:
        upper = max(C1[-1][1],C2[-1][1])
        C1 += [[upper+2,upper+2,0,[],0,0,{}],[upper+4,upper+4,0,[],0,0,{}]] #pad out the end of C1
        C2 += [[upper+2,upper+2,0,[],0,0,{}],[upper+4,upper+4,0,[],0,0,{}]] #pad out the end of C2
        while j1+j2 < n1+n2:  #pivioting dual ordinal indecies scan left to right on C1, C2
            a = <long>C1[j1][0]-<long>C2[j2][0]
            b = <long>C1[j1][0]-<long>C2[j2][1]
            c = <long>C1[j1][1]-<long>C2[j2][0]
            d = <long>C1[j1][1]-<long>C2[j2][1]
            if    C1[j1][0:2]==C2[j2][0:2]:    #[7] c1 and c2 are equal on x
                orientation_7_no_idx(C1,j1,C2,j2,I,U,D1,D2)
                j1 += 1
                j2 += 1 
            elif  c<0:               #[1] c1 disjoint of left of c2
                orientation_1_no_idx(C1,j1,C2,j2,I,U,D1,D2)                   
                j1 += 1    
            elif  b>0:               #[6] c1 disjoint right of c2
                orientation_6_no_idx(C1,j1,C2,j2,I,U,D1,D2)             
                j2 += 1 
            elif  a<0 and d<0:       #[2] c1 right overlap to c2 left no envelopment
                orientation_2_no_idx(C1,j1,C2,j2,I,U,D1,D2)           
                j1 += 1 
            elif  a>0 and d>0:       #[4] c1 left overlap of c2 right no envelopment
                orientation_4_no_idx(C1,j1,C2,j2,I,U,D1,D2) 
                j2 += 1 
            elif  a<=0 and d>=0:     #[3] c1 envelopment of c2
                orientation_3_no_idx(C1,j1,C2,j2,I,U,D1,D2)
                j2 += 1 
            elif  a>=0 and d<=0:     #[5] c1 enveloped by c2
                orientation_5_no_idx(C1,j1,C2,j2,I,U,D1,D2)
                j1 += 1 
            if j1>=n1: j1,j2 = n1,j2+1 #sticky indecies wait for eachother
            if j2>=n2: j2,j1 = n2,j1+1 #sticky indecies wait for eachother
        #pop off extras for each features (at most two at the end)
        while len(C1) > 0 and C1[-1][0] > upper:  C1.pop()    
        while len(C2) > 0 and C2[-1][0] > upper: C2.pop()
        while len(I)  > 0 and I[-1][0]>upper:      I.pop()
        if len(I) > 0 and I[-1][1]>upper:         I[-1][1] = upper
        while len(U) > 0 and U[-1][0]>upper:      U.pop()
        if len(U) > 0 and U[-1][1]>upper:         U[-1][1] = upper
        while len(D1) > 0 and D1[-1][0]>upper:    D1.pop()
        if len(D1) > 0 and D1[-1][1]>upper:       D1[-1][1] = min(C2[-1][0]-1,C1[-1][1])
        while len(D2) > 0 and D2[-1][0]>upper:    D2.pop()
        if len(D2) > 0 and D2[-1][1]>upper:       D2[-1][1] = min(C1[-1][0]-1,C2[-1][1])
    else:
        if   n1==1:  
            if n2>1: U,D2 = U+C2,D2+C2
        elif n2==1:
            if n1>1: U,D1 = U+C1,D1+C1 
    return I[1:],U[1:],D1[1:],D2[1:] #cut the disjoint padding at the start of the features 

#each feature is a svult
def features(dict C1,dict C2):
    cdef dict I,U,D1,D2
    cdef list S1,S2
    cdef int t
    cdef set types = set(C1.keys()).union(set(C2.keys())) #all types of both sets
    I,U,D1,D2 = {},{},{},{}
    for t in types: #y is the type with value t
        if t < 6:
            S1,S2 = [],[]
            if C1.has_key(t):    S1 = merge_regions(C1[t])
            if C2.has_key(t):    S2 = merge_regions(C2[t])            
            I[t],U[t],D1[t],D2[t] = LR(S1,S2)
        else:
            I[t],U[t],D1[t],D2[t]= [],[],[],[]
    return {'I':I,'U':U,'D1':D1,'D2':D2}  

def features_no_idx(dict C1,dict C2):
    cdef dict I,U,D1,D2
    cdef list S1,S2
    cdef int t
    cdef set types = set(C1.keys()).union(set(C2.keys())) #all types of both sets
    I,U,D1,D2 = {},{},{},{}
    for t in types: #y is the type with value t
        if t < 6:
            S1,S2 = [],[]
            if C1.has_key(t):    S1 = merge_regions(C1[t])
            if C2.has_key(t):    S2 = merge_regions(C2[t])            
            I[t],U[t],D1[t],D2[t] = LR_no_idx(S1,S2)
        else:
            I[t],U[t],D1[t],D2[t]= [],[],[],[]
    return {'I':I,'U':U,'D1':D1,'D2':D2}

#union of interections by number    
def cascading_pileup(dict M, list k):
    cdef dict U
    cdef int i
    if len(M) > 2 and len(k) > 2:
        U = M[k[0]]['I']
        for i in range(1,len(k)):
            U = features(U,M[k[i]]['I'])['U']
    elif len(M) > 1 and len(k) > 1:
        U  = features(M[k[0]]['I'],M[k[1]]['I'])['U']
    elif len(M) > 0 and len(k) > 0:
        U  = M[k[0]]['I']
    return U

@cython.boundscheck(False)
@cython.nonecheck(False)
def combinatoric_features(dict S, bint no_idx=True, bint ascending=False, bint singles=False):
    cdef dict M,P,I,F
    cdef int i,n,y
    cdef tuple x,k
    n = len(S)
    M = {}
    if no_idx:
        if not singles:
            M = {i:{} for i in range(2,n+1)}
            for i in M: #generate combinations
                P = {}
                for k in [tuple(sorted(x,reverse=not ascending)) for x in it.combinations(S.keys(),i)]:
                    P[k] = {}
                M[i] = P
            P = {}
            for i in S: #singletons 
                P[(i,)] = {'U':S[i],'I':S[i]} 
            M[1] = P
            #cascade on the new k[-1] value into the old
            for i in range(2,len(M)):
                for k in M[i]:
                    x,y = k[0:-1],k[-1]
                    I = features_no_idx(M[i-1][x]['I'],S[y])['I']
                    F = features_no_idx(M[i-1][x]['U'],S[y])
                    M[i][k] = {'I':I,'U':F['U'],'D1':F['D1'],'D2':F['D2']}
            for k in M[1]:
                M[1][k].pop('U') #take out the singles union
        else:
            M[1] = {(i,):{'I':S[i]} for i in S}
    else:
        if not singles:
            M = {i:{} for i in range(2,n+1)}
            for i in M: #generate combinations
                P = {}
                for k in [tuple(sorted(x,reverse=not ascending)) for x in it.combinations(S.keys(),i)]:
                    P[k] = {}
                M[i] = P
            P = {}
            for i in S: #singletons 
                P[(i,)] = {'U':S[i],'I':S[i]} 
            M[1] = P
            #cascade on the new k[-1] value into the old
            for i in range(2,len(M)):
                for k in M[i]:
                    x,y = k[0:-1],k[-1]
                    I = features(M[i-1][x]['I'],S[y])['I']
                    F = features(M[i-1][x]['U'],S[y])
                    M[i][k] = {'I':I,'U':F['U'],'D1':F['D1'],'D2':F['D2']}
            for k in M[1]:
                M[1][k].pop('U') #take out the singles union
        else:
            M[1] = {(i,):{'I':S[i]} for i in S}
        
    return M
    
#given the combinatroial features of a group
#for each value k, merge k-combination intersections
def combinatoric_weights(dict M):
    cdef dict W
    cdef int w,t,i
    W = {w:{} for w in M}    
    for w in M: #the weights
        W[w] = cascading_pileup(M[w],M[w].keys()) #union of all intersections
        for t in W[w]:
            for i in range(len(W[w][t])):
                W[w][t][i][4],W[w][t][i][5] = w*1.0,w*1.0
    return W

#given a weighted pileup and combination
#make for each weight/type a new comb key as: (w1, ),'UoI', (w2,),'UoI', etc...
def merge_pile_comb(dict W,dict M):
    cdef int c,t
    for c in M:
        k = tuple([-1*c for i in range(c)])
        M[c][k] = {'U':{}} #intersection of unions tag as union
        for t in W[c]:
            M[c][k]['U'][t] = W[c][t]
    return M

#return the subset of M with combination c
#with_comb means either the side with all entries that contain the query combination c
#or the other side that only has the containment of the query combination c
#leave_single means after with_comb is working, put back the singleton version of one c
#depp_copy means a new dict structure is constructed leaving M in tact
def features_subset(dict M,list c,with_comb=False,leave_single=False,deep_copy=False):
    #if type(c) is int or type(c) is long: c = [c]
    cdef dict N
    N = {}
    for w in M:
        for k in M[w]:
            if with_comb:
                if set(c).issubset(set(k)):
                    if deep_copy:
                        if N.has_key(w): N[w][k] = copy.deepcopy(M[w][k])
                        else:            N[w] = {k:copy.deepcopy(M[w][k])}
                    else:
                        if N.has_key(w): N[w][k] = M[w][k]
                        else:            N[w] = {k:M[w][k]}
            else:
                if not set(c).issubset(set(k)):
                    if deep_copy:
                        if N.has_key(w): N[w][k] = copy.deepcopy(M[w][k])
                        else:            N[w] = {k:copy.deepcopy(M[w][k])}
                    else:
                        if N.has_key(w): N[w][k] = M[w][k]
                        else:            N[w] = {k:M[w][k]}
        if leave_single and w <= 1: #put back single instances if flagged
            for k in M[w]:
                for i in range(len(c)):
                    if (c[i],)==k:
                        if deep_copy:
                            if N.has_key(w): N[w][k] = copy.deepcopy(M[w][k])
                            else:            N[w] = {k:copy.deepcopy(M[w][k])}
                        else:
                            if N.has_key(w): N[w][k] = M[w][k]
                            else:            N[w] = {k:M[w][k]}
                        
    return N
################################################################################   