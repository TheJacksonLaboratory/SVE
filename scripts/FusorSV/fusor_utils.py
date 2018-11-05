import copy
import hashlib
import random
import math
from bisect import bisect_left
import socket
import cPickle as pickle
import glob
import time
import itertools as it
import numpy as np
import scipy.stats as stats
#local libs
import fusion_utils as fu

#puts D2 into D1
def tb_merge(D1,D2):
    for t in D2:
        for b in D2[t]:
            if D1.has_key(t): D1[t][b] = D2[t][b]
            else: D1[t] = D2[t]
    
#use a shuffled md5-hashed hostname with a specified length
def get_identifier(length=1000):
    l = int(round(math.log(length),0))+1
    return ''.join(random.sample(hashlib.md5(socket.gethostname()).hexdigest(),min(l,hashlib.md5().digestsize)))

#given a partition P[t][b][c][s] write a josn for each
def write_partitions_by_sample(sname_partition_path,P):
    for t in P:
        for b in P[t]:
            for c in P[t][b]:
                for s in P[t][b][c]:
                    path = sname_partition_path+'_S%s_T%s_B%s.pickle'%(c,t,b)
                    S = {t:{b:{c:{s:P[t][b][c][s]}}}}
                    with open(path,'wb') as f:
                        pickle.dump(S,f)
    return True

#for a partition of type and bin, pool all samples and all callers
def read_partitions_by_caller(partition_path,callers,exclude,t,b,verbose=False):
    P = {t:{b:{c:{} for c in callers}}}
    for c in set(callers).difference(set(exclude)):
        samples = glob.glob(partition_path+'*_S%s_T%s_B%s.pickle'%(c,t,b)) #keys are in the pickler
        for sample in samples:
            with open(sample, 'rb') as f:
                start = time.time()
                if verbose: print('reading %s'%sample)
                S = pickle.load(f)
                s = S[t][b][c].keys()[0]
                stop = time.time()
                if verbose: print('finished loading %s in %s sec'%(sample,round(stop-start,2)))
            P[t][b][c][s] = copy.deepcopy(S[t][b][c][s])
    return P

#for a partition get all the call sets aka the partition
def assemble_partition(L,t,b,exclude=[]):
    P = {t:{b:{}}}
    for i in L:
        S = i[1] #i[0] is sname and should be s inside S[t][b][c].keys()
        for c in S[t][b]:
            P[t][b][c] = S[t][b][c]
    return P

#get one sample worth of data AKA all partitions
def read_partitions_by_sample(partition_path,sname):
    sample_data = glob.glob(partition_path+'/%s*.pickle'%sname)
    Q = {}
    for i in range(len(sample_data)):
        with open(sample_data[i],'rb') as f: S = pickle.load(f)
        for t in S:
            if not Q.has_key(t):
                Q[t] = {}
            for b in S[t]:
                if not Q[t].has_key(b):
                    Q[t][b] = {}
                for c in S[t][b]:
                    if not Q[t][b].has_key(c):
                        Q[t][b][c] = {}
                    for s in S[t][b][c]:
                        Q[t][b][c][s] = S[t][b][c][s]
    return Q

#[(t,b),J,D,E,alpha]
def assemble_model(L):
    J,D,E,alpha,K = {},{},{},{},{}
    for i in range(len(L)):
        (t,b),j,d,e,a,k = L[i]
        if J.has_key(t):
            J[t][b]     = j[t][b]
            D[t][b]     = d[t][b]
            E[t][b]     = e[t][b]
            alpha[t][b] = a[t][b]
            K[t][b]     = k[t][b]
        else:
            J[t]     = {b:j[t][b]}
            D[t]     = {b:d[t][b]}
            E[t]     = {b:e[t][b]}
            alpha[t] = {b:a[t][b]}
            K[t]     = {b:k[t][b]}
    return J,D,E,alpha,len(L),K
    
#save the main variables or the pooled target search into J pickled
def read_pickle_metric_search(pickle_path):
    J = {}
    with open(pickle_path+'metric_search.pickle', 'rb') as f:
        start = time.time()
        print('reading %s'%pickle_path+'metric_search.pickle')
        J = pickle.load(f)
        stop = time.time()
        print('finished loading %s in %s sec'%(pickle_path+'metric_search.pickle',round(stop-start,2)))
    return J

#this is a table based export for visualization
def export_distance_matrix(D,callers,types,bins,path,sim=True):
    header = ['C1','C2','type','bin','j']
    s = '\t'.join(header)+'\n'
    if sim: #simularity instead of distance
        for t in D:
            b = sorted(D[t].keys())
            for i in range(len(D[t])):
                for c1,c2 in sorted(D[t][b[i]]):
                    s += '\t'.join([callers[c1],callers[c2],types[t],bins[t][i],str(1.0-D[t][b[i]][(c1,c2)])])+'\n'
    else:#reagular 0.0 to 1.0 distance measure
        for t in D:
            b = sorted(D[t].keys())
            for i in range(len(D[t])):
                for c1,c2 in sorted(D[t][b[i]]):
                    s += '\t'.join([callers[c1],callers[c2],types[t],bins[t][i],str(D[t][b[i]][(c1,c2)])])+'\n'
    with open(path,'w') as f:
        f.write(s)
        return True
    return False

#this is a table based export for visualization
def export_caller_performance(cfs,callers,path):
    header = ['sname','caller','type','prec','rec','f1','j','n','m','l_mu','r_mu','l_sd','r_sd']
    s = '\t'.join(header)+'\n'
    for c in cfs:
        for t in cfs[c]:
            for sample in cfs[c][t]:#['sname','caller','type','prec','rec','f1','j','l_mu','r_mu']
                s += '\t'.join([str(i) for i in [sample[0],callers[c],sample[1]]+sample[2:6]+sample[7:]])+'\n' #metric you want to plot
    with open(path,'w') as f:
        f.write(s)
        return True
    return False

def export_detailed_performance(ds,callers,path):
    header = ['sname','caller','type','bin','prec','rec','f1','j','n','m','l_mu','l_sd','r_mu','r_sd']
    s = '\t'.join(header)+'\n'
    for c in ds:
        for t in ds[c]:
            for b in ds[c][t]:
                for sample in ds[c][t][b]:
                    s += '\t'.join([str(i) for i in [sample[0],callers[c],sample[1]]+sample[2:6]+sample[7:]])+'\n'
    with open(path,'w') as f:
        f.write(s)
        return True
    return False
    
def export_caller_by_type_and_bin(E,alpha,callers,types,bins,path):
    header = ['caller','type','bin','j']
    s = '\t'.join(header)+'\n'
    for t in E:
        b = sorted(E[t].keys())
        for i in range(len(E[t])):
            written = {}
            for g in E[t][b[i]]:
                if len(g)==1 and g[0] is not None:
                    s += '\t'.join([callers[g[0]],types[t],bins[t][i],str(E[t][b[i]][g])])+'\n'
                    written[g[0]] = True
            for c in callers:
                if c not in written:
                    s += '\t'.join([callers[c],types[t],bins[t][i],'0.0'])+'\n'
            s += '\t'.join(['alpha',types[t],bins[t][i],str(alpha[t][i])])+'\n'
    with open(path,'w') as f:
        f.write(s)
        return True
    return False
        
#given a path write out the model that was computed
def export_fusion_model(B,J,D,E,alpha,n,K,path):
    with open(path,'wb') as f:
        pickle.dump({'B':B,'J':J,'D':D,'E':E,'alpha':alpha,'n':n,'K':K},f)
        return True
    return False

#given a path write out the model that was computed
def import_fusion_model(path):
    B,J,E,alpha,K = {},{},{},{},{}
    with open(path,'rb') as f:
        P = pickle.load(f)
        B,J,D,E,alpha,n,K = P['B'],P['J'],P['D'],P['E'],P['alpha'],P['n'],P['K']
    return B,J,D,E,alpha,n,K

#filter out these regions
#Q[t][id][sname]
def filter_samples(Q,R,ids=[-1],r=0.0):
    for t in Q:
        print('filtering svtype %s '%t)
        for i in ids:
            for sname in Q[t][i]:
                print('filtering sample %s'%sname)
                Q[t][i][sname] = fu.filter_regions(Q[t][i][sname],R,r)

#filter any touch
def filter_samples2(Q,R,ids=[-1]):
    for t in Q:
        print('fitering svtype %s'%t)
        for i in ids:
            for sname in Q[t][i]:
                print('filtering sample %s'%sname)
                Q[t][i][sname] = fu.filter_regions2(Q[t][i][sname],R)

#new pooled optimal group selection
#slice out the samples for a [[svult_1,sname_1],[svult_2,sname_2],...[svult_s,snmae_s]]
def slice_samples(L,exclude=[]):
    #do empty for types here for downstream binning
    all_types = set([])
    for i in L:
        for c in i[1]:
            for t in i[1][c]:
                all_types.add(t)
    all_types    = list(all_types)      
    #target_types = list(set([j for x in [s[1][k].keys() for s in L] for j in x])) #target types
    Q = {}
    
    # DEBUG
#     index = 0
#     print "types list: "+str(all_types)
    
    for t in all_types:
        Q[t] = {}
        for i in range(len(L)): #L[i][0] is sname, L[i][1] is the data
        
            # DEBUG
#             if index < 10:
#                 print "Breakpoint 1: "+str(L[i][0])+" "+str(L[i][1])
        
            if not L[i][0] in exclude: 
                
                # DEBUG
#                 if index < 10:
#                     print "Breakpoint 2"
                
                for c in L[i][1]:
                    if Q[t].has_key(c):
                        if L[i][1][c].has_key(t):
                            Q[t][c][L[i][0]] = L[i][1][c][t]
                    else:
                        if L[i][1][c].has_key(t):
                            Q[t][c] = {L[i][0]:L[i][1][c][t]}
                
            # DEBUG            
            # index += 1
                            
        # DEBUG
#         print "Q:"
#         for key, value in Q.iteritems():
#             print str(key)+"; "+str(value)
#             index += 1
#             if index > 10:
#                 break
                
    return Q

#cluster the calls from several samples worth of SVUL
#for each type and caller make one large pooled and clustered sample set
def pre_cluster_samples(P,r):
    C = {}
    for t in P: #types
        C[t] = {}
        for b in P[t]: #callers
            C[t][b] = {}
            for c in P[t][b]:
                L,C[t][b][c] = [],{'CLUSTER':[]}
                for s in P[t][b][c]: #samples
                    L += idx_sample(P[t][b][c][s],s)
                if len(L) > 0:
                    M = {}
                    try:
                        M = coordinate_cluster_overlap(L,r)
                    except Exception:
                        print t,b,c,s
                        pass
                    for m in M: #for each cluster in M, merge it down
                        s1 = np.uint32(round(1.0*sum([x[0] for x in M[m]])/len(M[m]),0))
                        s2 = np.uint32(round(1.0*sum([x[1] for x in M[m]])/len(M[m]),0))
                        idx = idx_sample_merge([x[6] for x in M[m]])
                        C[t][b][c]['CLUSTER'] += [[s1,s2,M[m][0][2],M[m][0][3],M[m][0][4],M[m][0][5],idx]]
    return C

#given a list of I merge to {idx}={c_id:[sname_row,sname_row,ect...]}
#so we can still pull all the supporting data rows from all the callers
def idx_sample_merge(I):
    X = {}
    for idx in I:
        for c_id in idx:
            if not X.has_key(c_id): X[c_id] = idx[c_id]
            else:                   X[c_id] = X[c_id].union(idx[c_id])
    return X

#general idx string conversion
def idx_sample(L,sample):
    for i in range(len(L)):
        L[i] =  L[i][0:6]+[idx_to_str(L[i][6],sample)]
    return L

def idx_to_str(idx,sample):
    n_idx = {}
    for c in idx:
        n_idx[c] = set([sample+'_'+str(i) for i in idx[c]])
    return n_idx

#pretty fast and pretty good
#not testing as functional
def coordinate_cluster_overlap(L,r):
    C = {}
    if len(L)>1:
        d = sorted(L,key=lambda x:x[0])
        i,j,x = 0,0,0
        while i < len(d)-1:
            x,j = i,i+1
            while overlap(d[i][0],d[i][1],d[j][0],d[j][1]) >=r and j < len(d)-1: j += 1
            C[i] = d[i:j]
            i = j
        if overlap(d[x][0],d[x][1],d[j][0],d[j][1]) >=r : C[x] += [d[j]]
        else:                                             C[j]  = [d[j]]
    elif len(L)<=1: C[0] = L
    return C               

#give the clustering ratio P >> C
def cluster_ratio(P,C,c_id):
    R = {}
    for t in P:
        R[t] = {}
        for b in P[t]:
            if P[t][b].has_key(c_id):
                p = sum([len(P[t][b][c_id][k]) for k in P[t][b][c_id]])
                c = sum([len(C[t][b][c_id][k]) for k in C[t][b][c_id]])
                if c <= 0.0: R[t][b] = 0.0
                else:        R[t][b] = 1.0*p/c   
    return R
    
#compute the bin positions for each type is a type is not present in Q[k]
#if a bin is not present in the target, use equally spaced from start to stop
#Q[t][c][s] = B[t]
def distribute_bins(Q,k,n_b=10,m_b=2000,lower=None,upper=None,event=False):
    bins = {t:[] for t in Q}
    for t in Q:
        if Q[t].has_key(k):
            bins[t] = equal_bins(Q,t,k,n_b,m_b,lower,upper,event=event,target=True)
        else:               
            bins[t] = equal_bins(Q,t,k,n_b,m_b,lower,upper,event=event,target=False)
    return bins 
 
#given a pooled svul, histogram and set the bins sizes from start to stop
#using an equal spacing for each bin allocation trying to number of bins n_b
#event=True is for equal events, event=False is for equal bp distribution
#if a bin is not present in the target, use equally spaced from start to stop
def equal_bins(Q,t,k=0,n_b=10,m_b=2000,lower=None,upper=None,event=True,target=True):
    #do upper and lower limits if there is no avlues in the Histogram
    C = []
    if target:
        for s in Q[t][k]: C += Q[t][k][s]
        if len(C) < 1:
            for c in Q[t]:
                for s in Q[t][c]:
                    C += Q[t][c][s] 
    else:
        for c in Q[t]:
            for s in Q[t][c]:
                C += Q[t][c][s]
    #lower n_b if needed using the criteria of len(C)
    if len(C) > m_b:
        while n_b>0 and len(C)/n_b < m_b: n_b -= 1           
        n_b = max(1,n_b)
    else:
        n_b = 1
    H = fu.discrete_mag_hist(C)
    B =[] #take a lower at the lower boundry != 1
    if lower is not None:
        B += [max(np.uint32(1),min(sorted(H)[0]-np.uint32(1),np.uint32(lower)))]
    else:
        B += [max(np.uint32(1),sorted(H)[0]-np.uint32(1))]
    if event:   #do number of events per bin: b  
        i,b = np.uint32(0),np.uint32(round(1.0*len(C)/n_b,0)) #:::TO DO::: this can be more robust
        for j in sorted(H):
            i += np.uint32(H[j])  #accumulate sum
            if i >= b: #check if it is over the target
                if j>np.uint32(1): B += [np.uint32(j)]
                i = np.uint32(0)
    else: #try to place an equal amount of bp per bin
        i,n = np.uint64(0),np.uint64(0) #can't have more than the human genome size...
        for j in H: n += np.uint64(j)*np.uint64(H[j]) #size of bp
        b = np.uint64(round(1.0*n/n_b,0))          #:::TO DO::: this can be more robust
        for j in sorted(H):
            i += np.uint64(j)*np.uint64(H[j])
            if i >= b:
                if j>np.uint64(1): B += [np.uint32(j)]
                i = np.uint32(0)
    #update the value for the last bin
    if upper is not None: #extend the bin to include values not in max([H[i] for i in H])
        if i > 0: B += [max(sorted(H)[-1]+np.uint32(1E3),np.uint32(upper))]
        else:     B[-1] = np.uint32(upper)
    else:                 #use the values in max([H[i] for i in H])+1E3 wiggle
        if i > 0: B += [sorted(H)[-1]+np.uint32(1E3)]
        else:     B[-1] += np.uint32(1E3)
    return B

#apply precomputed partitions to each sample in ||
def partition_sample(S,B):
    P,types = {},list(set([i for j in [S[k].keys() for k in S] for i in j]))
    for t in types:     #union of all types in the sample
        P[t] = {}
        for c in S:
            if S[c].has_key(t):
                N = fu.partition_by_mag(S[c][t],B[t])           #non interpolating partitioning
                for i in range(len(B[t])-1):                    #precomputed bin sizes here
                    if P[t].has_key(i): P[t][i][(c,)] = N[i]    #use single tuples as keys
                    else:               P[t][i] = {(c,):N[i]}   #instead of the single ints
    return P  
    
#apply bin partitions to the sliced samples
#allow caller keys to be filtered out here with exclude=[k1,k2,...]
def partition_sliced_samples(Q,B,exclude=[]):
    P,N = {},{}
    for t in Q:
        P[t] = {}
        C = set(Q[t].keys()).difference(set(exclude))
        for c in C:
            for s in Q[t][c]:
                N = fu.partition_by_mag(Q[t][c][s],B[t]) #this will need to be distributed
                for i in range(len(B[t])-1):
                    if P[t].has_key(i):
                        if P[t][i].has_key(c): P[t][i][c][s] = N[i]
                        else:                  P[t][i][c] = {s:N[i]}
                    else:
                        P[t][i] = {c:{s:N[i]}}
    return P

#def partition_by_mag(C,B,self_merge=False,interpolate=False):
#    b = len(B)
#    M,P = [],{i:[] for i in range(b-1)}
#    if self_merge: C = fu.merge_regions(C)
#    if interpolate:
#        M = [(B[i]+B[i+1])/2.0 for i in range(b-1)] #bin midpoints
#        for i in range(len(C)):
#            v = fu.get_xy_mag(C[i])
#            if   v < M[0]: #check first and last
#                P[0]   += [[C[i][0],C[i][1],C[i][2],C[i][3],1.0*C[i][4],1.0*C[i][5],C[i][6]]]
#            elif v >= M[b-2]:
#                P[b-1] += [[C[i][0],C[i][1],C[i][2],C[i][3],1.0*C[i][4],1.0*C[i][5],C[i][6]]]
#            else:
#                for j in range(b-2):
#                    if v >= M[j] and v < M[j+1]: #bin center nearest neighbors
#                        w_a = 1.0-(v-M[j])/(M[j+1]-M[j]) #linear interpolation
#                        w_b = 1.0 - w_a                  #and complement
#                        P[j]  += [[C[i][0],C[i][1],C[i][2],C[i][3],w_a*C[i][4],w_a*C[i][5],C[i][6]]]
#                        P[j+1]+= [[C[i][0],C[i][1],C[i][2],C[i][3],w_b*C[i][4],w_b*C[i][5],C[i][6]]]
#    else:
#        for i in range(len(C)):
#            n = fu.get_xy_mag(C[i])
#            for j in range(b-1):
#                if n >= B[j] and n < B[j+1]:
#                    P[j] += [C[i]]
#    return P
    
#agregate all the bins back together again for Q in OOC|| version
def unpartition_sliced_samples(P):     
    Q = {}
    for t in P:
        Q[t] = {}
        for b in P[t]: #agregation point is here-----------------
            for c in P[t][b]:
                if not Q[t].has_key(c):
                    Q[t][c] = {}
                for s in P[t][b][c]:
                    if not Q[t][c].has_key(s):
                        Q[t][c][s] = copy.deepcopy(P[t][b][c][s])
                    else:
                        Q[t][c][s] += copy.deepcopy(P[t][b][c][s])
    for t in Q:
        for c in Q[t]:
            for s in Q[t][c]:
                Q[t][c][s] = sorted(Q[t][c][s],key=lambda x: x[0])
    return Q
        
#for each sample get pairwise feature magnitudes and pool together
#P[t][b][c][s] => M[t][b][(i,j)][|I|,|U|,|D1|,|D2|]
def all_samples_all_pairs_magnitudes(P,snames,self_merge=True):
    M,N,F = {},[],[]
    for t in P:
        M[t] = {}
        for b in P[t]:
            M[t][b] = {}
            for i,j in sorted(it.combinations(P[t][b].keys(),2)):
                N = [np.uint64(0),np.uint64(0),np.uint64(0),np.uint64(0)] #initial
                for s in snames:
                    C1,C2 = [],[]
                    if P[t][b][i].has_key(s): C1 = P[t][b][i][s]
                    if P[t][b][j].has_key(s): C2 = P[t][b][j][s]                   
                    F = fu.feature_magnitudes(C1,C2,self_merge)
                    N[0] += np.uint64(F[0])
                    N[1] += np.uint64(F[1])
                    N[2] += np.uint64(F[2])
                    N[3] += np.uint64(F[3])
                M[t][b][(i,j)] = N
    return M

#i <- is.na(A) | is.na(B)
#dist(rbind(A[!i], B[!i])) * sqrt(length(A) / length(A[!i]))
def group_brkpt_stats(X1,X2):
    X3 = X1+X2 #add to propigate nans
    X4 = X3[~np.isnan(X3)]
    x_n,x_sm,x_mu,x_md,x_sd,x_sk,x_ks = 0,0,0.0,0.0,0.0,0.0,0.0
    x_n = len(X4)
    x_sm = sum(X4)
    if x_n>0:
        x_mu = np.mean(X4)
        x_md = np.median(X4)
        x_sd = np.std(X4)
        x_sk = stats.skew(X4)
        x_ks = stats.kurtosis(X4)
    return [x_n,x_sm,x_mu,x_md,x_sd,x_sk,x_ks]    
    

def brkpt_stats(X):
    X1 = X[~np.isnan(X)]
    x_n,x_sm,x_mu,x_md,x_sd,x_sk,x_ks = 0,0,0.0,0.0,0.0,0.0,0.0
    x_n = len(X1)
    x_sm = sum(X1)
    if x_n>0:
        x_mu = np.mean(X1)
        x_md = np.median(X1)
        x_sd = np.std(X1)
        x_sk = stats.skew(X1)
        x_ks = stats.kurtosis(X1)
    return [x_n,x_sm,x_mu,x_md,x_sd,x_sk,x_ks]
    
#breakpoint differentials to the target
def all_samples_all_callers_bkpts(P,snames,k=0,r=0.5,self_merge=False,exclude=[1],verbose=False):
    M = {}
    for t in P:
        M[t] = {}
        for b in P[t]:
            M[t][b] = {}
            #marginal distributions
            for i in set(P[t][b].keys()).difference(set(exclude+[k])): #just this section when partitioning in OOC||
                M[t][b][(i,)] = {'L':[],'R':[],'L_stats':[],'R_stats':[]}
                for sname in snames:
                    C1,C2 = [],[]
                    if P[t][b].has_key(k) and P[t][b][k].has_key(sname): C1 = P[t][b][k][sname]
                    if P[t][b].has_key(i) and P[t][b][i].has_key(sname): C2 = P[t][b][i][sname]
                    s = fu.metric_score(C1,C2,r,self_merge=False,d1_ids=False,d2_ids=False)
                    for brkpt in s['b']: #leave in the None now....
                        if brkpt is not None:
                            M[t][b][(i,)]['L'] += [np.float64(brkpt[0])]
                            M[t][b][(i,)]['R'] += [np.float64(brkpt[1])]
                        else:
                            M[t][b][(i,)]['L'] += [np.float64(np.nan)]
                            M[t][b][(i,)]['R'] += [np.float64(np.nan)]
                M[t][b][(i,)]['L'] = np.array(M[t][b][(i,)]['L'])
                M[t][b][(i,)]['R'] = np.array(M[t][b][(i,)]['R'])
                M[t][b][(i,)]['L_stats'] = brkpt_stats(M[t][b][(i,)]['L'])
                M[t][b][(i,)]['R_stats'] = brkpt_stats(M[t][b][(i,)]['R'])
            #powerset of joint distributions
            C =  set([v for w in M[t][b] for v in w]).difference(set([k]))
            for i in range(2,len(C)+1): #generate combinations
                for g in [tuple(sorted(x,reverse=False)) for x in it.combinations(C,i)]:
                    M[t][b][g] = {'L':[],'R':[],'L_stats':[],'R_stats':[]}
                    #[1] average over the dimension
                    M[t][b][g]['L'] = np.sum(np.array([(1.0/len(g))*M[t][b][(j,)]['L'] for j in g]),axis=0)
                    M[t][b][g]['R'] = np.sum(np.array([(1.0/len(g))*M[t][b][(j,)]['R'] for j in g]),axis=0)
                    #[2] get stats on the averaged group
                    M[t][b][g]['L_stats'] = brkpt_stats(M[t][b][g]['L'])
                    M[t][b][g]['R_stats'] = brkpt_stats(M[t][b][g]['R'])
            for g in M[t][b]: #toss away intermediates
                M[t][b][g].pop('L')
                M[t][b][g].pop('R')
            #maybe keep one big distance matrix on this?
    return M

#given the all samples all callers breakpt differential map
#compute a normalized pairwise L,R brkpt distance matrix
def brkpt_pair_matrix(M):
    D = {}
    for t in M:
        D[t] = {}
        for b in M[t]:
            D[t][b] = {}
            for i,j in it.combinations(M[t][b],2):
                i,j = min(i,j),max(i,j)
                D[t][b][(i,j)] = []

#given two brkpt differential rows with None for NA values,
#clean out the rows that have None in either
def brkpt_clean(x1,x2):
    return []
    
                
#return the index in C1 where i appears
#:::TO DO binary search:::
def get_call_index(C1,c,i,IDX=6):
    j = -1
    for x in range(len(C1)):
        if C1[x][IDX][c]==i:
            j = x
            break
    return j

def get_all_call_index(P,IDX=6):
    I = {}
    for t in P:
        I[t] = {}
        for b in P[t]:
            I[t][b] = {}
            for c in P[t][b]:
                I[t][b][c] = {}
                for s in P[t][b][c]:
                    I[t][b][c][s] = {}
                    for i in range(len(P[t][b][c][s])):
                        for x in P[t][b][c][s][i][IDX][c]:
                            I[t][b][c][s][x] = i
    return I

#given the availble groupings in c, check for the
#tightest brkpt smoothing options looking at std in 0.0 < x <= std
def tightest_brkpt_comb(K,M,t,b,c,x=4):
    a_l,a_r = M[t][b]['L_stats'][x],M[t][b]['R_stats'][x]
    for i in range(1,len(c)+1):
        for g in it.combinations(c,i):
            g = tuple(sorted(list(g))) #sort keys
            if K[t][b][g]['L_stats'][x]>0.0 and K[t][b][g]['L_stats'][x]<K[t][b][a_l]['L_stats'][x]:
                a_l = g
            if K[t][b][g]['R_stats'][x]>0.0 and K[t][b][g]['R_stats'][x]<K[t][b][a_r]['R_stats'][x]:
                a_r = g
    if K[t][b][a_l]['L_stats'][x] >= K[t][b][M[t][b]['L_stats'][x]]['L_stats'][x]:
        a_l = None
    if K[t][b][a_r]['R_stats'][x] >= K[t][b][M[t][b]['R_stats'][x]]['R_stats'][x]:
        a_r = None    
    return [a_l,a_r]

#n = number of stats    
def max_brkpt_stats(K,n=7):
    M = {}
    for t in K:
        M[t] = {}
        for b in K[t]:
            M[t][b] = {'L_stats':[],'R_stats':[]}
            for j in range(n):
                M[t][b]['L_stats'] += [K[t][b].keys()[np.argmax([K[t][b][x]['L_stats'][j] for x in K[t][b]])]]
                M[t][b]['R_stats'] += [K[t][b].keys()[np.argmax([K[t][b][x]['R_stats'][j] for x in K[t][b]])]]
    return M                

#take the calls and stats and integrate the answer
#only on thesource coordinate space currently
def smooth_brkpt(f1,C,L_stats,R_stats):
    ls,rs,p = [],[],float(len(C))
    l_sd = np.sum([L_stats[c][4] for c in L_stats])
    r_sd = np.sum([R_stats[c][4] for c in R_stats])
    l_w  = list(np.array([l_sd/(L_stats[c][4]+p) for c in L_stats])/sum([l_sd/(L_stats[c][4]+p) for c in L_stats]))
    l_w  = {L_stats.keys()[i]:l_w[i]for i in range(len(l_w))}
    r_w  = list(np.array([r_sd/(R_stats[c][4]+p) for c in R_stats])/sum([r_sd/(R_stats[c][4]+p) for c in R_stats]))
    r_w  = {R_stats.keys()[i]:r_w[i]for i in range(len(r_w))}
    for c in C:
        for i in range(len(C[c])):
            ls += [l_w[c]*(int(f1[0])-int(C[c][i][0])+int(L_stats[c][2]))] #difference of the mean shift
            rs += [r_w[c]*(int(f1[1])-int(C[c][i][1])+int(R_stats[c][2]))] #difference of the mean shift
    if sum(ls) is np.nan or sum(rs) is np.nan: ls,rs = [0],[0]
    l = np.uint32(max(int(round(sum(ls)+f1[0],0)),0))
    r = np.uint32(max(int(round(sum(rs)+f1[1],0)),0))
    if(r<l): l,r = f1[0],f1[1]
    return [l,r]+f1[2:]

def best_smooth_brkpt(f1,C,K,M,t,b,x_stat=4):
    l,r,ls,rs,g = f1[0],f1[1],[],[],tuple(sorted(C.keys()))
    [l_g,r_g] = tightest_brkpt_comb(K,M,t,b,g,x_stat)
    if not l_g is None:
        for c in C:
            for i in range(len(C[c])):
                if c in l_g: ls += [C[c][i][0]]
        l = np.uint32(np.round(np.mean(ls)+K[t][b][g]['L_stats'][2],0)) #mean shift
    if not r_g is None:
        for c in C:
            for i in range(len(C[c])):
                if c in r_g: rs += [C[c][i][1]]
        r = np.uint32(np.round(np.mean(rs)+K[t][b][g]['R_stats'][2],0)) #mean shift
    return [l,r]+f1[2:]

def best_smooth_brkpt_samples_search(F,K,P,I,M,exclude=[1],IDX=6,x_stat=4):
    S = {}
    for t in F:
        S[t] = {}
        for b in F[t]:
            S[t][b] = {}
            for s in F[t][b]:
                S[t][b][s] = []
                for i in range(len(F[t][b][s])):       #each call now
                    idx   = F[t][b][s][i][IDX]         #pull the caller map
                    c_ids,C = sorted(idx.keys()),{}    #attach keys to pull calls from
                    #retrive the distributions from K
                    for c in set(c_ids).difference(set(exclude)):
                        C[c] = [P[t][b][c][s][I[t][b][c][s][x]] for x in idx[c]]
                    S[t][b][s] += [best_smooth_brkpt(F[t][b][s][i],C,K,M,t,b,x_stat)]
    return S

def best_smooth_brkpt_samples(F,K,P,exclude=[1],IDX=6,x_stat=4):
    I,M,S = get_all_call_index(P),max_brkpt_stats(K),{}
    for t in F:
        S[t] = {}
        for b in F[t]:
            S[t][b] = {}
            for s in F[t][b]:
                S[t][b][s] = []
                for i in range(len(F[t][b][s])):       #each call now
                    idx   = F[t][b][s][i][IDX]         #pull the caller map
                    c_ids,C = sorted(idx.keys()),{}    #attach keys to pull calls from
                    #retrive the distributions from K
                    for c in set(c_ids).difference(set(exclude)):
                        C[c] = [P[t][b][c][s][I[t][b][c][s][x]] for x in idx[c]]
                    S[t][b][s] += [best_smooth_brkpt(F[t][b][s][i],C,K,M,t,b,x_stat)]
    return S
   
#for each partition read out the brkpt distributions for L,R
#and integrate over them using the breakpoint smoothing algorithm
#F[t][b][s]
def smooth_brkpt_samples(F,K,P,exclude=[1],IDX=6):
    S,I = {},get_all_call_index(P)
    for t in F:
        S[t] = {}
        for b in F[t]:
            S[t][b] = {}
            for s in F[t][b]:
                S[t][b][s] = []
                for i in range(len(F[t][b][s])):       #each call now
                    idx   = F[t][b][s][i][IDX]         #pull the caller map
                    c_ids = sorted(idx.keys())         #attach keys to pull calls from
                    L_stats,R_stats,C = {},{},{}
                    #retrive the distributions from K
                    for c in set(c_ids).difference(set(exclude)):
                        L_stats[c] = K[t][b][(c,)]['L_stats']
                        R_stats[c] = K[t][b][(c,)]['R_stats']
                        if L_stats[c][0]>5 and R_stats[c][0]>5:
                            C[c] = [P[t][b][c][s][I[t][b][c][s][x]] for x in idx[c]]
                    S[t][b][s] += [smooth_brkpt(F[t][b][s][i],C,L_stats,R_stats)]
    return S

#uses the M[t][b][(i,j)][|I|,|U|,|D1|,|D2|]
def pooled_distance(M,mode='j'):
    D,NN = {},{}
    for t in M:
        D[t],NN[t] = {},{}
        for b in M[t]:
            D[t][b],NN[t][b] = {},{}
            #get the combintaion distances
            for g in M[t][b]:
                D[t][b][g] = np.float128(1.0)  #default is the  max distance
                if mode=='j' and np.float128(M[t][b][g][1]) > 0.0:
                    #1-(|I|/|U|)
                    D[t][b][g] = np.float128(1.0)-np.float128(M[t][b][g][0])/np.float128(M[t][b][g][1])
                if mode=='u' and np.float128(M[t][b][g][1]) > 0.0 and np.float128(M[t][b][g][2]) > 0.0:
                    #1-2*((|I|/|U|)*(|D1|/(|D1|+|D2|))/((|I|/|U|)+(|D1|/(|D1|+|D2|))
                    j = np.float128(M[t][b][g][0])/np.float128(M[t][b][g][1])
                    u = np.float128(M[t][b][g][2])/(np.float128(M[t][b][g][2])+np.float128(M[t][b][g][3]))
                    D[t][b][g] = np.float128(1.0)-np.float128(2.0)*(j*u)/(j+u)
            K = sorted(list(set([v for w in M[t][b] for v in w])))        
            for k in K:
                N = []
                for i,j in sorted(M[t][b].keys()): #distance then key
                    if k==i: N += [[D[t][b][(i,j)],j]]
                    if k==j: N += [[D[t][b][(i,j)],i]]
                NN[t][b][k] = sorted(N,key=lambda x: x[0])
    return D,NN

#given prior and new data smooth based on the magnitude
#of observations for the features for each pair
def additive_magnitude_smoothing(J,J_new,k=None):
    A = {}
    for t in J:
        A[t] = {}
        for b in J[t]:
            A[t][b],upper = {},np.uint64(0)
            for g in J[t][b]: #find largest observed union to penalize caller drop-out
                if J[t][b][g][1] > upper: upper = J[t][b][g][1]
            for g in J[t][b]:
                F,F_new = J[t][b][g],[np.uint64(1),upper,np.uint64(1),upper] #init default pseudo counts
                if J_new.has_key(t) and J_new[t].has_key(b) and J_new[t][b].has_key(g):
                    if k is None: F_new = J_new[t][b][g] #found key so eliminate the pseudo counts now
                    elif k in g:  F_new = J[t][b][g]     #swap old to new if needed for true
                else:
                    if k in g:    F_new = J[t][b][g]     #swap old to new if needed for true
                A[t][b][g] = [F[0]+F_new[0],F[1]+F_new[1],F[2]+F_new[2],F[3]+F_new[3]]
    return A

#agregation operations=======================================================================
#give it the all pairs magnitudes
#get back A[t][b][g][row][I,U,D1,D2]
def agregate_point_magitudes(J):
    G,S,A = set([]),{},{}
    for t in J:
        S[t] = {}
        for b in J[t]:
            S[t][b] = {}
            for g in J[t][b]:
                G.add(g)
                if not S[t][b].has_key(g[0]):
                    S[t][b][g[0]] = []
                if not S[t][b].has_key(g[1]):
                    S[t][b][g[1]] = []
    #collect the full matrix row from the lower triangle representation
    G = sorted(list(G)) #sort the groups to look at
    for t in J:
        for b in J[t]:
            for g in J[t][b]:
                S[t][b][g[0]] += [[g[1]]+J[t][b][g]]
                S[t][b][g[1]] += [[g[0]]+J[t][b][g]]
    for t in S:
        A[t] = {}
        for b in S[t]:
            A[t][b] = {}
            for g in sorted(S[t][b]):
                A[t][b][g] = np.zeros((len(S[t][b][g]),5),dtype='u8')
                for i in range(len(S[t][b][g])):
                    A[t][b][g][i] = S[t][b][g][i]
                #imposing ordering of caller keys
                x = np.argsort(A[t][b][g][:,0])  #0 is the caller id slot
                A[t][b][g] = A[t][b][g][x] #reordering in place
    return A

#  2*((D1/U)*(D2/U))/((D1/U)+(D2/U)) #harmonic difference
#  2*((D1/U)*(I/U))/((D1/U)+(I/U))   #harmonic uniqueness

#given the agregation A[t][b][g][row][I,U,D1,D2]
def pair_sum(A):
    X = {}
    for t in A:
        X[t] = {}
        for b in A[t]:
            X[t][b] = {}
            for g in A[t][b]:
                X[t][b][g] = np.array([np.sum(A[t][b][g][:,1]),np.sum(A[t][b][g][:,2]),
                                       np.sum(A[t][b][g][:,3]),np.sum(A[t][b][g][:,4])],dtype='f16')
    return X

#given the agregation A[t][b][g][row][I,U,D1,D2]
def pair_mean(A):
    X = {}
    for t in A:
        X[t] = {}
        for b in A[t]:
            X[t][b] = {}
            for g in A[t][b]:
                X[t][b][g] = np.array([np.mean(A[t][b][g][:,1]),np.mean(A[t][b][g][:,2]),
                                       np.mean(A[t][b][g][:,3]),np.mean(A[t][b][g][:,4])],dtype='f16')
    return X

#given the agregation A[t][b][g][row][I,U,D1,D2]
def pair_std(A):
    X = {}
    for t in A:
        X[t] = {}
        for b in A[t]:
            X[t][b] = {}
            for g in A[t][b]:
                X[t][b][g] = np.array([np.std(A[t][b][g][:,1]),np.std(A[t][b][g][:,2]),
                                       np.std(A[t][b][g][:,3]),np.std(A[t][b][g][:,4])],dtype='f16')
    return X

#given A[t][b][g][row][I,U,D1,D2]
#take the mean and std of each row with the true removed
#and compute its magnitude difference to use to modify the true
def impute_true(A,k):
    I = {}
    for t in A:
        I[t] = {}
        for b in A[t]:
            I[t][b] = {}
            for g in A[t][b]:
                I[t][b][g] = [np.uint64(0),np.uint64(0),np.uint64(0),np.uint64(0)]
                for i in range(4):
                    I[t][b][g][i] = np.float128(np.mean(A[t][b][g][np.where(A[t][b][g][:,0]!=k),i+1])+1)
    return I

def inject_imputed_true(J_new,I,k):
    for t in J_new:
        for b in J_new[t]:
            for g in sorted(I[t][b]):
                J_new[t][b][(k,g)] = [np.uint64(i) for i in I[t][b][g]]
    return J_new
                
#agregation operations=======================================================================

#get the pairwise distance matrix from the D:
#D[t][b][g]
def get_pwdm(D,G):
    P = {}
    for t in D:
        P[t] = {}
        for b in D[t]:
            P[t][b] = np.zeros((len(G),len(G)),dtype=np.float128)
            for i,j in it.combinations(range(len(G)),2):
                if D[t][b].has_key((G[i],G[j])):
                    P[t][b][i][j] = P[t][b][j][i] = D[t][b][(G[i],G[j])]
                elif D[t][b].has_key((G[j],G[i])):
                    P[t][b][i][j] = P[t][b][j][i] = D[t][b][(G[j],G[i])]
                else: #no key for (i,j) is missing
                    P[t][b][i][j] = P[t][b][j][i] = np.float128(1.0)
    return P

def sum_pairs(P):
    S = {}
    for t in P:
        S[t] = {}
        for b in P[t]:
            S[t][b] = []
            for row in P[t][b]:
                S[t][b] += [np.sum(row)/np.float128(len(row)-1)]
    return S

def exp_hist(e):
    h = {}
    for g in e:
        if h.has_key(e[g]): h[e[g]] += [g]
        else:               h[e[g]]  = [g]
    return np.array([[v,len(h[v])] for v in sorted(h)])
    
def exp_stats(e,e_post,a):
    if a < np.float128(1.0):
        H,H_post = exp_hist(e),exp_hist(e_post)
        upper,lower,value = 0,0,np.float128(1.0)
        for row in H:
            if row[0] >= a: upper += int(row[1])
        for row in H_post:
            lower += int(row[1])
            if lower >= upper:
                value = row[0]
                break
        h = np.array([e[g] for g in e])
        x = [a,value,len(e)-upper,upper,np.median(h),np.mean(h),np.std(h)]
    else:
        x = [a,a,len(e),0,np.float128(0.0),np.float128(0.0),np.float128(0.0)]
    return x
      
#given the old group expectation E and cutoff t,b value alpha
#estimates a alpha_post value that uses the histogram of the
#posterior estimate E_post to adjust the alpha value
def post_filter_cutoff(E,E_post,alpha):
    alpha_post = {}
    for t in E:
        alpha_post[t] = {}
        for b in E[t]:
            alpha_post[t][b] = exp_stats(E[t][b],E_post[t][b],alpha[t][b])[1]
    return alpha_post
                    
#given the target T and Group G with distance calculation D and Nearest Neibors NN
#compute the weight of the group which is its expected value given a pileup
def group_weight(g,k,D,NN,t,b):
    L,Q,w = [],[],np.float128(0.0)
    if NN[t][b].has_key(k):
        for i in NN[t][b][k]: #sort the group G
            if i[1] in g: Q += [i[1]]
        for i in range(len(Q)):
            c,j,n = np.float128(1.0),0,len(NN[t][b][Q[i]])
            while j<n and NN[t][b][Q[i]][j][1] == k: j+=1
            if i > 0: c = NN[t][b][Q[i]][j][0]        
            L += [np.float128(1.0)-(np.float128(1.0)-D[t][b][(k,Q[i])])*c]
        w = np.float128(1.0)-np.float128(np.prod(L))
    return w


#for all types and bins precompute all group weights
#M is the all pair all sample magnitudes (expected distance matrix input)
#k is the target key to use to calculate optimal groups in M[t][b] on
def all_group_weights(M,k,mode='j'):
    W = {}      
    D,NN = pooled_distance(M,mode)
    for t in M:
        W[t] = {}
        for b in M[t]:
            W[t][b] = {}
            C =  set([v for w in M[t][b] for v in w]).difference(set([k]))
            for i in range(1,len(C)+1): #generate combinations
                for g in [tuple(sorted(x,reverse=False)) for x in it.combinations(C,i)]:
                    W[t][b][g] = group_weight(g,k,D,NN,t,b)
    return W
#::::::::: optimal group selection ::::::::::::::::::::::::::::x

#given the type and bin, select an optimal group for pileup
#:::TO DO:::with each i in g: (i,)>alpha and CU{i}>C+alpha/2
#don't have to worry about the target key here as it has been applied
def select_groups(W,gamma=0.0):
    G = {}
    for t in W:
        G[t] = {}
        for b in W[t]:
            C,G[t][b] = [],{}
            for g in W[t][b]:
                if len(g)<=1 and W[t][b][g]>=gamma: C += [g[0]]
            if len(C)<=0:
                G[t][b] = {(None,):np.float128(0.0)}
            else:
                C = sorted(C)
                for i in range(1,len(C)+1):
                    for j in it.combinations(C,i):
                        G[t][b][j] = W[t][b][j]
    return G
                
#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp  
#given a svul, join idxs from list j
################################################################################
#given a svul, join idxs from list j
def join_idx(C,j,p):
    A = {}
    for idx in [C[i][p] for i in j]: #hard coded idx=6
        for k in idx:
            for i in idx[k]:
                if A.has_key(k): A[k].add(i)
                else:            A[k] =  {i}        
    return A

#given two idx merge them into one
def merge_idx(idx1,idx2):
    A = {}
    for k in idx1:
        for i in idx1[k]:
            if A.has_key(k): A[k].add(i)
            else:            A[k] =  {i}
    for k in idx2:
        for i in idx2[k]:
            if A.has_key(k): A[k].add(i)
            else:            A[k] =  {i}
    return A

#given a svul, join the y ranges from list [C[j] for j in x]
#if they overlap, otherwise leave them alone
def join_y(C,x,p):
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
    return A+[[B[m][1],B[m][0]] for m in range(len(B))]       #----->REV
  
#if the y values fron Y1 and Y2 overlap, merge them together
#otherwise add them as unique entries and return list of lists
def merge_y(Y1,Y2):
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
    return A+[[B[m][1],B[m][0]] for m in range(len(B))] 
 
#returns G[g][d][C[g]=>i, C[g][wx]] 
#as well as a flattened combination list for use
#with existing join_y,join_idx functions
#key f G are the coordinates x1,x2
def weight_graph(C):
    G,X = {},[]
    c_i,i_c = group_indecies(C)
    if len(i_c)>0:
        for g in sorted(C.keys()):
            X += C[g]
            for i in range(len(C[g])):
                for j in range(2): #just for x:wx for now
                    if G.has_key(C[g][i][j]):
                        if G[C[g][i][j]].has_key(g):
                            G[C[g][i][j]][g][j] = [i+c_i[g],C[g][i][4]]
                        else:
                            G[C[g][i][j]][g] = {j:[i+c_i[g],C[g][i][4]]} 
                    else:
                        G[C[g][i][j]] = {g:{j:[i+c_i[g],C[g][i][4]]}}
        last = sorted(G.keys())[-1]+2          #get last x2 value in G
        G[last] = {(None,):{None:[None,None]}} #to pad out a terminal vertex
    return X,G,c_i,i_c

#given the sorted keys of C
#calculate the offset you need to add to each groups
#index to resolve the correct row entry
def group_indecies(C):
    c_i,i_c,i = {},{},0
    for g in sorted(C.keys()):
        c_i[g]=i
        for j in range(i,i+len(C[g])): i_c[j]=g
        i+=len(C[g])
    return c_i,i_c

#return a list of indecies into X, given an edge in G
def get_x_i(A):
    return [A[g][A[g].keys()[0]][0] for g in A]

#A or active edges gets put into here, need to see what A is doing
#and be able to use the selected groups here instead
#::::TO DO... set pileup to avergae or don't use gamma value =0.0
def get_weight(A,E,average=0):
    if E.has_key((None,)): #default average weighting when no target is available
        if average>0:        #to do is to connect the independant expectations here
            c = np.float128(1.0)
            n = c/np.float128(average)
            w  = c-np.prod([c-n*np.float128((A[g][A[g].keys()[0]][1])) for g in A])
        else:
            w = np.float128(0.0) #0.0
    else:
        w = E[tuple(sorted([a[0] for a in A]))]
    return w
    
#active edges A, next vertext edges B
#clear off terminal edges of A and add B
def del_edges(A):
    k = A.keys()
    for i in range(len(k)): #take off the keys that have a 1
        if A[k[i]].has_key(1):
            #print('pop:%s'%k[i])
            A.pop(k[i])  

#active edges A, next vertext edges B
#clear off terminal edges of A and add B
def add_edges(A,B): 
    for k in B:
        if A.has_key(k):
            for i in B[k]: #put the new keys of E into A
                #print('push:%s,%s'%(k,i))
                A[k][i] = B[k][i]
        else:
            #print('push:%s'%(k))
            A[k] = {i:B[k][i] for i in B[k]}
    
#given a pileup graph g scan the sorted indecies
#added the dynamic table lookup W
def apply_weight_graph(X,G,c_i,i_c,EW,average=0):
    P,V,A,B,D = [],[],{},{},set([])
    if len(X)>0:
        t = X[0][2] #get the type
        V = sorted(G.keys())  #V are the sorted vertices (x1,x2 values)
        for i in range(0,len(V)-1):     #scan i-1,i,i+1 up to padding
            B = G[V[i]]                 #get edges for v[i+1]
            D = {m for l in [G[V[i]][k] for k in G[V[i]]] for m in l}  #check edge direction for i
            if len(A) <= 0:         #[1] len(a) <= 0 (f has 0 edge)   #section starting, start new  p+=[]
                #section is starting
                add_edges(A,B)
                x_i =  get_x_i(A)
                w   =  get_weight(A,EW,average)
                P += [[np.uint32(V[i]),np.uint32(V[i]),t,join_y(X,x_i,3),w,w,join_idx(X,x_i,6)]]
                #check for singles and new section together   
                if D == set([0,1]): #check for singles
                    #print('vertex i=%s'%i)
                    del_edges(A)    #clean the singles
                    if len(A)>0:    #with new sections that are not singles
                        x_i =  get_x_i(A)
                        w   =  get_weight(A,EW,average)
                        P += [[np.uint32(V[i]+1),np.uint32(V[i]+1),t,join_y(X,x_i,3),w,w,join_idx(X,x_i,6)]]
            else:
                if D == set([0]):   #[2] len(a) > 0 and f has 0 edge  #close subsection, start new  p[-1],p+=[]
                    #close the last open section and merge
                    x_i =  get_x_i(A)                                
                    P[-1][1] = np.uint32(V[i]-1)                     
                    P[-1][3] = merge_y(P[-1][3],join_y(X,x_i,3))     
                    P[-1][6] = merge_idx(P[-1][6],join_idx(X,x_i,6)) 
                    #clean the closed edges and start the new section
                    del_edges(A)
                    add_edges(A,B)                                
                    x_i =  get_x_i(A)
                    w   =  get_weight(A,EW,average)
                    P += [[np.uint32(V[i]),np.uint32(V[i]),t,join_y(X,x_i,3),w,w,join_idx(X,x_i,6)]]
                if D == set([0,1]): #[3] len(a) > 0 and f has 0 and 1 #close subsection, set single p[-1],p+=[]
                    #close the last open section
                    x_i =  get_x_i(A)                                
                    P[-1][1] = np.uint32(V[i]-1)                     
                    P[-1][3] = merge_y(P[-1][3],join_y(X,x_i,3))     
                    P[-1][6] = merge_idx(P[-1][6],join_idx(X,x_i,6))
                    #start the single section
                    add_edges(A,B)                                
                    x_i =  get_x_i(A)
                    w   =  get_weight(A,EW,average)
                    P += [[np.uint32(V[i]),np.uint32(V[i]),t,join_y(X,x_i,3),w,w,join_idx(X,x_i,6)]]
                    #add the new section past the single section
                    del_edges(A)
                    if len(A)>0:                                
                        x_i =  get_x_i(A)
                        w   =  get_weight(A,EW,average)
                        P += [[np.uint32(V[i]+1),np.uint32(V[i]+1),t,join_y(X,x_i,3),w,w,join_idx(X,x_i,6)]]
                if D == set([1]):   #[4] len(a) > 0 and f has 1       #section closing,  fix last   p[-1]
                    #close and clean the last section
                    x_i =  get_x_i(A)                                
                    P[-1][1] = np.uint32(V[i])                     
                    P[-1][3] = merge_y(P[-1][3],join_y(X,x_i,3))     
                    P[-1][6] = merge_idx(P[-1][6],join_idx(X,x_i,6))
                    add_edges(A,B)
                    del_edges(A)
                    #check to see if any open section remain
                    if len(A)>0:
                        x_i =  get_x_i(A)
                        w   =  get_weight(A,EW,average)
                        P += [[np.uint32(V[i]+1),np.uint32(V[i]+1),t,join_y(X,x_i,3),w,w,join_idx(X,x_i,6)]]
    for p in P: #clean up dangling sections
        if p[1]<p[0]: P.remove(p)    
    return P   
#pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp

#need to check average behavior...
#given the partition P[t][b][c][s]
#this will work without the key k 
def pileup_group_by_sample(P,E,true_key=(0,),average=0):
    A = {}
    for t in P:
        A[t] = {}
        for b in P[t]: #select the best groups for each type and bin or use average
        
            # DEBUG
#             print "T: "+str(t)
#             print "B: "+str(b)
            
            if b in E[t]:
                A[t][b] = {}
                G = set([i if len(i)<=1 else () for i in E[t][b]]).difference(set([()])) #best group using E[t][b]
                #swap out P[t][b][c][s] => P[t][b][s][c]
                #S = list(set([i for j in [P[t][b][c].keys() for c in P[t][b]] for i in j]))
                Q = {}
                for c in P[t][b]:
                    for s in P[t][b][c]:
                        if Q.has_key(s): Q[s][(c,)] = P[t][b][c][s] #the svult
                        else:            Q[s] = {(c,):P[t][b][c][s]}
                #swap out P[t][b][c][s] => P[t][b][s][c]       
                for s in Q.keys(): #now apply to each sample, if calls are avaible, otherwise use the average
                    GS = copy.deepcopy(G)             #don't use the target in this calculation      
                    if GS.issuperset(set([(None,)])): GS = set(Q[s]).difference((set([true_key,()])))
                    GS = list(GS)
                    C = {}
                    for i in GS:
                        if Q[s].has_key(i): C[i] = copy.deepcopy(Q[s][i])
                    X,W,c_i,i_c = weight_graph(C)
                    A[t][b][s] = apply_weight_graph(X,W,c_i,i_c,E[t][b],average)
            else:
                A[t][b][s] = []
    return A

#setup partitioned target by sample
def target_by_sample(P,k):
    T = {}
    for t in P:
        T[t] = {}
        for b in P[t]:
            T[t][b] = {}
            if P[t][b].has_key(k):
                for s in P[t][b][k]:
                    T[t][b][s] = P[t][b][k][s]
            #not sure what to do when a sample doesn't have a target key
    return T

#a good starting point for alpha is the maximum performing single caller cuttoff
def get_alpha(E):
    alpha = {}
    for t in E:
        alpha[t] = {}
        for b in E[t]:
            alpha[t][b] = max([1.0 if g[0] is None else E[t][b][g] if len(g)<=1 else 0.0 for g in E[t][b]])
    return alpha

#apply a static filter on everything
def target_filter_cuttoff_static(B,static_filter):
    alpha = {}
    for t in B:
        alpha[t] = {}
        for b in B[t]:
            alpha[t][b] = np.float128(static_filter)
    return alpha

#:::TO DO::: do the search but not inside each bin
def target_filter_cuttoff_static_exhaustive(A,E,T):
    alpha = {}
    for t in A:
        alpha[t] = {}
        for b in A[t]:
            alpha[t][b] = np.float128(1.0)
    return alpha

#when gamm = 0.0 will search all combinations
def target_filter_cutoff_exhaustive(A,E,T):
    alpha = {}
    for t in A:
        alpha[t] = {}
        for b in A[t]:
            alpha[t][b] = np.float128(1.0)
            J = {k:[np.uint64(0),np.uint64(0)] for k in sorted(list(set([E[t][b][z] for z in E[t][b]])))}
            if len(J)>1:
                for s in A[t][b]:
                    for k in J:
                        if T[t][b].has_key(s):
                            N = fu.feature_magnitudes(T[t][b][s],filter_single_bin(A[t][b][s],k))
                            J[k][0],J[k][1] = J[k][0]+np.uint64(N[0]),J[k][1]+np.uint64(N[1])
                for k in J:
                    if J[k][1]>0.0: J[k] = float(np.float128(J[k][0])/np.float128(J[k][1]))
                    else:           J[k] = 0.0 
                alpha[t][b] = sorted([[k,J[k]] for k in J],key=lambda x: x[1], reverse=True)[0][0]
            elif len(J)==1 and J.keys()[0] > 0.0:
                for s in A[t][b]:
                    for k in J:
                        if T[t][b].has_key(s):
                            N = fu.feature_magnitudes(T[t][b][s],filter_single_bin(A[t][b][s],k))
                            J[k][0],J[k][1] = J[k][0]+np.uint64(N[0]),J[k][1]+np.uint64(N[1])
                for k in J:
                    if J[k][1]>0.0: J[k] = float(np.float128(J[k][0])/np.float128(J[k][1]))
                    else:           J[k] = 0.0 
                alpha[t][b] = sorted([[k,J[k]] for k in J],key=lambda x: x[1], reverse=True)[0][0]
    return alpha

#when gamm = 0.0 will search all combinations
def target_filter_cutoff_exhaustive_brkpt(A,P,T,E,K):
    alpha = {}
    I,M = get_all_call_index(P),max_brkpt_stats(K)
    for t in A:
        alpha[t] = {}
        for b in A[t]:
            alpha[t][b] = np.float128(1.0)
            J = {p:[np.uint64(0),np.uint64(0)] for p in sorted(list(set([E[t][b][z] for z in E[t][b]])))}
            if len(J)>1:
                for s in A[t][b]:
                    for p in J:
                        if T[t][b].has_key(s):
                            F = {t:{b:{s:filter_single_bin(A[t][b][s],p)}}} #cut by expectation
                            S = best_smooth_brkpt_samples_search(F,K,P,I,M)
                            N = fu.feature_magnitudes(T[t][b][s],S[t][b][s])
                            J[p][0],J[p][1] = J[p][0]+np.uint64(N[0]),J[p][1]+np.uint64(N[1])
                for p in J:
                    if J[p][1]>0.0: J[p] = float(np.float128(J[p][0])/np.float128(J[p][1]))
                    else:           J[p] = 0.0 
                alpha[t][b] = sorted([[p,J[p]] for p in J],key=lambda x: x[1], reverse=True)[0][0]
            elif len(J)==1 and J.keys()[0] > 0.0:
                for s in A[t][b]:
                    for p in J:
                        if T[t][b].has_key(s):
                            F = {t:{b:filter_single_bin(A[t][b][s],p)}} #cut by expectation
                            S = best_smooth_brkpt_samples_search(F,K,P,I,M)
                            N = fu.feature_magnitudes(T[t][b][s],S[t][b][s])
                            J[p][0],J[p][1] = J[p][0]+np.uint64(N[0]),J[p][1]+np.uint64(N[1])
                for p in J:
                    if J[p][1]>0.0: J[p] = float(np.float128(J[p][0])/np.float128(J[p][1]))
                    else:           J[p] = 0.0 
                alpha[t][b] = sorted([[p,J[p]] for p in J],key=lambda x: x[1], reverse=True)[0][0]
    return alpha

def optimal_groups(E,alpha):
    G = {}
    for t in E:
        G[t] = {}
        for b in E[t]:
            o = set([])
            for g in E[t][b]:
                if E[t][b][g] >= alpha[t][b]:
                    for e in g: o.add(e)
            G[t][b] = tuple(sorted(list(o)))
    return G
    
#EM algorithm that finds a good filter cutoff value for each slice given the target
#set a step size automatically based on the statistics like mean, std, etc
def target_filter_cutoff_simple(A,E,T,step=0.01,decay=0.9,zeta=1E-9):
    alpha = get_alpha(E) #inital is the maximal single caller per slice
    for t in A:
        for b in A[t]:
            #do the EM filter cuttoff for each type/bin-----------------------------
            J = alpha_samples_step(A,T,t,b,alpha,step) #initial search
            while J['-']>J['0'] or J['+']>J['0'] and abs(J['+']-J['-'])>zeta:
                if J['-']>J['0']: alpha[t][b] = max(0.0,alpha[t][b]-step)
                if J['+']>J['0']: alpha[t][b] = min(1.0,alpha[t][b]+step)                
                if step > zeta: step *= decay #progress ending here
                J = alpha_samples_step(A,T,t,b,alpha,step)
            print('final t:%s b:%s alpha:%s'%(t,b,alpha[t][b]))
    return alpha

#EM step scoring calculation for upper lower and current value
def alpha_samples_step(A,T,t,b,alpha,step):
    V = {'-':{0:np.float64(0.0),1:np.float64(0.0)},
         '0':{0:np.float64(0.0),1:np.float64(0.0)},
         '+':{0:np.float64(0.0),1:np.float64(0.0)}}
    J = {'-':np.float64(0.0),'0':np.float64(0.0),'+':np.float64(0.0)}
    for s in A[t][b]:
        if T[t][b].has_key(s):
            N = fu.feature_magnitudes(filter_single_bin(A[t][b][s],min(1.0,alpha[t][b]+step)),T[t][b][s])
            V['+'][0],V['+'][1] = V['+'][0]+np.float64(N[0]),V['+'][1]+np.float64(N[1]) 
            N = fu.feature_magnitudes(filter_single_bin(A[t][b][s],alpha[t][b]),T[t][b][s])
            V['0'][0],V['0'][1] = V['0'][0]+np.float64(N[0]),V['0'][1]+np.float64(N[1])
            N = fu.feature_magnitudes(filter_single_bin(A[t][b][s],max(0.0,alpha[t][b]-step)),T[t][b][s])
            V['-'][0],V['-'][1] = V['-'][0]+np.float64(N[0]),V['-'][1]+np.float64(N[1])
    for k in V:
        if V[k][1]>0.0: J[k] = V[k][0]/V[k][1]
    return J

#apply the filter with given alpha values fuse step
def filter_pileup_by_sample(A,alpha,E_post=None,wx=4,wy=5,idx=6,leave_in=True):
    F = {}
    if E_post is None:
        for t in A:
            F[t] = {}
            for b in A[t]:
                F[t][b] = {}
                for s in A[t][b]:
                    #F[t][b][s] = filter_bin(A[t][b][s],alpha[t][b],E_post,leave_in)
                    F[t][b][s] = filter_bin(A,alpha,t,b,s,E_post,wx,wy,idx,leave_in)
    else:
        for t in A:
            F[t] = {}
            for b in A[t]:
                F[t][b] = {}
                for s in A[t][b]:
                    #F[t][b][s] = filter_bin(A[t][b][s],alpha[t][b],E_post[t][b],leave_in)
                    F[t][b][s] = filter_bin(A,alpha,t,b,s,E_post,wx,wy,idx,leave_in)
    return F
    
#macro that does the type and bin iteration putting the optimal group
#and expectation combinations intot th apply weight graph computation
#this can be done in || for speedup...
#S[t][b][(c,)], E[t][b][g], remove the target key if needed...
def pileup_groups(S,E,k,average=0):
    P = {}
    for t in S:
        P[t] = {}
        for b in S[t]: #select the best groups for each type and bin or use average
            G = set([i if len(i)<=1 else () for i in E[t][b]]).difference(set([()]))
            if G.issuperset(set([(None,)])): G = set(S[t][b]).difference((set([(k,),()])))
            G = list(G)
            C = {i:copy.deepcopy(S[t][b][i]) for i in G} 
            X,G,c_i,i_c = weight_graph(C)
            P[t][b] = apply_weight_graph(X,G,c_i,i_c,E[t][b],average)
    return P
    
#given a pileup with applied weights, slice out the regions as merged svultb
def filter_pileup(A,alpha,t,b,s=None,E_post=None,wx=4,wy=5,idx=6,leave_in=False):
    #alpha = get_alpha(E) already will have this for fusion
    F = {}
    for t in A:
        F[t] = {}
        for b in A[t]:
            F[t][b] = filter_bin(A,alpha,t,b,s,E_post,wx,wy,idx,leave_in)
    return F

def filter_single_bin(a,f,leave_in=False):
    C = []
    if leave_in:
        for i in range(len(a)):
                if a[i][4]>=f and a[i][5]>=f:
                    C += [a[i]+[1]]
                else:
                    C += [a[i]+[-1]]
    else: #stadard case is here--------------
        for i in range(len(a)):
                if a[i][4]>=f and a[i][5]>=f:
                    C += [a[i]+[1]]
    return C

#filter a pileup call set a using cuttoff f
#this merge will now average weight by svlen contribution
def filter_bin(A,alpha,t,b,s=None,E_post=None,wx=4,wy=5,idx=6,leave_in=False):
    C = []
    if s is not None: a = A[t][b][s]
    else:             a = A[t][b]
    if leave_in:
        if E_post is None: 
            for i in range(len(a)):
                if a[i][wx]>=alpha[t][b] and a[i][wy]>=alpha[t][b]:
                    C += [a[i]+[1]]
                else:
                    C += [a[i]+[-1]]
        else:
            for i in range(len(a)):
                g = tuple(sorted(a[i][idx]))
                if E_post[t][b].has_key(g):
                    xw = yw = E_post[t][b][g]
                    if a[i][wx]>=alpha[t][b] and a[i][wy]>=alpha[t][b]:
                        C += [a[i][0:wx]+[xw,yw]+a[i][wy+1:]+[1]]
                else:
                    C += [a[i]+[0.0,0.0]+a[i]+[-1]]
    else:
        if E_post is None:
            for i in range(len(a)):
                if a[i][wx]>=alpha[t][b] and a[i][wy]>=alpha[t][b]:
                    C += [a[i]]
            C = fu.merge_regions(C)
        else:
            for i in range(len(a)):
                g = tuple(sorted(a[i][idx]))
                if E_post[t][b].has_key(g):
                    xw = yw = E_post[t][b][g]
                    if a[i][wx]>=alpha[t][b] and a[i][wy]>=alpha[t][b]:
                        C += [a[i][0:wx]+[xw,yw]+a[i][wy+1:]+[1]]
                    # SL: print SVult
                    #print ["Filter_bin", ';', t, ';', b, ';', s, ';', g, ';', [a[i][0:wx]+[xw,yw]+a[i][wy+1:]+[1]], ';', a[i][wx], ';', a[i][wy], ';', alpha[t][b], ';', a[i][wx]>=alpha[t][b], ';', a[i][wy]>=alpha[t][b]]
            C = fu.merge_regions(C)    
    return C
#given a piled and filtered set, merge back into one call set
def merge_filtered(F,check_merge=True):
    M = {}
    if check_merge:
        for t in F:
            M[t] = []
            for b in F[t]:
                M[t] = append_bin(M[t],F[t][b])
    else:
        for t in F:
            M[t] = []
            for b in F[t]:
                M[t] += F[t][b]
            M[t] = sorted(M[t],key=lambda x: x[0])        
    return M

def overlap(c_i0,c_i1,c_j0,c_j1,check_x=True,check_y=False):
    if c_i0<=c_i1: a,b = np.float64(c_i0),np.float64(c_i1)
    else:          a,b = np.float64(c_i1),np.float64(c_i0)    
    if c_j0<=c_j1: c,d = np.float64(c_j0),np.float64(c_j1)
    else:          c,d = np.float64(c_j1),np.float64(c_j0)
    i = abs(a-c)+abs(b-d)                           
    u = min((b-a+1)+(d-c+1),max(b,d)-min(a,b)+1) 
    return max(np.float64(0.0),(u-i)/u)

#this could be done using a merging weight graph
#this is the entry point for updating the E[g] values?
def append_bin(C1,C2,cutoff=0.0):
    A,B,C3 = {},{},[]
    n,m = len(C1),len(C2)
    if m > 0:
        for i in range(n): #find the overlapping indexcies of C2 to C1
            for j in range(m):
                c = overlap(int(C1[i][0]),int(C1[i][1]),int(C2[j][0]),int(C2[j][1]))
                if c > cutoff:
                    if A.has_key(j): A[j] += [i]
                    else:            A[j]  = [i]                    
                    if B.has_key(i): B[i] += [j]
                    else:            B[i]  = [j]
        S = sorted(list(set(range(m)).difference(set(A.keys()))))#no overlap of C2  
        T = sorted(list(set(range(n)).difference(set(B.keys()))))#no overlap of C1
        C3 = [C1[i] for i in T]+[C2[j] for j in S]
        for j in A: #now resolve the overlap by highest weight
            if all([True if C2[j][4]>C1[i][4] else False for i in A[j]]):  C3 += [C2[j]]
        for i in B:
            if all([True if C1[i][4]>=C2[j][4] else False for j in B[i]]): C3 += [C1[i]]
        C3 = fu.merge_regions(sorted(C3,key=lambda x: x[0])) #sort to use merge_regions
    else:
        C3 = C1
    return C3   
    
#given a piled and filtered set, merge back into one call set
def merge_filtered_sample(F,s,cutoff=0.0):
    M = {}
    if cutoff>0.0:
        for t in F:
            M[t] = []
            for b in F[t]:
                if F[t][b].has_key(s):
                    M[t] = append_bin(M[t],F[t][b][s],cutoff)
    else:
        for t in F:
            M[t] = []
            for b in F[t]:
                if F[t][b].has_key(s):
                    M[t] += F[t][b][s]
            M[t] = sorted(M[t],key=lambda x: x[0])        
    return M

#put a call set into Q:  Q[t][c][s] <= C[t][s]
def merge_caller_samples(Q,C,c,snames):
    for t in Q:
        Q[t][c] = {}
        for s in snames:
            Q[t][c][s] = C[t][s]

#given the original Q[t][c][s] insert a merged F[t][b][s]
def merge_filtered_samples(Q,F,c,snames,exclude,cutoff):
    N = {}
    for s in snames:
        if s not in exclude:
            N[s] = merge_filtered_sample(F,s,cutoff)
    for t in Q:
        Q[t][c] = {}
        for s in snames:
            if s not in exclude:
                Q[t][c][s] = N[s][t]

#given a Q[t][c][s], do a pooled overal score (no bins...)
def score_merged(Q,k,snames):
    S = {}
    for t in Q:
        for i in set(Q[t].keys()).difference(set([k])):
            if not S.has_key(i): S[i] = {}           
            I,U = 0.0,0.0
            S[i][t] = 0.0
            if Q[t].has_key(k):
                for s in snames:
                    C1,C2 = [],[]
                    if Q[t][i].has_key(s): C1 = Q[t][i][s]
                    if Q[t][k].has_key(s): C2 = Q[t][k][s]                
                    N = fu.feature_magnitudes(C1,C2,self_merge=True)
                    I += N[0]
                    U += N[1]
                if U > 0.0: S[i][t] = I/U
            else:
                for s in snames:
                    if Q[t].has_key(i):
                        S[i][t] = 0.0
                        break
    return S

#given the original unbinned call set S and target key k
#with the binned,optimized,filtered,merged fusion call set F
#calculate the overal score for each type
def score_sample(S,k,F,self_merge=True):
    C = set(S.keys()).difference(set([k]))
    T = {}
    for c in C: #each caller here
        T[c] = {}
        for t in S[k]:
            T[c][t] = 0.0
            if S[c].has_key(t):
                T[c][t] = fu.jaccard_score(S[k][t],S[c][t],self_merge)
    T[-1] = {}
    for t in S[k]:
        T[-1][t] = 0.0
        if F.has_key(t):
            T[-1][t] = fu.jaccard_score(S[k][t],F[t],self_merge)
    return T

#temperary helper function for viewing
def pretty_stats(Q,types,t,k,c_id,sname,r=0.5,verbose=True):
    C1,C2 = [],[]
    if Q[t].has_key(k) and Q[t][k].has_key(sname):       C1 = Q[t][k][sname]
    if Q[t].has_key(c_id) and Q[t][c_id].has_key(sname): C2 = Q[t][c_id][sname]            
    s = fu.metric_score(C1,C2,r,self_merge=False,d1_ids=True,d2_ids=True)
    q = fu.metric_score(C1,C2,r,self_merge=True,d1_ids=False,d2_ids=False)
    d2_ids,ids = s['d2_ids'],[]               #now get the idx indecies from C1
    for i in d2_ids: ids += [i] #and make a list of idx entries
    prec,rec,f1,j = s['n:m'],s['m:n'],0.0,0.0
    L,R = [],[]
    for brkpt in s['b']:
        if not brkpt is None:
            L,R = L+[brkpt[0]],R+[brkpt[1]]
    l_mu,l_sd,r_mu,r_sd = 0.0,0.0,0.0,0.0
    if len(L)>0: l_mu,l_sd = np.mean(L),np.std(L)
    if len(R)>0: r_mu,r_sd = np.mean(R),np.std(R)
    if s['n:m']+s['n:m']>0.0: f1 = 2.0*prec*rec/(prec+rec)
    if q['j'][1]>0.0: j = q['j'][0]/q['j'][1]            
    if verbose: print('sample=%s\tt=%s\tprec=%s\trec=%s\tf1=%s\tj=%s'%(sname,types[t],prec,rec,f1,j))
    return [sname,types[t],prec,rec,f1,j,ids,len(C1),len(C2),l_mu,r_mu,l_sd,r_sd]
    #0=sname, 1=prec, 2=rec, 3=f1, 4=j, 5=ids, 6=n, 7=m, 8=l_mu, 9=r_mu, 10=l_sd, 11=r_sd

#repartitioning is easiest and then return all the bins
def pretty_bin_stats(Q,types,t,B,bins,k,c_id,sname,r=0.5,verbose=True):
    P = partition_sliced_samples(Q,B,exclude=[])
    D = {}
    for i in range(len(B[t])-1):
        C1,C2 = [],[]
        if P[t][i].has_key(k) and P[t][i][k].has_key(sname):       C1 = P[t][i][k][sname]
        if P[t][i].has_key(c_id) and P[t][i][c_id].has_key(sname): C2 = P[t][i][c_id][sname]
        s = fu.metric_score(C1,C2,r,self_merge=False,d1_ids=True,d2_ids=True)
        q = fu.metric_score(C1,C2,r,self_merge=True,d1_ids=False,d2_ids=False)
        d1_ids,d2_ids = s['d1_ids'],s['d2_ids']   #now get the idx indecies from C1
        prec,rec,f1,j = s['n:m'],s['m:n'],0.0,0.0
        L,R = [],[]
        for brkpt in s['b']:
            if not brkpt is None:
                L,R = L+[brkpt[0]],R+[brkpt[1]]
        if s['n:m']+s['n:m']>0.0: f1 = 2.0*prec*rec/(prec+rec)
        if q['j'][1]>0.0: j = q['j'][0]/q['j'][1]
        if verbose: print('sample=%s\tt=%s\tbin=%s\tprec=%s\trec=%s\tf1=%s\tj=%s'%(sname,types[t],bins[t][i],prec,rec,f1,j))
        D[bins[t][i]] = [sname,types[t],bins[t][i],prec,rec,f1,j,len(C1),len(C2),d1_ids,d2_ids,L,R]
    return D
    
def assemble_stats(L):
    cross_fold_stats,detailed_stats,hist = {},{},{}
    for i in L:
        for c in i[1]:
            if not cross_fold_stats.has_key(c):
                cross_fold_stats[c] = {}
            for t in i[1][c]:
                if not cross_fold_stats[c].has_key(t):
                    cross_fold_stats[c][t] = []
                cross_fold_stats[c][t] += i[1][c][t]
    for i in L:
        for s in i[2]:
            hist[s] = copy.deepcopy(i[2][s])
    for i in L:
        for c in i[3]:
            if not detailed_stats.has_key(c):
                detailed_stats[c] = {}
            for t in i[3][c]:
                if not detailed_stats[c].has_key(t):
                    detailed_stats[c][t] = {}
                for b in i[3][c][t]:
                    if not detailed_stats[c][t].has_key(b):
                        detailed_stats[c][t][b] = []
                    detailed_stats[c][t][b] += i[3][c][t][b]
    return cross_fold_stats,hist,detailed_stats
    
