import os
import sys
import pysam
import subprocess32 as subprocess
import itertools as it
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class cram2bam_split_all(stage_wrapper.Stage_Wrapper):
    #path will be where a node should process the data using the in_ext, out_ext
    #stage_id should be pre-registered with db, set to None will require getting
    #a new stage_id from the  db by writing and registering it in the stages table
    def __init__(self,wrapper,dbc,retrieve,upload,params):
        #inheritance of base class stage_wrapper    
        stage_wrapper.Stage_Wrapper.__init__(self,wrapper,dbc,retrieve,upload,params)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        return 0  
    
    #override this function in each wrapper...
    def run(self,run_id,inputs):
        #workflow is to run through the stage correctly and then check for error handles
    
        #[1a]get input names and output names setup
        in_names  = {'.fa':inputs['.fa'][0],'.cram':inputs['.cram'][0],'chroms':inputs['chroms']}
        out_ext = self.split_out_exts()[0]
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['.cram'],'.cram'))
            out_names = {'.cram.bam' : out_dir+stripped_name}
        else:
            cascade = self.strip_in_ext(in_names['.cram'],'.cram')
            out_names = {'.cram.bam' :cascade}  

        #load and print params
        defaults,params = self.params,[]
        params = [k+'='+str(defaults[k]['value']) for k in defaults]        
        print(params)
        p = defaults['p']['value']
        
        #-----------------------------------
#        p = 1
#        chroms = [str(i) for i in range(1,6)]
#        directory = '/home/tbecker/data/test/'
#        in_names  = {'chroms':chroms,'.fa':'/home/tbecker/data/human_g1k_v37_decoy.fa','.cram':'/home/tbecker/data/crams/NA12878.low_cov.cram'}
#        if not os.path.exists(directory): os.makedirs(directory)
#        out_names = {'.cram.bam':'/home/tbecker/data/test/NA12878'}         
        #---------------------------------------------------------
        
        #formulate chr-chr pairs and partitions into p bam file streams lower bound is 2
        #lower bound for p is 2, upper bound is c + choose(c+1,2) for c chroms
        #take each chrA->chrA and then add in some chrA->chrB to chrA partition
        #these are stored in sorted pair ordering, but have to be checked in bth directions
        single = [(str(i),'=') for i in in_names['chroms']]
        splits = [tuple(sorted([i,j])) for i,j in sorted(it.combinations(in_names['chroms']+['*'],2),key=lambda x: x[0])]
        print('single chroms:') 
        print(single)
        print('split pair chroms:')
        print(splits)
        #formulate the partitions
        handles = [{'chr'+s[0]+'-chr'+s[0]:[s]} for s in single] #this is the first set,key is the name
        n = len(splits)/p #number of split files
        m = len(splits)%p
        for i in range(p):
            split = splits[i*n:(i+1)*n]
            if i >= n-1 and m > 0: split += [splits[(i+1)*n:]]
            handles += [{'split'+str(i):split}]
        #reverse the keys for line by line lookups back into the handle
        K,L = {},[]
        for h in handles:
            P = h.keys()[0]          #single handle here
            L += [P]
            for k in h[P]: K[k] = out_names['.cram.bam']+'.'+P+'.cram.bam'
        #K[tuple(sorted([s,t])) for each s=$3 and t=$7 from each line
        H = [out_names['.cram.bam']+'.'+l+'.sam' for l in L]
        G = [out_names['.cram.bam']+'.'+l+'.cram.bam' for l in L]
        print('total split files:')
        print(H)
        #[2a]build command args
#        samtools = '/home/tbecker/software/samtools-1.2/samtools'
        samtools = self.software_path+'/samtools-1.2/samtools'       
        cram2sam  = [samtools,'view','-T',in_names['.fa'],'-Sh', in_names['.cram']] #read from in
        sam2bam   = [samtools,'view','-bh','%s','-o','%s']

        self.db_start(run_id,in_names['.cram'])
        output,err = '',{}
        #[3a]execute the command here----------------------------------------------------
        try:
            #read each line of one in file, write to many output handles
            samtoolsin  = subprocess.Popen(cram2sam,stdout=subprocess.PIPE)
            savesams    = {h:open(h,'w') for h in H}
            bams    = {H[i]:G[i] for i in range(len(H))}
            #read each line and write out to the mapped handle
            for line in samtoolsin.stdout:
                if(line.startswith("@")):
                    for h in H:
                        savesams[h].write(line+'\n')
                else:
                    linesplit = line.split('\t') #alignment reacords are length 11
                    i = tuple(sorted([linesplit[2],linesplit[6]])) #get sorted tuple refs for pairs
                    if savesams.has_key(K[i]):
                        savesams[K[i]].write(line+'\n')
            for h in savesams:
                savesams[h].close()
            for h in bams:
                output += subprocess.check_output(' '.join(sam2bam)%(h,bams[h]),
                                                  stderr=subprocess.STDOUT, shell=True)
            #put back the headers for each bam file
#            for h in H:
#                output += subprocess.check_output(reheader+[h],stderr=subprocess.STDOUT)
        #catch all errors that arise under normal call behavior
        except subprocess.CalledProcessError as E:
            print('call error: '+E.output)        #what you would see in the term
            err['output'] = E.output
            #the python exception issues (shouldn't have any...
            print('message: '+E.message)          #?? empty
            err['message'] = E.message
            #return codes used for failure....
            print('code: '+str(E.returncode))     #return 1 for a fail in art?
            err['code'] = E.returncode
        except OSError as E:
            print('os error: '+E.strerror)        #what you would see in the term
            err['output'] = E.strerror
            #the python exception issues (shouldn't have any...
            print('message: '+E.message)          #?? empty
            err['message'] = E.message
            #the error num
            print('code: '+str(E.errno))
            err['code'] = E.errno
        except Exception as E:
            print('some other error:'+' '.join(E.message))
#        finally:
#            print('closing files')
#            for h in savesams:
#                savesams[h].close()
                
        print('output:\n'+output)
        self.command = ' '.join(cram2sam)+' '.join(sam2bam)
        #print(self.get_command_str())
        #[3b]check results--------------------------------------------------
        if err == {}:
            self.db_stop(run_id,{'output':output},'',True)
            results = H #this needs to have all .bam files
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("sucessfull........")
                return results   #return a list of names
            else:
                print("failure...........")
                return False
        else:
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None