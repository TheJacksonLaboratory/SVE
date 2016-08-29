import os
import sys
import pysam
import subprocess32 as subprocess
import itertools as it
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class bam_split_all(stage_wrapper.Stage_Wrapper):
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
        #[1a]get input names and output names setup
        in_names  = {'.bam':inputs['.bam'][0],'chroms':inputs['chroms']}
        out_ext = self.split_out_exts()[0]
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['.bam'],'.bam'))
            out_names = {'.bam' : out_dir+stripped_name}
        else:
            cascade = self.strip_in_ext(in_names['.bam'],'.bam')
            out_names = {'.bam' :cascade}  
        
        H = [out_names['.bam']+'.chr'+str(i)+'.bam' for i in in_names['chroms']+['UM']] #UM for unmapped
        K = {}
        for i in in_names['chroms']+['*']:
            a = i
            if i == '*': a = 'UM'
            else:
                K[(str(i),str(i))] = [out_names['.bam']+'.chr'+str(a)+'.bam'] #use UM instead of *
        
        for i,j in sorted(it.combinations(in_names['chroms']+['*'],2),key=lambda x: x[0]):
            a,b = i,j
            if i == '*': a = 'UM'
            if j == '*': b = 'UM'
            K[(i,j)] = [out_names['.bam']+'.chr'+str(a)+'.bam',out_names['.bam']+'.chr'+str(b)+'.bam']
        print('processing bam files %s'%(H,))
        self.db_start(run_id,in_names['.bam'])
        output,err = '',{}
        #[3a]execute the command here----------------------------------------------------
        try:
            inbam = pysam.AlignmentFile(in_names['.bam'], "rb")
            header = inbam.header #copy out the header to each
            refs   = list(inbam.references)
            rmap = {k:refs[k] for k in range(len(refs))}
            rmap[-1] = '*' #unmapped bin
            #open output bams for writing
            outbams = {h:pysam.AlignmentFile(h,'wb',header=header) for h in H}
            bam_iter = inbam.fetch(until_eof=True,multiple_iterators=False) 
            for read in bam_iter:   
                i = tuple(sorted([rmap[read.rname],rmap[read.rnext]]))
                if K.has_key(i):
                    for j in K[i]: #duplicate the splits to both bams outputs
                        outbams[j].write(read)
            inbam.close()
            for h in outbams:
                outbams[h].close()
            print('finished writing split bams')
            samtools = self.software_path+'/samtools-1.3/samtools' #should still be in sorted order
#            sort  = [samtools,'sort','-T','sorting','%s','-O',"'bam'",'-o','%s']
            index = [samtools,'index','%s']
            for h in outbams:
#                print('sorting %s'%h)
#                output += subprocess.check_output(' '.join(sort)%(h,h),
#                                                  stderr=subprocess.STDOUT,shell=True)
                print('indexing %s'%h)
                output += subprocess.check_output(' '.join(index)%h,
                                                  stderr=subprocess.STDOUT,shell=True)
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
        finally: #have to work this out with some tests....
            print('closing files')
            inbam.close()
            for h in outbams:
                outbams[h].close()  
        print('output:\n'+output)
        self.command = ' '
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