import os
import sys
import time
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class mrfast_index(stage_wrapper.Stage_Wrapper):
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
    #~/software/svtoolkit/lib/...
    def run(self,run_id,inputs):
        #workflow is to run through the stage correctly and then check for error handles
             
        #[1a]get input names and output names setup
        in_name = {'.fa':inputs['.fa'][0]}
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        out_name = {'.fa.fai' :self.strip_in_ext(in_name['.fa'],'.fa')+out_exts[0],
                    '.fa.index' :self.strip_in_ext(in_name['.fa'],'.fa')+out_exts[1]}
        
        defaults,params = self.params,[]
        params = [k+' '+str(defaults[k]['value']) for k in defaults]#wspace delimited
        
        #[a]use to run several sub scripts via command line/seperate process
        mrfast = self.software_path+'/mrfast-2.6.1.0/mrfast'
        mrfast_i = [mrfast,'--index',in_name['.fa']]+params
        samtools = self.software_path+'/samtools-1.0/bin/samtools'
        samtools_i = [samtools,'faidx',in_name['.fa']]
        
        self.db_start(run_id,in_name['.fa']) 
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print(' '.join(mrfast_i))
            output += subprocess.check_output(' '.join(mrfast_i),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
            print(' '.join(samtools_i))
            output += subprocess.check_output(' '.join(samtools_i),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
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
            print('vcf write os/file IO error')
            err['output'] = 'vcf write os/file IO error'
            err['message'] = 'vcf write os/file IO error'
            err['code'] = 1
        print('output:\n'+output)
                                                
        #[3b]check results--------------------------------------------------
        if err == {}:
            self.db_stop(run_id,{'output':output},'',True)
            results = [out_name['.fa.index']]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("sucessfull........")
                return results   #return a list of names
            else:
                print("failure...........")
                return False
        else:
            print("failure...........")
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None