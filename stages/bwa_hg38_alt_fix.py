import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class bwa_index(stage_wrapper.Stage_Wrapper):
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
        #[1b]get some metadata for I/O names
        
        #[2]build command args
        samtools = self.software_path+'/samtools-1.3/samtools'
        view = [samtools, 'view', '-Sh', inputs['.bam']]

        path = self.software_path+'/bwa.kit'
        k8 = path + '/k8'
        postalt   = path + '/bwa-postalt.js'
        alt_files = path + '/resource-GRCh38/hs38DH-extra.fa ' 
        alt_files += path + '/resource-GRCh38/hs38DH.fa.alt'
        alt_fix = [k8, postalt, '-p', alt_files]

        out_file = inputs['out_file']
        if out_file is None:
            out_file = self.strip_in_ext(inputs['.bam'],'.bam') + '.alt.bam'
        view2 = [samtools, 'view', '-1', '-', '-o', out_file]
        
        #[1a]make start entry which is a new staged_run row
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(view + ['|'] + alt_fix + ['|'] + view2))
            output = subprocess.check_output(' '.join(view + ['|'] + alt_fix + ['|'] + view2),stderr=subprocess.STDOUT, shell=True)
        #catch all errors that arise under normal behavior
        except subprocess.CalledProcessError as E:
            print('call error: '+E.output)             #what you would see in the term
            err['output'] = E.output
            #the python exception issues (shouldn't have any...
            print('message: '+E.message)          #?? empty
            err['message'] = E.message
            #return codes used for failure....
            print('code: '+str(E.returncode))     #return 1 for a fail in art?
            err['code'] = E.returncode
        except OSError as E:
            print('os error: '+E.strerror)             #what you would see in the term
            err['output'] = E.strerror
            #the python exception issues (shouldn't have any...
            print('message: '+E.message)          #?? empty
            err['message'] = E.message
            #the error num
            print('code: '+str(E.errno))
            err['code'] = E.errno
        print('output:\n'+output)
        
        #[3b]check results--------------------------------------------------
        if err == {}:
            #self.db_stop(run_id,{'output':output},'',True)
            results = [out_file]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("<<<<<<<<<<<<<bwa index sucessfull>>>>>>>>>>>>>>>\n")
                return results
            else:
                print("<<<<<<<<<<<<<bwa index failure>>>>>>>>>>>>>>>\n")
                return None
        else:
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
