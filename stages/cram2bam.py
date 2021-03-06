import os
import sys
import subprocess32 as subprocess
import itertools as it
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class bam2cram(stage_wrapper.Stage_Wrapper):
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
    #bwa sampe ref.fa L.sai R.sai L.fq R.fq -f out.sam
    def run(self,run_id,inputs):
        #workflow is to run through the stage correctly and then check for error handles
    
        #[1a]get input names and output names setup
        in_names  = {'.fa':inputs['.fa'][0],'.cram':inputs['.cram'][0]}
        out_ext = self.split_out_exts()[0]
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['.cram'],'.cram'))
            out_names = {'.cram.bam' : out_dir+stripped_name+'.cram.bam'}
        else:
            cascade = self.strip_in_ext(in_names['.cram'],'.cram')
            out_names = {'.cram.bam' :cascade+'.cram.bam'}  
        
        #[2a]build command args
        #single bam2cram here, extend with multiple version later
        samtools = self.software_path+'/samtools-1.3/samtools'
        #cut out a headerless sam file that includes anything to do with chr A
        bam  =  [samtools,'view','-T', in_names['.fa'],'-bh',in_names['.cram'],'-o',out_names['.cram.bam']]
        index = [samtools,'index',out_names['.cram.bam']]
        #self.db_start(run_id,in_names['.cram'])
        output,err = '',{}
            #[3a]execute the command here----------------------------------------------------
        try:
            output += subprocess.check_output(' '.join(bam),stderr=subprocess.STDOUT,shell=True)
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
        print('output:\n'+output)
        self.command = ' '.join(bam)
        #print(self.get_command_str())
        #[3b]check results--------------------------------------------------
        if err == {}:
            #self.db_stop(run_id,{'output':output},'',True)
            results = [out_names['.cram.bam']] #this needs to have all .bam files
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("<<<<<<<<<<<<<cram2bamsucessfull>>>>>>>>>>>>>>>\n")
                return results   #return a list of names
            else:
                print("<<<<<<<<<<<<<cram2bam failure>>>>>>>>>>>>>>>\n")
                return False
        else:
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
