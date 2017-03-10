import os
import sys
import time
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper
import stage_utils as su

#function for auto-making svedb stage entries and returning the stage_id
class bwa_mem(stage_wrapper.Stage_Wrapper):
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
    
        #[1b]
        in_names  = {'.fa':inputs['.fa'][0],'.fq':inputs['.fq']}
        out_ext = self.split_out_exts()[0]
        #add tigra.ctg.fa bam by allowing a single '.fa' input key
        
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(su.get_common_string_left(in_names['.fq']))
            out_names = {'.sam' : out_dir+stripped_name+out_ext}
        else: #untested...
            right = in_names['.fa'][0].rsplit('/')[-1]
            out_dir = in_names['.fq'][0].replace(right,'')
            cascade = self.strip_in_ext(in_names['.fq'][0],'.fq')
            out_names = {'.sam' :cascade+out_ext}
        #[2]build command args
        if not os.path.exists(out_dir): os.makedirs(out_dir)
        if inputs.has_key('SM'):
            SM = inputs['SM'][0]
        else:
            SM = stripped_name
        threads = str(self.get_params()['-t']['value'])
        bwa = self.software_path+'/bwa-master/bwa'
        sample = stripped_name+'RG'
        #'@RG\tID:H7AGF.2\tLB:Solexa-206008\tPL:illumina\tPU:H7AGFADXX131213.2\tSM:HG00096\tCN:BI'
        RG = r'\t'.join(["'@RG",'ID:'+sample,'LB:'+'Solexa'+sample,'PL:'+inputs['platform_id'][0],
                            'PU:'+sample,'SM:'+SM+"'"])
        command = [bwa,'mem','-M','-t',threads,'-R',RG,in_names['.fa']]+in_names['.fq']+['>',out_names['.sam']]
        #[1a]make start entry which is a new staged_run row  
        self.command = command
        print(self.get_command_str())
        self.db_start(run_id,in_names['.fq'][0])
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            output = subprocess.check_output(' '.join(command),stderr=subprocess.STDOUT,shell=True)
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
            self.db_stop(run_id,{'output':output},'',True)
            results = [out_names['.sam']]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("bwa mem sucessfull........")
                return [out_names['.sam']]
            else:
                print("bwa mem failure...........")
                return False
        else:
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None
