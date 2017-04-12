import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper
import stage_utils as su
from stages.utils.CheckGenerateRG import CheckRG
from stages.utils.CheckGenerateRG import GenerateRG

#function for auto-making svedb stage entries and returning the stage_id
class speedseq_realign(stage_wrapper.Stage_Wrapper):
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
    
    def run(self,run_id,inputs):
        #workflow is to run through the stage correctly and then check for error handles
        #[1b]
        stripped_name = self.strip_path(inputs['.bam'])
        stripped_name = self.strip_in_ext(stripped_name,'.bam')
        out_name = inputs['out_dir']+'/'+stripped_name
        #[2]build command args
        #:::TO DO::: ALLOW THE USER TO GENERATE AND ENTER RG
        threads = str(self.get_params()['-t']['value'])
        mem = str(self.get_params()['-m']['value'])
        speedseq = self.software_path+'/speedseq/bin/speedseq'

        #'@RG\tID:H7AGF.2\tLB:Solexa-206008\tPL:illumina\tPU:H7AGFADXX131213.2\tSM:HG00096\tCN:BI'
        
        realign = [speedseq,'realign','-t',str(inputs['threads']),'-M',str(inputs['mem']),'-o',out_name]
        if inputs['RG'] != '':
            realign += ['-R "'+inputs['RG']+'"']
        else:
           result = []
           result = CheckRG(self.software_path+'/samtools-1.3/samtools',inputs['.bam'], out_name, result)
           if len(result) == 0:
               rg = GenerateRG(stripped_name)
               print "ERROR: " + inputs['.bam'] + "doesn't have RG. " + rg + " is generated."
               realign += ['-R "'+rg+'"']
               
        realign += [inputs['.fa'],inputs['.bam']]

        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print('starting speedseq re-alignment')
            print(' '.join(realign))
            output += subprocess.check_output(' '.join(realign),stderr=subprocess.STDOUT,shell=True) #clean up the inputs now
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
            print(E)
        print('output:\n'+output)
        
        #[3b]check results--------------------------------------------------
        if err == {}:
            #self.db_stop(run_id,{'output':output},'',True)
            results = [out_name+'.bam']
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("speedseq realign sucessfull........")
                return results   #return a list of names
            else:
                print("speedseq realign failure...........")
                return False
        else:
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
