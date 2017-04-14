import os
import sys
import subprocess32 as subprocess
import stage_wrapper
import stage_utils as su
from stages.utils.CheckGenerateRG import GenerateRG

#function for auto-making svedb stage entries and returning the stage_id
class speedseq_align(stage_wrapper.Stage_Wrapper):
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
        #[1b]
        out_ext = self.split_out_exts()[0]
        #add tigra.ctg.fa bam by allowing a single '.fa' input key
        if len(inputs['.fq']) == 2:   csl = su.get_common_string_left(inputs['.fq'])
        else:                         csl = inputs['.fq'][0]
        stripped_name = self.strip_path(csl)
        if (stripped_name[-1:] == '_' or stripped_name[-1:] == '.'): stripped_name = stripped_name[:-1]
        out_dir = inputs['out_dir']
        out_name = out_dir + '/' + stripped_name
        RG = ''
        if (not 'RG' in inputs) or (inputs['RG'] == ''): # RG is not defined
            RG = GenerateRG(stripped_name)
	else:
            RG = inputs['RG']

        #[2]build command args
        threads = str(inputs['threads'])
        mem = str(inputs['mem'])
	speedseq = self.tools['SPEEDSEQ']
        #'@RG\tID:H7AGF.2\tLB:Solexa-206008\tPL:illumina\tPU:H7AGFADXX131213.2\tSM:HG00096\tCN:BI'
        align = [speedseq,'align','-t',threads,'-R','"'+RG+'"','-M',mem,'-T',out_dir,
                 '-o',out_name,inputs['.fa']]+inputs['.fq'] #out_name: speedseq -o is prefix (without.bam)
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print(' '.join(align))
            output += subprocess.check_output(' '.join(align),
                                              stderr=subprocess.STDOUT,
                                              shell=True) #clean up the inputs now
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
            err['output']= str(E)
        print('output:\n'+output)
        
        #[3b]check results--------------------------------------------------
        if err == {}:
            #self.db_stop(run_id,{'output':output},'',True)
            results = [out_name+'.bam']
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("<<<<<<<<<<<<<speedseq align sucessfull>>>>>>>>>>>>>>>\n")
                return results   #return a list of names
            else:
                print("<<<<<<<<<<<<<speedseq align failure>>>>>>>>>>>>>>>\n")
                return False
        else:
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
