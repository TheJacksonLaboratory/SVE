import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class picard_sam_convert(stage_wrapper.Stage_Wrapper):
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
        in_name  = {'.sam':inputs['.sam'][0]}
        #out_ext = self.split_out_exts()[0] #change this to '.bam' to be destructive??
        out_name = self.strip_in_ext(in_name['.sam'],'.sam').rstrip('_')

        #[2a]build command args
        rg = '.'.join(in_name['.sam'].rsplit('/')[-1].rsplit('.')[0:-1])+'_RG'#may need to do addreplacerg
        java   = self.tools['JAVA-1.8']
        mem    = '-Xmx32g'
        picard = self.tools['PICARD']
        sort   =  [java,mem,'-jar',picard,'SortSam','I=',in_name['.sam'],
                   'O=',out_name+'.sorted.bam','SORT_ORDER=coordinate']                  #can delete .sam after this step
        mark   =  [java,mem,'-jar',picard,'MarkDuplicates','I=',out_name+'.sorted.bam',
                   'O=',out_name+'.bam','METRICS_FILE=',out_name+'.picard.metrics.txt']  #delete .sorted.bam after this steps
        index  = [java,mem,'-jar',picard,'BuildBamIndex','I=',out_name+'.bam']  #no .bam.bai here ?
        rename = ['mv',out_name+'.bai',out_name+'.bam.bai']        
        clean = ['rm',in_name['.sam'],out_name+'.sorted.sam',out_name+'.sorted.bam']
        #[2b]make start entry which is a new staged_run row
        self.command = sort
        print(self.get_command_str())
        #self.db_start(run_id,in_name['.sam'])
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            output += subprocess.check_output(sort,stderr=subprocess.STDOUT)
            output += subprocess.check_output(mark,stderr=subprocess.STDOUT)
            output += subprocess.check_output(index,stderr=subprocess.STDOUT)
            output += subprocess.check_output(rename,stderr=subprocess.STDOUT)
            output += subprocess.check_output(clean,stderr=subprocess.STDOUT) #clean up the inputs now
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
        
        #[3b]check results--------------------------------------------------
        if err == {}:
            #self.db_stop(run_id,{'output':output},'',True)
            results = [out_name+'.bam']
            print results
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("<<<<<<<<<<<<<picard_sam_convert sucessfull>>>>>>>>>>>>>>>\n")
                return results   #return a list of names
            else:
                print("<<<<<<<<<<<<<picard_sam_convert failure>>>>>>>>>>>>>>>\n")
                return False
        else:
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
