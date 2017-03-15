import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class bwa_sampe(stage_wrapper.Stage_Wrapper):
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
        in_names  = {'.fa':inputs['.fa'][0],'.fq':inputs['.fq'],'.sai':inputs['.sai']}
        out_ext = self.split_out_exts()[0]
        out_name = self.strip_in_ext(in_names['.fq'][0],'.fq')[:-1]+'_S'+str(self.stage_id)

        threads = str(self.get_params()['-t']['value'])
        bwa = self.software_path+'/bwa-master/bwa' #latest release
        samtools = self.software_path+'/samtools-1.3/samtools'
        java   = self.software_path+'/jre1.8.0_51/bin/java'
        sambamba = self.software_path+'/sambamba_v0.6.6'
        mem    = '-Xmx%sg'%str(self.get_params()['-m']['value'])
        picard = self.software_path+'/picard-tools-2.5.0/picard.jar' #latest release here
     

        #[2]build command args
        sampe = [bwa,'sampe',in_names['.fa']]
        sampe += in_names['.sai']+in_names['.fq']
        for k in self.params:
            param = self.params[k]
            if param['type']=='bool': command += [k]
            else:                     command += [k, str(param['value'])]  
        
        view = ['|',samtools,'view','-Sb','-','-o',out_name+'.bam']
   
        sort   =  [sambamba,'sort','-o',out_name+'.sorted.bam','-l','5','-t',threads,out_name+'.bam']

        mark   =  [java,mem,'-jar',picard,'MarkDuplicates','I='+out_name+'.sorted.bam',
                  'O='+out_name+'.bam','METRICS_FILE='+out_name+'.picard.metrics.txt',
                  'MAX_RECORDS_IN_RAM='+str(250000*16)] #delete .sorted.bam after this steps
        
        index  =  [sambamba,'index','-t',threads,out_name+'.bam'] 

        clean = ['rm','-rf',out_name+'.sorted.bam',
                out_name+'.sorted.bam.bai',
                out_name+'.picard.metrics.txt'] #clean just the sorted bam file when done
        #[1a]make start entry which is a new staged_run row
        self.command = sampe + view
        print(self.get_command_str())
        self.db_start(run_id,in_names['.fq'])
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            output = subprocess.check_output(sampe+view,stderr=subprocess.STDOUT)
            output += subprocess.check_output(' '.join(sort),stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(mark),stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(clean),stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(index),stderr=subprocess.STDOUT,shell=True)
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
            results = [out_name]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("bwa sampe sucessfull........")
                return [out_name+out_ext]
            else:
                print("bwa sampe failure...........")
                return False
        else:
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None
