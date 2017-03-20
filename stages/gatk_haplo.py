import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class gatk_haplo(stage_wrapper.Stage_Wrapper):
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
    #./software/jre1.8.0_25/bin/java -jar ./software/GATK_3.3/GenomeAnalysisTK.jar 
    #-T HaplotypeCaller -R ./node_data/rg1.fa -I ./node_data/node02/rg1_R1_S1_S8.bam 
    #--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 
    #-o ./node_data/node02/rg1_S1_S8_S13.vcf"""
    def run(self,run_id,inputs):
        #workflow is to run through the stage correctly and then check for error handles
    
        #[1a]get input names and output names setup
        in_names  = {'.fa':inputs['.fa'][0],'.bam':inputs['.bam']}
        
        out_exts = self.split_out_exts()
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['.bam'][0],'.bam'))
            out_names = {'.vcf'  : out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0]}
        else:
            cascade = self.strip_in_ext(in_names['.bam'][0],'.bam')
            out_names = {'.vcf'  :cascade+'_S'+str(self.stage_id)+out_exts[0]} 
        
        #[2a]build command args
        java = self.software_path+'/jre1.8.0_51/bin/java'
        gatk = self.software_path+'/GATK_3.7/GenomeAnalysisTK.jar'
        command = [java,'-Xmx12g','-jar',gatk,'-T','HaplotypeCaller',
                   '-R',in_names['.fa'],'-I'] + in_names['.bam']
        #add this param function
        for k in self.params:
            param = self.params[k]
            if param['type']=='bool':
                if param['value']: command += [k]
            else: command += [k, str(param['value'])]  
        command += ['-o',out_names['.vcf']]
        #[2b]make start entry which is a new staged_run row
        self.command = command
        print(self.get_command_str())
        self.db_start(run_id,in_names['.bam'][0])
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            output = subprocess.check_output(command,stderr=subprocess.STDOUT)
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
            results = [out_names['.vcf']]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("GATK sucessfull........")
                self.db_stop(run_id,self.vcf_to_vca(out_names['.vcf']),'',True)
                return results   #return a list of names
            else:
                print("GATK failure...........")
                self.db_stop(run_id,{'output':output},'',False)
                return False
        else:
            print("GATK failure...........")
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None
