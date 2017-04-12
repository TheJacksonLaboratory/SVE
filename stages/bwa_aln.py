import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper
import stage_utils as su
from stages.utils.CheckGenerateRG import GenerateRG

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
        #add tigra.ctg.fa bam by allowing a single '.fa' input key
        if len(inputs['.fq']) == 2: csl = su.get_common_string_left(inputs['.fq'])
        else:                       csl = inputs['.fq'][0]
        stripped_name = self.strip_path(csl)
        self.strip_in_ext(stripped_name,'.fq')
        if (stripped_name[-1:] == '_' or stripped_name[-1:] == '.'): stripped_name = stripped_name[:-1]

        out_dir = inputs['out_dir']
        out_name = out_dir + '/' + stripped_name

        threads = str(inputs['threads'])
        bwa = self.software_path+'/bwa-master/bwa' #latest release
        samtools = self.software_path+'/samtools-1.3/samtools'
        sambamba = self.software_path+'/sambamba_v0.6.6'

        #[2]build command args
        aln1 = [bwa,'aln','-t',threads,inputs['.fa'],inputs['.fq'][0],'-f',out_name+'_1.sai']
        aln2 = [bwa,'aln','-t',threads,inputs['.fa'],inputs['.fq'][1],'-f',out_name+'_2.sai']
        #'@RG\tID:H7AGF.2\tLB:Solexa-206008\tPL:illumina\tPU:H7AGFADXX131213.2\tSM:HG00096\tCN:BI'
        RG = inputs['RG']
        if RG == '': # RG is not defined
            RG = GenerateRG(stripped_name)

        sampe = [bwa,'sampe','-r',"'"+RG+"'",inputs['.fa'],out_name+'_1.sai',out_name+'_2.sai',inputs['.fq'][0],inputs['.fq'][1]]+['|']
        view  = [samtools,'view','-Sb','-','-o',out_name+'.bam']
        sort  = [sambamba,'sort','-o',out_name+'.sorted.bam','-l','5','-t',threads,out_name+'.bam']
        
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(aln1))
            output += subprocess.check_output(' '.join(aln1),stderr=subprocess.STDOUT,shell=True)
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(aln2))
            output += subprocess.check_output(' '.join(aln2),stderr=subprocess.STDOUT,shell=True)
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(sampe+view))
            output += subprocess.check_output(' '.join(sampe+view),stderr=subprocess.STDOUT,shell=True)
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(sort))
            output += subprocess.check_output(' '.join(sort),stderr=subprocess.STDOUT,shell=True)
            move = ['mv',out_name+'.sorted.bam',out_name+'.bam']
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(move))
            output += subprocess.check_output(' '.join(move),stderr=subprocess.STDOUT,shell=True)
            move = ['mv',out_name+'.sorted.bam.bai',out_name+'.bam.bai']
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(move))
            #output += subprocess.check_output(' '.join(move),stderr=subprocess.STDOUT,shell=True)
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
            results = [out_name+'.bam']
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("<<<<<<<<<<<<<bwa sampe sucessfull>>>>>>>>>>>>>>>\n")
                return out_name+'.bam'
            else:
                print("<<<<<<<<<<<<<bwa sampe failure>>>>>>>>>>>>>>>\n")
                return False
        else:
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
