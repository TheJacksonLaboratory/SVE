import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper
import stage_utils as su

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
        in_names  = {'.fa':inputs['.fa'][0],'.fq':inputs['.fq'],}
        out_ext = self.split_out_exts()[0]

        #add tigra.ctg.fa bam by allowing a single '.fa' input key
        if len(in_names['.fq']) == 2: csl = su.get_common_string_left(in_names['.fq'])
        else:                         csl = in_names['.fq'][0]
        stripped_name = self.strip_path(csl)
        self.strip_in_ext(stripped_name,'.fq')

        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            out_name = out_dir+stripped_name
            out_name = out_name.rstrip('_')
        else: #untested...
            right = in_names['.fa'][0].rsplit('/')[-1]
            out_dir = in_names['.fq'][0].replace(right,'')
            cascade = self.strip_in_ext(in_names['.fq'][0],'.fq')
            out_name = cascade
            out_name = out_name.rstrip('_')

        if not os.path.exists(out_dir): os.makedirs(out_dir)


        threads = str(self.get_params()['-t']['value'])
        bwa = self.tools['BWA'] #latest release
        samtools = self.tools['SAMTOOLS']
        java   = self.tools['JAVA-1.8']
        sambamba = self.tools['SAMBAMBA']
        mem    = '-Xmx%sg'%str(self.get_params()['-m']['value'])
        picard = self.tools['PICARD'] #latest release here
     

        #[2]build command args
        aln1 = [bwa,'aln','-t',threads,in_names['.fa'],in_names['.fq'][0],'-f',out_name+'_1.sai']
        aln2 = [bwa,'aln','-t',threads,in_names['.fa'],in_names['.fq'][1],'-f',out_name+'_2.sai']
        sampe = [bwa,'sampe',in_names['.fa'],out_name+'_1.sai',out_name+'_2.sai',in_names['.fq'][0],in_names['.fq'][1]]+['|']
        #for k in self.params:
        #    param = self.params[k]
        #    if param['type']=='bool': command += [k]
        #    else:                     command += [k, str(param['value'])]  
        
        view = [samtools,'view','-Sb','-','-o',out_name+'.bam']
   
        sort   =  [sambamba,'sort','-o',out_name+'.sorted.bam','-l','5','-t',threads,out_name+'.bam']

        mark   =  [java,mem,'-jar',picard,'MarkDuplicates','I='+out_name+'.sorted.bam',
                  'O='+out_name+'.bam','METRICS_FILE='+out_name+'.picard.metrics.txt',
                  'MAX_RECORDS_IN_RAM='+str(250000*16)] #delete .sorted.bam after this steps
        
        index  =  [sambamba,'index','-t',threads,out_name+'.bam'] 

        clean = ['rm','-rf',out_name+'.sorted.bam',
                out_name+'.sorted.bam.bai',
                out_name+'.picard.metrics.txt'] #clean just the sorted bam file when done
        #[1a]make start entry which is a new staged_run row
        self.command = aln1
        print(self.get_command_str())
        #self.db_start(run_id,in_names['.fq'][0])
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            output += subprocess.check_output(' '.join(aln1),stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(aln2),stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(sampe+view),stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(sort),stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(mark),stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(clean,stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(index,stderr=subprocess.STDOUT,shell=True)
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
            results = [out_name]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("<<<<<<<<<<<<<bwa sampe sucessfull>>>>>>>>>>>>>>>\n")
                return [out_name+out_ext]
            else:
                print("<<<<<<<<<<<<<bwa sampe failure>>>>>>>>>>>>>>>\n")
                return False
        else:
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
