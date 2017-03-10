import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper
import stage_utils as su

#function for auto-making svedb stage entries and returning the stage_id
class bam_clean(stage_wrapper.Stage_Wrapper):
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
        in_names  = {'.header':inputs['.header'][0],'.valid':inputs['.valid'][0],'.bam':inputs['.bam'][0]}
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['.header'],'.header'))
            out_name = {'.clean.bam' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0]}
        else:
            cascade = self.strip_in_ext(in_names['.header'],'.header') #default path
            out_name = {'.clean.bam' :cascade+'_S'+str(self.stage_id)+out_exts[0]}
            
        #[2]build command args
        if not os.path.exists(out_dir): os.makedirs(out_dir)
        
        #valid check: No errors found
        valid_string = ''
        with open(in_names['.valid'],'r') as f: valid_string = f.readlines()
        valid = (len(valid_string) == 1 and valid_string[0].startswith('No errors found'))
        #conditional execution starts with this parsed string
        
        mem      = '-Xmx4g'
        samtools = self.software_path+'/samtools-1.3/samtools'
        java   = self.software_path+'/jre1.8.0_51/bin/java'
        picard = self.software_path+'/picard-tools-2.5.0/picard.jar' #latest release here
        phred64to33 = self.software_path+'/SVE/stages/utils/phred64to33.py'
        
        reheader  = [samtools,'reheader','-i',in_names['.header']+'.rg',in_names['.bam']]
        cleansam  = [java,mem,'-jar',picard,'CleanSam','I=%s'%in_names['.bam'],'O=%s'%in_names['.bam']]
        fixmate   = [java,mem,'-jar',picard,'FixMateInformation','I=%s'%in_names['.bam'],'O=%s'%in_names['.bam']]
        #samtools view -Sh old.bam | SVE/stages/utils/phred64to33.py | samtools view -Sb - > ./phred33.bam
        fixphred  = [samtools,'view','-Sh',in_names['.bam'],'|',phred64to33,'|',samtools,'view','-Sb','-','>',out_name['.clean.bam']]

        #[2b]make start entry which is a new staged_run row
        #[1a]make start entry which is a new staged_run row  
        self.command = ''
        print(self.get_command_str())
        self.db_start(run_id,in_names['.bam'][0])
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print(valid_string)
            print('validated bam file = %s'%valid)
#            output += subprocess.check_output(' '.join(reheader),stderr=subprocess.STDOUT,shell=True)
#            output += subprocess.check_output(' '.join(cleansam),stderr=subprocess.STDOUT,shell=True)
#            output += subprocess.check_output(' '.join(fixmate),stderr=subprocess.STDOUT,shell=True)
#            output += subprocess.check_output(' '.join(fixphred),stderr=subprocess.STDOUT,shell=True)
        #catch all errors that arise under normal cleaning behavior
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
            self.db_stop(run_id,{'output':output},'',True)
            results = [out_name+'.bam']
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("sucessfull........")
                return results   #return a list of names
            else:
                print("failure...........")
                return False
        else:
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None
