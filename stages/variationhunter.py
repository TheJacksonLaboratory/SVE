import os
import sys
import time
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class variationhunter(stage_wrapper.Stage_Wrapper):
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
             
        #[1a]get input names and output names setup
        in_names = {'.fa':inputs['.fa'][0],
                    'L.fq':inputs['L.fq'][0],
                    'R.fq':inputs['R.fq'][0]}
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        out_dir = self.strip_name(in_names['L.fq']) #default output directory
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['L.fq'],'L.fq'))
            out_names = {'.calls' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0]}
        else:
            cascade = self.strip_in_ext(in_names['.bam']['case'],'.bam')
            out_names = {'.calls' :cascade+'_S'+str(self.stage_id)+out_exts[0]} 
            
        #load params    
        defaults,params = self.params,[]
        params = [k+' '+str(defaults[k]['value']) for k in defaults]
        
        
        #[a]use to run several sub scripts via command line/seperate process
        vh     = self.software_path+'/variationhunter_0.4/'
        cl  = vh+'clustering/VH'
        sc  = vh+'selection/multiInd_SetCover'

        #clustering step
        #genome.satellite.bed: Create this file from the reference genome you have. It's basically RepeatMasker output, with "satellite" only, and converted to a BED file.
        cluster = []
        
        #selection step
        #This is a BED file that shows locations of the assembly gaps
        select  = []
        
        vcf     = [] #convert from vh format to VCF
        
        self.db_start(run_id,in_names['L.fq']) 
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print(' '.join(cluster))
            output += subprocess.check_output(' '.join(cluster),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
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
            print('vcf write os/file IO error')
            err['output'] = 'vcf write os/file IO error'
            err['message'] = 'vcf write os/file IO error'
            err['code'] = 1
        print('output:\n'+output)
                                                
        #[3b]check results--------------------------------------------------
        if err == {}:
            if all([os.path.exists(r) for r in results]):
                print("sucessfull........")
                self.db_stop(run_id,self.vcf_to_vca(out_names['.vcf']),'',True)
                return results   #return a list of names
            else:
                print("failure...........")
                self.db_stop(run_id,{'output':output},err['message'],False)
                return False
        else:
            print("failure...........")
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None