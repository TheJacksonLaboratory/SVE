import os
import sys
import subprocess32 as subprocess
import utils.breakdancer2vcf as bd
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class breakdancer(stage_wrapper.Stage_Wrapper):
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
                    '.bam':inputs['.bam']}
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        out_dir = self.strip_name(in_names['.fa']) #default output directory
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['.bam'][0],'.bam'))
            out_names = {'.calls' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0],
                         '.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[1]}
        else:
            cascade = self.strip_in_ext(in_names['.bam'][0],'.bam')
            out_names = {'.calls' :cascade+'_S'+str(self.stage_id)+out_exts[0],
                         '.vcf' :cascade+'_S'+str(self.stage_id)+out_exts[1]}  
        #[2a]build command args
        #PATH = os.environ['PATH'] 
#        PATH = os.path.dirname(os.path.abspath('~'))+'/software/perl/bin:'+os.environ['PATH']
#        PERL5LIB = os.path.dirname(os.path.abspath('~'))+'/software/perl/lib/site_perl/5.20.1:'+\
#                   os.path.dirname(os.path.abspath('~'))+'/software/perl/lib/5.20.1'#+os.environ['PERL5LIB']
        cfg    = out_dir+"bd_confg.txt" #new version 1.1.2 working!
#        perl   = os.path.dirname(os.path.abspath('~'))+'/software/perl/bin/perl'
        config = self.software_path+'/breakdancer-1.4.5/perl/bam2cfg.pl'        
        breakd = self.software_path+'/breakdancer-1.4.5/bin/breakdancer-max'
        configure = ['perl',config,'-q','30','-n','10000'] + in_names['.bam']+['>',cfg]
        sv_call   = [breakd, cfg,'>', out_names['.calls']]
        
#        breakdancer_max $OUT$BAM2CFG > $OUT$BD
#        python ./breakdancer2vcf.py $REF $OUT         
#                
        self.db_start(run_id,in_names['.bam'][0])        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            output += subprocess.check_output(' '.join(configure),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
            output += subprocess.check_output(' '.join(sv_call),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
            table = bd.read_breakdancer(out_names['.calls'])
            bd.write_vcf(out_names['.vcf'],bd.vcf_header(in_names['.fa']),bd.build_vcf(table))
#            output += subprocess.check_output(' '.join(['rm',out_names['.calls']]),
#                                              stderr=subprocess.STDOUT,shell=True)+'\n'
#            output += subprocess.check_output(' '.join(['rm',cfg]),
#                                              stderr=subprocess.STDOUT,shell=True)+'\n'
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
            results = [out_names['.vcf']]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("sucessfull........")
                self.db_stop(run_id,self.vcf_to_vca(out_names['.vcf']),'',True)
                return results   #return a list of names
            else:
                print("failure...........")
                self.db_stop(run_id,{'output':output},'',False)
                return False
        else:
            print("failure...........")
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None