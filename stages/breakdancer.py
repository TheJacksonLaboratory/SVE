import os
import sys
import subprocess32 as subprocess
import utils.breakdancer2vcf as bd
sys.path.append('../') #go up one in the modules
import stage_wrapper
from stages.utils.CheckVcf import GetCallCount

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
        if ('.fa' not in inputs) or ('.bam' not in inputs) or ('out_dir' not in inputs):
            print "ERROR: .fa, .bam, and out_dir are required for genome_strip.py"
            return None
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        out_dir = inputs['out_dir'] + '/'
        stripped_name = ''
        if len(inputs['.bam']) == 1: stripped_name = self.strip_path(self.strip_in_ext(inputs['.bam'][0],'.bam'))
        else: stripped_name = 'joint'
        #[2a]build command args
        
        #build temp directory to work in
        sub_dir = out_dir+stripped_name+'_S'+str(self.stage_id)+'/'
        if not os.path.exists(sub_dir): os.makedirs(sub_dir)

        out_names = {'.calls' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0],
                     '.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[1]}
        #[2a]build command args
        PERL5LIB = self.tools['PERL_LIB_PATH'] + '/lib/perl5'
        if os.environ.has_key('PERL5LIB'): PERL5LIB += ':' + os.environ['PERL5LIB']
        PATH = self.tools['SAMTOOLS_PATH'] + ':' + os.environ['PATH']
        cfg    = sub_dir+"bd_confg.txt" #new version 1.1.2 working!
        config = self.tools['BREAKDANCER_PATH'] + '/perl/bam2cfg.pl'        
        breakd = self.tools['BREAKDANCER_PATH'] + '/build/bin/breakdancer-max'
        configure = ['perl',config,'-q','30','-n','10000'] + inputs['.bam']+['>',cfg]
        sv_call   = [breakd, cfg,'>', out_names['.calls']]
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(configure))
            output += subprocess.check_output(' '.join(configure),
                                              stderr=subprocess.STDOUT,shell=True, env={'PERL5LIB':PERL5LIB, 'PATH':PATH})+'\n'
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(sv_call))
            output += subprocess.check_output(' '.join(sv_call),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
            table = bd.read_breakdancer(out_names['.calls'])
            bd.write_vcf(out_names['.vcf'],bd.vcf_header(inputs['.fa']),bd.build_vcf(table))
            os.remove(out_names['.calls'])
            os.remove(cfg)
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
        if err != {}:
            print err
        if GetCallCount(out_names['.vcf']) > 0:
            print("<<<<<<<<<<<<<breakdancer sucessfull>>>>>>>>>>>>>>>\n")
            return out_names['.vcf']   #return a list of names
        else:
            print("<<<<<<<<<<<<<breakdancer failure>>>>>>>>>>>>>>>\n")
            return None
