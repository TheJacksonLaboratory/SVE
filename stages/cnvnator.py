import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper
import read_utils as sr
from stages.utils.CheckVcf import GetCallCount

#function for auto-making svedb stage entries and returning the stage_id
class cnvnator(stage_wrapper.Stage_Wrapper):
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
    def run(self,run_id,inputs):
        #workflow is to run through the stage correctly and then check for error handles
    
        #[1a]get input names and output names setup
        if ('.fa' not in inputs) or ('.bam' not in inputs) or ('out_dir' not in inputs):
            print "ERROR: .fa, .bam, and out_dir are required for genome_strip.py"
            return None
        
        out_exts = self.split_out_exts()
        out_dir = inputs['out_dir'] + '/'
        stripped_name = ''
        if len(inputs['.bam']) == 1: stripped_name = self.strip_path(self.strip_in_ext(inputs['.bam'][0],'.bam'))
        else: stripped_name = 'joint'
        sub_dir = out_dir+stripped_name+'_S'+str(self.stage_id)+'/'
        if not os.path.exists(sub_dir): os.makedirs(sub_dir)
        out_names = {'.root' : sub_dir+'temp',
                     '.calls': sub_dir+'temp'+out_exts[1],
                     '.vcf'  : out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[2]}
        #[2a]build command args
        
        #split the ref seq into seperate chroms...
        cnvnator = self.tools['CNVNATOR']
        cnv2vcf  = self.tools['CNVNATOR2VCF']

        #[self.strip_in_ext(self.strip_path(i),'.fa') for i in in_names['.fa']] #by using the input list
        bin_size = str(150)
        #refd = self.strip_name(inputs['.fa']) #this is a bit hackish
        
        extr     = [cnvnator, '-unique','-root', out_names['.root']+'.tree.root','-tree']+ inputs['.bam']
        hist     = [cnvnator, '-root', out_names['.root']+'.tree.root',
                    '-outroot',out_names['.root']+'.his.root','-his', bin_size, '-d', sub_dir]
        stats    = [cnvnator, '-root', out_names['.root']+'.his.root', '-stat', bin_size]
        sig      = [cnvnator, '-root', out_names['.root']+'.his.root', '-partition', bin_size]
        call     = [cnvnator, '-root', out_names['.root']+'.his.root','-call', bin_size, '>', out_names['.calls']]
        conv     = ['perl', cnv2vcf, out_names['.calls'], '>', out_names['.vcf']]
        
        #[2b]make start entry which is a new staged_run row

        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (" ".join(extr))
            output += subprocess.check_output(' '.join(extr),stderr=subprocess.STDOUT, shell=True)+'\n'
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (" ".join(hist))
            output += subprocess.check_output(' '.join(hist),stderr=subprocess.STDOUT, shell=True)+'\n'
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (" ".join(stats))
            output += subprocess.check_output(' '.join(stats),stderr=subprocess.STDOUT, shell=True)+'\n'
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (" ".join(sig))
            output += subprocess.check_output(' '.join(sig),stderr=subprocess.STDOUT, shell=True)+'\n'
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (" ".join(call))
            output += subprocess.check_output(' '.join(call),stderr=subprocess.STDOUT, shell=True)+'\n'
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (" ".join(conv))
            output += subprocess.check_output(' '.join(conv),stderr=subprocess.STDOUT, shell=True)+'\n'
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print ('rm -rf %s'%sub_dir)
            output += subprocess.check_output('rm -rf %s'%sub_dir,
                                              stderr=subprocess.STDOUT, shell=True)
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
        if err != {}:
            print err
        if GetCallCount(out_names['.vcf']) > 0:
            print("<<<<<<<<<<<<<cnvnator sucessfull>>>>>>>>>>>>>>>\n")
            return out_names['.vcf']   #return a list of names
        else:
            print("<<<<<<<<<<<<<cnvnator failure>>>>>>>>>>>>>>>\n")
            return None
