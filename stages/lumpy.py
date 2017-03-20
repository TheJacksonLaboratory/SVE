import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class lumpy(stage_wrapper.Stage_Wrapper):
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
        in_names = {'.bam':inputs['.bam']}
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        out_dir = self.strip_name(in_names['.bam'][0]) #default output directory
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            if inputs.has_key('multisample'):
                stripped_name = 'multisample'
            else:
                stripped_name = self.strip_path(self.strip_in_ext(in_names['.bam'][0],'.bam'))
            out_names = {'.calls' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0],
                         '.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[1]}
        else:
            if inputs.has_key('multisample'):
                cascade = 'multisample'
            else:
                cascade = self.strip_in_ext(in_names['.bam'][0],'.bam')
            out_names = {'.calls' :out_dir+cascade+'_S'+str(self.stage_id)+out_exts[0],
                         '.vcf' :out_dir+cascade+'_S'+str(self.stage_id)+out_exts[1]} 
            
        defaults,params = self.params,[]
        params = [k+'='+str(defaults[k]['value']) for k in defaults]
        print(params)
        #[2a]build command args       
        lumpy   = self.software_path+'/lumpy-sv/bin/lumpyexpress'
        vcfsort = self.software_path+'/vcftools_0.1.12b/bin/vcf-sort'
        PERL = self.software_path+'/vcftools_0.1.12b/perl'
        if os.environ.has_key('PERL5LIB'):
            PERL += ':'+os.environ['PERL5LIB']
        sv_call   = [lumpy,'-B']+ [','.join(in_names['.bam'])]+\
                     ['-m 4','-T',out_dir+'temp','-P','-o',out_names['.calls']] #more work on params
        #sv_fast can do a version that checks for matching .bam, .split.bam and .disc.bam triples (prior samblasted)
        sort_vcf  = [vcfsort,out_names['.calls'],'>',out_names['.vcf']]        
        self.db_start(run_id,in_names['.bam'][0])        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            output += subprocess.check_output(' '.join(sv_call),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
            output += subprocess.check_output(' '.join(sort_vcf),
                                              stderr=subprocess.STDOUT,shell=True,
                                              env={'PERL5LIB':PERL})+'\n'
            #output += subprocess.check_output(' '.join(['rm',out_names['.calls']]),
            #                                  stderr=subprocess.STDOUT,shell=True)+'\n'
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
                print("lumpy sucessfull........")
                self.db_stop(run_id,self.vcf_to_vca(out_names['.vcf']),'',True)
                return results   #return a list of names
            else:
                print("lumpy failure...........")
                self.db_stop(run_id,{'output':output},'',False)
                return False
        else:
            print("lumpy failure...........")
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None
