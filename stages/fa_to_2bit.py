import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class fa_to_2bit(stage_wrapper.Stage_Wrapper):
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
        in_names = {'.fa':inputs['.fa']}
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        out_dir = self.strip_name(in_names['.fa'][0]) #default output directory
        out_names = {'.2bit':[]}        
        for i in in_names['.fa']:
            out_names['.2bit'] += [self.strip_in_ext(i,'.fa')+out_exts[0]]
        #[2a]build command args

        faToTwoBit = self.software_path+'/faToTwoBit'
        #self.db_start(run_id,','.join(in_names['.fa']))        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            for i in range(len(in_names['.fa'])):
                call = [faToTwoBit, in_names['.fa'][i], out_names['.2bit'][i]]
                output += subprocess.check_output(' '.join(call),
                                                  stderr=subprocess.STDOUT,shell=True)
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
            results = out_names['.2bit']
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("<<<<<<<<<<<<<fa_to_2bit sucessfull>>>>>>>>>>>>>>>\n")
                #self.db_stop(run_id,{'output':output},'',True)
                return results   #return a list of names
            else:
                print("<<<<<<<<<<<<<fa_to_2bit failure>>>>>>>>>>>>>>>\n")
                #self.db_stop(run_id,{'output':output},'',False)
                return False
        else:
            print("<<<<<<<<<<<<<fa_to_2bit failure>>>>>>>>>>>>>>>\n")
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
