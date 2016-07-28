import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class bam_stats(stage_wrapper.Stage_Wrapper):
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
    
    def summary_as_list(self,summary):
        out = []
        lines = summary.split('\n')
        for line in lines:
            line = line.replace(':','')            #get rid of the:
            expand = line.split('\t')              #samtools v1.0 uses \t
            if len(expand)>1: out += [expand[0:2]] #chop the #...
        return out
        
        
    #override this function in each wrapper...
    #bwa sampe ref.fa L.sai R.sai L.fq R.fq -f out.sam
    def run(self,run_id,inputs):
        #workflow is to run through the stage correctly and then check for error handles
    
        #[1a]get input names and output names setup
        in_names  = {'.bam':inputs['.bam'][0]}
        
        out_ext = self.split_out_exts()[0]
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['.bam'],'.bam'))
            out_name = out_dir+stripped_name+'_S'+str(self.stage_id)+out_ext
        else:
            out_name = self.strip_in_ext(in_names['.bam'],'.bam')+out_ext
        
        #[2a]build command args
        samtools = self.software_path+'/samtools-1.0/bin/samtools'
        summary =   [samtools,'stats',in_names['.bam'],'| grep ^SN | cut -f 2-']
        size    =   [samtools,'view','-H', in_names['.bam'],
                     "| grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'"]
        average =   [samtools,'depth',in_names['.bam'], "| awk '{sum+=$3} END {print sum/%s}'"]
        #write   =   ['echo',' > ',out_name]             
        #[2b]make start entry which is a new staged_run row
        self.command = summary+size+average
        print(self.get_command_str())
        self.db_start(run_id,in_names['.bam'])
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print('calculating total ref size by sum of all contigs')
            L = []
            output = subprocess.check_output(' '.join(size),stderr=subprocess.STDOUT,shell=True)
            output = output[0:-1]
            if type(output) is str and len(output)>0: L += [['ref size',output]]
            else: L += [['ref size',0]]
            average[-1] = average[-1]%output
            print('calculating the average by reading bam depth')
            output = subprocess.check_output(' '.join(average),stderr=subprocess.STDOUT,shell=True)
            output = output[0:-1]
            if type(output) is str and len(output)>0: L += [['average depth',output]]
            else: L += [['average depth',0.0]]
            print('calculating summaries of read statistics')
            output = subprocess.check_output(' '.join(summary),stderr=subprocess.STDOUT,shell=True)
            if type(output) is str and len(output)>0: L += self.summary_as_list(output)
            output = 'calculating summaries of read statistics\n'            
            for line in L:
                if len(line)>1:
                    output += '%s = %s\n'%(line[0],line[1])
            with open(out_name, 'w') as f:
                f.write(output)
            #output += subprocess.check_output(' '.join(['echo',output,'>',out_name]),stderr=subprocess.STDOUT,shell=True)
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
        except IOError as E:
            print('IO error: '+E.strerror)        #what you would see in the term
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
            self.db_stop(run_id,{'output':output},'',True)
            results = [out_name]
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