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
            out_name = out_dir+stripped_name+'_S'+str(self.stage_id)
        else:
            out_name = self.strip_in_ext(in_names['.bam'],'.bam')

        if inputs.has_key('chroms'): #subset of the chroms in the bam header
            chroms = inputs['chroms'].split(',')
        else:
            chroms = [str(i) for i in range(1,23)] #none selected will do a full stats run

        #[2a]build command args
        samtools = self.software_path+'/samtools-1.3/samtools'
        summary =   [samtools,'stats',in_names['.bam'],'| grep ^SN | cut -f 2-']
        header  =   [samtools, 'view', '-SH', in_names['.bam']]
        #some routines here for X:Y analysis for gender estimation

        #write   =   ['echo',' > ',out_name]             
        #[2b]make start entry which is a new staged_run row
        self.command = summary
        print(self.get_command_str())
        self.db_start(run_id,in_names['.bam'])
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            h = subprocess.check_output(' '.join(header),stderr=subprocess.STDOUT,shell=True)
            with open(out_name+'.header','w') as f:
                f.write(h)
            #get sequence names and lengths
            seqs = {}
            for line in h.split('\n'):
                l = line.split('\t')
                if l[0].startswith('@SQ'):
                    seqs[l[1].split(':')[-1]] = [int(l[2].split(':')[-1])]
            #get converage over each sequence
            for k in sorted(seqs):
                seq_cov = [samtools, 'depth', '-r %s'%k, in_names['.bam'], "| awk '{sum+=$3} END {print sum}'"]
                c = subprocess.check_output(' '.join(seq_cov),stderr=subprocess.STDOUT,shell=True)
                x = 0
                try: x = int(c)
                except ValueError: pass
                seqs[k] += [x]
            c,x,y = '',0,0
            for k in sorted(seqs):
                if k in chroms or k in ['chr'+i for i in chroms]:
                    c += k+'='+str(seqs[k][0])+':'+str(seqs[k][1])+'\n'
                    x += seqs[k][0]
                    y += seqs[k][1]
            c += 'average coverage = %s\n'%(int(round(1.0*x/y,0)))
            c += 'over total length of %s\n'%x
            with open(out_name+'.cov','w') as f:
                f.write(c)
            #get the summary
            s = subprocess.check_output(' '.join(summary), stderr=subprocess.STDOUT, shell=True)
            with open(out_name+'.summary','w') as f:
                f.write(s)
            output = '' #put together the pass back for variant_processor.py here
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