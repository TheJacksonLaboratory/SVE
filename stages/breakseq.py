import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class breakseq(stage_wrapper.Stage_Wrapper):
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
        in_names = {'.fa':inputs['.fa'][0], 
                    '.gff':inputs['.gff'][0], #base gff and brkpt lib files here
                    '.bam':inputs['.bam']}
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        out_dir = self.strip_name(in_names['.fa']) #default output directory
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['.bam'][0],'.bam'))
            out_names = {'.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0]}
        else:
            cascade = self.strip_in_ext(in_names['.bam'][0],'.bam')
            out_names = {'.vcf' :cascade+'_S'+str(self.stage_id)+out_exts[0]}  
        #[2a]build command args

        #build temp directory to work in
        sub_dir = out_dir+'/'+'S'+str(self.stage_id)+'/'
        if not os.path.exists(sub_dir): os.makedirs(sub_dir)
            
        python    = 'python'
        samtools  = self.software_path+'/samtools-0.1.19/samtools'
        bwa       = self.software_path+'/bwa-master/bwa'
        breakseq  = self.software_path+'/breakseq2-2.2/scripts/run_breakseq2.py'
        #brkptlib = human_g1k_v37_decoy_S35.brkptlib.gff etc files ...
        #run_breakseq2.py --reference b37.fasta --bams bwamem.bam --work work --bwa bwa-0.7.12/bwa --samtools samtools-0.1.19/samtools --bplib_gff bplib.gff --nthreads 4 --sample NA12878
        w = str(self.get_params()['window']['value']) 
        j = str(self.get_params()['junction']['value'])
        call      = [python,breakseq,'--bwa',bwa,'--samtools',samtools,
                     '--reference',in_names['.fa'],'--bplib_gff',in_names['.gff'],
                     '--work',sub_dir,'--bams']+in_names['.bam']+\
                    ['--nthreads',str(4),'--min_span',str(2),'--window',max(100,w),
                     '--min_overlap',str(2),'--junction_length',max(200,j)] #junctio =2x lead length
        #decompress the .vcf.gz
        decomp    = ['gzip','-d',sub_dir+'breakseq.vcf.gz']
        #copy up to ../
        copy      = ['cp',sub_dir+'breakseq.vcf',out_names['.vcf']]
        #delete the SID folder and all contents
        clean     = ['rm','-rf',sub_dir]
                      
        self.db_start(run_id,','.join(in_names['.bam']))        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print(" ".join(call))
            output += subprocess.check_output(' '.join(call),
                                              stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(decomp),
                                              stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(copy),
                                              stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(clean),
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
