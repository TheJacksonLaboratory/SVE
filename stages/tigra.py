import os
import sys
import subprocess32 as subprocess
import utils.tigra2vcf as tv
sys.path.append('../') #go up one in the modules
import stage_wrapper
import read_utils as ru

#function for auto-making svedb stage entries and returning the stage_id
class tigra(stage_wrapper.Stage_Wrapper):
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
        #[1a]get input names and output names setup
        if ('.fa' not in inputs) or ('.bam' not in inputs) or ('out_dir' not in inputs) or ('.vcf' not in inputs):
            print "ERROR: .fa, .bam, .vcf and out_dir are required for genome_strip.py"
            return None
                    
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        out_dir = inputs['out_dir'] + '/'
        stripped_name = ''
        if len(inputs['.bam']) == 1: stripped_name = self.strip_path(self.strip_in_ext(inputs['.bam'][0],'.bam'))
        else: stripped_name = 'joint'
        out_names = {'.fa' : out_dir+'/tigra.ctg.fa',
                     '.bed': [out_dir+'/a.bed',out_dir+'/b.bed',out_dir+'/intersect.bed'],
                     '.vcf': out_dir+'/'+stripped_name+'_S'+str(self.stage_id)+out_exts[2]} #will overwrite the vcf

        sub_dir = out_dir+stripped_name+'S'+str(self.stage_id)+'/'
        if not os.path.exists(sub_dir): os.makedirs(sub_dir)
        if not os.path.exists(out_dir): os.makedirs(out_dir)

        #[2a]build command args
        tigra_ext = self.tools['TIGRA-EXT']
        tigra_dir = self.tools['TIGRA_PATH'] + '/'
        #tigra = soft+'/tigra/tigra-sv'
        bwa = self.tools['BWA']
        samtools = self.tools['SAMTOOLS']
        bedtools = self.tools['BEDTOOLS']
        sample  = inputs['.vcf'].replace('.vcf','').replace('.calls','').rsplit('/')[-1].rsplit('_')[0].split('.')[0]
        bamlist = sub_dir+'/'+sample+'.bamlist'
        with open(bamlist,'w') as f:
            for bam in inputs['.bam']:
                f.write(sample+':'+bam+'\n')
        #tigra -b -R NCBI36.example.fa -o output1.fa example.breakdancer.sv example1.bam 2> output1.log
        #perl TIGRA-ext.pl -k 15,25,35, -T /home/tbecker/software/tigra/ -r fasta -f bamlist [-b] .sv(bd or vcf) file...
        fullp = ['perl', tigra_ext, '-r', inputs['.fa'],'-T', tigra_dir, '-I',bedtools,'-B',bwa, #exe paths
                 '-d',sub_dir,'-o',out_names['.vcf'], inputs['.vcf']] + inputs['.bam']   #ref and out paths
        sam2bam  = [samtools, 'view', '-Sb',sub_dir+'/tigra.sam']
        bamsort  = [samtools, 'sort', '-o', sub_dir+'/tigra.sorted.bam', '-']
        bamindex = [samtools, 'index', sub_dir+'/tigra.sorted.bam']
        LD_LIB = self.tools['HTSLIB_PATH'] + ':' + os.environ['LD_LIBRARY_PATH']
        
#        if inputs['.calls'].endswith('.calls'):
#            fullp += ['-b'] #set breakdancer format flag
        output,err = '',{}
        
        #step [1] run tigra-ext
        try: #currently set for one sample targeted assembly
            print(' '.join(fullp))
            #check for tigra-ext files and skip if possible
            output += subprocess.check_output(' '.join(fullp),stderr=subprocess.STDOUT,shell=True, env={'LD_LIBRARY_PATH':LD_LIB})
            output += subprocess.check_output(' '.join(sam2bam, ['|'], bamsort),stderr=subprocess.STDOUT,shell=True)
            output += subprocess.check_output(' '.join(bamindex),stderr=subprocess.STDOUT,shell=True)
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
        print('tigra-ext call: \n'+output)
        
        #[3b]check results--------------------------------------------------
        if err == {}: #check all the tigra expected output files
            results = out_names['.bed']+[out_names['.fa'],out_names['.vcf']]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("sucessfull........")
                return results   #return a list of names
            else:
                print("failure...........")
                return False
        else:
            print("failure...........")
            print(err['message'])
            return None
