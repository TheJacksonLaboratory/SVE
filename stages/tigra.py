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
        #if self.db_get_ref_name(run_id): ref_name = self.ref_name        
        in_names = {'.fa':inputs['.fa'][0],'.calls':inputs['.calls'][0],
                    '.vcf':inputs['.vcf'][0],'.bam':inputs['.bam'][0]}
        out_exts = self.split_out_exts()
        if inputs.has_key('out_dir'):
            out_dir = inputs['out_dir'][0]
            stripped_name = self.strip_path(self.strip_in_ext(in_names['.bam'],'.bam'))
            out_names = {'.fa' : out_dir+'/tigra.ctg.fa',
                         '.bed': [out_dir+'/a.bed',out_dir+'/b.bed',out_dir+'/intersect.bed'],
                         '.vcf': out_dir+'/'+stripped_name+'_S'+str(self.stage_id)+out_exts[2]} #will overwrite the vcf
        else:
            right = in_names['.bam'].rsplit('/')[-1]
            out_dir = in_names['.bam'].replace(right,'')
            cascade = self.strip_in_ext(in_names['.bam'],'.bam')
            out_names = {'.fa' :out_dir+'/tigra.ctg.fa',
                         '.bed':[out_dir+'/a.bed',out_dir+'/b.bed',out_dir+'/intersect.bed'],
                         '.vcf':out_dir+'/'+cascade+'_S'+str(self.stage_id)+out_exts[2]}  #will overwrite the vcf               
        #[2a]build command args
        soft = self.software_path
        tigra_ext = soft+'/tigra-ext/TIGRA-ext.pl'
        tigra_dir = soft+'/tigra/'
        #tigra = soft+'/tigra/tigra-sv'
        bwa = soft+'/bwa-0.7.10/bwa'
        samtools = soft+'/samtools-0.1.19/samtools' #old samtools is needed     
        bedtools = soft+'/bedtools2.24/bin/intersectBed'
        sample  = in_names['.vcf'].replace('.vcf','').replace('.calls','').rsplit('/')[-1].rsplit('_')[0].split('.')[0]
        bamlist = out_dir+'/'+sample+'.bamlist'
        if not os.path.exists(out_dir): os.makedirs(out_dir)
        with open(bamlist,'w') as f:
            f.write(sample+':'+in_names['.bam']+'\n')
        F = str(self.get_params()['F']['value'])
        L = str(self.get_params()['L']['value'])
        A = str(self.get_params()['A']['value'])
        Q = str(self.get_params()['Q']['value'])
        P = str(self.get_params()['P']['value'])
        H = str(self.get_params()['H']['value'])
        threads = str(self.get_params()['p']['value'])
        #tigra -b -R NCBI36.example.fa -o output1.fa example.breakdancer.sv example1.bam 2> output1.log
        #perl TIGRA-ext.pl -k 15,25,35, -T /home/tbecker/software/tigra/ -r fasta -f bamlist [-b] .sv(bd or vcf) file...
        fullp = ['perl', tigra_ext, '-r', in_names['.fa'],'-T', tigra_dir, '-I',bedtools,'-B',bwa, #exe paths
                 '-k','15,25','-d',out_dir,'-o',out_names['.vcf'],'-F',F,'-p',P,'-h',H]  #ref and out paths
        sam2bam  = [samtools, 'view', '-Sb',out_dir+'/tigra.sam', '-@', threads,'-o',out_dir+'/tigra.bam']
        bamsort  = [samtools, 'sort', '-@', threads, '-f',out_dir+'/tigra.bam' ,out_dir+'/tigra.sorted.bam']
        bamindex = [samtools, 'index', out_dir+'/tigra.sorted.bam']
        clean    = ['rm', out_dir+'/tigra.sam',out_dir+'/tigra.bam']
        
        if in_names['.calls'].endswith('.calls'):
            fullp += ['-b'] #set breakdancer format flag
        fullp += [in_names['.vcf'],in_names['.bam']]                          #input vcf and bam file 
#        assem = [tigra, '-b', '-R', in_names['.fa'], '-o', out_names['.fa'],'-c','20',#only chrom one
#                 in_names['.calls'],in_names['.bam']]
        self.db_start(run_id,in_names['.bam'])        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        #step [1] run tigra-ext
        try: #currently set for one sample targeted assembly
            print(' '.join(fullp))
            #check for tigra-ext files and skip if possible
            if not os.path.exists(out_dir+'/tigra.sam'):
                output += subprocess.check_output(' '.join(fullp),stderr=subprocess.STDOUT,shell=True)
            else:
                output += '\n::: '+out_dir+'/tigra.sam was found, skipping tigra-ext step\n'
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
        #step [2] samtools sam to bam to sorted and indexed bam file for IGV
        try:
            #if the tigra-ext files are present skip to this section and process the contig alignment
            #if the contig alignment is complete skip to the final vcf fixing stage...
            if not os.path.exists(out_dir+'/tigra.bam'):
                output += subprocess.check_output(' '.join(sam2bam),stderr=subprocess.STDOUT,shell=True)
            else:
                output += '\n::: '+out_dir+'/tigra.bam was found, skipping samtools view step\n'
            if not os.path.exists(out_dir+'/tigra.sorted.bam'):
                output += subprocess.check_output(' '.join(bamsort),stderr=subprocess.STDOUT,shell=True)
            else:
                output += '\n::: '+out_dir+'/tigra.sorted.bam was found, skipping samtools sort step\n'
            if not os.path.exists(out_dir+'/tigra.sorted.bam.bai'):
                output += subprocess.check_output(' '.join(bamindex),stderr=subprocess.STDOUT,shell=True)
            else:
                output += '\n::: '+out_dir+'/tigra.sorted.bam.bai was found, skipping samtools index step\n'
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
        #step [3] repair/convert the VCF file
        try:
            #fix the tigra_output and overwrite the partial vcf-----------------------------------
            #produce a quick report of how many SV calls were confirmed...
            ref_seq = {in_names['.fa'].rsplit('/')[-1].rsplit('.fa')[-1]:ru.read_fasta(in_names['.fa'],True)}
            header_path = os.path.dirname(os.path.abspath(__file__))+'/../data/header_template.vcf' 
            output += str(tv.tigra_ext_bed_to_vcf(out_names['.bed'][2],sample,ref_seq,out_names['.vcf'],header_path))
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
        
        #[3b]check results--------------------------------------------------
        if err == {}: #check all the tigra expected output files
            results = out_names['.bed']+[out_names['.fa'],out_names['.vcf']]
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
            print(err['message'])
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None