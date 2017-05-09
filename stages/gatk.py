import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class gatk(stage_wrapper.Stage_Wrapper):
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
    #./software/jre1.8.0_25/bin/java -jar ./software/GATK_3.3/GenomeAnalysisTK.jar 
    #-T HaplotypeCaller -R ./node_data/rg1.fa -I ./node_data/node02/rg1_R1_S1_S8.bam 
    #--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 
    #-o ./node_data/node02/rg1_S1_S8_S13.vcf"""
    def run(self,run_id,inputs):
        #workflow is to run through the stage correctly and then check for error handles
    
        #[1a]get input names and output names setup

        #[1a]get input names and output names setup
        if ('.fa' not in inputs) or ('.bam' not in inputs) or ('out_dir' not in inputs):
            print "ERROR: .fa, .bam, and out_dir are required for genome_strip.py"
            return None
                    
        #will have to figure out output file name handling
        out_exts = self.split_out_exts()
        out_dir = inputs['out_dir'] + '/'
        stripped_name = self.strip_path(self.strip_in_ext(inputs['.bam'][0],'.bam'))
        out_names = {'.vcf'   : out_dir+stripped_name+'_S'+str(self.stage_id) + '.vcf',
                     '.g.vcf' : out_dir+stripped_name+'_S'+str(self.stage_id) + '.g.vcf'}
        
        #[2a]build command args
        java = self.tools['JAVA-1.8']
        gatk = self.tools['GATK']
        call    = [java, '-jar', gatk, '-T', 'HaplotypeCaller',
                   '-R', inputs['.fa'], '-I', inputs['.bam'][0], '-o', out_names['.g.vcf'],
                   '-ERC', 'GVCF', '-variant_index_type', 'LINEAR', '-variant_index_parameter', str(128000)]
        combine = [java, '-jar', gatk, '-T', 'CombineGVCFs',
                   '-R', inputs['.fa'], '-V', out_names['.g.vcf'], '-o', out_names['.vcf']]
        
        #[3a]execute the command here----------------------------------------------------
        print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
        print (' '.join(call))
        subprocess.check_output(' '.join(call), stderr=subprocess.STDOUT, shell=True)
        print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
        print (' '.join(combine))
        subprocess.check_output(' '.join(combine), stderr=subprocess.STDOUT, shell=True)
        
        #[3b]check results--------------------------------------------------
        results = [out_names['.vcf']]
        #for i in results: print i
        if all([os.path.exists(r) for r in results]):
            print("GATK sucessfull........")
            return results   #return a list of names
        else:
            print("GATK failure...........")
            return False
