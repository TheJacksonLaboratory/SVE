import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class genome_strip(stage_wrapper.Stage_Wrapper):
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
        out_names = {'.del.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0],
                     '.del.genotype.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[1],
                     '.cnv.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[2],}

        #[2a]build command args
        #environment variable passing here
        SV_DIR = self.tools['GENOME_STRIP_PATH']
        SV_TMPDIR = out_dir+stripped_name+'_S'+str(self.stage_id)+'/temp'
        if not os.path.exists(SV_TMPDIR): os.makedirs(SV_TMPDIR)

        #reused paths and files...
        sv = self.tools['GENOME_STRIP_PATH']
        classpath = sv + '/lib/SVToolkit.jar:' + sv + '/lib/gatk/GenomeAnalysisTK.jar:' + sv + '/lib/gatk/Queue.jar'
        java  = self.tools['JAVA-1.8'] + ' -Xmx4g'
        qcmd  = 'org.broadinstitute.gatk.queue.QCommandLine'
        qs    = sv + '/qscript/SVQScript.q'
        gatk  = sv + '/lib/gatk/GenomeAnalysisTK.jar'
        conf  = sv + '/conf/genstrip_parameters.txt'
        pp    = sv+'/qscript/SVPreprocess.q'
        gmask, ploidy, rdmask = '', '', ''
#        cnmask = in_names['.gcmask.fasta']
        
        sub_dir = out_dir+stripped_name+'_S'+str(self.stage_id)+'/'
        if not os.path.exists(sub_dir): os.makedirs(sub_dir)
        
        #[2] bam file list is needed        
        bams = sub_dir+'/bam_files.list'
        bam_names = '\n'.join(inputs['.bam'])
        with open(bams,'w') as f: f.write(bam_names) #try comma, tab and newline

        #[3] gender_map this is for each sample...?
        #for each sample do a samtools view -SH pull out the RG sample name
        gender_map = sub_dir+'/sample_gender.map' #1 is Paternal 2 is Maternal?
        s = ''
        for bam in inputs['.bam']:
            s += bam.split('/')[-1].split('.')[0]+' '+'1\n' #check to see if this will work
        with open(gender_map,'w') as f: #write the file
            f.write(s)
        
        md   = sub_dir+'/md'
        if not os.path.exists(md): os.makedirs(md)
        rd   = sub_dir+'/run'
        if not os.path.exists(rd): os.makedirs(rd)
        logs = sub_dir+'/logs'
        if not os.path.exists(logs): os.makedirs(logs)
        scheduler = []
            
        #[0] Preprocess The Bam Data and Generate MetaData
        preprocess = [java+' -cp ' + classpath,
                      qcmd,
                      '-S %s' %pp,
                      '-S %s' %qs,
                      '-md %s' %md,
                      '-jobLogDir %s' %logs,
                      '-cp %s' %classpath,
                      '-gatk %s' %gatk,
                      '-R %s' %inputs['.fa'],
                      '-configFile %s' %conf,
                      '-reduceInsertSizeDistributions true',
                      '-computeGCProfiles true',
                      '-bamFilesAreDisjoint true',
                      '-rmd /home/leew/SVE/data/Homo_sapiens_assembly19/',
                      '-I %s' %bams,
                      '-run']
        """
        if '.svmask.fasta' in inputs:   preprocess += ['-genomeMaskFile ' + inputs['.svmask.fasta']]
        if '.ploidymap.txt' in inputs:  preprocess += ['-ploidyMapFile ' + inputs['.ploidymap.txt']]
        if '.rdmask.bed' in inputs:     preprocess += ['-readDepthMaskFile ' + inputs['.rdmask.bed']]
        if '.gendermask.bed' in inputs: preprocess += ['-genderMaskBedFile ' + inputs['.gendermask.bed']]
        if '.gcmask.fasta' in inputs:   preprocess += ['-copyNumberMaskFile ' + inputs['.gcmask.fasta']]
        """
        #script = self.generate_bash_script() + '\n' + (' '.join(preprocess))
        #with open(rd+'/preprocess.sh','w') as f: f.write(script)
        
        gs_bwa_path = self.tools['GENOME_STRIP_PATH'] + '/bwa'
        PATH = self.tools['JAVA-1.8_PATH'] + ':' + self.tools['SAMTOOLS_PATH'] + ':' + self.tools['BCFTOOLS_PATH'] + ':' + self.tools['HTSLIB_PATH']
        if os.environ.has_key('PATH'): PATH += ':' + os.environ['PATH']
        LD_LIB = gs_bwa_path 
        if os.environ.has_key('LD_LIBRARY_PATH'): LD_LIB += ':' + os.environ['LD_LIBRARY_PATH']
        output, err = '', {}
        
        try:
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(preprocess))
            output = subprocess.check_output(' '.join(preprocess), stderr=subprocess.STDOUT, shell=True,
                                             env={'SV_DIR': self.tools['GENOME_STRIP_PATH'], 'LD_LIBRARY_PATH': LD_LIB, 'PATH': PATH})
        # catch all errors that arise under normal call behavior
        except subprocess.CalledProcessError as E:
            print('call error: ' + E.output)  # what you would see in the term
            err['output'] = E.output
            # the python exception issues (shouldn't have any...
            print('message: ' + E.message)  # ?? empty
            err['message'] = E.message
            # return codes used for failure....
            print('code: ' + str(E.returncode))  # return 1 for a fail in art?
            err['code'] = E.returncode
        except OSError as E:
            print('os error: ' + E.strerror)  # what you would see in the term
            err['output'] = E.strerror
            # the python exception issues (shouldn't have any...
            print('message: ' + E.message)  # ?? empty
            err['message'] = E.message
            # the error num
            print('code: ' + str(E.errno))
            err['code'] = E.errno
        print('output:\n' + output)
        
        #[1] Initial Pooled Deletion Discovery
        dd = sv+'/qscript/SVDiscovery.q'
        del_discovery = [java,'-cp %s'%classpath,
                         qcmd,
                         '-S %s' %dd,
                         '-S %s' %qs,
                         '-gatk %s' %gatk,
                         '-cp %s' %classpath,
                         '-configFile %s' %conf,
                         '-runDirectory %s' %rd,
                         '-R %s' %inputs['.fa'],
                         '-md %s' %md,
                         '-jobLogDir %s' %logs,
                         '-minimumSize %s' %100,
                         '-maximumSize %s' %1000000,
                         '-genomeMaskFile %s' %gmask,
                         #'-genderMapFile %s' %gender_map,
                         '-suppressVCFCommandLines',
                         '-P select.validateReadPairs:false',
                         '-I %s'%bams,
                         '-O %s'%out_names['.del.vcf'],
                         '-run']

        if '.svmask.fasta' in inputs:  del_discovery += ['-genomeMaskFile ' + inputs['.svmask.fasta']]
        if '.ploidymap.txt' in inputs: del_discovery += ['-ploidyMapFile ' + inputs['.ploidymap.txt']]

        """
        try:
            output = subprocess.check_output(' '.join(del_discovery), stderr=subprocess.STDOUT, shell=True,
                                             env={'SV_DIR': self.tools['GENOME_STRIP_PATH'], 'LD_LIBRARY_PATH': LD_LIB, 'PATH': PATH})
        # catch all errors that arise under normal call behavior
        except subprocess.CalledProcessError as E:
            print('call error: ' + E.output)  # what you would see in the term
            err['output'] = E.output
            # the python exception issues (shouldn't have any...
            print('message: ' + E.message)  # ?? empty
            err['message'] = E.message
            # return codes used for failure....
            print('code: ' + str(E.returncode))  # return 1 for a fail in art?
            err['code'] = E.returncode
        except OSError as E:
            print('os error: ' + E.strerror)  # what you would see in the term
            err['output'] = E.strerror
            # the python exception issues (shouldn't have any...
            print('message: ' + E.message)  # ?? empty
            err['message'] = E.message
            # the error num
            print('code: ' + str(E.errno))
            err['code'] = E.errno
        print('output:\n' + output)
        """
        
        #[2] Genotype Individual Deleteions (this needs the GS_DEL_VCF_splitter.py)
        dg = sv+'/qscript/SVGenotyper.q'
        del_genotyping = [java,'-cp %s'%classpath,
                          qcmd,
                          '-S %s' %dg,
                          '-S %s' %qs,
                          '-gatk %s' %gatk,
                          '-cp %s' %classpath,
                          '-configFile %s' %conf,
                          '-tempDir %s' %SV_TMPDIR,
                          '-R %s' %inputs['.fa'],
                          '-runDirectory %s' %rd,
                          '-md %s' %md,
                          '-jobLogDir %s' %logs,
                          '-genomeMaskFile %s' %gmask,
                          '-genderMapFile %s' %gender_map,
                          '-I %s' %bams,
                          '-vcf %s' %out_names['.del.vcf'],
                          '-O %s' %out_names['.del.genotype.vcf'],
                          '-parallelRecords 100',
                          '-run']

        if '.svmask.fasta' in inputs:  del_genotyping += ['-genomeMaskFile ' + inputs['.svmask.fasta']]
        if '.ploidymap.txt' in inputs: del_genotyping += ['-ploidyMapFile ' + inputs['.ploidymap.txt']]
        """
        try:
            output = subprocess.check_output(' '.join(del_genotyping), stderr=subprocess.STDOUT, shell=True,
                                             env={'SV_DIR': self.tools['GENOME_STRIP_PATH'], 'LD_LIBRARY_PATH': LD_LIB, 'PATH': PATH})
        # catch all errors that arise under normal call behavior
        except subprocess.CalledProcessError as E:
            print('call error: ' + E.output)  # what you would see in the term
            err['output'] = E.output
            # the python exception issues (shouldn't have any...
            print('message: ' + E.message)  # ?? empty
            err['message'] = E.message
            # return codes used for failure....
            print('code: ' + str(E.returncode))  # return 1 for a fail in art?
            err['code'] = E.returncode
        except OSError as E:
            print('os error: ' + E.strerror)  # what you would see in the term
            err['output'] = E.strerror
            # the python exception issues (shouldn't have any...
            print('message: ' + E.message)  # ?? empty
            err['message'] = E.message
            # the error num
            print('code: ' + str(E.errno))
            err['code'] = E.errno
        print('output:\n' + output)
        """
        #[3] GenomeSTRiP2.0 CNV algorithm (this needs the gs_slpit_merge.py)
        cnv = sv+'/qscript/discovery/cnv/CNVDiscoveryPipeline.q'
        cnv_discovery = [java,'-cp %s'%classpath,
                         qcmd,
                         '-S %s' %cnv,
                         '-S %s' %qs,
                         '-gatk %s' %gatk,
                         '-cp %s' %classpath,
                         '-configFile %s' %conf,
                         '-R %s' %inputs['.fa'],
                         '-runDirectory %s' %rd,
                         '-md %s' %md,
                         '-jobLogDir %s' %logs,
                         '-genderMapFile %s' %gender_map,
#                          '-ploidyMapFile %s'%ploidy,
                         '-I %s' %bams,
                         '-tilingWindowSize %s' %5000,
                         '-tilingWindowOverlap %s' %2500,
                         '-maximumReferenceGapLength %s' %25000,
                         '-boundaryPrecision %s' %200,
                         '-minimumRefinedLength %s' %2500,
                         '-run']

        if '.svmask.fasta' in inputs:  cnv_discovery += ['-genomeMaskFile ' + inputs['.svmask.fasta']]
        if '.ploidymap.txt' in inputs: cnv_discovery += ['-ploidyMapFile ' + inputs['.ploidymap.txt']]
        """
        try:
            output = subprocess.check_output(' '.join(cnv_discovery), stderr=subprocess.STDOUT, shell=True,
                                             env={'SV_DIR': self.tools['GENOME_STRIP_PATH'], 'LD_LIBRARY_PATH': LD_LIB, 'PATH': PATH})
        # catch all errors that arise under normal call behavior
        except subprocess.CalledProcessError as E:
            print('call error: ' + E.output)  # what you would see in the term
            err['output'] = E.output
            # the python exception issues (shouldn't have any...
            print('message: ' + E.message)  # ?? empty
            err['message'] = E.message
            # return codes used for failure....
            print('code: ' + str(E.returncode))  # return 1 for a fail in art?
            err['code'] = E.returncode
        except OSError as E:
            print('os error: ' + E.strerror)  # what you would see in the term
            err['output'] = E.strerror
            # the python exception issues (shouldn't have any...
            print('message: ' + E.message)  # ?? empty
            err['message'] = E.message
            # the error num
            print('code: ' + str(E.errno))
            err['code'] = E.errno
        print('output:\n' + output)
        """
        
        # #then some renaming, conversion and clean up
        # move_vcf = ['mv',sub_dir+'/results/gs_cnv.genotypes.vcf.gz'] #this seems to be hard coded
        # convert_vcf = ['python','gs_slpit_merge.py']
        # clean = ['rm','-rf',SV_TMPDIR, sub_dir] #delete temp, stage_id folder after vcfs are copied
        
        final_vcf = rd+'/S14.vcf'
        #[3b]check results--------------------------------------------------
        if err == {}:
            #self.db_stop(run_id,{'output':output},'',True)
            results = [final_vcf]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("sucessfull........")
                return results   #return a list of names
            else:
                print("failure...........")
                return False
        else:
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
