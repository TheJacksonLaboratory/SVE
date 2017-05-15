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
                     '.cnv.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[2]}

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
        
        #[2] bam file list is needed        
        bams = SV_TMPDIR+'/bam_files.list'
        bam_names = '\n'.join(inputs['.bam'])
        with open(bams,'w') as f: f.write(bam_names) #try comma, tab and newline

        md   = SV_TMPDIR+'/md'
        if not os.path.exists(md): os.makedirs(md)
        rd   = SV_TMPDIR+'/run'
        if not os.path.exists(rd): os.makedirs(rd)
        logs = SV_TMPDIR+'/logs'
        if not os.path.exists(logs): os.makedirs(logs)
        scheduler = []

        gender_map = SV_TMPDIR + '/gender.map'
        if os.path.isfile(gender_map): os.remove(gender_map)  
        for b in inputs['.bam']:
          print b
          sample_name_command = [self.tools['SAMTOOLS'], 'view', '-H', b]
          sample_name_command += ['| grep "RG" | awk -F\'\\t\' \'{for (i=1;i<NF;++i) if ($i ~ \"SM:\") print $i}\' | awk -F\':\' \'{print $2}\' | sort | uniq']
          sample_name_command += ['| awk \'{print $1\"\\t\"\"M\"}\'']
          sample_name_command += ['>>', gender_map]
          subprocess.check_output(' '.join(sample_name_command), stderr=subprocess.STDOUT, shell=True)
        
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
                      '-I %s' %bams]
        if inputs['bundle_dir'] != None: preprocess += ['-rmd', inputs['bundle_dir']]
        else:
            if inputs['.svmask.fasta'] != None:   preprocess += ['-genomeMaskFile ' + inputs['.svmask.fasta']]
            if inputs['.ploidymap.txt'] != None:  preprocess += ['-ploidyMapFile ' + inputs['.ploidymap.txt']]
            if inputs['.rdmask.bed'] != None:     preprocess += ['-readDepthMaskFile ' + inputs['.rdmask.bed']]
            if inputs['.gendermask.bed'] != None: preprocess += ['-genderMaskBedFile ' + inputs['.gendermask.bed']]
            if inputs['.gcmask.fasta'] != None:   preprocess += ['-copyNumberMaskFile ' + inputs['.gcmask.fasta']]

        preprocess += ['-run']
        
        gs_bwa_path = self.tools['GENOME_STRIP_PATH'] + '/bwa'
        PATH = self.tools['JAVA-1.8_PATH'] + ':' + self.tools['SAMTOOLS_PATH'] + ':' + self.tools['BCFTOOLS_PATH'] + ':' + self.tools['HTSLIB_PATH'] + ':' + self.tools['R_PATH'] + '/bin'
        if os.environ.has_key('PATH'): PATH += ':' + os.environ['PATH']
        LD_LIB = gs_bwa_path 
        if os.environ.has_key('LD_LIBRARY_PATH'): LD_LIB += ':' + os.environ['LD_LIBRARY_PATH']
       
        output,err = '', {}
        try: 
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(preprocess))
            output = subprocess.check_output(' '.join(preprocess), stderr=subprocess.STDOUT, shell=True,
                                            env={'SV_DIR': self.tools['GENOME_STRIP_PATH'], 'LD_LIBRARY_PATH': LD_LIB, 'PATH': PATH})

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

        #[1] Initial Pooled Deletion Discovery
        dd = sv+'/qscript/SVDiscovery.q'
        del_discovery = [java,'-cp %s'%classpath,
                         qcmd,
                         '-S %s' %dd,
                         '-S %s' %qs,
                         '-cp %s' %classpath,
                         '-gatk %s' %gatk,
                         '-configFile %s' %conf,
                         '-R %s' %inputs['.fa'],
                         '-I %s'%bams,
                         '-md %s' %md,
                         '-runDirectory %s' %rd,
                         '-jobLogDir %s' %logs,
                         #'-minimumSize %s' %100,
                         #'-maximumSize %s' %1000000,
                         #'-genderMapFile %s' %gender_map,
                         '-suppressVCFCommandLines',
                         #'-P select.validateReadPairs:false',
                         '-O %s'%out_names['.del.vcf']]

        if inputs['bundle_dir'] != None: del_discovery += ['-rmd', inputs['bundle_dir']]
        else:
            if inputs['.svmask.fasta'] != None:  del_discovery += ['-genomeMaskFile ' + inputs['.svmask.fasta']]
            if inputs['.ploidymap.txt'] != None: del_discovery += ['-ploidyMapFile ' + inputs['.ploidymap.txt']]
        del_discovery += ['-run']
        try: 
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(del_discovery))
            output = subprocess.check_output(' '.join(del_discovery), stderr=subprocess.STDOUT, shell=True,
                                             env={'SV_DIR': self.tools['GENOME_STRIP_PATH'], 'LD_LIBRARY_PATH': LD_LIB, 'PATH': PATH})
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
                          #'-genomeMaskFile %s' %gmask,
                          #'-genderMapFile %s' %gender_map,
                          '-I %s' %bams,
                          '-vcf %s' %out_names['.del.vcf'],
                          '-O %s' %out_names['.del.genotype.vcf']]
                          #'-parallelRecords 100',
        
        if inputs['bundle_dir'] != None: del_genotyping += ['-rmd', inputs['bundle_dir']]
        else:
            if inputs['.svmask.fasta'] != None:  del_genotyping += ['-genomeMaskFile ' + inputs['.svmask.fasta']]
            if inputs['.ploidymap.txt'] != None: del_genotyping += ['-ploidyMapFile ' + inputs['.ploidymap.txt']]
        del_genotyping += ['-run']

        try:
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(del_genotyping))
            output = subprocess.check_output(' '.join(del_genotyping), stderr=subprocess.STDOUT, shell=True,
                                    env={'SV_DIR': self.tools['GENOME_STRIP_PATH'], 'LD_LIBRARY_PATH': LD_LIB, 'PATH': PATH})
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
                         '-jobRunner Shell',
                         '-gatkJobRunner Shell']
                         #'-rmd /home/leew/SVE/data/Homo_sapiens_assembly19/',
        if inputs['bundle_dir'] != None: cnv_discovery += ['-rmd', inputs['bundle_dir']]
        else:
            if inputs['.svmask.fasta'] != None:  cnv_discovery += ['-genomeMaskFile ' + inputs['.svmask.fasta']]
            if inputs['.ploidymap.txt'] != None: cnv_discovery += ['-ploidyMapFile ' + inputs['.ploidymap.txt']]
        cnv_discovery += ['-run']
        try:
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print (' '.join(cnv_discovery))
            output += subprocess.check_output(' '.join(cnv_discovery), stderr=subprocess.STDOUT, shell=True,
                                   env={'SV_DIR': self.tools['GENOME_STRIP_PATH'], 'LD_LIBRARY_PATH': LD_LIB, 'PATH': PATH})
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
        
        final_vcf = rd+'/S14.vcf'
        #[3b]check results--------------------------------------------------
        #self.db_stop(run_id,{'output':output},'',True)
        results = [final_vcf]
        #for i in results: print i
        if all([os.path.exists(r) for r in results]):
            print("GenomeSTRiP sucessfull........")
            return results   #return a list of names
        else:
            print("GenomeSTRiP failure...........")
            return False
