import os
import sys
import time
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class hydra(stage_wrapper.Stage_Wrapper):
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
        #if self.db_get_ref_name(run_id): ref_name = self.ref_name        
        ref_name = inputs['.fa'].rsplit('/')[-1].rsplit('.')[0]
        
        out_exts = self.split_out_exts()
        out_dir = inputs['out_dir'] + '/'
        stripped_name = ''
        if len(inputs['.bam']) == 1: stripped_name = self.strip_path(self.strip_in_ext(inputs['.bam'][0],'.bam'))
        else: stripped_name = 'joint'
        sub_dir = out_dir+stripped_name+'_S'+str(self.stage_id)+'/'
        if not os.path.exists(sub_dir): os.makedirs(sub_dir)
        out_names = {'.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0]}
            
        #[a]use to run several sub scripts via command line/seperate process
        python = sys.executable
        hydra  = self.tools['HYDRA_PATH'] + '/'
        hydra_to_vcf = self.tools['SVE_HOME'] + '/stages/utils/hydra_to_vcf.py'
        #ENV
        PATH = hydra+'bin:'+hydra+'scripts:'+\
               self.tools['SAMTOOLS-0.1.19_PATH']
        if os.environ.has_key('PATH'): PATH += ':' + os.environ['PATH']

        #[0] stub file generation        
        bams = 'bam.stub'
        bam_names = '\n'.join(['sample%s'%i+'\t'+inputs['.bam'][i] for i in range(len(inputs['.bam']))])
        with open(bams,'w') as f: f.write(bam_names) #follow readme.md tenplate
        
        #[1] make a config file
        cfg = sub_dir+'bam.stub.config'
        #s is number of sample pairs, n is the max unit of variation
        #python scripts/make_hydra_config.py -i config.stub.txt > config.hydra.txt
        make_cfg = [python,hydra+'scripts/make_hydra_config.py','-i',bams,
                    '-s',str(int(1E5)),'-n',str(16),'>',cfg]
                    
        #[2] extract discordant alignments for each sample .bam file
        #python scripts/extract_discordants.py -c config.hydra.txt -d <sample_name>
        #--min_mapq=INT,--allow_dups=FLAG,--mem=INT=2E9?2GB-4GB?
        extract =  [python,hydra+'scripts/extract_discordants.py','-c',cfg,'-d']
        #[3] run hydra router
        #hydra-router -config config.hydra.txt -routedList routed-files.txt
        routed_bams = sub_dir+'bam.routed'
        route   =  ['hydra-router',
                    '-config',cfg,'-routedList',routed_bams]
        
        #[4] assemble SV breakpoint clusters
        #sh scripts/assemble-routed-files.sh routed-files-test.txt config.hydra.txt 1
        #assemble-routed-files.sh <config file> <routed file list file> <number of processes> <punt parameter>
        #punt should be 5x the average read depth over all samples
        assemble_command = hydra + 'scripts/assemble-routed-files.sh'
        assemble =[assemble_command,
                   cfg,routed_bams,str(1),str(60)]
        
        #[5] merge SV assembly files
        #sh scripts/combine-assembled-files.sh /full/path/to/assembled/files/ all.assembled
        asm = sub_dir+'all.assembled'
        merge_command = hydra + 'scripts/combine-assembled-files.sh'
        merge = [merge_command,'.', asm]
        
        #[6] finalize SV breakpoints
        #scripts/forceOneClusterPerPairMem.py -i all.assembled -o all.sv-calls
        svs = sub_dir+'all-sv.calls'
        cluster = [python,hydra+'scripts/forceOneClusterPerPairMem.py','-i',asm,'-o',svs]
        
        #[7] annotate SV breakpoints on samples
        #scripts/frequency.py -f all.sv-calls.final -d all.sv-calls.detail > all.sv-calls.freq
        freq_name = svs+'.freq'
        freqs = [python,hydra+'scripts/frequency.py',
                 '-c',cfg,'-f',svs+'.final','-d',svs+'.detail','>',freq_name]
        
        #[8] change footprint intervals into breakpoint intervals
        final_name = svs+'.final' #not sure if this is correct in general case...
        vcf_name   = svs+'.vcf'
        bkpts_name = svs+'.bkpts'
        #grep -v "#" all.hydra.sv.freq | python ~/bin/hydraToBreakpoint.py -i stdin > all.hydra.sv.bkpts
        bkpts = ['grep','-v','"#"',freq_name,'|',
                 python,hydra+'scripts/hydraToBreakpoint.py','-i','stdin','>',bkpts_name]
        #bkpts = [python,hydra+'scripts/hydraToBreakpoint.py','-i',freq_name,'>',bkpts_name]
        
        #[9] convert to VCF using the utils/hydra_to_vcf.py tool
        fasta_2bit = inputs['.fa']+'.2bit'
        bkpt2vcf = [python,hydra_to_vcf,final_name,fasta_2bit]#assumes a .2bit for ref is there...
        
        #[10] Copy out and Clean up files
        copy  = ['cp',vcf_name,out_names['.vcf']]
        clean = ['rm','-rf',sub_dir]
        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print('making the hydra configuration')
            print(' '.join(make_cfg))
            output += subprocess.check_output(' '.join(make_cfg),
                                              stderr=subprocess.STDOUT,shell=True,
                                              env={'PATH':PATH})+'\n'

                                              
            for k in ['sample%s'%i for i in range(len(inputs['.bam']))]:
                print('extracting discordants for %s'%k)
                print(' '.join(extract+[k]))
                output += subprocess.check_output(' '.join(extract+[k]),
                                                  stderr=subprocess.STDOUT,shell=True,
                                                  env={'PATH':PATH})+'\n'
                                                  
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print('routing all samples into hydra router')
            print(' '.join(route))
            output += subprocess.check_output(' '.join(route),
                                              stderr=subprocess.STDOUT,shell=True,
                                              env={'PATH':PATH})+'\n'
                                              
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print('combining hydra assembly files')
            print(' '.join(assemble))
            output += subprocess.check_output(' '.join(assemble),
                                              stderr=subprocess.STDOUT,shell=True,
                                              env={'PATH':PATH})+'\n'
                                              
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print('merging results')
            print(' '.join(merge))
            output += subprocess.check_output(' '.join(merge),
                                              stderr=subprocess.STDOUT,shell=True,
                                              env={'PATH':PATH})+'\n'
                                              
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print('starting hydra clustering')
            print(' '.join(cluster))
            output += subprocess.check_output(' '.join(cluster),
                                              stderr=subprocess.STDOUT,shell=True,
                                              env={'PATH':PATH})+'\n'
                                              
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print('computing hydra frequencies')
            print(' '.join(freqs))
            output += subprocess.check_output(' '.join(freqs),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
                                              
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print('converting hydra to vcf format')
            print(' '.join(bkpt2vcf))
            if not os.path.isfile(fasta_2bit):
                generate_fasta_2bit = [self.tools['FATO2BIT'], inputs['.fa'], fasta_2bit]
                output += subprocess.check_output(' '.join(generate_fasta_2bit), stderr=subprocess.STDOUT,shell=True)+'\n'
            output += subprocess.check_output(' '.join(bkpt2vcf),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
                                              
            print ("<<<<<<<<<<<<<SVE command>>>>>>>>>>>>>>>\n")
            print('copying files and cleaning sub directory')
            output += subprocess.check_output(' '.join(copy),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
            output += subprocess.check_output(' '.join(clean),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
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
        
        print('vcf file %s exists=%s'%(out_names['.vcf'],os.path.exists(out_names['.vcf'])))
        print('computing hydra breakpoints')
        print(' '.join(bkpts))
        try:
            output = subprocess.check_output(' '.join(bkpts),
                                              stderr=subprocess.STDOUT,shell=True)+'\n'
            #if os.path.exists(out_names['.vcf']):
            #    output += subprocess.check_output(' '.join(clean),
            #                                      stderr=subprocess.STDOUT,shell=True)
            print('all hydra stages completed')
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
        #print('output:\n'+output)                                        
        #[3b]check results--------------------------------------------------
        if err == {}:
            results = [out_names['.vcf']]
            #for i in results: print i
            if all([os.path.exists(r) for r in results]):
                print("<<<<<<<<<<<<<hydra sucessfull>>>>>>>>>>>>>>>\n")
                #self.db_stop(run_id,self.vcf_to_vca(out_names['.vcf']),'',True)
                return results   #return a list of names
            else:
                print("<<<<<<<<<<<<<hydra failure>>>>>>>>>>>>>>>\n")
                #self.db_stop(run_id,{'output':output},'',False)
                return False
        else:
            print("<<<<<<<<<<<<<hydra failure>>>>>>>>>>>>>>>\n")
            #self.db_stop(run_id,{'output':output},err['message'],False)
            return None
