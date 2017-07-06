import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper
from stages.utils.CheckVcf import GetCallCount

#function for auto-making svedb stage entries and returning the stage_id
class delly(stage_wrapper.Stage_Wrapper):
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

    def call(command, output):
      output = subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT,shell=True)+'\n'
      return output

    #override this function in each wrapper...
    #this is an older pipeline and not for somatic, will do a new somatic
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
        stripped_name = self.strip_path(self.strip_in_ext(inputs['.bam'][0],'.bam'))
        out_names = {'.vcf' :out_dir+stripped_name+'_S'+str(self.stage_id)+out_exts[0]}
        #[2a]build command args
        sub_dir = out_dir+stripped_name+'_S'+str(self.stage_id)+'/'
        if not os.path.exists(sub_dir): os.makedirs(sub_dir)
        
        #add load libs parameters for OPEN_MP to do || processing        
        #will have to make some connection changes here
        delly = self.tools['DELLY']
        excl = ''
        if inputs['genome'] == 'hg19': excl  = self.files['DELLY-HG19']
        elif inputs['genome'] == 'hg38': excl  = self.files['DELLY-HG38']

        bcfs = {}
        type_list = ['del', 'dup', 'inv', 'bnd', 'ins']
        for type in type_list:
            bcfs[type] = sub_dir+type+'.bcf'

        #self.db_start(run_id,in_names['.bam'][0])        
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try: #should split these up for better robustness...
            count = 0
            # Delly call
            #if threads in inputs: p1 = mp.Pool(processes = inputs['threads'])
            #p1 = mp.Pool(processes = 1)
            for bam in inputs['.bam']:
                delly_call = [delly, 'call', '-g', inputs['.fa'], '-n']
                if excl != '': delly_call += ['-x', excl]
                for type in type_list:
                    type_call = delly_call + ['-t', type.upper(), '-o', sub_dir+str(count)+'.'+type+'.bcf'] + [bam]
                    print (" ".join(type_call))
                    #p1.apply_async(call,args=(type_call, output),callback=collect_results)
                    output += subprocess.check_output(' '.join(type_call), stderr=subprocess.STDOUT,shell=True)+'\n'
                count += 1

            # Delly merge
            if count > 1:
                delly_merge = [delly, 'merge', '-r', str(0.5), '-b', str(500), '-n', str(1000000), '-m', str(500)]
                for type in type_list:
                    type_merge = delly_merge + ['-t', type.upper(), '-o', 'b_geno_'+bcfs[type]]
                    for i in range (count): 
                        type_merge += [sub_dir+str(count)+'.'+type+'.bcf']
                    output += subprocess.check_output(' '.join(type_merge), stderr=subprocess.STDOUT,shell=True)+'\n'

                # Delly renotype
                for bam in inputs['.bam']:
                    delly_geno = [delly, 'call', '-g',inputs['.fa']]
                    if excl != '': delly_geno += ['-x', excl]
                    for type in type_list:
                        type_geno = delly_geno + ['-v', 'b_geno_'+bcfs[type], '-t', type.upper(), '-o', sub_dir+str(count)+'.'+type+'.geno.bcf'] + [bam]
                        print (" ".join(type_geno))
                        output += subprocess.check_output(' '.join(type_geno), stderr=subprocess.STDOUT,shell=True)+'\n'

                # Merge regeno bcf
                delly_geno_merge = [delly, 'merge', '-m', 'id', '-O', 'b']
                for type in type_list:
                    type_geno_merge += ['-o', bcfs[type]]
                    for i in range (count):
                        type_geno_merge += [sub_dir+str(count)+'.'+type+'.geno.bcf']
                    output += subprocess.check_output(' '.join(type_geno_merge), stderr=subprocess.STDOUT,shell=True)+'\n'

            elif count == 1:
                for type in type_list:
                    bcfs[type] = sub_dir+str(count - 1)+'.'+type+'.bcf'

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
        print('output:\n'+output)

        #merge/filter all the calls into one .vcf with vcftools
        bcftools = self.tools['BCFTOOLS']
        concat = [bcftools,'concat','-a','-o',out_names['.vcf'],'-O','v']
        for type in type_list: concat += [bcfs[type]]

        try:
            print(' '.join(concat))
            output += subprocess.check_output(' '.join(concat),
                                   stderr=subprocess.STDOUT,shell=True)
            print('rm -rf %s'%sub_dir)
            subprocess.check_output('rm -rf %s'%sub_dir, stderr=subprocess.STDOUT,shell=True)                        
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
        print('output:\n'+output)
                                                
        #[3b]check results--------------------------------------------------
        if err != {}:
            print err
        if GetCallCount(out_names['.vcf']) > 0:
            print("<<<<<<<<<<<<<delly sucessfull>>>>>>>>>>>>>>>\n")
            return out_names['.vcf']   #return a list of names
        else:
            print("<<<<<<<<<<<<<delly failure>>>>>>>>>>>>>>>\n")
            return None
