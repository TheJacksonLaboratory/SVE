import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class sambamba_index(stage_wrapper.Stage_Wrapper):
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

    def run(self,run_id,inputs):
        sambamba = self.software_path+'/sambamba_v0.6.6'
        threads = str(self.get_params()['-t']['value'])
        input_bam = inputs['.bam'][0]
        command = [sambamba,'index','-t',threads,input_bam]

        output = subprocess.check_output(command,stderr=subprocess.STDOUT)

        if os.path.exists(input_bam+'.bai'):
            print("sambamba index sucessfull........")
            return input_bam+'.bai'   #return output file
        else:
            print("sambamba index failure...........")
            return False
