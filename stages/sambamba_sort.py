import os
import sys
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import stage_wrapper

#function for auto-making svedb stage entries and returning the stage_id
class sambamba_sort(stage_wrapper.Stage_Wrapper):
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
    #bwa sampe ref.fa L.sai R.sai L.fq R.fq -f out.sam

    def run(input,output,threads):
        sambamba = self.tools['SAMBAMBA']
        command = [sambamba,sort,'-o',outout,'-l','5','-t',str(threads),inout]

        #[3a]execute the command here----------------------------------------------------
        print (' '.join(command))
        output = subprocess.check_output(command,stderr=subprocess.STDOUT)
        #[3b]check results--------------------------------------------------
        if os.path.exists(output):
            print("<<<<<<<<<<<<<sambamba sort sucessfull>>>>>>>>>>>>>>>\n")
            return output   #return output file
        else:
            print("<<<<<<<<<<<<<sambamba sort failure>>>>>>>>>>>>>>>\n")
            return False

