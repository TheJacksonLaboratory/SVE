#write a generic stage template to help with wrapper
#workflow is to pass arguments, parse with argparse and then
#develop a basic API that can be used to bind calls via polymorphism....
#could keep the db / schema free to allow use on multiple systems..

import os
import json
import sys
import socket
import HTSeq as ht
import subprocess32 as subprocess
sys.path.append('../') #go up one in the modules
import svedb

#function for auto-making svedb stage entries and returning the stage_id
class Stage_Wrapper(object):
    #path will be where a node should process the data using the in_ext, out_ext
    #stage_id should be pre-registered with db, set to None will require getting
    #a new stage_id from the  db by writing and registering it in the stages table
    def __init__(self,wrapper,dbc,retrieve,upload,params):
        self.srv = dbc['srv'] #MYSQL server hostname
        self.db  = dbc['db']  #db name  this is sve
        self.uid = dbc['uid'] #user id  this is tied to the sve user
        self.pwd = dbc['pwd'] #password this is tied to the sve user
        if retrieve:
            try:
                self.db_get_stage_info(wrapper) #use the db to gather info (one workflow)
            except IndexError:
                print('stage_information not availble via SVEDB, trying local configuration')
                #try the json instead
                self.json_get_stage_info(upload)
        else:
            self.json_get_stage_info(upload) #read meta data from JSON and upload          
            
        if params is None:
            #use the minimum values with all optional flags set to off...
            self.params = self.default_params()
        else:
            self.params = params
        self.software_path = os.path.dirname(os.path.abspath(__file__))+'/../../'
        
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        return 0
    
    #returns the default params in the param_map range
    def default_params(self):
        defaults,i = {},1
        for k in self.param_map:
            if self.param_map[k].has_key('rank'):
                defaults[k] = {'type':self.param_map[k]['type'],
                               'value':self.param_map[k]['min'],
                               'rank':self.param_map[k]['rank']}
            else:
                defaults[k] = {'type':self.param_map[k]['type'],
                               'value':self.param_map[k]['min'],
                               'rank':i}
                i+=1
        return defaults
        
    #will load new_stage_data of the form:
    #wrapper_name.json into wrapper_name.WrapperName:self...
    def json_get_stage_info(self,upload=False):
        info = {}
        path = os.path.dirname(os.path.abspath(__file__))+'/'
        json_name = self.__str__().split('.')[1]+'.json' #could be unsafe?
        #check to see if a JSON file is there? TO DO!!!
        with open(path+json_name,'r') as f:
            info = json.load(f)
            self.stage_id  = info['stage_id']
            self.typ       = info['type']
            self.name      = info['name']
            self.version   = info['version']
            self.wrapper   = info['wrapper']
            self.in_ext    = info['in_ext']    #update these to take CSV
            self.out_ext   = info['out_ext']   #update these to take CSV
            self.param_map = info['param_map']
            print('loaded param_map from: %s'%json_name)
            if upload: return self.db_set_stage_info(True)
        return False
        
    def get_host(self):
        return socket.gethostname()
        
    def split_in_exts(self):
        return self.in_ext.split(',')
        
    def split_out_exts(self):
        return self.out_ext.split(',')
    
    #returns '_S[stage_id].out_ext' for every out in out_ext
    def get_stage_exts(self):
        return ['_S'+str(self.stage_id)+'.'+ext for ext in self.split_out_exts()]

    def strip_in_ext(self,name,ext):
        i = name.rfind(ext)
        if i > 0: return name[0:i]
        else: return name
    
    def strip_name(self,full_path):
        i = full_path.rfind('/')
        if i > 0: return full_path[0:i]
        else: return full_path
        
    def strip_path(self,full_path):
        i = full_path.rfind('/')
        if i > 0: return full_path[i+1:]
        else: return full_path
        
    def strip_all_stages(self,name):
        i = name.find('_S')
        if i > 0: return name[0:i]
        else: return name
        
    def get_common_string_left(self,L):
        S = ''
        if len(L)>1:
            j,n = 0,min([len(i) for i in L])
            for i in range(n):
                if not all([L[0][i]==m[i] for m in L[1:]]): break
                else:                                       j+=1
            S = L[0][0:j]
        if len(L)==1: S = L
        return S
        
    def db_get_run_info(self,run_id):
        with svedb.SVEDB(self.srv, self.db, self.uid, self.pwd) as dbo:
            info = dbo.get_run_info(run_id)
            self.ref_id   = info['ref_id']
            self.mut_mag  = info['mut_mag']
            self.mut_len  = info['mut_len']
            self.mut_type = info['mut_type']
            return True
        return False
    
    def db_get_ref_name(self,run_id):
        with svedb.SVEDB(self.srv, self.db, self.uid, self.pwd) as dbo:
            try:
                self.ref_name = dbo.get_ref_name(run_id)
                return True
            except IndexError:
                self.ref_name = 'unkown_genome'
                print('reference name not avaible from SVEDB using unkown_genome')
                return False


    #get from the svedb given the unique stage_id allowing auto-retrieval
    #of all fields including the param_map to explore param optimizations
    def db_get_stage_info(self,wrapper):
        t,k,info = 'stages',{'stage_id':-1},{} #fail state stage_id = -1
        with svedb.SVEDB(self.srv, self.db, self.uid, self.pwd) as dbo:
            k['stage_id'] = dbo.get_stage_id(wrapper)
            info = dbo.select_row(t,k)
            self.stage_id  = info['stage_id']
            self.typ       = info['type']
            self.name      = info['name']
            self.version   = info['version']
            self.wrapper   = info['wrapper']
            self.in_ext    = info['in_ext']    #update these to take CSV
            self.out_ext   = info['out_ext']   #update these to take CSV
            self.param_map = info['param_map']
            return True
        return False

    #write a new entry into the svedb using the supplied param_map etc...
    def db_set_stage_info(self,new_stage=False):
        with svedb.SVEDB(self.srv, self.db, self.uid, self.pwd) as dbo:
            if new_stage:
                dbo.new_stage(self.stage_id,self.typ,self.name,
                              self.version,self.wrapper,
                              self.in_ext,self.out_ext,self.param_map)
            else:
                dbo.update_stage(self.stage_id,self.typ,self.name,
                                 self.version,self.wrapper,
                                 self.in_ext,self.out_ext,self.param_map)
            return True
        return False
    
    #this is the start of a staged_run which executes this stage
    def db_start(self,run_id,in_files):
        if type(in_files) is not list: in_files = [in_files] #make a list if not
        self.in_files_size = self.total_files_size(in_files)
        self.in_files = self.trim_in_file(in_files)
        with svedb.SVEDB(self.srv, self.db, self.uid, self.pwd) as dbo:
            dbo.new_staged_run(run_id,self.stage_id,self.in_files,self.in_files_size,self.params)
            return True
        return False
        
    #this is the stop of a stagged_run which should load results and errors and update the done flag    
    def db_stop(self,run_id,results,errors,done):
        with svedb.SVEDB(self.srv, self.db, self.uid, self.pwd) as dbo:
            try:
                debug = dbo.get_run_info(run_id)['debug']
            except Exception:
                debug = False
            if not debug: results = {} #can turn off results now
            dbo.update_staged_run(run_id,self.stage_id,self.in_files,results,errors,done)
            return True
        return False
        #search through the input names to see which have the correct in_ext...
    
    def vcf_to_vca(self,vcf_path):
        vca = []
        try:
            vcfr = ht.VCF_Reader(vcf_path)
            for vc in vcfr:
                vca += [vc]
        except Exception as E:
            pass
        return vca
    
    def total_files_size(self,in_files):
        if type(in_files) is not list: in_files = [in_files]
        byte_size = 0        
        for f in in_files:
            if os.path.exists(f):
                byte_size += os.path.getsize(f)
        return self.bytes_to_str(byte_size)
                
    def bytes_to_str(self,num_bytes):
        x,m = [num_bytes,0],0
        while x[0]/1024>0:
            x[1] = x[0]%1024 
            x[0],m = x[0]/1024,m+1
        size_map = {0:'',1:'K',2:'M',3:'G',4:'T',5:'P',6:'E'}
        return str(x[0])+'.'+str(x[1])+' '+size_map[m]+'B'
        
    #take up to 255 chars of the file name and trim the path info --> '/*/'
    def trim_in_file(self,in_file):
        trimmed = []
        for f in in_file:
            i = f.rfind('/')
            if i > 0 : trimmed += [f[i+1:]]
            else:      trimmed += [f]
        s = ', '.join(trimmed)
        l = min(255,len(s))
        return s[0:l]
        
    #look at the dict named inputs, each key is a in_ext...
    #inputs = {'.fa':['/home/tjb09009/node_data/rg1.fa'],
    #          '1.fq':['/home/tjb09009/node_data/node02/rg1_R11.fq']}
    #return a dict of the input extensions that are valid
    def search_inputs(self,inputs):
        valid = {}
        in_exts = self.split_in_exts()
        for k in inputs:
            if k in in_exts:
                valid[k]=inputs[k]
        return valid
    
    #search input bam names for grouping stems? or can pass a sample label mapping to inputs...    
    #def get_input_groups(self,inputs):
        
        
    #basic accessor for command string include paths, etc    
    def get_command_str(self):
        try:
            return ' '.join(self.command)
        except AttributeError:
            return ''
    
    def get_stage_id(self):
        try:
            return ' '.join(self.stage_id)
        except AttributeError:
            return ''
        
    #basic param accessor for params used for the command without paths
    #this is comparable between runs and invariant of paths used...
    def get_params(self):
        return self.params
    
    def set_params(self,params):
        self.params = params
        
    #override this function in each wrapper...
    def run(self,run_id,inputs={}):
        in_file = ''
        
        #retrieve all the parameters and map them to the series of calls to be executed...
        command = ['ls','-als']
        self.db_start(run_id,in_file,self.params)
        
        #workflow is to run through the stage correctly and then check for error handles
        #[3a]execute the command here----------------------------------------------------
        output,err = '',{}
        try:
            output = subprocess.check_output(command,stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as E:
            print('ouput: '+E.output)             #what you would see in the term
            err['output'] = E.output
            #the python exception issues (shouldn't have any...
            print('message: '+E.message)          #?? empty
            err['message'] = E.message
            #return codes used for failure....
            print('code: '+str(E.returncode))     #return 1 for a fail in art?
            err['code'] = E.returncode
        except OSError as E:
            print('ouput: '+E.strerror)             #what you would see in the term
            err['output'] = E.strerror
            #the python exception issues (shouldn't have any...
            print('message: '+E.message)          #?? empty
            err['message'] = E.message
            #the error num
            print('code: '+str(E.errno))
            err['code'] = E.errno
        print('output:\n'+output)
        #[3a]execute the command here----------------------------------------------------
        #[3b]do a os directory/data check or a ls type command to check the size
        #[3b]of the produced data to ensure that everything is fine...        
        if err == {}:
            self.db_stop(run_id,{'output':output},'',True)
            return output
        else:
            self.db_stop(run_id,{'output':output},err['message'],False)
            return None
