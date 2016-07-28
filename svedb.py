#SVE DB Atomic Wrappers v 1.3, 10/25/2014-11/10/2014
#Timothy Becker, UCONN/SOE/CSE Phd Student
#includes common utilities such as pickling,gziping,byte calculations, current system time etc...

import StringIO
import datetime
import gzip
import cPickle as pickle
import mysql_connector as mysql

class SVEDB:
    #constructor
    def __init__(self,srv,db,uid,pwd):
        self.srv = srv           #MYSQL server hostname
        self.db = db             #db name  this is sve
        self.uid = uid           #user id  this is tied to the sve user
        self.pwd = pwd           #password this is tied to the sve user
        self.schema = {}         #master in-python schema, table-names and all field info
        
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        return 0
        
    #wrappers for all highest-level functions::::::::::::::::::::::::::::::::::::::::::::::::::::
    #query the max ref_id and ++ to generate a new one
    #pickle and coompress the seq object which is a string
    def new_ref(self,name,ref_len,seq_names,seq_lens,seqs='',url=''):
        #check the name for its ref_id first...
        try:
            response = self.select_fields_row(table='refs',primary_keys={'name':name},
                                              fields=['ref_id','ref_len'])
            if len(response)>0: ref_id = response[0]['ref_id']
        except IndexError:
            registered = self.get_max_key('refs')
            if registered is None: ref_id = 1
            else: ref_id = registered+1
        print(ref_id)
        #use the newly assigned ref_id or the existing one to try an insert op    
        seq_blob = self.obj_to_blob(seqs,True)
        pk = {'ref_id':ref_id}
        v = {'ref_id':ref_id,'name':name,'ref_len':ref_len,'seq_names':seq_names,
             'seq_lens':seq_lens,'seqs':seq_blob,'url':''}
        try:
            response = self.select_row('refs',pk)
            if type(response) is dict and len(response)>0:
                print('ref_id = %s already registered as a PK'%ref_id)
        except IndexError:
            return self.insert('refs',v) #True if things work, False otherwise...
        
    #query the max run_id and ++ to genearte a new one
    def new_run(self,platform_id,node_id,ref_id,calibrating=False,
                mut_mag=0,mut_len=0,mut_type='',mut_true_vc='',mut_ens_vc='',
                stage_id_list='',stage_depth=0,curr_stage_id=-1):
        run_id = self.get_max_key('runs')+1 #0 will be an empty table
        v = {'run_id':run_id,'platform_id':platform_id,'node_id':node_id,'ref_id':ref_id,
             'calibrating':calibrating,'mut_mag':mut_mag,'mut_len':mut_len,'mut_type':mut_type,
             'mut_true_vc':mut_true_vc,'mut_ens_vc':mut_ens_vc,'stage_id_list':stage_id_list,
             'stage_depth':stage_depth,'curr_stage_id':curr_stage_id,
             'start':self.time()} #stop left off until curr_stage == -1 => finsihed
        return (self.insert('runs',v),run_id) #True if things work, False otherwise...
        
    #query the max stage_id and ++ to generate a new one
    #typ is one of: 'dna synth', 'simulator', 'aligner', 'variant caller' and 'utility'
    #name is the registered name used in the UNIX/Terminal environment bound to the workers
    #in_ext, out_ext is the naming extension that defines the input output behavior of the stage
    #param_ranges is a map of the form: {'-p1':{'min':1.0,'max':2.0,'step':0.1}, ...}
    #this can be refrenced when exploring which params were used to get high scores etc...
    #to use you pass in the dict and this function will pickle/compress to a blob
    def new_stage(self,stage_id,typ,name,version,wrapper,in_ext,out_ext,param_map):
        #stage_id = self.get_max_key('stages')+1
        param_blob = self.obj_to_blob(param_map,True)
        v = {'stage_id':stage_id,'type':typ,'name':name,
             'version':version,'wrapper':wrapper,
             'in_ext':in_ext,'out_ext':out_ext,'param_map':param_blob}
        return self.insert('stages',v) #True if things work, False otherwise...
    
    #need to have a valid run_id & stage_id in mind to start
    #results/errors are optional and used to store any important data < 4GB inside the DB...
    #the main plan is to store a .vcf.gz blob for each variant caller to construct the mut_ens_vc array
    def new_staged_run(self,run_id,stage_id,in_files,in_files_size,params):
        param_blob = self.obj_to_blob(params,True)
        v = {'run_id':run_id,'stage_id':stage_id,'in_files':in_files,
             'in_files_size':in_files_size,'params':param_blob,
             'start':self.time(),'done':False}
        return self.insert('staged_runs',v) #True if things work, False otherwise...
        
    #have to have a run_id and stage_id in mind to update
    #common operations are updating the curr_stage pointer for
    #easy run monitoring of the full pipeline and the stop which
    #verifies that the run was finished (gives a full time estimate as well...)
    #:::TO DO:::
    def update_run(self,run_id,stage_id):
        pk = {'run_id':run_id,'stage_id':stage_id}
        v = {}
        return self.update('runs',pk,v)
    
    def update_stage(self,stage_id,typ,name,version,wrapper,in_ext,out_ext,param_map):
        pk = {'stage_id':stage_id}
        param_blob = self.obj_to_blob(param_map,True)
        v = {'stage_id':stage_id,'type':typ,'name':name,
             'version':version,'wrapper':wrapper,
             'in_ext':in_ext,'out_ext':out_ext,'param_map':param_blob}
        return self.update('stages',pk,v)
        
    #have to have a pre existing (run_id,stage_id) key to update
    #common operations involve updating the stop time to not be null
    #signifying that the stage is completed work and the pipeline can continue
    #its work by use of the stage_id_list(which the both distpatch and worker nodes need to have...)
    def update_staged_run(self,run_id,stage_id,in_files,results,errors,done):
        pk = {'run_id':run_id,'stage_id':stage_id,'in_files':in_files}
        results_blob = self.obj_to_blob(results,True)
        v = {'results':results_blob,'errors':errors,'stop':self.time(),'done':done} 
        return self.update('staged_runs',pk,v) #True if things work, False otherwise...
    
    #used to retrieve the last key for a given table
    #assuming that the key is monotonically increasing 2^31-1 max id key given int(11) MYSQL type used
    def get_max_key(self,table):
        if self.schema.has_key(table):
            pk = self.schema[table]['pk'][0] #assume a single key here
            sql = """SELECT max(%s) AS %s FROM %s.%s;"""%(pk,pk,self.db,table)
            v,res = (),{}
            with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
                res = ms.query(sql,v,True)
                return max(0,res[0][pk]) #assume one max of course
            return None
            
    def get_stage_ids_names(self):
        sql = """SELECT stage_id,name FROM %s.stages;"""%(self.db)
        v,res = (),{}
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(sql,v,True)
            return res
        return None
        
    def get_stage_id(self,wrapper):
        ss,data,table = [],[],'stages'
        select = 'SELECT stage_id FROM %s.%s'%(self.db,table)
        where = """\nWHERE wrapper = "%s";"""%wrapper #build this as k1 = v1 and k2 = v
        select += where
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(select,ss,True)
            for i in range(0,len(res)): #this is the number of rows
                row = res[i] #this is a dict return type
                for f in row: #this is the number of fields in each row f is now a field key for the dict
                    if type(row[f]) is bytearray: #is this only (long)blobs ?
                        row[f] = self.blob_to_obj(row[f])
                    else: row[f] = row[f]
                data.append(row) #insert uncompressed objects into data
        if len(data)==1: return data[0]['stage_id']
        else: raise IndexError

    def get_run_info(self,run_id):
        data = []
        run_info_sql = """
        SELECT ref_id,mut_mag,mut_len,mut_type
        FROM %s.runs WHERE run_id = %s"""%(self.db,run_id)
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(run_info_sql,(),True)
            for r in res: data.append(r)
        if len(data)==1: return data[0]
        else: raise IndexError
        
    def get_run_true_vc(self,run_id):
        data = []
        run_info_sql = """
        SELECT mut_true_vc
        FROM %s.runs WHERE run_id = %s"""%(self.db,run_id)
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(run_info_sql,(),True)
            for r in res: data.append(r)
        obj = None
        try:
            obj = self.blob_to_obj(data[0]['mut_true_vc'],status=True)
        except Exception:
            obj = data
        return obj

    #get the variant calls from one stagged run
    def get_staged_run_vc(self,run_id,stage_id):
        data = [] #don't worry about errors for now
        staged_run_info_sql = """
        SELECT params,results
        FROM %s.staged_runs WHERE run_id = %s AND stage_id = %s;
        """%(self.db,run_id,stage_id)
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(staged_run_info_sql,(),True)
            for r in res: data.append(r)
        if len(data)==1: return [self.blob_to_obj(data[0],status=True),
                                 self.blob_to_obj(data[1],status=True)]
        else: raise IndexError
    
    #get all the variant calls from a run (done or not)   
    def get_run_vc(self,run_id):
        data = []
        sql = """
        SELECT SR.stage_id, SR.run_id, S.name, SR.results, SR.done
        FROM %s.staged_runs as SR
            JOIN %s.stages as S
                ON (SR.stage_id = S.stage_id)
        WHERE SR.run_id = %s and S.out_ext like '%.vcf%';"""%(self.db,self.db,run_id)
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
                res = ms.query(sql,(),True)
                for r in res: data.append(r)
        return data
        
    def get_ref_id(self,ref_name):
        data = []
        sql = """
        SELECT ref_id FROM %s.refs WHERE name like '%s';"""%(self.db,ref_name)
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(sql,(),True)
            for r in res: data.append(r)
        if len(data)>=1: return data[0]['ref_id']
        else: raise IndexError
        
    def get_ref_name(self,run_id):
        data = []
        ref_name_sql = """
        SELECT name FROM %s.refs WHERE ref_id IN
        ( SELECT ref_id FROM %s.runs WHERE run_id = %s );"""%(self.db,self.db,run_id)
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(ref_name_sql,(),True)
            for r in res: data.append(r)
        if len(data)==1: return data[0]['name']
        else: raise IndexError
    
    def get_ref_info(self,ref_name):
        data = []
        sql = """
        SELECT * FROM %s.refs WHERE name like '%s';"""%(self.db,ref_name)
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(sql,(),True)
            for r in res: data.append(r)
        if len(data)==1: return data[0]
        else: raise IndexError
        
    #wrappers for all functions::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    def get_debug(self):
        debug = {'srv':self.srv,'db':self.db,'uid':self.uid,'pwd':self.pwd,'schema':self.schema}
        return debug
        
    #queries and retrives the full schema saving to self.schema{}
    def embed_schema(self):
        ts = self.select_tables()
        for t in ts['names']:
            fs = self.select_fields(t)
            self.schema[t]=fs
            
    def new(self):
        drop = "DROP TABLE IF EXISTS %s;"
        create = """
        CREATE TABLE %s (%s)
        ROW_FORMAT=DYNAMIC;
        """
        #new connection is safe and closes automatically
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            #DB workflow is to get the max id for each table and increment if a unique key
            #is obtainible from the checks instead of a new table for each simulator and aligner
            
            #Generated and Real References
            #ref_len is the total bp such as human ~3E7bp
            #seq is the fa.gz compressed string
            ms.query(drop%'refs')
            fields = """
            ref_id int primary key not null,
            name varchar(255) not null,
            ref_len bigint,
            seq_names text,
            seq_lens text,
            seqs longblob,
            url text
            """
            ms.query(create%('refs',fields))
            print('Created sve.refs')
            
            #Runs-----------------------------
            #a run is a set of simulations on parameters from a given ref-seq
            #the ref-seq is generated first and then a run is generated which is
            #a unique set of parameters (that induces randomized edits)
            #variant caller stages can be re-run to optimize parameters for a given platform and alignment
            ms.query(drop%'runs')
            fields = """
            run_id int primary key not null,
            platform_id varchar(50) not null,
            node_id varchar(50) not null,
            ref_id int not null,
            calibrating bit(1) not null,
            mut_mag int,
            mut_len bigint,
            mut_type char(3),
            mut_true_vc longblob,
            mut_ens_vc longblob,
            stage_id_list blob,
            stage_depth int,
            curr_stage_id int,
            start datetime,
            stop datetime
            """
            ms.query(create%('runs',fields))
            print('Created sve.runs')
            
            #Stages
            #stage_id a unique identifier for a new stage (a run is a (||)sequence of stages)
            #name is the text for the bound name in the enc, IE 'bwa'
            #param_ranges is a JSON string that gives {param:{start:1,stop:10,step:1}} for fixed step or
            #{param:{v1,v2,...,v3}} for fixed value iterations
            ms.query(drop%'stages')
            fields = """
            stage_id int primary key not null,
            type varchar(25) not null,
            name varchar(50) not null,
            version varchar(50),
            exe_path varchar(255),
            wrapper varchar(50),
            in_ext varchar(255),
            out_ext varchar(255),
            param_map  blob
            """
            ms.query(create%('stages',fields))
            print('Created sve.stages')
            
            #Staged_Runs links the stage+parameters to the run_id
            #so the stage_list and progress can be observed on a large
            #scale with full persistance and the possibility of multiple
            #distpatching centers building out unique runs into the DB
            #results can be used to store the output of each stage 4GB max compressed size...
            ms.query(drop%'staged_runs')
            fields = """
            run_id int not null,
            stage_id int not null,
            in_files varchar(255) not null,
            in_files_size varchar(255),
            params blob,
            results longblob,
            errors text,
            start datetime,
            stop datetime,
            done bit(1),
            primary key(run_id,stage_id,in_files)
            """
            ms.query(create%('staged_runs',fields))
            print('Created sve.staged_runs')
    
    #general table insertion wrapper that pickles and compresses
    #objects like dict,VariantCall, etc.. as blob/long blob fields   
    def insert(self,table,field_values):
        ss = [] #use tuple for mysql.connector compatibilty
        ins  = """INSERT INTO %s.%s ("""%(self.db,table)
        for f in field_values: ins += f+", "
        ins = ins[:-2]+") VALUES("
        for f in field_values:
            v = field_values[f]
            if type(v) is str  or type(v) is unicode or v is None:
                ins += "%s,"
                ss.append(v)
            else:
                ins += str(v)+","
        ins = ins[:-1]+");"
        ss = tuple(ss)
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            ms.query(ins,ss,False)
            return True #success?
        return False    #no success!

    #table name = table, primary_key value pairs to lookup = primary_key
    #field_values is the field names and values to update in the row
    def update(self,table,primary_keys,field_values):
        ss = []
        upd = """UPDATE %s.%s\nSET """%(self.db,table)
        for f in field_values:
            v = field_values[f]
            if type(v) is str or type(v) is unicode or v is None:
                upd += "%s = "%f
                upd += "%s,\n"
                ss.append(v)
            else:
                upd += """%s = %s,\n"""%(f,v)
        upd = upd[:-2] #seek back and eat the last ,\n
        where = """\nWHERE """ #build this as k1 = v1 and k2 = v2
        for k in primary_keys:
            v = primary_keys[k]
            if type(v) is str or type(v) is unicode or v is None:
                where += "%s = "%k
                where += "%s and "
                ss.append(v)
            else:
                where += """%s = %s and """%(k,primary_keys[k])
        where = where[:-4]+";"
        upd += where
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            ms.query(upd,ss,False)
            return True #success?
        return False    #no success!
        
    #returns a list of table bound to the given schema db = sve
    def select_tables(self):
        tbls = {'names':[]}
        tables = """
        SELECT TABLE_NAME
        FROM INFORMATION_SCHEMA.TABLES
        WHERE TABLE_SCHEMA = '%s';"""
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            sql,v = tables%(self.db),()
            res = ms.query(sql,v,True)
            for r in res: tbls['names']+=[str(r[u'TABLE_NAME'])]
        return tbls            
        
    #returns a ditc of field names as key
    def select_fields(self,table):
        flds,pk = {},[]
        fields = """
        SELECT COLUMN_NAME, COLUMN_TYPE, ORDINAL_POSITION, COLUMN_KEY 
        FROM INFORMATION_SCHEMA.COLUMNS 
        WHERE TABLE_SCHEMA = '%s' AND TABLE_NAME = '%s';"""
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            sql,v = fields%(self.db,table),()
            res = ms.query(sql,v,True)
            for r in res:
                flds[str(r[u'COLUMN_NAME'])]={'pos':r[u'ORDINAL_POSITION'],
                                              'type':str(r[u'COLUMN_TYPE'])}
                if r[u'COLUMN_KEY'] == u'PRI': pk+=[str(r[u'COLUMN_NAME'])]
            flds['pk']=pk
        return flds
            
    #general select all wrapper that finds and decompresses blobs
    #returns data as a list of dictionaries with embedded Python objects
    def select_all(self,table):
        data = [] #final uncompressed dictionary
        select,v = 'SELECT * FROM %s.%s;'%(self.db,table),()
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(select,v,True)
            for i in range(0,len(res)): #this is the number of rows
                row = res[i] #this is a dict return type
                for f in row: #this is the number of fields in each row f is now a field key for the dict
                    if type(row[f]) is bytearray: #is this only (long)blobs ?
                        row[f] = self.blob_to_obj(row[f])
                    else: row[f] = row[f]
                    #maybe also do a check to convert u'' to '' via type(row[f]) is unicode
                data.append(row) #insert uncompressed objects into data
        return data

    #sleect a row that has a PKm returns as a single dict
    def select_row(self,table,primary_keys):
        ss,data = [],[]
        select = 'SELECT * FROM %s.%s'%(self.db,table)
        where = """\nWHERE """ #build this as k1 = v1 and k2 = v2
        for k in primary_keys:
            v = primary_keys[k]
            if type(v) is str or type(v) is unicode or v is None:
                where += "%s = "%k
                where += "%s and "
                ss.append(v)
            else:
                where += """%s = %s and """%(k,primary_keys[k])
        where = where[:-4]+";"
        select += where
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(select,ss,True)
            for i in range(0,len(res)): #this is the number of rows
                row = res[i] #this is a dict return type
                for f in row: #this is the number of fields in each row f is now a field key for the dict
                    if type(row[f]) is bytearray: #is this only (long)blobs ?
                        row[f] = self.blob_to_obj(row[f])
                    else: row[f] = row[f]
                data.append(row) #insert uncompressed objects into data
        if len(data)==1: return data[0]
        else: raise IndexError
    
    #sleect a row that has a PKm returns as a single dict
    def select_fields_row(self,table,primary_keys,fields):
        ss,data = [],[]
        select = 'SELECT '
        for f in fields: select += f+', '
        select = select[:-2] #seek back one ', '
        select += ' FROM %s.%s'%(self.db,table)
        where = """\nWHERE """ #build this as k1 = v1 and k2 = v2
        for k in primary_keys:
            v = primary_keys[k]
            if type(v) is str or type(v) is unicode or v is None:
                where += "%s = "%k
                where += "%s and "
                ss.append(v)
            else:
                where += """%s = %s and """%(k,primary_keys[k])
        where = where[:-4]+";"
        select += where
        with mysql.MYSQL(self.srv,self.db,self.uid,self.pwd) as ms:
            res = ms.query(select,ss,True)
            for i in range(0,len(res)): #this is the number of rows
                row = res[i] #this is a dict return type
                for f in row: #this is the number of fields in each row f is now a field key for the dict
                    if type(row[f]) is bytearray: #is this only (long)blobs ?
                        row[f] = self.blob_to_obj(row[f])
                    else: row[f] = row[f]
                data.append(row) #insert uncompressed objects into data
        if len(data)>=1: return data
        else: raise IndexError    
    
    #a general select with a where clause may need to write more specific queries here
    #def select_rows(self,table,where):
    
    def obj_to_blob(self,O,status=False):
        results = pickle.dumps(O,protocol=pickle.HIGHEST_PROTOCOL)             #pickle and object in-memory
        if status:
            print("pickle length is %s mega bytes"%(round(len(results)/pow(1024.0,2),2)))
        gzip_out = StringIO.StringIO()                                          #and then gzip compress it
        with gzip.GzipFile(fileobj=gzip_out, mode="wb") as f: f.write(results) #to save space and prepare
        results  = gzip_out.getvalue()                                         #for storage as a blob
        if status:        
            print("compressed pickle length is %s mega bytes"%(round(len(results)/pow(1024.0,2),2)))  
            #with open(vca_path, 'wb') as pkl:
            #    pickle.dump(results,pkl,pickle.HIGHEST_PROTOCOL)              #save as a file
        return results
        
    def blob_to_obj(self,blob,status=False):
        s,buff = '',StringIO.StringIO(blob)                 #buffer for longblobs <= 4GB
        with gzip.GzipFile(fileobj=buff) as f: s = f.read() #read and attach full <= 4GB
        s = pickle.loads(s)                                 #unpickle the VariantCall objects
        return s
    
    #gets the current time as a MYSQL format datetime
    def time(self):
        return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    #returns the  number of bytes given number of M or G
    def toB(self,size):
        m,e = size*1.0,0
        while m/10>=1.0: m,e = m/10,e+1
        return int(m*pow(10,e%3)*pow(1024,e/3))