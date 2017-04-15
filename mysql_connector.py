#MYSQL v 0.1, 10/25/2014
#Timothy Becker, UCONN/SOE/CSE Phd Student
#MYSQL Connection Factory, wraps up sophisticated functionality in
#an easy to use extensible class...

import os
import re
import sys
import getpass
import mysql.connector as msc #pyodbc not easy to configure on mac, pypyodbc not encoding/decoding

class MYSQL:
    
    #constructor
    def __init__(self,srv,db,uid=False,pwd=False):
        #The MSSQL variables for injection safe connection strings
        self.srv = srv           #MYSQL server hostname to connect to
        self.db = db             #db name
        self.uid = uid
        self.pwd = pwd
        self.errors = ''
        self.conn = self.start() #keeps the connection object
        self.SQL =  []           #list of Qs
        self.V   =  []           #list of values for SQL
        
    def __enter__(self):
        return self

    #hype is that the DB error is generating stange files
    def __exit__(self, type, value, traceback):
        #saves sve_errors.txt to DATA folder that sits at the same level as this class.py file
        #directory = os.path.dirname(os.path.abspath(__file__))+'/'
        #if not os.path.exists(directory): os.makedirs(directory)
        #with open(directory+self.db+'_errors.txt', 'a') as f:
        #    f.write(self.errors)
        #close up the connection
        try: self.conn.close()
        except RuntimeError:
            print('__exit__():ER1.ODBC')
            self.errors += '__exit__():ER1.ODBC' + '\n'
        except Exception as err:
            print('__exit__():ER2.Unknown_Error: {}'.format(err))
            self.errors += '__exit__():ER2. Unknown Error: {}'.format(err)+'\n'
            pass
            
    def start(self):
        conn = None
        if (not self.uid) and (not self.pwd):
            print('uid: '),
            self.uid = sys.stdin.readline().replace('\n','')                   
            self.pwd = getpass.getpass(prompt='pwd: ',stream=None)#was stream=sys.sdin
        try:#connection start
            conn = msc.connect(host=self.srv,database=self.db,user=self.uid,password=self.pwd)
        except RuntimeError:
            print('start():ER3.ODBC')
            self.errors += 'start():ER3.ODBC' + '\n'
        except msc.errors.ProgrammingError:
            print('start():ER4.Connection')
            self.errors += 'start():ER4.Connection' + '\n'
        except Exception as err:
            print('start():ER5.Unknown_Error: {}'.format(err))
            self.errors += 'start():ER5.Unknown_Error: {}'.format(err)+'\n'
            pass
        return conn
        
    def query(self,sql,v=[],r=False):
        res = {}
        try:  #execute one sql and v list
            if r:
                cursor = self.conn.cursor(dictionary=True)
                cursor.execute(sql,v)
                #for row in cursor: res.append(row)
                res = cursor.fetchall()
                cursor.close()
            else: #this could be an insert command
                cursor = self.conn.cursor()
                cursor.execute(sql,v)
                cursor.close()
            self.conn.commit()
        except msc.errors.ProgrammingError as err:
            print('query():ER6.SQL_Malformed Error: {}'.format(err))
            #print('SQL GIVEN: '+sql+str(v))
            self.errors += 'query():ER6.SQL_Malformed\n' +sql+str(v)+ '\n'
        except msc.errors.DataError:
            print('query():ER7.Data_Not_Matching_Template')
            self.errors += 'query():ER7.Data_Not_Matching_Template' + '\n'
            print('SQL code:\n'+sql+'\nValue List:\n')
            self.errors += ('SQL code:\n'+sql+'\nValue List:\n')
            self.errors += str(v) + '\n'
        except msc.errors.IntegrityError:
            print('query():ER8.SQL_Constraint_Violation')
            self.errors += 'query():ER8.SQL_Constraint_Violation' + '\n'
            #print('SQL code:\n'+sql+'\nValue List:\n')
            self.errors += ('SQL code:\n'+sql+str(v)+'\nValue List:\n')
            self.errors += str(v) + '\n'
        except UnicodeDecodeError:
           print('query():ER9.Unicode_Decoding_Error')
           self.errors += 'query():ER9.unicode decoding issues\n'
        except Exception as err:
           print('query():ER10.Unknown_Error: {}'.format(err))
           self.errors += 'query():ER10.Unknown_Error: {}'.format(err)+'\n'
           pass
        return res
