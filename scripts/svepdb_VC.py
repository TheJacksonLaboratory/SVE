#!/usr/bin/env python
#get the staged runs for a given run id and extract the VC
#svcpdb_test.py
import os
import sys
import HTSeq as ht #need this to reconstitute the pickled VC object
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import svedb

srv = 'arc-gis.ad.engr.uconn.edu'
db  = 'svcp'
uid = 'sv_calibrator'
pwd = 'sv_calibrator'

data = []
with svedb.SVEDB(srv, db, uid, pwd) as dbo:
    vca = dbo.get_run_true_vc(run_id=0) #this now works!
