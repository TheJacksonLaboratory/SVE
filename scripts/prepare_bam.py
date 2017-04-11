#!/usr/bin/env python
#BAM input metacalling pipeline
#given a reference directory (same name as the refrence) that already contains all of the
#preprossesing steps such as indexing, masking, chrom spliting, etc...

#simple_sve.py ref_path bam_path out_dir
import os
import sys
import glob
import socket
import time
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import stage_utils as su
import svedb
import stage
import subprocess32 as subprocess
from stages.utils.ParseParameters import ParseParameters
from stages.utils.ParseParameters import para_dict


def Dedup_Sort(in_bam, mem, threads):
    # Mark duplications
    st = stage.Stage('picard_mark_duplicates',dbc)
    params = st.get_params()
    params['-t']['value'] = threads
    params['-m']['value'] = mem
    st.set_params(params)
    dedup_bam = st.run(run_id,{'.bam':[sorted_bam]})
    if (dedup_bam == False):
        print "ERROR: picard_mark_duplicates fails"
    else:
        subprocess.call(["mv",dedup_bam,sorted_bam])
        # Sort bam
        st = stage.Stage('sambamba_index',dbc)
        st.set_params(params)
        st.run(run_id,{'.bam':[sorted_bam]})

if __name__ == '__main__':
    paras = para_dict
    ParseParameters(paras)

paras['machine'] = socket.gethostname
dbc = {'srv':'','db':'','uid':'','pwd':''}
run_id = 0
   
if paras['command'] == "realign":
    print "realign"
    """
    if args.mark_duplicates:
        d_start = time.time()
        st = stage.Stage('picard_mark_duplicates',dbc)
        outs = st.run(run_id,{'.bam':[bam]})
        d_stop = time.time()
        print('SVE:picard_mark_duplicates time was % hours'%round((d_stop-d_start)/(60**2),1))
    if args.replace_rg: #set to do one at a time only for now...
        r_start = time.time()
        base = bam.rsplit('/')[-1].rsplit('.')[0].rsplit('_')[0].rsplit('-')
        if SM is None: SM = base
        st = stage.Stage('picard_replace_rg',dbc)
        outs = st.run(run_id,{'.bam':[bam],
                              'platform_id':['illumina'],
                              'SM':[SM]})
        r_stop = time.time()
        print('SVE:picard_replace_rg time was %s sec'%round((r_stop-r_start)/(60**2),1))
    """
    if True:
        a_start = time.time()
        print paras['out_dir']
        st = stage.Stage('speedseq_realign',dbc)
        outs = st.run(run_id, {'.fa':paras['ref'],'.bam':paras['BAM'],'out_dir':paras['out_dir'],'threads':paras['threads'],'mem':paras['mem'],'RG':paras['RG']})
        a_stop = time.time()
        print('SVE:picard_mark_duplicates time was % hours'%round((a_stop-a_start)/(60**2),1))
else:
    print "align"
    a_start = time.time()
    """
    base = su.get_common_string_left(reads).rsplit('/')[-1].rsplit('.')[0]
    if SM is None: SM = base
    aligner_params = {'.fa':[ref_fa_path],'.fq':reads,'platform_id':['illumina'],'SM':[SM],'out_dir':[directory]}
    # Set the stage's parameters
    st = stage.Stage('fq_to_bam_piped',dbc)
    aligner_stage_params = st.get_params()
    aligner_stage_params['-t']['value'] = threads
    aligner_stage_params['-m']['value'] = mem
    if algorithm == 'mem':
        st = stage.Stage('fq_to_bam_piped',dbc)
        st.set_params(aligner_stage_params)
        # outs will receive ".sorted.bam"
        sorted_bam = st.run(run_id,aligner_params)
        Dedup_Sort(sorted_bam, mem, threads)
    elif algorithm == 'speed_seq':
        st = stage.Stage('speedseq_align',dbc)
        st.set_params(aligner_stage_params)
        outs = st.run(run_id,aligner_params)
    elif algorithm == 'aln':
        st = stage.Stage('bwa_aln',dbc)
        st.set_params(aligner_stage_params)
        # outs will receive ".sorted.bam"
        sorted_bam = st.run(run_id,aligner_params)
        Dedup_Sort(sorted_bam, mem, threads)
    """     
    a_stop = time.time()
    #print('SVE:BAM:%s was completed in %s hours'%(algorithm,round((a_stop-a_start)/(60.0**2),4)))
