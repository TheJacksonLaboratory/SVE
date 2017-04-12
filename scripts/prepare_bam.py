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
    dedup_bam = st.run(run_id,{'.bam':in_bam,'mem':mem})
    if (dedup_bam == False):
        print "ERROR: picard_mark_duplicates fails"
    else:
        subprocess.call(["mv",dedup_bam,in_bam])
        # Sort bam
        st = stage.Stage('sambamba_index',dbc)
        st.run(run_id,{'.bam':[in_bam],'threads':threads})

if __name__ == '__main__':
    paras = para_dict
    ParseParameters(paras)

paras['machine'] = socket.gethostname
dbc = {'srv':'','db':'','uid':'','pwd':''}
run_id = 0

# Index FASTA if they are not there
if not all ([os.path.isfile(paras['ref'] + '.' + suffix) for suffix in ['amb','ann','bwt','pac','sa']]):
    st = stage.Stage('bwa_index',dbc)
    st.run(run_id, {'.fa':[paras['ref']]})

if paras['command'] == "realign":
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
    a_start = time.time()
    st = stage.Stage('speedseq_realign',dbc)
    outs = st.run(run_id, {'.fa':paras['ref'],'.bam':paras['BAM'],'out_dir':paras['out_dir'],'threads':paras['threads'],'mem':paras['mem'],'RG':paras['RG']})
    a_stop = time.time()
    print('SVE:picard_mark_duplicates time was % hours'%round((a_stop-a_start)/(60**2),1))
elif paras['command'] == "align":
    a_start = time.time()
    aligner_params = {'.fa':paras['ref'],'.fq':paras['FASTQ'],'out_dir':paras['out_dir'],'threads':paras['threads'],'mem':paras['mem'],'RG':paras['RG']}
    if paras['algorithm'] == 'bwa_mem':
        st = stage.Stage('fq_to_bam_piped',dbc)
        # outs will receive ".sorted.bam"
        sorted_bam = st.run(run_id,aligner_params)
        Dedup_Sort(sorted_bam, paras['mem'], paras['threads'])
    elif paras['algorithm'] == 'speed_seq':
        st = stage.Stage('speedseq_align',dbc)
        outs = st.run(run_id,aligner_params)
    elif paras['algorithm'] == 'bwa_aln':
        st = stage.Stage('bwa_aln',dbc)
        # outs will receive ".sorted.bam"
        sorted_bam = st.run(run_id,aligner_params)
        Dedup_Sort(sorted_bam, paras['mem'], paras['threads'])
    a_stop = time.time()
    print('SVE:BAM:%s was completed in %s hours'%(algorithm,round((a_stop-a_start)/(60.0**2),4)))
