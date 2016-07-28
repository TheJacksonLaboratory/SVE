import numpy as np
import scipy.stats as st
import ctypes
import multiprocessing as mp
import pysam
import pysamstats

bam_file = '/Users/tbecker/Documents/CourseWork/12_2014_Fall/GRAD/SVE/data/rg1_R1_CASE_S1_S4_S8.bam'
bam = pysam.AlignmentFile(bam_file,'rb')
L = {s['SN']:s['LN'] for s in bam.header['SQ']}

cvg = mp.Array(ctypes.c_int32,2*L['chr1'],lock=False)

cov = pysamstats.load_coverage_strand(bam,chrom='chr1',
                                      truncate=True,pad=True,
                                      fields=['reads_pp_fwd','reads_pp_rev'])


