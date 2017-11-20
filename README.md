Structural Variation Engine (SVE)
=================================

SVE is a python script based execution engine for Structural Variation (SV) detection and can be used for any levels of data inputs, raw FASTQs, aligned BAMs, or variant call format (VCFs), and generates a unified VCF as its output.
By design, SVE consists of alignment, realignment and the ensemble of eight state-of-the-art SV-calling algorithms by default. 
They are BreakDancer, BreakSeq, cnMOPS, CNVnator, DELLY, GenomeSTRiP, Hydra and LUMPY.
FusorSV is also embedded that is a data mining approach to assess performance and merge callsets from an ensemble of SV-calling algorithms.
![Alt text](overview.jpg?raw=true "SVE")

Requirements
------------
SVE requires the following to run.
	- python 2.7, numpy, scipy, subprocess32, scp, HTSeq
	- gcc 4.8 or greater
	- cmake(https://cmake.org/)
	- Root(https://root.cern.ch/)
	- R
Please set ROOT enviorment.
	- export ROOTSYS=/ROOT_Build_Path
	- export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib

FusorSV requires the following to run.
	- python 2.7, numpy, scipy, subprocess32, scp, HTSeq

Installation
------------
For SVE
```bash
git clone --recursive https://github.com/wanpinglee/SVE.git
cd SVE
make
```

Please check python2.7 header files and modify "CFLAGS_FUSOR_SV" in Makefile.
The header files may be on "/usr/include/python2.7" and use "CFLAGS_FUSOR_SV=-I /usr/include/python2.7" instead.
For FusorSVE
```bash
make FusorSV
```

Usage
------------
align
```
align [options] <-r FASTA> <FASTQ1 [FASTQ2]>
```

realign
```
realign <-r FASTA> <BAM>
```

SV call
```
call <-r FASTA> <-g hg19|hg38|others> <-a breakdancer|breakseq|cnvnator|hydra|delly|lumpy|cnmops> <BAM [BAM ...]>
```
