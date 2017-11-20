<<<<<<< HEAD
Structural Variation Engine (SVE)
=================================

![Alt text](fusorSVlogo.jpg?raw=true "Logo")<br>

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
=======

![Alt text](fusorSVlogo.jpg?raw=true "Logo")<br>

##Structural Variation Engine<br>
(c) 2017 Timothy Becker & Wan-Ping Lee<br><br>
A python script based execution engine for SV calling that abstracts separate SV calling pipelines into a stage.
Each stage has a set of configurations for runtime which is stored as a JSON format parameter map.
Each SV caller stage has access to a set of standard inputs as well as reference specific and SV caller
specific files and ecosystems that are needed for execution.  Additional metadata files directories are
automatically generated at runtime into user specified directories. Each SV calling stage produces
a VCF formated file for use as input to the FusorSV data fusion and arbitration method. Additional features include process spawning, output checking, file conversion and database integration.  Adapted for use on single systems, clusters or
cloud systems via docker images.  Easily extensible for addition of new SV calling algorithms
and data sources.  Several common pre and post processing stages are included.<br>

###Requirements (docker) full SVE
docker toolbox (or engine) version 1.13.0+<br>
a full docker image can be obtained by:<br>
>>>>>>> master
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
bin/sve align [options] <-r FASTA> <FASTQ1 [FASTQ2]>
```

realign
```
bin/sve realign <-r FASTA> <BAM>
```

SV call
```
bin/sve call <-r FASTA> <-g hg19|hg38|others> <-a breakdancer|breakseq|cnvnator|hydra|delly|lumpy|cnmops> <BAM [BAM ...]>
```

After calling, each sample may have mulitple VCFs depending on how many callers used.
Please collect VCFs of a sample in a folder.

Example input vcf files can be organized as follows. Please note that vcfFiles is the argument for -i for FusorSV.
vcfFiles/sample1/sample1_S11.vcf  ## S11 is for delly, please check the table for tool and ID matching
vcfFiles/sample1/sample1_S10.vcf
vcfFiles/sample1/sample1_S4.vcf
vcfFiles/sample1/sample1_S0.vcf  ## optional, for truth
vcfFiles/sample2/sample2_S11.vcf
vcfFiles/sample2/sample2_S10.vcf
vcfFiles/sample2/sample2_S4.vcf
vcfFiles/sample2/sample2_S0.vcf  ## optional, for truth

Merge VCFs
```
python scripts/FusorSV/FusorSV.py -f scripts/FusorSV/data/models/human_g1k_v37_decoy.P3.INV2.11and17.pickle -L DEFAULT <-r FASTA> -i <vcfFiles> -p <THREADS> -o <OUT_DIR>
```
