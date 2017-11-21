# Structural Variation Engine (SVE)

![Alt text](fusorSVlogo.jpg?raw=true "Logo")<br>
(c) 2017 Timothy Becker & Wan-Ping Lee<br><br>

SVE is a python script based execution engine for Structural Variation (SV) detection and can be used for any levels of data inputs, raw FASTQs, aligned BAMs, or variant call format (VCFs), and generates a unified VCF as its output.
By design, SVE consists of alignment, realignment and the ensemble of state-of-the-art SV-calling algorithms by default. 
They are BreakDancer, BreakSeq, cnMOPS, CNVnator, DELLY, Hydra and LUMPY.
FusorSV is also embedded that is a data mining approach to assess performance and merge callsets from an ensemble of SV-calling algorithms.
![Alt text](overview.jpg?raw=true "SVE")

## Requirements
* python 2.7, HTSeq, numpy, scipy, subprocess32, bx-python, CrossMap and mygene
* gcc 4.8 or greater
* cmake(https://cmake.org/)
* Root(https://root.cern.ch/)
* R 3.2 or greater
	
Please set ROOT enviorment.
```
export ROOTSYS=/ROOT_Build_Path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
```

## Installation
### For SVE

```bash
git clone --recursive https://github.com/wanpinglee/SVE.git
cd SVE
make
```
### For FusorSV
Please check python2.7 header files and modify "CFLAGS_FUSOR_SV" in Makefile.
The header files may be on "/usr/include/python2.7" and use "CFLAGS_FUSOR_SV=-I /usr/include/python2.7" instead.
```
make FusorSV
```
Or, you can install FusorSV by setup.py
```
cd SVE/scripts/FusorSV/
python setup.py build_ext --inplace
tar -zxvf data.tar.gz
```

## Usage
### Align
Short reads in FASTQ will be mapped against the given FASTA and a sorted BAM will be generated.
```
bin/sve align [options] <-r FASTA> <FASTQ1 [FASTQ2]>
```
### Realign
If the reads are given by BAM format, realign will remap reads against FASTA and generate a sorted BAM.
We use SpeedSeq to accomplish realign.
```
bin/sve realign <-r FASTA> <BAM>
```
### Call SVs
There are seven SV calling algorithms that can be used for SV calling. VCF will be generated.
```
bin/sve call <-r FASTA> <-g hg19|hg38|others> <-a breakdancer|breakseq|cnvnator|hydra|delly|lumpy|cnmops> <BAM [BAM ...]>
```

### Merge VCFs
After calling, each sample may have mulitple VCFs depending on how many callers used.
Please collect VCFs of a sample in a folder.

Example input vcf files can be organized as follows. Please note that vcfFiles is the argument for -i for FusorSV.
* vcfFiles/sample1/sample1_S11.vcf
* vcfFiles/sample1/sample1_S10.vcf
* vcfFiles/sample1/sample1_S4.vcf
* vcfFiles/sample1/sample1_S0.vcf
* vcfFiles/sample2/sample2_S11.vcf
* vcfFiles/sample2/sample2_S10.vcf
* vcfFiles/sample2/sample2_S4.vcf

```
python scripts/FusorSV/FusorSV.py -f scripts/FusorSV/data/models/human_g1k_v37_decoy.P3.INV2.11and17.pickle -L DEFAULT -r <FASTA> -i <vcfFiles> -p <THREADS> -o <OUT_DIR>
```

## Docker
Dockerfile is provided for docker users.
```
cd SVE
docker build .
```

## License
The project is licensed under the GPL-3.0 License.
