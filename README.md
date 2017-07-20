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
*[python 2.7, numpy, scipy, subprocess32, scp, HTSeq]
*[gcc 4.8 or greater]
*[Root](https://root.cern.ch/)


FusorSV requires the following to run.
*[python 2.7, numpy, scipy, subprocess32, scp, HTSeq]

Installation
------------
For SVE
```bash
git clone --recursive https://github.com/wanpinglee/SVE.git
cd SVE
make
```

For FusorSVE
```bash
make FusorSV
```


