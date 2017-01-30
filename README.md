
##Structural Variation Engine<br>
(c) Timothy Becker, January 30 2016<br><br>
A script based execution engine for SV calling that abstracts seperate SV calling pipelines into a stage.
Each stage has a set of configurations for runtime which is stored as a JSON format parameter map.
Each SV caller stage has access to a set of standard inputs as well as reference specific and SV caller
specific files and ecosystems that are needed for execution.  Additional metadata files directories are
automatically generated at runtime into user specified directories. Each SV calling stage produces
a VCF formated file for use as input to fusorSV. Additional features include process spawning, output checking
file conversion and database integration.  Adapted for use on single systems, clusters or
cloud systems via VM images such as docker.  Easily extensible for addition of new SV calling algorithms
and data sources.  Several comon pre and post processing stages are included.<br>
###Requirements
python 2.7.10+, numpy, scipy, subprocess32, paramiko, scp, HTSeq, mysql.connector<br>
automated bash configuration of requirements is included for docker or container use

###Alternative docker image of full SVE
a full docker image can be obtained by:<br>
```bash
docker pull timothyjamesbecker/sve
```
###Current SV Callers
1.  BreakDancer (with VCF converter)<br>
2.  BreakSeq2<br>
3.  cnMOPS (with VCF converter)<br>
4.  CNVnator<br>
5.  Delly2<br>
6.  Hydra-Multi (with VCF converter)<br>
7.  GATK Haplotype Caller 3.6(with SV size VCF filter)<br>
8.  GenomeSTRiP2.0 (both SVDiscovery and CNVdiscovery) (with DEL/DUP VCF converter)<br>
9.  Lumpy-SV<br>
10. Tigra-SV (and EXT pipeline)<br>

###Planned Future SV Callers
1. SVelter<br>
2. MindTheGap<br>
3. TakeABreak<br>

###Current Metacalling Methods:
1.  FusorSV (with optional crossmap liftover)

###Current Pre and Post Processing Tools
1.  art_illumina<br>
2.  samtools<br>
3.  picard_tools<br>
4.  bwa aln, bwa mem<br>
5.  fa_to_2bit<br>
6.  vcf_tools<br>
7.  sambamba<br>
8.  samblaster<br>
9.  phred base quality encoding in BAM files<br>
10. bam_stats tool (samtools flagstat, coverage by sequence, BAM validation, phred sensing, read-group checking,ect)<br>
11. bam_clean (broken/problematic BAM file cleaning conditional routines)<br>

##Usage

####Alignment of FASTQ and generation of BAM files

####Structural Variation Calling on Bam files and generation of VCF files

####Merging Structural Variation Call Sets (VCF files)

##Core Frameworks and Extending


