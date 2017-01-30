
##Structural Variation Engine<bd>
(c) Timothy Becker, January 30 2016<bd>
A script based execution engine for SV calling that abstracts seperate SV calling pipelines into a stage.
Each stage has a set of configurations for runtime which is stored as a JSON format parameter map.
Each SV caller stage has access to a set of standard inputs as well as reference specific and SV caller
specific files and ecosystems that are needed for execution.  Additional metadata files directories are
automatically generated at runtime into user specified directories. Each SV calling stage produces
a VCF formated file for use as input to fusorSV. Additional features include process spawning, output checking
file conversion and database integration.  Adapted for use on single systems, clusters or
cloud systems via VM images such as docker.  Easily extensible for addition of new SV calling algorithms
and data sources.  Several comon pre and post processing stages are included.<bd>
###Requirements
python 2.7.10+, numpy, scipy, subprocess32, paramiko, scp, HTSeq, mysql.connector<bd>
automated bash configuration of requirements is included for docker or container use

###Alternative docker image of full SVE
a full docker image can be obtained by:<br>
```bash
docker pull timothyjamesbecker/sve
```
###Current SV Callers
1.  BreakDancer (with VCF converter)<bd>
2.  BreakSeq2<bd>
3.  cnMOPS (with VCF converter)<bd>
4.  CNVnator<bd>
5.  Delly2<bd>
6.  Hydra-Multi (with VCF converter)<bd>
7.  GATK Haplotype Caller 3.6(with SV size VCF filter)<bd>
8.  GenomeSTRiP2.0 (both SVDiscovery and CNVdiscovery) (with DEL/DUP VCF converter)<bd>
9.  Lumpy-SV<bd>
10. Tigra-SV (and EXT pipeline)<bd>

###Planned Future SV Callers
1. SVelter<bd>
2. MindTheGap<bd>
3. TakeABreak<bd>

###Current Metacalling Methods:
1.  FusorSV (with optional crossmap liftover)

###Current Pre and Post Processing Tools
1.  art_illumina<bd>
2.  samtools<bd>
3.  picard_tools<bd>
4.  bwa aln, bwa mem<bd>
5.  fa_to_2bit<bd>
6.  vcf_tools<bd>
7.  sambamba<bd>
8.  samblaster<bd>
9.  phred BAM file sensing and
10. bam_stats tool (samtools flagstat, coverage by sequence, BAM validation, phred sensing, read-group checking,ect)
11. bam_clean (broken/problematic BAM file cleaning conditional routines)

##Usage
