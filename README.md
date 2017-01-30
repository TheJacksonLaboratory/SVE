
Structural Variation Engine
(c) Timothy Becker, November 17 2016

requirements: python 2.7.10+, numpy, scipy, subprocess32, paramiko, scp, HTSeq, mysql.connector
automated bash configuration of requirements is included for docker or container use

A script based execution engine for SV calling that abstracts seperate SV calling pipelines into a stage.
Each stage has a set of configurations for runtime which is stored as a JSON format parameter map.
Each SV caller stage has access to a set of standard inputs as well as reference specific and SV caller
specific files and ecosystems that are needed for execution.  Additional metadata files directories are
automatically generated at runtime into user specified directories. Each SV calling stage produces
a VCF formated file for use as input to fusorSV. Additional features include process spawning, output checking
file conversion and database integration.  Adapted for use on single systems, clusters or
cloud systems via VM images such as docker.  Easily extensible for addition of new SV calling algorithms
and data sources.  Several comon pre and post processing stages are included.

Current SV Callers:
1.  BreakDancer (with VCF converter)
2.  BreakSeq2
3.  cnMOPS (with VCF converter)
4.  CNVnator
5.  Delly
6.  Hydra-Multi (with VCF converter)
7.  GATK Haplotype Caller (with SV size VCF filter)
8.  GenomeSTRiP (both SVDiscovery and CNVdiscovery) (with DEL/DUP VCF converter)
9.  Lumpy-SV
10. Tigra-SV (and EXT pipeline)

Future SV Callers:
1. SVelter
2. MindTheGap
3. TakeABreak

Current Metacalling Methods:
1.  FusorSV (with optional crossmap liftover)

Pre and Post Processing:
1.  art_illumina
2.  samtools
3.  picard_tools
4.  bwa
5.  fa_to_2bit
6.  vcf_tools
7.  sambamba
8.  samblaster
9.  phred64to33 base quality conversion
