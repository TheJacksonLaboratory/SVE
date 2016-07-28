Structural Variation Engine
(c) Timothy Becker, July 28 2016

current requirements: python 2.7.10+, subprocess32, paramiko, scp, HTSeq, mysql.connector
future requirements: python 2.7.10+, subprocess32

An engine for SV calling that abstracts seperate SV calling pipelines into a stage.
Each stage has a set of configurations for runtime which is stored as a JSON format parameter map.
Each SV caller stage has access to a set of standard inputs as well as reference specific and SV caller
specific and ecosystems that are needed for execution.  Additional metadata files directories are
automatically generated at runtime into user specified directories. Each SV calling stage produces
a VCF file format for use as input to fusorSV. Aditional features include process spawning, output checking
file conversion, database integration.  Adapted for use on single systems, clusters or
cloud systems via VM images such as docker.  Easily extensible for addition of new SV calling algorithms
and data sources.

Current SV Callers:
1.  BreakDancer
2.  BreakSeq2
3.  cnMOPS
4.  CNVnator
5.  Delly
6.  Hydra-Multi
7.  GATK Haplotype Caller
8.  GenomeSTRiP (both SVDiscovery and CNVdiscovery)
9.  Lumpy-SV
10. Tigra-SV

Current Metacalling Method:
1.  fusorSV
