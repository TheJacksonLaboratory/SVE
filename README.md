
##Structural Variation Engine<br>
(c) 2016 Timothy Becker<br><br>
A script based execution engine for SV calling that abstracts seperate SV calling pipelines into a stage.
Each stage has a set of configurations for runtime which is stored as a JSON format parameter map.
Each SV caller stage has access to a set of standard inputs as well as reference specific and SV caller
specific files and ecosystems that are needed for execution.  Additional metadata files directories are
automatically generated at runtime into user specified directories. Each SV calling stage produces
a VCF formated file for use as input to the FusorSV data fusion and arbitration method. Additional features include process spawning, output checking
file conversion and database integration.  Adapted for use on single systems, clusters or
cloud systems viadocker images.  Easily extensible for addition of new SV calling algorithms
and data sources.  Several comon pre and post processing stages are included.<br>
###Requirements (Non-docker)
python 2.7.10+, numpy, scipy, subprocess32, paramiko, scp, HTSeq, mysql.connector<br>
automated bash configuration of requirements is included for docker or container use

###Requirements (docker) full SVE
docker toolbox (or engine) version (:::) +
a full docker image can be obtained by:<br>
```bash
docker pull timothyjamesbecker/sve
```
this alternative requires a docker toolbox or docker engine to be installed on your system<br>

###Current SV Callers
1.  BreakDancer (with VCF converter)<br>https://github.com/genome/breakdancer2.<br>
2.  BreakSeq2<br>https://github.com/bioinform/breakseq2<br>
3.  cnMOPS (with VCF converter)<br>http://bioconductor.org/packages/release/bioc/html/cn.mops.html<br>
4.  CNVnator<br>https://github.com/abyzovlab/CNVnator<br>
5.  Delly2<br>https://github.com/dellytools/delly<br>
6.  Hydra-Multi (with VCF converter)<br>https://github.com/arq5x/Hydra<br>
7.  GATK Haplotype Caller 3.6(with SV size VCF filter)<br>https://software.broadinstitute.org/gatk/download<br>
8.  GenomeSTRiP2.0 (both SVDiscovery and CNVdiscovery) (with DEL/DUP VCF converter)<br> http://software.broadinstitute.org/software/genomestrip/download-genome-strip<br>
9.  Lumpy-SV<br> https://github.com/arq5x/lumpy-sv<br>
10. Tigra-SV (and EXT pipeline)<br> https://bitbucket.org/xianfan/tigra<br> https://bitbucket.org/xianfan/tigra-ext<br>

###Planned Future SV Callers
1. SVelter<br> https://github.com/mills-lab/svelter<br>
2. MindTheGap<br> https://gatb.inria.fr/software/mind-the-gap/<br>
3. TakeABreak<br> https://gatb.inria.fr/software/takeabreak/<br>

###Current Metacalling Methods:
1.  FusorSV (with optional crossmap liftover) <br> https://github.com/timothyjamesbecker/FusorSV<br>

###Current Pre and Post Processing Tools
1.  art_illumina<br> https://www.niehs.nih.gov/research/resources/software/biostatistics/art/<br>
2.  samtools (also bcftools, tabix, ect)<br> https://github.com/samtools/samtools<br>
3.  picard_tools<br> https://broadinstitute.github.io/picard/<br>
4.  bedtools2 <br> https://github.com/arq5x/bedtools2<br>
5.  bwa aln, bwa mem<br> https://github.com/lh3/bwa<br>
6.  fa_to_2bit<br>http://hgdownload.soe.ucsc.edu/admin/exe/<br>
7.  vcf_tools<br> https://github.com/vcftools/vcftools<br>
8.  sambamba<br>https://github.com/lomereiter/sambamba<br>
9.  samblaster<br>https://github.com/GregoryFaust/samblaster<br>
10. phred base quality encoding in BAM files<br>
11. bam_stats tool (samtools flagstat, coverage by sequence, BAM validation, phred sensing, read-group checking,ect)<br>
12. bam_clean (broken/problematic BAM file cleaning conditional routines)<br>

##Core Frameworks and Extending

##Usage
docker usage requires docker toolbox or other docker engine be installed, otherwise all executables should be build and eplaced inside a ...some_path/software/ folder where the SVE scritps should reside at: ...some_path/software/SVE.  Alternativly sym links can be used to redirect the script paths.

####(1) Alignment of FASTQ and generation of BAM files
The first step is to align FASTQ paired end reads to a reference genome.  The 1000 Genomes phase 3 reference fasta is currently sugested and tested against: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz  (Hg38 and mm10 are planned)
```bash
```

####(2) Structural Variation Calling on Bam files and generation of VCF files
```bash
```
####(3) Merging Structural Variation Call Sets (VCF files)
```bash
```


