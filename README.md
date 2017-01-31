
##Structural Variation Engine<br>
(c) 2016 Timothy Becker<br><br>
A python script based execution engine for SV calling that abstracts separate SV calling pipelines into a stage.
Each stage has a set of configurations for runtime which is stored as a JSON format parameter map.
Each SV caller stage has access to a set of standard inputs as well as reference specific and SV caller
specific files and ecosystems that are needed for execution.  Additional metadata files directories are
automatically generated at runtime into user specified directories. Each SV calling stage produces
a VCF formated file for use as input to the FusorSV data fusion and arbitration method. Additional features include process spawning, output checking
file conversion and database integration.  Adapted for use on single systems, clusters or
cloud systems via docker images.  Easily extensible for addition of new SV calling algorithms
and data sources.  Several common pre and post processing stages are included.<br>

###Requirements (docker) full SVE
docker toolbox (or engine) version 1.13.0+<br>
a full docker image can be obtained by:<br>
```bash
docker pull timothyjamesbecker/sve
```

###Requirements (non-docker)
python 2.7.10+, numpy, scipy, subprocess32, paramiko, scp, HTSeq, mysql.connector<br>
All callers and pre post processing executables will have to be built and tested. See the links provided for each of these.
An automated bash configuration of the python requirements will setup the python distribution for you, but not individual callers or algorithms.<br>

###Current SV Callers
(see each of these links for liscensing information and citations)
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

##Core Frameworks and Extension
![Alt text](overview.jpg?raw=true "SVE")<br>
SVE uses a Object Orientated Design (OOD) pattern to decrease the perceived complexity and configuration of arbitrary pipelines.  The heart of the system invokes the registration of a stage which is extended by inheriting from the stage_wrapper class.  A stage is written by this extension and then its run method must be defined (overloaded).  Execution start by collecting and passing arguments into the run runing stage, but before this occurs the API can be configured with a database that can access frequent information such as reference identifies (using MD5 hashes of the sequence lengths) or other mechanisms.  Each stage has a unique runtimesruntimes identifier that allows the databases or SVEDB to keep track of runtimes for each file which can be used to gather information such as the average runtime for SV calling on a sample ect.  If no DB is configured, the SVE will continue in a autonomous mode (with fewer features).  Once a stage has been loaded at runtime, the name is used to load the configurations and executable paths which includes a stage param_map (JSON file) that passes valuable information such as the window size to use, ect.  Once the run method is invoked, the SVE spawns a new process on the host system and puts in a robust block on the parent SVE script.  Each line is run inside this child process so that the output is fully checked meaning that if the stages fails, the SVE does not and can rerun stages or continue to the next stage in a list ect...An example of this workflow can be seen inside the prepare_bam.py and variant_processor.py scripts which dynamically load stages as needed. Complex pipelines should be written in a way that each part of the pipeline is run and checked before continuing, so that meaningful error messages can be obtained from the resulting standard error stream<br>
##Usage
docker usage requires docker toolbox or other docker engine be installed, otherwise all executables should be build and eplaced inside a ...some_path/software/ folder where the SVE scritps should reside at: ...some_path/software/SVE.  Alternatively sym links can be used to redirect the script paths.  Additionally if you are using docker, you can update to the latest SVE scripts and pre-built executables (will just do a file diff and will be much faster than the first pull):
```bash
docker pull timothyjamesbecker/sve
```
####(1) Alignment of FASTQ and generation of BAM files
The first step is to align FASTQ paired end reads to a reference genome.  The 1000 Genomes phase 3 reference fasta is currently sugested and tested against: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz  (Hg38 and mm10 are planned)
```bash
docker run -v /data:/data timothyjamesbecker/sve /software/SVE/scripts/prepare_bam.py\
-a piped_mem\
-r /data/human_g1k_v37_decoy/human_g1k_v37_decoy.fa\
-f /data/sample1/sample1_1.fq,/data/sample1/sample1_2.fq\
-o /data/bams/\
-p 4
```
-a or --algorithm selects between the ```bash bwa aln, bwa mem and bwa mem | samtools view -Sb -``` workflows. The default is the highest performing that also produces the compressed BAM directly: piped_mem.<br>
-r or --ref is a fasta reference path.  The files should already have indecies produced.  You can perfrom these steps by using the optional script /software/SVE/scripts/prepare_ref.py described below/<br>
-f or --fqs is a comma seperated list of the FASTQ files you wish to align and map to the reference FASTA file listed with the -r argument.<br>
-o or --out_dir is the output directory that you wish to put intermediary files such as the unsorted (by coordinate) BAM files.  If this directory does not exist, a new directory will be generated in the files system and then files will be written to this location.<br>
-p or --cpus is the number of threads/cpus allocated to the alignment and mapping step which is passed to any section that is multithread or multiprocessor aware.<br>
-h or --help provides the latest build options.

####(2) Structural Variation Calling on Bam files and generation of VCF files
```bash
docker run -v /data:/data timothyjamesbecker/sve /software/SVE/scripts/variant_processor.py\
-r /data/human_g1k_v37_decoy/human_g1k_v37_decoy.fa\
-b /data/bams/sample1.bam\
-o /data/vcfs/\
-s breakdancer,breakseq,cnmops,cnvnator,delly,hydra,lumpy
```
-r or --ref is a fasta reference path.  The files should already have indecies produced.  You can perfrom these steps by using the optional script /software/SVE/scripts/prepare_ref.py described below/<br>
-b or --bams is the coordinate sorted and index BAM files that was generated in step (1) by the prepare_bam.py script<br>
-o or --out_dir is the output directory that you wish to put intermediary files such as FASTA assembly and partitioned BAM and BED files.  If this directory does not exist, a new directory will be generated in the files system and then files will be written to this location.<br>
-s a comma seperated listing of Structural Varaition calling pipelines you wish to invoke in series.  This workflow will continue through the list even if one caller fails until each caller pipeline has either finished returning a success or failed returning the failure status. If the -s argument list is left out, the currently supported names of caller will be displayed to assit the user.<br>
```bash
unknown processor name used as input argument: ham
availble stages are:
------------------------------
art_illumina
bam2cram
bam_split_all
bam_split_simple
bam_stats
breakdancer
breakseq
bwa_aln
bwa_index
bwa_mem
bwa_sampe
cnmops
cnvnator
cram2bam
cram2bam_split_all
delly
exomedepth
fa_to_2bit
fq_to_bam_piped
gatk_haplo
genome_strip
genome_strip_prepare_ref
hydra
lumpy
mrfast_divet
mrfast_index
picard_dict
picard_index
picard_mark_duplicates
picard_merge
picard_replace_rg
picard_sam_convert
pindel
samtools_fasta_index
samtools_index
samtools_merge
samtools_snp
samtools_sort
samtools_view
svseq
tigra
variationhunter
vcftools_filter
```
Optional Arguments:<br>
-D or --read_depth will privide average read depth guidance to RD callers that will assist with auto-setting internal paramters such as window size and junction length.<br><br>
-L or --read_length will provide guidnce on the read length used in the BAM file which will assist with auto-setting internal parameters such as window size and junction length.<br><br>
If -D or -L are left out, the bam_stats information gathering stage with determine this information for you and pass it on to the serially executed SV callers attached to the -s argument above<br><br>
-t or -targets option is used for target assembly for in silico SV callng validation or SV calling breakpoint refinement.  This step is recommened after all calling is completed and the step (3) FusorSV data fusion and arbitration step has produce SV calls for each sample which will be passed to this -t argument.  See the FusorSV manual for more information.<br><br>
-c or --chroms will attemp to run SV callers on a subset of the sequences present in the BAM file, effectively skipping alignments that fall on the undesired sequences.<br>
<br>

####(3) Merging Structural Variation Call Sets (Using a default fusion model with FusorSV)
```bash
docker run -v /data:/data timothyjamesbecker/sve /software/FusorSV/FusorSV.py\
-r /data/human_g1k_v37_decoy/human_g1k_v37_decoy.fa\
-c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y\
-i /data/vcfs/\
-o /data/fused/\
-f /software/FusorSV/data/models/human_g1k_v37_decoy.current.pickle\
-p 4\
-M 0.5\
-L
```
<br>
-r or --ref is a fasta reference path.  The files should already have indecies produced.<br><br>
-c or --chroms is an optional chrom list if you only want calls that were made on specific sequences.The default is 1-22,X,Y,MT with auto-sensing of (chr prefix)<br><br>
-i or --in_dir takes the VCF directory that was produced from step (2) the variant_processor.py script.  This will auto-append the stag or SV caller name to each samples VCF file in addition to creating a per sample sirectory with all callers inside which will look like this if you had samle1 and sample2 already processed:<br>
```bash
ls /data/vcfs/*/*
/data/vcfs/sample1/
  bd_config.txt
  sample1_S3.cov
  sample1_S3.header
  sample1_S3.header.rg
  sample1_S3.valid
  sample1_S4.vcf
  sample1_S9.vcf
  sample1_S10.vcf
  sample1_S11.vcf
  sample1_S13.vcf
  sample1_S14,vcf
  sample1_S17.vcf
  sample1_S18.vcf
  sample1_S35.vcf
  sample1_S36.vcf
  sample1_S38.vcf
/data/vcfs/sample2/
  bd_config.txt
  sample2_S3.cov
  sample2_S3.header
  sample2_S3.header.rg
  sample2_S3.valid
  sample2_S4.vcf
  sample2_S9.vcf
  sample2_S10.vcf
  sample2_S11.vcf
  sample2_S13.vcf
  sample2_S14,vcf
  sample2_S17.vcf
  sample2_S18.vcf
  sample2_S35.vcf
  sample2_S36.vcf
  sample2_S38.vcf
```
<br>
-o or --out_dir is the resulting directory where FusorSV will write to.  Output for FUsorSV includes a single VC file for each sample as well as a single merged VCF file across all samples attached to the -i arguments search.  Additionally the -L command will run the crossmap liftover tool using one of the internal UCSC chain files and partition the resulting VCF files into either [mapped]: meaning that there is a one to one mapping in the chain file on the coordinates or [unmapped]: meaning there was either a one to zero or a one to many mapping from the source into the destination sections of the chain file used.<br><br>
-f or --apply_fusion_model_path is used to apply a default model or to apply a new model you have created using the Training commands.  For more information on how to generate new models or update/append into models see the full FusorSV documentation.<br><br>
-p or --cpus sets the number of processor cores that will be used.  This will increase the amount of RAM but will speed up the processing by p as all major parts of the FusorSV processing are out-of-core and optimally || due to complete independance in the data stream by partitions.<br><br>
-M or cluster_overlap sets the amount of overlap permitted in the resulting VCF file for all samples attached to the -i argument.  This currently uses a (non-optimal) reciprocal overlap scanning proceedure to join together toaching calls, where the final resulting breakpoints with be averaged across all calls.  Individual VCF files will still remain, but this single master file will provide approximated genotype information with the support f every caller and the expectation under the fusion model giving you a clear an concise file to move forward with the biological relavence of the high-accuracy FusorSV calls.<br>
-L or --lift_over argument is used to pass a valid chain file to the FusorSV machinery to apply crossmap lift over<br>

