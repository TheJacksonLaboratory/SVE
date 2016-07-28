library('Rsamtools');
library('cn.mops');

start   <- proc.time()
ref     <- c('chr1','chr2','chr3','chr4','chr5')
f_names <- '~/Documents/CourseWork/12_2014_Fall/GRAD/SVCP/data/rg1_R1_CASE_DEL.bam'
s_names <- 'CASE'
paired  <- 'paired'
window  <- 1
cores   <- 1

header  <- scanBamHeader(f_names)
seqs <- header[[1]]$targets #names(seqs)=>c('chr1','chr2','chr3',...)

#cnmops version seems a bit low on read counts 23946 for 123521 mapped and 27119 mapped pairs ~ 50 bp long
#should give a per base average of around 6 -> read count sum should be ~ 5E5*6 = 3E6 ~= 2.7E3 ??? 
#data    <- getReadCountsFromBAM(BAMFiles=f_names, sampleNames=s_names, 
#                             	refSeqName = ref, WL=window, mode = paired, 
#                             	parallel=(cores-1));

stop    <- proc.time()
time    <- (stop-start)[3] #in seconds
time
counts <- mcols(data)$CASE

