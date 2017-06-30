# Timothy James Becker, UCONN CSE, 2015
# Fully automated, non-interactive cn.mops R wrapper V1.0
# Rscript cnmops.R ref_n=chr1 ctrl_bam= case_bam= ...

#cat('cnmops now running:',args['-v'][[1]],'\n'); #check out the list of list object named args
#is.null(args['k'][[1]]) can check the parser for existance of key = 'k'

# #START of preliminaries-----------------------------------------------------
options(warn=-1);
# in_bams <- c('~/Documents/CourseWork/12_2014_Fall/GRAD/SVCP/data/rg1_R1_CTRL_DEL.bam',
             # '~/Documents/CourseWork/12_2014_Fall/GRAD/SVCP/data/rg1_R1_CASE_DEL.bam')
# labels  <- c('CTRL','CASE')
# out_vcf <- '~/Documents/CourseWork/12_2014_Fall/GRAD/SVCP/data/rg1_R1_S16.vcf'
# cores   <- 1
# paired  <- 'paired'
# normal  <- 'poisson'
# cutoff  <- 0.25
# upper   <- 0.5
# lower   <- -0.9
# prior   <- 0
# min_seg <- 2
# min_cnt <- 0
# cir_seg <- 'DNAcopy'
# window  <- 100
# mode    <- 0

#Rscript cmd_parser.R commands...
cmd_args <- commandArgs();
scriptpath <- gsub('cnmops.R','',gsub('--file=','',cmd_args[4]));
# change library path
#print(scriptpath)
localLib <- paste(scriptpath,"../../src/R-package",sep="")
.libPaths(localLib)
#print(.libPaths())
source(paste(scriptpath,'cmd_parser.R',sep='')); #after this runs => args = a list of list object
in_bams  <- args['in_bams'][[1]]
in_chroms<- args['in_chroms'][[1]]
labels   <- args['labels'][[1]]
out_vcf  <- args['out_vcf'][[1]]
cores    <- args['cores'][[1]]
paired   <- args['paired'][[1]] #coded True/False
if(paired){ paired <- 'paired'; } else { paired <- 'unpaired'; }
normal   <- args['normal'][[1]] #coded 0,1,2,3,4
if (normal == 0){normal  <- 'mean'; }
if (normal == 1){normal  <- 'median'; }
if (normal == 2){normal  <- 'quant'; }
if (normal == 3){ normal <- 'poisson'; }
if (normal == 4){ normal <- 'mode'; }
cutoff   <- args['cutoff'][[1]]
upper    <- args['upper'][[1]]
lower    <- args['lower'][[1]]
prior    <- args['prior'][[1]]
min_seg  <- args['min_seg'][[1]]
min_cnt  <- args['min_cnt'][[1]]
cir_seg  <- args['cir_seg'][[1]] #coded True/False
if(cir_seg){ cir_seg <- 'DNAcopy'; } else { cir_seg <- 'fast'; }
window   <- args['window'][[1]]
mode     <- args['mode'][[1]]

#start up the dynamically loaded libraries and run cn.mops
source("http://bioconductor.org/biocLite.R")
#biocLite("Rsamtools");
#biocLite("cn.mops");
library('Rsamtools');
library('cn.mops');
f_names <- in_bams;
cat('bam input files are: ',f_names,'\n');

s_names <- as.vector(matrix('',nrow=length(in_bams),ncol=1))
for(i in 1:length(in_bams)){ 
	s_in_bam   <- strsplit(in_bams[i],'/')[[1]]
	s_in_bam   <- s_in_bam[length(s_in_bam)]
	s_names[i] <- strsplit(s_in_bam,'.bam')[[1]]
}
cat('sample labels are: ');
cat(paste(s_names));

#foreach chrom in the union of bam file set...
ref <- c();
for(i in 1:length(f_names)){
	header <- scanBamHeader(f_names[i]);
	ref <- union(ref,names(header[[1]]$targets));
}

res <- GRanges(); #empty place holder for error catching
#calculate copy variation number regions   
if(mode==0){ #mode = 0 -> cn.mops(regular multi sample WGS) 
    #try to make this more robust or alter input ref seq list = ref	
	tryCatch({
		cat('----------using cn.mops->multi sample WGS mode\n');
		data <- getReadCountsFromBAM(BAMFiles=f_names, sampleNames=s_names, 
                            	 	 refSeqName = ref, WL=window, 
                             	 	 parallel=(cores-1));              
		cat(paste('----------finished reading ',dim(mcols(data))[1],' bins\n',sep=''));
		#cleaning routine to suppress zero count issues.....cnmopsBUG
		for(contig in ref){ #remove from ref before runing anlaysis...
			if(sum(mcols(data[seqnames(data)%in%contig])[,1])==0){
				ref <- setdiff(ref,contig);
			}
		}
		cat('cleaned sequences are:')
		cat(paste(ref,sep=','))
		cat('\n')
		data <- data[seqnames(data)%in%ref]; #adjust data using cleaned ref seqs
		seqlevels(data) <- ref;
		seqnames(data)  <- droplevels(seqnames(data));
		res  <- suppressWarnings(cn.mops(data,
                		  	     I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                	      	     classes = paste('CN',seq(0,8,1),sep=''),
                	      	     priorImpact = prior, parallel = (cores-1),
                	      	     normType = normal, normQu = cutoff, norm = 1,
                	      	     upperThreshold = upper, lowerThreshold = lower,
                	      	     minWidth = min_seg, minReadCount = min_cnt)); #took out DNAcopy
    	cat(paste('----------found ',dim(mcols(cnvr(res)))[1],' variations\n',sep=''));
	}, error = function(err){
		cat(paste(err));
		cat('\n-------------------------------error----------------------------\n')
		#save and debug the image here...
		image_m <- strsplit(out_vcf,'.vcf')[[1]]
		image_n <- sub('.bam','',image_m[length(image_m)])
		outimage <- paste(image_n,'.RData',sep='');
		save.image(file=outimage); #saves objects for remote debugging/ect..
	}, finally = {});
    cat(warnings());            	
}
if(mode==1){ #mode = 1 -> referencecn.mops(case ctrl WGS)
	cat('----------using referencecn.mops->case/control WGS mode\n'); #future for cancer calling
	#TO DO split the s_names and f_names appropriatly into case_f_names,...
	ctrl_f_names <- as.vector(f_names[1]);
	ctrl_s_names <- as.vector(f_names[1]);
	case_f_names <- as.vector(f_names[2]);
	case_s_names <- as.vector(s_names[2]);
	case_data <- getReadCountsFromBAM(BAMFiles=case_f_names, sampleNames=case_s_names, 
                             	      refSeqName = ref, WL=window, 
                             	      parallel=(cores-1));
    cat(paste('----------finished reading ',dim(mcols(data))[1],' case bins\n',sep='')); 
    ctrl_data <- getReadCountsFromBAM(BAMFiles=ctrl_f_names, sampleNames=ctrl_s_names, 
                             	      refSeqName = ref, WL=window, 
                             	      parallel=(cores-1));
    cat(paste('----------finished reading ',dim(mcols(data))[1],' ctrl bins\n',sep=''));                          	              
	res  <- referencecn.mops(cases=case_data,controls=ctrl_data,
                	I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                	classes = paste('CN',seq(0,8,1),sep=''),
                	priorImpact = prior, parallel = (cores-1),
                	normType = normal, normQu = cutoff, norm = 1,
                	upperThreshold = upper, lowerThreshold = lower,
                	minWidth = min_seg, minReadCount = min_cnt,segAlgorithm = cir_seg);
    cat(paste('----------found ',dim(mcols(cnvr(res)))[1],' variations\n',sep=''));
             	
}
if(mode==2){ #mode = 2 -> exomecn.mops(whole exome sequencing WES)
    #try to make this more robust or alter input ref seq list = ref	
	tryCatch({
		cat('----------using exomecn.mops->WES mode\n');
		data <- getReadCountsFromBAM(BAMFiles=f_names, sampleNames=s_names, 
                            	 	 refSeqName = ref, WL=window, 
                             	 	 parallel=(cores-1));              
		cat(paste('----------finished reading ',dim(mcols(data))[1],' bins\n',sep=''));
		#cleaning routine to suppress zero count issues.....cnmopsBUG
		for(contig in ref){ #remove from ref before runing anlaysis...
			if(sum(mcols(data[seqnames(data)%in%contig])[,1])==0){
				ref <- setdiff(ref,contig);
			}
		}
		cat('cleaned sequences are:')
		cat(paste(ref,sep=','))
		cat('\n')
		data <- data[seqnames(data)%in%ref]; #adjust data using cleaned ref seqs
		seqlevels(data) <- ref;
		seqnames(data)  <- droplevels(seqnames(data));
		res  <- suppressWarnings(exomecn.mops(data,
                		  	     I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                	      	     classes = paste('CN',seq(0,8,1),sep=''),
                	      	     priorImpact = prior, parallel = (cores-1),
                	      	     normType = normal, normQu = cutoff, norm = 1,
                	      	     upperThreshold = upper, lowerThreshold = lower,
                	      	     minWidth = min_seg, minReadCount = min_cnt)); #took out DNAcopy
    	cat(paste('----------found ',dim(mcols(cnvr(res)))[1],' variations\n',sep=''));
	}, error = function(err){
		cat(paste(err));
		cat('\n-------------------------------error----------------------------\n')
		#save and debug the image here...
		image_m <- strsplit(out_vcf,'.vcf')[[1]]
		image_n <- sub('.bam','',image_m[length(image_m)])
		outimage <- paste(image_n,'.RData',sep='');
		save.image(file=outimage); #saves objects for remote debugging/ect..
	}, finally = {});
    cat(warnings());           	
}
if(mode==3){ #mode = 3 -> singlecn.mops(single sample WGS)
	#try to make this more robust or alter input ref seq list = ref	
        cat("BAM file: ", f_names, " ", s_names, " ", ref, " ", window, " ", paired)
        cat("\n")
	tryCatch({
		cat('----------using singlecn.mops->single sample WGS mode\n');
		data <- getReadCountsFromBAM(BAMFiles=f_names, sampleNames=s_names, 
                            	 	 refSeqName = ref, WL=window, 
                             	 	 parallel=(cores-1));              
		cat(paste('----------finished reading ',dim(mcols(data))[1],' bins\n',sep=''));
		#cleaning routine to suppress zero count issues.....cnmopsBUG
		for(contig in ref){ #remove from ref before runing anlaysis...
			if(sum(mcols(data[seqnames(data)%in%contig])[,1])==0){
				ref <- setdiff(ref,contig);
			}
		}
		cat('cleaned sequences are:')
		cat(paste(ref,sep=','))
		cat('\n')
		data <- data[seqnames(data)%in%ref]; #adjust data using cleaned ref seqs
		seqlevels(data) <- ref;
		seqnames(data)  <- droplevels(seqnames(data));
		res  <- suppressWarnings(singlecn.mops(data,
                		  	     I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4),
                	      	     classes = paste('CN',seq(0,8,1),sep=''),
                	      	     priorImpact = prior, parallel = (cores-1),
                	      	     normType = normal, normQu = cutoff, norm = 1,
                	      	     upperThreshold = upper, lowerThreshold = lower,
                	      	     minWidth = min_seg, minReadCount = min_cnt)); #took out DNAcopy
    	cat(paste('----------found ',dim(mcols(cnvr(res)))[1],' variations\n',sep=''));
	}, error = function(err){
		cat(paste(err));
		cat('\n-------------------------------error----------------------------\n')
		#save and debug the image here...
		image_m <- strsplit(out_vcf,'.vcf')[[1]]
		image_n <- sub('.bam','',image_m[length(image_m)])
		outimage <- paste(image_n,'.RData',sep='');
		save.image(file=outimage); #saves objects for remote debugging/ect..
	}, finally = {});
    cat(warnings());         	
}
                
#library(cn.mops);
#load('~/Desktop/TEMP/SVCP/data/NA19239_1_S9.RData');
#s_names <- c('NA19239_1')
#this needs to be redone to handle variable sample sizes and validated...

#TO DO ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#update for multiple sample analysis mode=0 and 1---------------------------------------------
#check that CNV were found...
#image_m <- strsplit(out_vcf,'.vcf')[[1]]
#image_n <- sub('.bam','',image_m[length(image_m)])
#outimage <- paste(image_n,'.RData',sep='');
#save.image(file=outimage); #saves objects for remote debugging/ect..
len <- 0;
if(class(res)[1]=="CNVDetectionResult"){
	len <- dim(mcols(cnvr(res)))[1][1]; #number of CNV reported by cn.mops
}
if(len > 0){
	#estimate the integer copy number over the CNVRs and store in the CNVDetectionResult->res               
	res_icn <- calcIntegerCopyNumbers(res);                
	cnvrs <- cnvr(res_icn); #GRanges object
	chr <- as.vector(seqnames(cnvrs));
	cnr <- as.data.frame(ranges(cnvrs)); #get the start, end, width of CNV regions
	ics <- as.data.frame(mcols(cnvrs));
	icn <- as.data.frame(matrix(0,nrow=len,ncol=length(s_names)));
	colnames(ics) <- s_names;
	colnames(icn) <- colnames(ics);
	#convert 'CN0'->0, 'CN1'->1, etc...s_names<-labels has the correct info here
	for(i in 1:len){
		for(j in s_names){
			icn[i,j] <- as.numeric(strsplit(toString(ics[i,j]),'CN')[[1]][2]);
		}
	}
	ctype<- matrix('',nrow=len,ncol=1);
	colnames(icn) <- s_names; #sample names
	d_r <- dim(cnr)[1]; #number of cnv regions
	s_n <- dim(icn)[2]; #number of samples read
	
	#count up duplications, deletions or uncertain calls...IE genotyping...
	if(mode==3){ #each mode have to do del, dup or cnv
		dupn <- icn[,1] > 2;
		deln <- icn[,1] < 2;
		eqn  <- icn[,1] == 2;
	}
	if(mode==1){ #each mode have to do del, dup or cnv
		dupn <- icn[,1] < icn[,2];
		deln <- icn[,1] > icn[,2];
		eqn  <- icn[,1] == icn[,2];
	}
	
	for(i in 1:len){
		if(dupn[i]){ ctype[i] <- 'DUP'; }
		if(deln[i]){ ctype[i] <- 'DEL'; }
		if(eqn[i]) { ctype[i] <- 'CNV'; }
	}
	#convert 'CN0'->0, 'CN1'->1, etc...s_names<-labels has the correct info here
	
	#build the default vectors to be inserted into the final table
	CHROM <- matrix(chr,nrow=len, ncol=1);
	POS   <- cnr[,'start']; #copy over start values
	ID    <- matrix('',nrow=len, ncol=1);
	ID    <- as.matrix(paste('cnmops_',seq(1,len,1),ID,sep=''));
	REF   <- matrix('N',nrow=len, ncol=1); #N for a DEL, full seq ACGT...AAACT for ref
	ALT   <- matrix(paste('<',ctype,'>',sep=''),nrow=len, ncol=1);
	QUAL  <- matrix('.',nrow=len, ncol=1);
	FILTER<- matrix('PASS',nrow=len, ncol=1);
	INFO  <- matrix('IMPRECISE',nrow=len, ncol=1);
	FORMAT<- matrix('CN',nrow=len, ncol=1);
	#SAMPLES already have in result matrix
	SVTYPE<- gsub('<','',ALT);
	SVTYPE<- gsub('>','',SVTYPE);
	#make the INFO field with the meta tags and END position
	for(i in 1:len){
		INFO[i] <-paste('END=',cnr[i,'end'],';SVTYPE=',SVTYPE[i],';SVLEN=',
	                    cnr[i,'width'],';IMPRECISE',sep=''); #no seperation on values
	}

	#bind the colums into a new table
	table <- data.frame(cbind(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO),
                        row.names=seq(1,len,1));
	d_t   <- dim(table)[2]; #default table size ~9
	#table[,(d_t+1):(d_t+s_n-1)] <- icn[2]; #insert integer copy numbers into the table
	#rename the colums for proper VCF compatibility
	c_n   <- c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'); 
	colnames(table) <- c_n; #attach the column names to the table
	#sort by chrom,pos the table if multi chrom...
	#table <- table[order(table[,'CHROM'],-table[,'POS']),];
	#redo the cmops_ids...
	#table[,'ID'] <- ID;
} else { table <- data.frame(); } #empty data frame -> no calls...

#VCF 4.0 Header construction-------------------------------------------------------
# '.' => don't know the value like <NA> or unknown
l01 <- '##fileformat=VCFv4.1\n';
l02 <- paste('##fileDate=',format(Sys.time(),'%Y%m%d'),'\n',sep='');
l03 <- '##source=cn.mops_CNV_calling\n';
l04 <- paste('##reference=',s_names[1],'\n',sep='');
l05 <- '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n';
l06 <- '##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n';
l07 <- '##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n';
l08 <- '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n';
l09 <- '##ALT=<ID=CNV,Description=\"Copy number variable region\">\n';
l10 <- '##ALT=<ID=DEL,Description=\"Deletion\">\n';
l11 <- '##ALT=<ID=DUP,Description=\"Duplication\">'; #put back \n if adding format
#l12 <- '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">';
header <- paste(l01,l02,l03,l04,l05,l06,l07,l08,l09,l10,l11,sep='');
#VCF 4.0 Header construction-------------------------------------------------------

#write the header and then append the finished table
#outvcf <- paste(outputdir,outputprefix,'.vcf',sep='');
write(header,file=out_vcf,append=F);
write.table(table,file=out_vcf,quote=F,append=T,col.names=T,row.names=F,sep='\t');

#save the plots...
#seg_plot <- paste(outputdir,outputprefix,'_segplot.pdf',sep='');
#pdf.options(width=8, height=10.5,onefile=T,paper='letter')
#pdf(file=seg_plot)
#segplot(res, seqname = chrom, sampleIdx = (1:n_sample))
#dev.off()

#read_plot <- paste(outputdir,outputprefix,'_readplot.pdf',sep='');	
#plot(res, which = (1:n_sample),toFile=T,filename=read_plot)

cat('\ncn.mops cnv calling complete........')
#END of cn.mops() cnv calling and visualization ===================================


