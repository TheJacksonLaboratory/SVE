#V0.1 Timothy James Becker © 2015
#Rscript Command Parser for basic Param passing
options(warn=-1); #suppress NA warnings among others

#args = c('-v=poisson','-h=1,2.45,1e-3','-k=true');    #for non-command-line testing
cmd_args = commandArgs(trailingOnly=T);       #white-space delimited...

L  <- vector('list',length(cmd_args));                    #parsed param list
ks <- as.vector(matrix('',ncol=1,nrow=length(cmd_args))); #param names
for (i in 1:length(cmd_args)){
	kv <- strsplit(cmd_args[i],'=')[[1]]  #use the = for params...
	if(length(kv)>1){ #matched the '=' char, continue
		ks[i]  <- kv[1];                 #grab the key for its name...
		vs <- strsplit(kv[2],',')[[1]];  #use the , to seperate lists...
		l <- vector(length=length(vs));  #for one param data row
		for(j in 1:length(vs)){          #for each value
			logical <- as.logical(vs[j]);
			numeric <- as.numeric(vs[j]);
			if(!is.na(logical)){         
				l[[j]] <- logical;
			} else {
				if(!is.na(numeric)){   
					l[[j]] <- numeric;
				} else {
					l[[j]] <- vs[j];    
				}
			}
		}
	} else {} #no matching '='
	L[[i]] <- l; #pack up the parsed values
	names(L) <- ks;
}
#L['-p'][[1]] #get out the data this way...
for(i in names(L)){
	#cat('class is: ',class(L[i][[1]]),'\n');
	#cat('length is: ',length(L[i][[1]]),'\n');
	cat(i,'=', L[i][[1]], ': class =', class(L[i][[1]]),
	    ': len =',length(L[i][[1]]),'\n', sep=' '); #print for diagnostics
}
cat('\n')
old <- c(ls(),'old'); #clean up the namespace
args <- L;
rm('list'= old);  #should only have args now
#cat('env=',ls(),'\n'); #check it out...
