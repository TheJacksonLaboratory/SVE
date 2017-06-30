source("http://bioconductor.org/biocLite.R")
print(getwd())
.libPaths(getwd())
biocLite("cn.mops",lib="./")
biocLite("Rsamtools",lib="./")

