#!/bin/sh

infile=$1
bam1=$2
json1=$3
outfile=$4
fastafile=$5
# read_group=$2
# base=`basename $infile`
svtyper -i ${infile} -B ${bam1} -l ${json1} -o ${outfile} -T ${fastafile}

#   --preemptible \
