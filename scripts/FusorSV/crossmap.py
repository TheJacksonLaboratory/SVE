#!/Applications/anaconda/bin/python
"""
---------------------------------------------------------------------------------------
CrossMap: lift over genomics coordinates between assemblies.
Support BED, GFF/GTF, BAM/SAM, BigWig/Wig, etc.
------------------------------------------------------------------------------------------
"""

import os,sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
        print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " CrossMap needs python2.7.*!\n"
        sys.exit()

import optparse
import collections
import sets
import subprocess
import string
from textwrap import wrap
from time import strftime

import pysam
from bx.bbi.bigwig_file import BigWigFile
from bx.intervals import *
import numpy as np
import datetime
from cmmodule  import ireader
from cmmodule  import BED
from cmmodule  import annoGene
from cmmodule  import bam_cigar
from cmmodule  import sam_header
from cmmodule  import myutils
from cmmodule  import fasta
from cmmodule  import bgrMerge

__author__ = "Liguo Wang, Hao Zhao"
__contributor__="Liguo Wang, Hao Zhao"
__copyright__ = "Copyright 2013, Mayo Clinic"
__credits__ = []
__license__ = "GPLv2"
__version__="0.2.3"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Production"

def printlog (mesg_lst):
	"""
	print progress into stderr
	"""
	if len(mesg_lst)==1:
		msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " +  mesg_lst[0]
	else:
		msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + ' '.join(mesg_lst)
	print >>sys.stderr,msg

def parse_header( line ):
        return dict( [ field.split( '=' ) for field in line.split()[1:] ] )
        
def wiggleReader( f ):
	"""
    Iterator yielding chrom, start, end, strand, value.
    Values are zero-based, half-open.
    Regions which lack a score are ignored.
	"""
	current_chrom = None
	current_pos = None
	current_step = None

    # always for wiggle data
	strand = '+'

	mode = "bed"
	for line in ireader.reader(f):
		if line.isspace() or line.startswith( "track" ) or line.startswith( "#" ) or line.startswith( "browser" ):
			continue
		elif line.startswith( "variableStep" ):
			header = parse_header( line )
			current_chrom = header['chrom']
			current_pos = None
			current_step = None
			if 'span' in header: current_span = int( header['span'] )
			else: current_span = 1
			mode = "variableStep"
		elif line.startswith( "fixedStep" ):
			header = parse_header( line )
			current_chrom = header['chrom']
			current_pos = int( header['start'] ) - 1
			current_step = int( header['step'] )
			if 'span' in header: current_span = int( header['span'] )
			else: current_span = 1
			mode = "fixedStep"
		elif mode == "bed":
			fields = line.split()
			if len( fields ) > 3:
				if len( fields ) > 5:
					yield fields[0], int( fields[1] ), int( fields[2] ), fields[5], float( fields[3] )
				else:
					yield fields[0], int( fields[1] ), int( fields[2] ), strand, float( fields[3] )
		elif mode == "variableStep": 
			fields = line.split()
			pos = int( fields[0] ) - 1
			yield current_chrom, pos, pos + current_span, strand, float( fields[1] )
		elif mode == "fixedStep":
			yield current_chrom, current_pos, current_pos + current_span, strand, float( line.split()[0] )
			current_pos += current_step
		else:
			raise "Unexpected input line: %s" % line.strip()

def bigwigReader(infile, chrom_sizes=None, bin_size = 2000):
	'''
	infile: bigwig format file
	chromsize: chrom_name: size.
	return: chrom, position (0-based), value
	'''
	bw_obj = BigWigFile(file=open(infile))
	for chr_name, chr_size in chrom_sizes.items():
		for chrom, st, end in BED.tillingBed(chrName = chr_name,chrSize = chr_size,stepSize = bin_size):
			sig_list = bw_obj.get_as_array(chrom,st,end)
			if sig_list is None:
				continue
			sig_list = np.nan_to_num(sig_list)
			if np.sum(sig_list)==0:
				continue	
			low_bound = st
			point_num = 1
			score = sig_list[0]
			for value in (sig_list[1:]):
				if value == score:
					point_num += 1
				else:
					yield(( chrom, low_bound,  low_bound + point_num, score))
					score = value
					low_bound = low_bound + point_num
					point_num = 1

			
def check_bed12(bedline):
	'''
	check if bed12 format is correct or not
	'''
	fields = bedline.strip().split()
	if len(fields) !=12:
		return False
	if fields[5] not in ['+','-','.']:
		return False
	try:
		chromStart = int(fields[1])
		chromEnd = int(fields[2])
		thickStart = int(fields[6])
		thickEnd = int(fields[7])
		blockCount = int(fields[9])
		blockSizes =  [int(i) for i in fields[10].rstrip(',').split(',')]
		blockStarts = [int(i) for i in fields[11].rstrip(',').split(',')]
	except:
		return False
	if chromStart > chromEnd or thickStart > thickEnd:
		return False
	if thickStart < chromStart or thickEnd > chromEnd:
		return False
	if len(blockSizes) != blockCount:
		return False
	if len(blockStarts) != blockCount:
		return False
	if blockCount <1:
		return False
	for i in blockSizes:
		if i < 0: return False
	for i in blockStarts:
		if i < 0: return False
	return True

def intersectBed((chr1, st1, end1), (chr2, st2, end2)):
	'''
	return intersection of two bed regions
	'''
	
	if int(st1) > int(end1) or int(st2) > int(end2):
		raise Exception ("Start cannot be larger than end")
	if chr1 != chr2:
		return None
	if int(st1) > int(end2) or int(end1) < int(st2):
		return None
	return (chr1, max(st1, st2), min(end1,end2))

def read_chain_file (chain_file, print_table=False):
	'''
	input chain_file could be either plain text, compressed file (".gz", ".Z", ".z", ".bz", ".bz2", ".bzip2"),
	or a URL pointing to the chain file ("http://", "https://", "ftp://"). If url was used, chain file must be plain text
	'''
	
	printlog(["Read chain_file: ", chain_file]),
	maps={}
	target_chromSize={}
	source_chromSize={}
	if print_table:
		blocks=[]
	
	for line in ireader.reader(chain_file):
		# Example: chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
		if not line.strip():
			continue
		line=line.strip()
		if line.startswith(('#',' ')):continue
		fields = line.split()
		
		if fields[0] == 'chain' and len(fields) in [12, 13]: 
  			score = int(fields[1])        # Alignment score
			source_name = fields[2]       # E.g. chrY
			source_size = int(fields[3])  # Full length of the chromosome
			source_strand = fields[4]     # Must be +
			if source_strand != '+':
				raise Exception("Source strand in an .over.chain file must be +. (%s)" % line)
			source_start = int(fields[5]) # Start of source region
			source_end = int(fields[6])   # End of source region
        
			target_name = fields[7]       # E.g. chr5
			target_size = int(fields[8])  # Full length of the chromosome
			target_strand = fields[9]     # + or -
			target_start = int(fields[10])
			target_end = int(fields[11])
			target_chromSize[target_name]= target_size
			source_chromSize[source_name] = source_size
			
			if target_strand not in ['+', '-']:
				raise Exception("Target strand must be - or +. (%s)" % line)
			#if target_strand == '+':
			#	target_start = int(fields[10])
			#	target_end = int(fields[11])
			#if target_strand == '-':
			#	target_start = target_size - target_end
			#	target_end = target_size - target_start
			chain_id = None if len(fields) == 12 else fields[12]
			if source_name not in maps:
				maps[source_name] = Intersecter()			
			
			sfrom, tfrom = source_start, target_start
       		
        # Now read the alignment chain from the file and store it as a list (source_from, source_to) -> (target_from, target_to) 		
		elif fields[0] != 'chain' and len(fields) == 3: 
			size, sgap, tgap = int(fields[0]), int(fields[1]), int(fields[2])
			if print_table:
				if target_strand == '+': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, tfrom, tfrom+size, target_strand))
				elif  target_strand == '-': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))
				
			if target_strand == '+':
				maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,tfrom, tfrom+size,target_strand)))
			elif  target_strand == '-':
				maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,target_size - (tfrom+size), target_size - tfrom, target_strand)))
				
			sfrom += size + sgap
			tfrom += size + tgap
		
		elif fields[0] != 'chain' and len(fields) == 1: 
			size = int(fields[0])
			if print_table:
				if target_strand == '+': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, tfrom, tfrom+size, target_strand))
				elif  target_strand == '-': blocks.append((source_name,sfrom, sfrom+size, source_strand, target_name, target_size - (tfrom+size), target_size - tfrom, target_strand))
					
			if target_strand == '+':		
				maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,tfrom, tfrom+size,target_strand)))
			elif target_strand == '-':
				maps[source_name].add_interval( Interval(sfrom, sfrom+size,(target_name,target_size - (tfrom+size), target_size - tfrom, target_strand)))
		else:
				raise Exception("Invalid chain format. (%s)" % line)
	if (sfrom + size) != source_end  or (tfrom + size) != target_end:
		raise Exception("Alignment blocks do not match specified block sizes. (%s)" % header)		
	
	if print_table:
		for i in blocks:
			print '\t'.join([str(n) for n in i])
	
	return (maps,target_chromSize, source_chromSize)	
	

def map_coordinates(mapping, q_chr, q_start, q_end, q_strand='+', print_match=False):
	'''
	Map coordinates from source assembly to target assembly
	'''
	
	matches=[]
	complement ={'+':'-','-':'+'}
	if q_chr not in mapping:
		return None
	else:
		targets = mapping[q_chr].find(q_start, q_end)
		if len(targets)==0:
			return None
		elif len(targets)==1:
			s_start = targets[0].start
			s_end = targets[0].end
			t_chrom = targets[0].value[0]
			t_start = targets[0].value[1]
			t_end = targets[0].value[2]
			t_strand = targets[0].value[3]
			
			(chr, real_start, real_end)  = intersectBed((q_chr,q_start,q_end),(q_chr,s_start,s_end))
			l_offset = abs(real_start - s_start)
			#r_offset = real_end - s_end
			size = abs(real_end - real_start)
			
			matches.append( (chr, real_start, real_end,q_strand))
			if t_strand == '+':
				i_start = t_start + l_offset
				if q_strand == '+':
					matches.append( (t_chrom, i_start, i_start + size, t_strand))
				else:
					matches.append( (t_chrom, i_start, i_start + size, complement[t_strand]))
			elif t_strand == '-':
				i_start = t_end - l_offset - size
				if q_strand == '+':
					matches.append( (t_chrom, i_start,  i_start + size, t_strand))
				else:
					matches.append( (t_chrom, i_start,  i_start + size, complement[t_strand]))
			else:
				raise Exception("Unknown strand: %s. Can only be '+' or '-'." % q_strand)
			
		elif len(targets) > 1:
			for t in targets:
				s_start = t.start
				s_end = t.end
				t_chrom = t.value[0]
				t_start = t.value[1]
				t_end = t.value[2]
				t_strand = t.value[3]
				
				(chr, real_start, real_end)  = intersectBed((q_chr,q_start,q_end),(q_chr,s_start,s_end))

				l_offset = abs(real_start - s_start)
				#r_offset = abs(real_end - s_end)
				size = abs(real_end - real_start)

				matches.append( (chr, real_start, real_end,q_strand) )
				if t_strand == '+':
					i_start = t_start + l_offset
					if q_strand == '+':
						matches.append( (t_chrom, i_start, i_start + size, t_strand))
					else:
						matches.append( (t_chrom, i_start, i_start + size, complement[t_strand]))
				elif t_strand == '-':
					i_start = t_end - l_offset - size
					if q_strand == '+':
						matches.append( (t_chrom, i_start,  i_start + size, t_strand))
					else:
						matches.append( (t_chrom, i_start,  i_start + size, complement[t_strand]))
				else:
					raise Exception("Unknown strand: %s. Can only be '+' or '-'." % q_strand)
	
	if print_match:
		print matches
		# input: 'chr1',246974830,247024835
		# output: [('chr1', 246974830, 246974833, '+' ), ('chr1', 248908207, 248908210, '+' ), ('chr1', 247024833, 247024835, '+'), ('chr1', 249058210, 249058212,'+')]
		# [('chr1', 246974830, 246974833), ('chr1', 248908207, 248908210)]
	
	return matches

def info_to_idx(info):
    svmethod = info.split('SVMETHOD=')[-1].split(';')[0]
    raw = svmethod.split(',')
    idx,prev = {},0    
    for i in range(len(raw)):
        if raw[i].find(':')>0:
            split = raw[i].split(':')
            prev = int(split[0])
            idx[prev] = set([int(split[1])])
        else:
            idx[prev].add(int(raw[i]))
    return idx

def vcf_row_to_coord(fields):
    start = int(fields[1])+1
    end = start + len(fields[3])
    try:
        end = int(fields[7].split('END=')[-1].split(';')[0])+1
    except Exception:
        pass
    return (start,end)

def update_vcf_info(fields,start,end,strand):
    info = fields[7].split(';')
    for i in range(len(info)):
        if info[i].startswith('END='):
            info[i] = 'END='+str(end)
        if info[i].startswith('SVLEN='):
            info[i] = 'SVLEN='+str(np.abs(end-start))
    return ';'.join(info)
	
def crossmap_vcf_file(mapping,infile,liftoverfile,refgenome):
	#index refegenome file if it hasn't been done
	if not os.path.exists(refgenome + '.fai'):
		printlog(["Creating index for", refgenome])
		pysam.faidx(refgenome)
	
	FILE_OUT = open(infile.rstrip('.vcf')+'_liftover.mapped.vcf' ,'w')
	UNMAP = open(infile.rstrip('.vcf')+'_liftover.unmapped.vcf','w')
	
	total = 0
	fail = 0
	withChr = False
	
	for line in ireader.reader(infile):
		if not line.strip():
			continue
		line=line.strip()
		if line.startswith('##'):
			print >>FILE_OUT,line
			print >>UNMAP, line
			continue
		fields = line.split()
		#fields = string.split(line,maxsplit=7)
		if fields[0].startswith('#'):
			print >>FILE_OUT, "##liftOverProgram=CrossMap(https://sourceforge.net/projects/crossmap/)"
			print >>FILE_OUT, "##liftOverFile=" + liftoverfile
			print >>FILE_OUT, "##new_reference_genome=" + refgenome
			print >>FILE_OUT, "##liftOverTime=" + datetime.date.today().strftime("%B%d,%Y")
			print >>FILE_OUT,line
			print >>UNMAP, line
		else:
			total += 1
			if fields[0].startswith('chr'):
				withChr = True
				chrom = fields[0]
			else:
				chrom = 'chr' + fields[0]
			start,end = vcf_row_to_coord(fields)
			a = map_coordinates(mapping, chrom, start, end,'+')
			if a is None:
				print >>UNMAP, line
				fail += 1
				continue
			if len(a) == 2 and a[0][0]==a[1][0]:
				# update chrom
				if withChr is False:
					fields[0] = str(a[1][0]).replace('chr','')
				# update start coordinate
				fields[1] = a[1][1] + 1
				fields[7] = update_vcf_info(fields,a[1][1]+1,a[1][2]+1,a[1][3])
				# update ref allele
				tmp = pysam.faidx(refgenome,str(a[1][0]) + ':' + str(a[1][1]+1) + '-' + str(a[1][2]))
				#fields[3] =tmp[-1].rstrip('\n\r').upper()
				
				if tmp[-1].rstrip('\n\r').upper() != fields[4]:
					print >>FILE_OUT, '\t'.join(map(str, fields))
				else:
					print >>UNMAP, line
					fail += 1
			else:
				print >>UNMAP, line
				fail += 1
				continue
	FILE_OUT.close()
	UNMAP.close()
	printlog (["Total entries:", str(total)])
	printlog (["Failed to map:", str(fail)])
				
			
def crossmap_bed_file(mapping, inbed,outfile=None):
	'''
	Convert genome coordinates (in bed format) between assemblies.
	BED format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
	'''
	
	# check if 'outfile' was set. If not set, print to screen, if set, print to file
	if outfile is not None:
		FILE_OUT = open(outfile,'w')
		UNMAP = open(outfile + '.unmap','w')
	else:
		pass
	
	for line in ireader.reader(inbed):
		if line.startswith('#'):continue
		if line.startswith('track'):continue
		if line.startswith('browser'):continue
		if not line.strip():continue
		line=line.strip()
		fields=line.split()
		strand = '+'
		
		# filter out line less than 3 columns
		if len(fields)<3:
			print >>sys.stderr, "Less than 3 fields. skip " + line
			if outfile:
				print >>UNMAP, line
			continue
		try:
			int(fields[1])
		except:
			print >>sys.stderr, "Start corrdinate is not an integer. skip " + line
			if outfile:
				print >>UNMAP, line
			continue
		try:
			int(fields[2])
		except:
			print >>sys.stderr, "End corrdinate is not an integer. skip " + line
			if outfile:
				print >>UNMAP, line
			continue
		if int(fields[1]) > int(fields[2]):
			print >>sys.stderr, "\"Start\" is larger than \"End\" corrdinate is not an integer. skip " + line
			if outfile:
				print >>UNMAP, line
			continue
			
		# deal with bed less than 12 columns
		if len(fields)<12:
			
			# try to reset strand
			try:
				for f in fields:
					if f in ['+','-']:
						strand = f
			except:
				pass
				
			a = map_coordinates(mapping, fields[0], int(fields[1]),int(fields[2]),strand)
			# example of a: [('chr1', 246974830, 246974833,'+'), ('chr1', 248908207, 248908210,'+')]

			try:
				if len(a) % 2 != 0:
					if outfile is None:
						print line + '\tFail'
					else:
						print >>UNMAP, line
					continue
				if len(a) == 2:
					
					#reset fields
					fields[0] = a[1][0]
					fields[1] = a[1][1]
					fields[2] = a[1][2]
					for i in range(0,len(fields)):	#update the strand information
						if fields[i] in ['+','-']:
							fields[i] = a[1][3]
					
					if outfile is None:
						print line + '\t->\t' + '\t'.join([str(i) for i in fields])
					else:
						print >>FILE_OUT, '\t'.join([str(i) for i in fields])
				if len(a) >2 :
					count=0
					for j in range(1,len(a),2):
						count += 1
						fields[0] = a[j][0]
						fields[1] = a[j][1]
						fields[2] = a[j][2]
						for i in range(0,len(fields)):	#update the strand information
							if fields[i] in ['+','-']:
								fields[i] = a[j][3]
						
						if outfile is None:
							print line + '\t'+ '(split.' + str(count) + ':' + ':'.join([str(i) for i in a[j-1]]) + ')\t' + '\t'.join([str(i) for i in fields])
						else:
							print >>FILE_OUT, '\t'.join([str(i) for i in fields])
			except:
				if outfile is None:
					print line + '\tFail'
				else:
					print >>UNMAP, line
				continue
		
		# deal with bed12 and bed12+8 (genePred format)
		if len(fields)==12 or len(fields)==20:
			strand = fields[5]
			if strand not in ['+','-']:
				raise Exception("Unknown strand: %s. Can only be '+' or '-'." % strand)
			fail_flag=False
			exons_old_pos = annoGene.getExonFromLine(line)	#[[chr,st,end],[chr,st,end],...]
			#print exons_old_pos
			exons_new_pos = []
			for e_chr, e_start, e_end in exons_old_pos:
				a = map_coordinates(mapping, e_chr, e_start, e_end, strand)
				if a is None:
					fail_flag =True
					break

				if len(a) == 2:
					exons_new_pos.append(a[1]) 
					# a has two elements, first is query, 2nd is target. # [('chr1', 246974830, 246974833,'+'), ('chr1', 248908207, 248908210,'+')]
				else:
					fail_flag =True
					break
			
			if not fail_flag:
				# check if all exons were mapped to the same chromosome the same strand
				chr_id = set()
				exon_strand = set()
			
				for e_chr, e_start, e_end, e_strand in exons_new_pos:
					chr_id.add(e_chr)
					exon_strand.add(e_strand)
				if len(chr_id) != 1 or len(exon_strand) != 1:
					fail_flag = True
				
				if not fail_flag:
					# build new bed 
					cds_start_offset = int(fields[6]) - int(fields[1])
					cds_end_offset = int(fields[2]) - int(fields[7])
					new_chrom = exons_new_pos[0][0]
					new_chrom_st = exons_new_pos[0][1]
					new_chrom_end = exons_new_pos[-1][2]
					new_name = fields[3]
					new_score = fields[4]
					new_strand = exons_new_pos[0][3]	
					new_thickStart = new_chrom_st + cds_start_offset
					new_thickEnd = new_chrom_end - cds_end_offset
					new_ittemRgb = fields[8]
					new_blockCount = len(exons_new_pos)
					new_blockSizes = ','.join([str(o - n) for m,n,o,p in exons_new_pos])
					new_blockStarts = ','.join([str(n - new_chrom_st) for m,n,o,p in exons_new_pos])
				
					new_bedline = '\t'.join(str(i) for i in (new_chrom,new_chrom_st,new_chrom_end,new_name,new_score,new_strand,new_thickStart,new_thickEnd,new_ittemRgb,new_blockCount,new_blockSizes,new_blockStarts))
					if check_bed12(new_bedline) is False:
						fail_flag = True
					else:
						if outfile is None:
							print line + '\t->\t' + new_bedline
						else:
							print >>FILE_OUT, new_bedline
			
			if fail_flag:
				if outfile is None:
					print line + '\tFail'
				else:
					print >>UNMAP, line
				
def crossmap_gff_file(mapping, ingff,outfile=None):			
	'''
	Convert genome coordinates (in GFF/GTF format) between assemblies.
	GFF (General Feature Format) lines have nine required fields that must be Tab-separated:
	
    1. seqname - The name of the sequence. Must be a chromosome or scaffold.
    2. source - The program that generated this feature.
    3. feature - The name of this type of feature. Some examples of standard feature types
       are "CDS", "start_codon", "stop_codon", and "exon".
    4. start - The starting position of the feature in the sequence. The first base is numbered 1.
    5. end - The ending position of the feature (inclusive).
    6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1
       for this annotation data set, the score value will determine the level of gray in
       which this feature is displayed (higher numbers = darker gray). If there is no score
       value, enter ".".
    7. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
    8. frame - If the feature is a coding exon, frame should be a number between 0-2 that
       represents the reading frame of the first base. If the feature is not a coding exon,
       the value should be '.'.
    9. group - All lines with the same group are linked together into a single item.
    
    GFF format: http://genome.ucsc.edu/FAQ/FAQformat.html#format3
    
    GTF (Gene Transfer Format) is a refinement to GFF that tightens the specification. The
    first eight GTF fields are the same as GFF. The group field has been expanded into a
    list of attributes. Each attribute consists of a type/value pair. Attributes must end
    in a semi-colon, and be separated from any following attribute by exactly one space. 
    
    GTF format: http://genome.ucsc.edu/FAQ/FAQformat.html#format4
    
    We do NOT check if features (exon, CDS, etc) belonging to the same gene originally were
    converted into the same chromosome/strand. 
    '''	
    
	if outfile is not None:
		FILE_OUT = open(outfile,'w')
		UNMAP = open(outfile + '.unmap', 'w')
    
	for line in ireader.reader(ingff):
		if line.startswith('#'):continue
		if line.startswith('track'):continue
		if line.startswith('browser'):continue
		if line.startswith('visibility'):continue
		if not line.strip():continue
		
		line=line.strip()
		fields=line.split('\t')
		# check GFF file
		#if len(fields) != 9:
		#	print >>sys.stderr, 'GFF file must has 9 columns. Skip ' + line
		#	continue
		try:
			start = int(fields[3]) - 1	#0-based
			end =  int(fields[4])/1
			feature_size = end - start
		except:
			print >>sys.stderr, 'Cannot recognize \"start\" and \"end\" coordinates. Skip ' + line
			if outfile:
				print >>UNMAP, line
			continue
		if fields[6] not in ['+','-','.']:
			print >>sys.stderr, 'Cannot recognize \"strand\". Skip ' + line
			if outfile:
				print >>UNMAP, line
			continue
		
		strand = '-' if fields[6] == '-' else '+'

		chrom = fields[0]
		a = map_coordinates(mapping, chrom,start,end,strand)
		if a is None:
			if outfile is None:
				print line + '\tfail (no match to target assembly)'
			else:
				print >>UNMAP, line
			continue			
		if len(a) !=2:
			if outfile is None:
				print line + '\tfail (multpile match to target assembly)'
			else:
				print >>UNMAP, line
		else:
			if (int(a[1][2]) - int(a[1][1])) != feature_size:	# check if it is exact match
				if outfile is None:
					print line + '\tfail (not exact match)'
				else:
					print >>UNMAP, line
			fields[0] = a[1][0]			# chrom
			fields[3] = a[1][1] + 1		# start, 1-based 
			fields[4] = a[1][2]
			fields[6] = a[1][3]
			
			if outfile is None:
				print line + '\t->\t' + '\t'.join([str(i) for i in fields])
			else:
				print >>FILE_OUT, '\t'.join([str(i) for i in fields])
			
def crossmap_bam_file(mapping, chainfile, infile,  outfile_prefix, chrom_size, IS_size=200, IS_std=30, fold=3):			
	'''
	Convert genome coordinates (in BAM/SAM format) between assemblies.
	BAM/SAM format: http://samtools.sourceforge.net/
	chrom_size is target chromosome size
	'''
	
	# bam or sam
	try:
		samfile = pysam.Samfile(infile,'rb')
		if len(samfile.header) ==0:
			print >>sys.stderr, "BAM file has no header section. Exit!"
			sys.exit(1)
		bam_format = True
	except:
		samfile = pysam.Samfile(infile,'r')
		if len(samfile.header) ==0:
			print >>sys.stderr, "SAM file has no header section. Exit!"
			sys.exit(1)
		bam_format = False
	
	sam_ori_header = samfile.header
	comments=['Liftover from original BAM/SAM file: ' + infile]
	comments.append('Liftover is based on the chain file: ' + chainfile)
	
	(new_header, name_to_id) = sam_header.bam_header_generator(orig_header = sam_ori_header, chrom_size = chrom_size, prog_name="CrossMap",prog_ver = __version__, format_ver=1.0,sort_type = 'coordinate',co=comments)
	
	if outfile_prefix is not None:
		if bam_format:
			OUT_FILE = pysam.Samfile( outfile_prefix + '.bam', "wb", header = new_header )
			OUT_FILE_UNMAP = pysam.Samfile( outfile_prefix + '.unmap.bam', "wb", template=samfile )
			printlog (["Liftover BAM file:", infile, '==>', outfile_prefix + '.bam'])
		else:
			OUT_FILE = pysam.Samfile( outfile_prefix + '.sam', "wh", header = new_header )
			OUT_FILE_UNMAP = pysam.Samfile( outfile_prefix + '.unmap.sam', "wh", template=samfile )
			printlog (["Liftover SAM file:", infile, '==>',  outfile_prefix + '.sam'])
	else:
		if bam_format:
			OUT_FILE = pysam.Samfile( '-', "wb", header = new_header )
			OUT_FILE_UNMAP = pysam.Samfile( infile.replace('.bam','') + '.unmap.bam', "wb", template=samfile )
			printlog (["Liftover BAM file:", infile])
		else:
			OUT_FILE = pysam.Samfile( '-', "w", header = new_header )
			OUT_FILE_UNMAP = pysam.Samfile( infile.replace('.sam','') + '.unmap.sam', "wh", template=samfile )
			printlog (["Liftover BAM file:", infile])
		#print >>sys.stderr, "aa"
		#sys.exit(2)
	
	total_item = 0
	failed = 0
	try:
		while(1):
			total_item += 1
			new_alignment = pysam.AlignedRead()	# create AlignedRead object
			old_alignment = samfile.next()
			
			#if old_alignment.is_unmapped:
			#	OUT_FILE_UNMAP.write(old_alignment)
			#	failed += 1
			#	continue
			if old_alignment.is_qcfail:
				OUT_FILE_UNMAP.write(old_alignment)
				failed += 1
				continue
			# new_alignment.qqual = old_alignment.qqual	# quality string of read, exclude soft clipped part
			# new_alignment.qlen = old_alignment.qlen	# length of the "aligned part of read"
			# new_alignment. aligned_pairs = old_alignment. aligned_pairs	#read coordinates <=> ref coordinates
			# new_alignment_aend = old_alignment.aend		# reference End position of the aligned part (of read)
			
			#new_alignment.qname = old_alignment.qname	# 1st column. read name.
			#new_alignment.flag = old_alignment.flag	# 2nd column. subject to change. flag value 
			#new_alignment.tid = old_alignment.tid		# 3rd column. samfile.getrname(tid) == chrom name
			#new_alignment.pos = old_alignment.pos		# 4th column. reference Start position of the aligned part (of read) [0-based]
			#new_alignment.mapq = old_alignment.mapq	# 5th column. mapping quality
			#new_alignment.cigar= old_alignment.cigar	# 6th column. subject to change. 
			#new_alignment.rnext = old_alignment.rnext	# 7th column. tid of the reference (mate read mapped to)
			#new_alignment.pnext = old_alignment.pnext	# 8th column. position of the reference (0 based, mate read mapped to)
			#new_alignment.tlen = old_alignment.tlen	# 9th column. insert size
			#new_alignment.seq = old_alignment.seq		# 10th column. read sequence. all bases.
			#new_alignment.qual = old_alignment.qual		# 11th column. read sequence quality. all bases.
			#new_alignment.tags = old_alignment.tags		# 12 - columns
			
			new_alignment.qname = old_alignment.qname	# 1st column. read name.
			new_alignment.seq = old_alignment.seq		# 10th column. read sequence. all bases.
			new_alignment.qual = old_alignment.qual		# 11th column. read sequence quality. all bases.
			new_alignment.tags = old_alignment.tags		# 12 - columns
			
			if old_alignment.is_paired:
				new_alignment.flag = 0x0001
				if old_alignment.is_read1: new_alignment.flag = new_alignment.flag | 0x40
				if old_alignment.is_read2: new_alignment.flag = new_alignment.flag | 0x80
				
				if old_alignment.is_unmapped:
					OUT_FILE_UNMAP.write(old_alignment)
					failed += 1
					continue

				else:	#old alignment is mapped
					try:
						read1_chr = samfile.getrname(old_alignment.tid)
					except:
						OUT_FILE_UNMAP.write(old_alignment)
						failed += 1
						continue
					read1_strand = '-' if old_alignment.is_reverse else '+'
					read1_start = old_alignment.pos
					read1_end = old_alignment.aend
					read1_maps = map_coordinates(mapping, read1_chr, read1_start, read1_end, read1_strand)
				
					if not old_alignment.mate_is_unmapped:
						try:
							read2_chr = samfile.getrname(old_alignment.rnext)									
							read2_strand = '-' if old_alignment.mate_is_reverse else '+'
							read2_start = old_alignment.pnext
							read2_end = read2_start + 1
							read2_maps = map_coordinates(mapping, read2_chr, read2_start, read2_end, read2_strand)
							# [('chr1', 246974830, 246974833, '+' ), ('chr1', 248908207, 248908210, '+' )]
						except:
							read2_maps = None
				
				
				# read 1 is unmapped after conversion	
				if read1_maps is None:														 
					new_alignment.flag = new_alignment.flag | 0x4
					
					# read2 is unmapped before conversion
					if old_alignment.mate_is_unmapped:
						OUT_FILE_UNMAP.write(old_alignment)
						failed += 1
						continue
					
					# read2 is unmapped after conversion	
					elif read2_maps is None:													
						#new_alignment.flag = new_alignment.flag | 0x8	
						#old_alignment.tid = name_to_id[samfile.getrname(old_alignment.tid)]				
						OUT_FILE_UNMAP.write(old_alignment)
						failed += 1
						continue
					
					# read2 is unique mapped
					elif len(read2_maps) == 2:	
						# 2											
						if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20		
						# 3
						new_alignment.tid = name_to_id[read2_maps[1][0]]
						# 4
						new_alignment.pos = 0
						# 5
						new_alignment.mapq = old_alignment.mapq
						# 6
						new_alignment.cigar = old_alignment.cigar
						# 7
						new_alignment.rnext = name_to_id[read2_maps[1][0]]
						# 8
						new_alignment.pnext =  read2_maps[1][1]	#start
						# 9
						new_alignment.tlen = 0
						OUT_FILE.write(new_alignment)

						continue
					
					# read2 is multiple mapped
					else:
						# 2
						if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20		
						# 2
						new_alignment.flag = new_alignment.flag | 0x100
						# 3
						new_alignment.tid = name_to_id[read2_maps[1][0]]
						# 4
						new_alignment.pos = 0
						# 5
						new_alignment.mapq = 255
						# 6
						new_alignment.cigar = old_alignment.cigar
						# 7
						new_alignment.rnext = name_to_id[read2_maps[1][0]]
						# 8
						new_alignment.pnext =  read2_maps[1][1]	#start
						# 9
						new_alignment.tlen = 0
						OUT_FILE.write(new_alignment)
						continue
				
				# read1 is unique mapped
				elif len(read1_maps) == 2:
					
					# 2
					if read1_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x10
					# 3
					new_alignment.tid = name_to_id[read1_maps[1][0]]	#chrom
					# 4
					new_alignment.pos = read1_maps[1][1]	#start
					# 5
					new_alignment.mapq = old_alignment.mapq	
					# 6
					new_alignment.cigar = old_alignment.cigar
					
					# read2 unmapped before or after conversion
					if (old_alignment.mate_is_unmapped) or (read2_maps is None):
						# 2
						new_alignment.flag = new_alignment.flag | 0x8
						# 7
						new_alignment.rnext = name_to_id[read1_maps[1][0]]
						# 8
						new_alignment.pnext =  0
						# 9
						new_alignment.tlen = 0
					
						OUT_FILE.write(new_alignment)
						continue
					
					# read2 is unique mapped
					elif len(read2_maps)==2:
						# 2
						if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
						# 7
						new_alignment.rnext = name_to_id[read2_maps[1][0]]	#chrom
						# 8
						new_alignment.pnext =  read2_maps[1][1]
						# 9						
						new_alignment.tlen = abs(new_alignment.pos - new_alignment.pnext) -1 - old_alignment.alen
						# 2
						if (read2_maps[1][3] != read1_maps[1][3]) and (new_alignment.tlen <= IS_size + fold * IS_std) and (new_alignment.tlen >= IS_size - fold * IS_std):
							new_alignment.flag = new_alignment.flag | 0x2
						OUT_FILE.write(new_alignment)
						continue
					
					# read2 is multiple mapped
					else:
						# 2 (strand)
						if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
						# 2 (secondary alignment)
						new_alignment.flag = new_alignment.flag | 0x100
						# 7
						new_alignment.rnext = name_to_id[read2_maps[1][0]]	#chrom
						# 8
						new_alignment.pnext =  read2_maps[1][1]	#start
						# 9
						new_alignment.tlen = abs(new_alignment.pos - new_alignment.pnext) -1 - old_alignment.alen
						
						OUT_FILE.write(new_alignment)
						continue
				
				# read1 is multiple mapped
				elif len(read1_maps) > 2 and len(read1_maps) % 2 ==0:
					# 2
					new_alignment.flag = new_alignment.flag | 0x100						
					# 2
					if read1_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x10
					# 3
					new_alignment.tid = name_to_id[read1_maps[1][0]]	#chrom
					# 4
					new_alignment.pos = read1_maps[1][1]	#start	
					# 5
					new_alignment.mapq = 255
					# 6
					new_alignment.cigar = old_alignment.cigar

					# read2 is unmapped
					if (old_alignment.mate_is_unmapped) or (read2_maps is None):
						# 2
						new_alignment.flag = new_alignment.flag | 0x8
						# 7
						new_alignment.rnext = name_to_id[read1_maps[1][0]]
						# 8
						new_alignment.pnext =  0
						# 9
						new_alignment.tlen = 0
					
						OUT_FILE.write(new_alignment)
						continue
					
					# read2 is unique mapped
					elif len(read2_maps)==2:
						# 2
						if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
						# 7
						new_alignment.rnext = name_to_id[read2_maps[1][0]]	#chrom
						# 8
						new_alignment.pnext =  read2_maps[1][1]
						# 9						
						new_alignment.tlen = abs(new_alignment.pos - new_alignment.pnext) -1 - old_alignment.alen
						
						OUT_FILE.write(new_alignment)
						continue		
					
					# read2 is multiple mapped
					else:
						# 2 (strand)
						if read2_maps[1][3] == '-': new_alignment.flag = new_alignment.flag | 0x20
						# 2 (secondary alignment)
						new_alignment.flag = new_alignment.flag | 0x100
						# 7
						new_alignment.rnext = name_to_id[read2_maps[1][0]]	#chrom
						# 8
						new_alignment.pnext =  read2_maps[1][1]	#start
						# 9
						new_alignment.tlen = abs(new_alignment.pos - new_alignment.pnext) -1 - old_alignment.alen
						
						OUT_FILE.write(new_alignment)
						continue				
				else:
					#old_alignment.tid = name_to_id[samfile.getrname(old_alignment.tid)]	
					OUT_FILE_UNMAP.write(old_alignment)	
					failed += 1
			# Singel end sequencing
			else:
				if old_alignment.is_unmapped:
					OUT_FILE_UNMAP.write(old_alignment)
					failed +=1
					continue
				new_alignment.flag = 0
				read_chr = samfile.getrname(old_alignment.tid)
				read_strand = '-' if old_alignment.is_reverse else '+'
				
				read_start = old_alignment.pos
				read_end = old_alignment.aend
				
				read_maps = map_coordinates(mapping, read_chr, read_start, read_end, read_strand)
				
				# unmapped
				if read_maps is None:
					#new_alignment.flag = new_alignment.flag | 0x4
					#old_alignment.tid = name_to_id[samfile.getrname(old_alignment.tid)]	
					OUT_FILE_UNMAP.write(old_alignment)
					failed += 1
					continue
				if len(read_maps)==2:
					# 2
					if read_maps[1][3] == '-': 	new_alignment.flag = new_alignment.flag | 0x10
					# 3
					new_alignment.tid = name_to_id[read_maps[1][0]]
					# 4
					new_alignment.pos = read_maps[1][1]
					# 5
					new_alignment.mapq = old_alignment.mapq
					# 6
					new_alignment.cigar = old_alignment.cigar
					OUT_FILE.write(new_alignment)
					continue
				if len(read_maps)==2 and len(read_maps) % 2 ==0:
					# 2
					new_alignment.flag = new_alignment.flag | 0x100
					# 2
					if read_maps[1][3] == '-': 	new_alignment.flag = new_alignment.flag | 0x10
					# 3
					new_alignment.tid = name_to_id[read_maps[1][0]]
					# 4
					new_alignment.pos = read_maps[1][1]
					# 5
					new_alignment.mapq = old_alignment.mapq
					# 6
					new_alignment.cigar = old_alignment.cigar
					OUT_FILE.write(new_alignment)
					continue
					
								
	except StopIteration:
		printlog(["Done!"])	
	printlog (["Total entries:", str(total_item)])
	printlog (["Failed to map:", str(failed)])
	
	if outfile_prefix is not None:
		if bam_format:
			try:
				printlog (['Sort "%s" ...' % (outfile_prefix + '.bam')])
				pysam.sort(outfile_prefix + '.bam', outfile_prefix + '.sorted')
				printlog (['Index "%s" ...' % (outfile_prefix + '.sorted.bam')])
				pysam.index(outfile_prefix + '.sorted.bam',outfile_prefix + '.sorted.bam.bai')
			except:
				printlog(["Warning: ","output BAM file was NOT sorted and indexed."]) 
	
	

def crossmap_wig_file(mapping, in_file, out_prefix, source_chrom_size, taget_chrom_size, in_format, binSize=100000):
	'''Convert genome coordinates (in wiggle/bigwig format) between assemblies.
	wiggle format: http://genome.ucsc.edu/goldenPath/help/wiggle.html
	bigwig format: http://genome.ucsc.edu/goldenPath/help/bigWig.html
	'''

	OUT_FILE1 = open(out_prefix + '.bgr','w')	# original bgr file
	OUT_FILE2 = open(out_prefix + '.sorted.bgr','w')	# sorted bgr file

	if in_format.upper() == "WIGGLE":
		printlog (["Liftover wiggle file:", in_file, '==>', out_prefix + '.bgr'])
		#print >>sys.stderr, 'liftover wiggle format file is strongly discouraged.\nIt is slow and genrated largefile.\nConvert it into bigwig will increase speed and reduce file size'
		
		for chr, start, end, strand, score in wiggleReader (in_file):
			#print chr,start,end,score
			# map_coordinates(mapping, q_chr, q_start, q_end, q_strand = '+', print_match=False):
			maps = map_coordinates(mapping, chr, start, end, '+')
			if maps is None:
				continue
			if len(maps) == 2:
				print >>OUT_FILE1, '\t'.join([str(i) for i in [maps[1][0],maps[1][1],maps[1][2], score]])
			else:
				continue
			maps[:]=[]
		OUT_FILE1.close()
		
		printlog (["Merging overlapped entries in bedGraph file ..."])
		for (chr, start, end, score) in bgrMerge.merge(out_prefix + '.bgr'):
			print >>OUT_FILE2, '\t'.join([str(i) for i in (chr, start, end, score )])
		OUT_FILE2.close()
		
		os.remove(out_prefix + '.bgr')
			
		wigToBigWig_cmd = myutils.which('wigToBigWig') 
		if wigToBigWig_cmd is not None:
			printlog (["Convert wiggle to bigwig ..."])
			tmp_file = bgrMerge.randomword(10) + '.chromsize'
			CS = open(tmp_file,'w')
			for k in sorted(taget_chrom_size):
				print >>CS, str(k) + '\t' + str(taget_chrom_size[k])
			CS.close()
			try:
				subprocess.call([wigToBigWig_cmd, "-clip",  out_prefix + '.sorted.bgr', tmp_file, out_prefix + '.bw'])
				#print ' '.join([wigToBigWig_cmd, "-clip",  out_prefix + '.3.bgr', 'tmp.chromSizes', out_prefix + '.bw'])
			except:
				printlog (["Failed!"])
			try:
				os.remove(tmp_file)
			except:
				pass
		else:
			printlog('Could not find wigToBigWig executable. Try to download it from "http://hgdownload.cse.ucsc.edu/admin/exe/" and convert wig to bigwig by yourself.')

		
	
	elif in_format.upper() == "BIGWIG":
		
		printlog (["Liftover bigwig file:", in_file, '==>', out_prefix + '.bgr'])
		for chr, start, end, score in bigwigReader(in_file, chrom_sizes =source_chrom_size, bin_size = binSize ):
			if score == 0:	#skip range with score of 0
				continue
			maps = map_coordinates(mapping, chr, start, end, '+')
			try:
				if maps is None: continue
				if len(maps) == 2:
					print >>OUT_FILE1, '\t'.join([str(i) for i in [maps[1][0],maps[1][1],maps[1][2], score]])
				else:
					continue
			except:
				continue
			maps[:]=[]
		OUT_FILE1.close()
		
		printlog (["Merging overlapped entries in bedGraph file ..."])
		for (chr, start, end, score) in bgrMerge.merge(out_prefix + '.bgr'):
			print >>OUT_FILE2, '\t'.join([str(i) for i in (chr, start, end, score )])
		OUT_FILE2.close()	
			
		wigToBigWig_cmd = myutils.which('wigToBigWig') 
		if wigToBigWig_cmd is not None:
			printlog (["Convert wiggle to bigwig ..."])
			tmp_file = bgrMerge.randomword(10) + '.chromsize'
			CS = open(tmp_file,'w')
			for k in sorted(taget_chrom_size):
				print >>CS, str(k) + '\t' + str(taget_chrom_size[k])
			CS.close()
			try:
				subprocess.call([wigToBigWig_cmd, "-clip",  out_prefix + '.sorted.bgr', tmp_file, out_prefix + '.bw'])
				#print ' '.join([wigToBigWig_cmd, "-clip",  out_prefix + '.3.bgr', 'tmp.chromSizes', out_prefix + '.bw'])
			except:
				printlog (["Failed!"])
			try:
				os.remove(tmp_file)
			except:
				pass				
		else:
			printlog(['Could not find wigToBigWig executable. Try to download it from "http://hgdownload.cse.ucsc.edu/admin/exe/" and convert wig to bigwig by yourself.'])
				
	else:
		raise Exception("Unknown foramt. Must be 'wiggle' or 'bigwig'")
	