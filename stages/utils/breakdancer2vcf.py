import sys
import time
import glob
import fileinput
import itertools
from time import strftime

#input arguments: sys.argv[n]
#1-script.py, 2-refname 3-bd_calls.txt 4-outputname

# opens a valid file breakdancer output file
def read_breakdancer(path):
   file = open(path)  #open file connection
   h,raw,table = False,[],[] #init variables
   #scan file for header, read raw rows
   for line in file:
      if h:
         s = line.replace('\n','')
         raw.append(s) 
      elif line.rfind('#Chr') != -1:
         s = line.replace('\n','')
         h = True
   file.close()
   #split each row by the \t and store in table variable
   for row in raw:
      table.append(row.split('\t'))     
   #CTX type issue makes the number of columns 11 instead of 12....
#   for i in range(1,len(table)):
#      if (len(table[i-1]) != len(table[i])):
#         return "Break Dancer Format Error:\nMissmatched Dimensions"
   return table

# writes a new .vcf file from a formatted table
def write_vcf(path,header,vcf_table):
   file = open(path, 'w')
   s = header
   t = ''
   for i in vcf_table:
      for j in i: t+=str(j)+'\t'
      t+='\n'
   file.write(s+t)
   file.close()

#Build the VCF4.0 Header
def vcf_header(ref):
   form = '##fileformat=VCFv4.1\n'
   date = '##fileDate=' + strftime('%Y%m%d',time.localtime()) + '\n' 
   src  = '##source=BreakDancer_Max-1.4.5\n'
   ref  = '##reference=' + ref + '\n'
   inf1 = '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n'
   inf2 = '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n'
   inf3 = '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n'
   inf4 = '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
   alt1 = '##ALT=<ID=CNV,Description="Copy number variable region">\n'
   alt2 = '##ALT=<ID=DEL,Description="Deletion">\n'
   alt3 = '##ALT=<ID=INS,Description="Insertion">\n'
   alt4 = '##ALT=<ID=TRA,Description="Translocation Event">\n' #this should be updated to a BND event?
   alt5 = '##ALT=<ID=DUP,Description="Duplication Event">\n'
   t_hd = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
   return (form+date+src+ref+inf1+inf2+inf3+inf4+alt1+alt2+alt3+alt4+alt5+t_hd)

#selects needed fields and converts what is needed to make a .vcf
#:::TO DO::: include additional CHR2 VCF tags in the INFO field
def build_vcf(table):
   vcf_table = []
   #main loop to work through the rows...
   for i in range(0,len(table)):
      CHR = table[i][0]
      POS = table[i][1]
      ID = 'breakdancer_' + str(i)
      REF = '.'
      #convert breakdancer type to vcf type
      t = table[i][6]
      SVTYPE  = 'CNV' #default type
      if   t == 'DEL': SVTYPE = 'DEL'
      elif t == 'INS': SVTYPE = 'INS'
      elif t == 'INV': SVTYPE = 'INV'
      elif t == 'ITX': SVTYPE = 'DUP' #intrachromasomal translocation chrx->chrx
      elif t == 'CTX': SVTYPE = 'TRA' #interchromasomal translocation chrx->chry
      elif t == 'DUP': SVTYPE = 'DUP'
      elif t == 'CNV': SVTYPE = 'CNV'
      ALT = '<' + SVTYPE + '>'
      QUAL  = table[i][8]
      FILTER = 'PASS'
      END = table[i][4]
      SVLEN = table[i][7]
      INFO = 'END='+END+';SVTYPE='+SVTYPE+';SVLEN='+SVLEN+';IMPRECISE;'
      vcf_table += [[CHR,POS,ID,REF,ALT,QUAL,FILTER,INFO]]
   max_seq = max([len(i[0]) for i in vcf_table])
   vcf_table = sorted(vcf_table,key=lambda x: (x[0].zfill(max_seq),int(x[1])))
   return vcf_table

#Test Code Here
#input arguments: sys.argv[n]
#0-script.py, 1-refname 2-bd_calls_dir
#glob_path = '/Users/tbecker/Documents/CourseWork/15_2015_Fall/test/'
#calls = glob.glob(glob_path+'*/*_S4.calls')
#for call in calls:
#    table = read_breakdancer(call)
#    write_vcf(call[0:-6]+'.vcf',vcf_header('human_g1k_v37_decoy'),build_vcf(table))

















