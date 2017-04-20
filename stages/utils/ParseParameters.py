#!/usr/bin/env python
import argparse
import os
import sys
import gzip

para_dict={
'machine':'',
'command':'',
'ref':'',
'refbase':'',
'genome':'',
'out_dir':'',
'out_file':'',
'threads':4,
'mem':8,
'algorithm':'speed_seq',
'RG':'',
'FASTQ':[],
'BAM':[]
}

class ParseParameters(object):
    def __init__(self, paras):
        parser = argparse.ArgumentParser(
            add_help=False,
            usage = '''sve <command> [options]\n
Command:\talign\tFASTQ->BAM
        \trealign\tBAM->BAM
        \thg38fix\tBAM->BAM
	\tcall\tBAM(s)->VCF
''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        sub_commands = ['align', 'realign', 'hg38fix', 'call']
        args = parser.parse_args(sys.argv[1:2])
        if not args.command in sub_commands:
            print args.command + ' is not a recognized command.\n'
            parser.print_usage()
            exit(1)

        # use dispatch pattern to invoke method with same name
        parser1 = getattr(self, args.command)()

        if args.command in sub_commands[0:2]:
            getattr(self, 'aln_parse')(args.command, parser1, paras)
        elif args.command in sub_commands[2]:
            getattr(self, 'fix_parse')(args.command, parser1, paras)
        elif args.command in sub_commands[3]:
            getattr(self, 'aln_parse')(args.command, parser1, paras)

    ##### sub commands #####
    def align(self):
        parser = argparse.ArgumentParser(usage = "prepare_bam align [options] <-r FILE> <FASTQ1 [FASTQ2]>")
        self.aln_common(parser)
        parser.add_argument('-a', dest='algorithm', type=str, metavar='STR', choices=['bwa_aln', 'bwa_mem', 'speed_seq'], default='speed_seq', help='the method used for alignment\t[speed_seq]')
        parser.add_argument('-R', dest='RG',type=str, metavar='STR', help='read group header line such as "@RG\\tID:id\\tSM:sampleName\\tLB:lib\\tPL:ILLUMINA"[null]')
        parser.add_argument('FASTQ',nargs='+', help='input FASTQs')
        return parser

    def realign(self):
        parser = argparse.ArgumentParser(usage = "prepare_bam realign [options] <-r FILE> <BAM>")
        self.aln_common(parser)
        parser.add_argument('-R', dest='RG',type=str, metavar='STR', help='read group header line such as "@RG\\tID:id\\tSM:sampleName\\tLB:lib\\tPL:ILLUMINA"\t[RG_in_input_BAM]')
        parser.add_argument('BAM', nargs=1, help='input BAM [null]')
        return parser

    def hg38fix(self):
        parser = argparse.ArgumentParser(usage = "prepare_bam hg38fix [options] <BAM>")
        #parser.add_argument('-p', dest='ref_alt', nargs='+', type=str, metavar='FILE', help='\t[hs38DH-extra hs38DH.fa.alt]')
        parser.add_argument('-o', dest='out_file', type=str, metavar='STR', help='output BAM\t[in_prefix.alt.bam]')
        parser.add_argument('BAM', nargs=1, help='input BAM [null]')
        return parser

    def call(self):
        parser = argparse.ArgumentParser(usage = "prepare_bam call [options] <-r FILE> <-g hg19|hg38|others> <BAM [BAM ...]>")
        self.aln_common(parser)
        parser.add_argument('-a', dest='algorithm', type=str, metavar='STR', choices=['breakdancer', 'breakseq', 'cnvnator', 'hydra', 'delly', 'lumpy', 'genome_strip', 'cnmops', 'gatk', 'tigra'], help='the method used for SV calling\t[NULL]')
        parser.add_argument('-g', dest='genome', type=str, metavar='STR', choices=['hg19', 'hg38', 'others'], help='tell us the input reference\t[NULL]')
        parser.add_argument('BAM', nargs='+', help='input BAM [null]')
        return parser
    ##### end sub commands #####

    def aln_common(self, parser):
        parser.add_argument('-r', dest='ref', type=str, metavar='FILE', help='FASTA reference file\t[null]')
        parser.add_argument('-o', dest='out_dir', type=str, metavar='STR', default='./output', help='output directory\t[./output]')
        parser.add_argument('-t', dest='threads',type=int, metavar='INT', help='number of threads per CPU\t[4]')
        parser.add_argument('-M', dest='mem',type=int, metavar='INT', help='ram in GB to use for per cpu/thread unit\t[8]')

    def load_args(self, parser):
        if len(sys.argv[2:]) == 0:
            parser.print_help()
            exit()
        global args 
        args = parser.parse_args(sys.argv[2:])
        return args

    def fix_parse(self, command, parser, paras):
        self.load_args(parser)
        paras['command'] = command

        ### BAM
        paras['BAM'] = args.BAM
        if not os.path.isfile(paras['BAM'][0]): 
            print "ERROR: Cannot open BAM file: " + paras['BAM'][0]
            exit()

        if args.out_file is not None: paras['out_file'] = args.out_file

    def aln_parse(self, command, parser, paras):
        self.load_args(parser)
        paras['command'] = command

        ### Output
        if args.out_dir is not None:    #optional reroute
            paras['out_dir'] = args.out_dir
        paras['out_dir'] = os.path.abspath(paras['out_dir'] )
        if paras['out_dir'] [:-1] == '/': paras['out_dir']  = paras['out_dir'][0:-1]

        ### Reference
        if args.ref is not None:
            paras['ref'] = args.ref
            if not os.path.isfile(paras['ref']):
                print "ERROR: Cannot open reference file: " + paras['ref']
                exit()
        else:
            print parser.print_help()
            exit()
        if (paras['ref'].endswith(('.fasta.gz', '.fa.gz', '.fasta', '.fa')) == False):
            print "ERROR: reference file should end with .fasta.gz, .fa.gz, .fasta or .fa."
        # Decompress FASTA
        if paras['ref'][-3:] == ".gz":
            decompress_fname = paras['ref'][:-3]
            if not os.path.isfile(decompress_fname):
                print "Decompressing " + paras['ref'] + " as " + decompress_fname
                with gzip.open(paras['ref'],'rb') as in_file:
                    s = in_file.read()
                with open(decompress_fname,'w') as f:
                    f.write(s)
            paras['ref'] = decompress_fname
        # Either one of them
        paras['refbase'] = paras['ref'].rsplit('/')[-1].split('.fasta')[0]
        paras['refbase'] = paras['ref'].rsplit('/')[-1].split('.fa')[0]

        ### Align
        if paras['command'] in ['align']:
            paras['FASTQ'] = args.FASTQ
            if len(paras['FASTQ']) > 2:
                print "ERROR: At most two FASTQs"
                exit()

            for r in paras['FASTQ']:
                if not os.path.isfile(r):
                    print "ERROR: Cannot open FASTQ file: " + r
                    exit()

        ### Realign and Call: for bam input
        if paras['command'] in ['realign']:
            paras['BAM'] = args.BAM
            if not os.path.isfile(paras['BAM'][0]): 
                print "ERROR: Cannot open BAM file: " + paras['BAM'][0]
                exit()

        ### Call: genome
        if paras['command'] in ['call']:
            if args.algorithm in ['genome_strip', 'delly'] and args.genome is None:
                print "ERROR: Please specify the input genome: -g hg19|hg38|others."
                exit()
            else:
                paras['genome'] = args.genome
            paras['BAM'] = args.BAM
            for bam in paras['BAM']:
                if not os.path.isfile(bam):
                    print "ERROR: Cannot open BAM file: " + bam
                    exit()
            

        ### Align and Call: for algorithm
        if paras['command'] in ['align', 'call']:
            if args.algorithm is not None: paras['algorithm'] = args.algorithm

        if paras['command'] in ['align', 'realign']:
            if args.RG is not None:        paras['RG'] = args.RG

        if args.threads is not None:   paras['threads'] = int(args.threads)
        if args.mem is not None:       paras['mem'] = int(args.mem)

        if not os.path.exists(paras['out_dir'] ): os.makedirs(paras['out_dir'] )

        
