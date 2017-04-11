#!/usr/bin/env python
import argparse
import os
import sys

para_dict={
'command':'',
'ref':'',
'refbase':'',
'out_dir':'',
'sample':'',
'threads':4,
'mem':8,
'algorithm':'speed_seq',
'FASTQ':['',''],
'BAM':''
}

class ParseParameters(object):
    def __init__(self):
        parser = argparse.ArgumentParser(
            add_help=False,
            usage = '''prepare_bam.py <command> [options]\n
Command:\talign\tFASTQ->BAM
        \trealign\tBAM->BAM
''')
        parser.add_argument('command', help='Subcommand to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print args.command + ' is not a recognized command.\n'
            parser.print_usage()
            exit(1)
        # use dispatch pattern to invoke method with same name
        parser1 = getattr(self, args.command)()
        getattr(self, 'common')(parser1, args.command)
        getattr(self, 'parse')(parser1)

    def align(self):
        parser = argparse.ArgumentParser(usage = "prepare_bam align [options] <-r FILE> <FASTQ1 [FASTQ2]>")
        parser.add_argument('-a', dest='algorithm', type=str, metavar='STR', choices=['aln|mem|speed_seq'], default='speed_seq', help='the method used for alignment\t[speed_seq]')
        parser.add_argument('FASTQ',nargs='+', help='input FASTQs')
        return parser

    def realign(self):
        parser = argparse.ArgumentParser(usage = "prepare_bam realign [options] <-r FILE> <BAM>")
        parser.add_argument('BAM', help='input BAM [null]')
        return parser

    def common(self, parser, command):
        parser.add_argument('-r', dest='ref', type=str, metavar='FILE', help='FASTA reference file\t[null]')
        parser.add_argument('-o', dest='out_dir', type=str, metavar='STR', default='./output', help='output directory\t[./output]')
        parser.add_argument('-s', dest='sample',type=str, metavar='STR', help='sample name\t[input]')
        parser.add_argument('-t', dest='threads',type=int, metavar='INT', help='number of threads per CPU\t[4]')
        parser.add_argument('-M', dest='mem',type=int, metavar='INT', help='ram in GB to use for per cpu/thread unit\t[8]')
        global args 
        args = parser.parse_args(sys.argv[2:])
        args.command = command

    def parse(self, parser):
        ### Output
        if args.out_dir is not None:    #optional reroute
            para_dict['out_dir'] = args.out_dir
        para_dict['out_dir'] = os.path.abspath(para_dict['out_dir'] )
        if para_dict['out_dir'] [:-1] == '/': para_dict['out_dir']  = para_dict['out_dir'][0:-1]

        ### Reference
        if args.ref is not None:
            para_dict['ref'] = args.ref
            if not os.path.isfile(para_dict['ref']):
                print "ERROR: Cannot open reference file: " + para_dict['ref']
                exit()
        else:
            print parser.print_help()
            exit()
        if (para_dict['ref'].endswith(('.fasta.gz', '.fa.gz', '.fasta', '.fa')) == False):
            print "ERROR: reference file should end with .fasta.gz, .fa.gz, .fasta or .fa."
        # Either one of them
        para_dict['refbase'] = para_dict['ref'].rsplit('/')[-1].split('.fasta.gz')[0]
        para_dict['refbase'] = para_dict['ref'].rsplit('/')[-1].split('.fa.gz')[0]
        para_dict['refbase'] = para_dict['ref'].rsplit('/')[-1].split('.fasta')[0]
        para_dict['refbase'] = para_dict['ref'].rsplit('/')[-1].split('.fa')[0]

        ### FASTQs
        if args.command == 'align':
            para_dict['FASTQ'] = args.FASTQ
            for r in para_dict['FASTQ']:
                if not os.path.isfile(r):
                    print "ERROR: Cannot open FASTQ file: " + r
                    exit()
            if args.algorithm is not None: para_dict['algorithm'] = args.algorithm
        ### BAM
        elif args.command == 'realign':
            para_dict['BAM'] = args.BAM
            if not os.path.isfile(para_dict['BAM']): 
                print "ERROR: Cannot open BAM file: " + para_dict['BAM']
                exit()

        if args.sample is not None:    para_dict['sample'] = args.sample
        if args.threads is not None:   para_dict['threads'] = args.threads
        if args.mem is not None:       para_dict['mem'] = int(args.mem)

        if not os.path.exists(para_dict['out_dir'] ): os.makedirs(para_dict['out_dir'] )

        
