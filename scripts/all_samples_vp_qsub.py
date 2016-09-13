#!/usr/bin/env python
import argparse
import os
import glob
import subprocess32 as subprocess  # to call qsub a bunch of times

des = """
script that auto generates PBS scripts for runing a single caller (or caller series) on a folder of bams via qsub
saving the results to a single root directory with sample-name subfolders for configuration/intermediate files
and temporary fasta, bed, bedpe, sam, bam, calls, vcf, root, histograms, jobs summaries and stats generated."""
parser = argparse.ArgumentParser(description=des)
parser.add_argument('-r', '--ref_path', type=str, help='reference fasta')
parser.add_argument('-b', '--bam_dir', type=str, help='bam file or directory')
parser.add_argument('-t', '--tar_dir', type=str, help='target file or directory')
parser.add_argument('-D', '--read_depth', type=int, help='Average Depth of Coverage in ROI')
parser.add_argument('-L', '--read_length', type=int, help='Read Length')
parser.add_argument('-s', '--stage_list', type=str, help='stage list to apply to all bam files')
parser.add_argument('-o', '--output_dir', type=str, help='outputdirectory to save into')
parser.add_argument('-a', '--all_at_once', action='store_true', help='run all bam files together')
parser.add_argument('-w', '--wall_time', type=str, help='wall time requested from cluster')
parser.add_argument('-m', '--memory', type=str, help='radom access memory needed')
parser.add_argument('-p', '--processors', type=int, help='processors needed')
parser.add_argument('-d', '--debug', action='store_true', help='save results/errors to db')
parser.add_argument('-e', '--email_address', type=str, help='cluster email results to this email address')
args = parser.parse_args()

if args.stage_list is not None:
    stages = args.stage_list.split(',')  # don't have to split, just pass in...
else:
    print('stage not found')
    raise IOError
if args.ref_path is not None:
    ref_path = args.ref_path
elif len(stages)==1 and stages[0]=='bam_stats':
    print('not using a reference for bam_stats')
    ref_path = ''
else:
    print('reference fasta not found')
    raise IOError

if args.bam_dir is not None:
    if args.bam_dir.endswith('.bam'):
        bam_list = [args.bam_dir]
    else:
        bam_list = glob.glob(args.bam_dir + '/*.bam')
    if len(bam_list) < 1:
        print('bam files were not found')
        raise IOError
    samples = {}
    samples = {i.rsplit('/')[-1].split('.')[0]: i for i in bam_list}  # file name from path
else:
    print('bam file or directory not found')
    raise IOError

if args.tar_dir is not None:
    if args.tar_dir.endswith('.vcf') or args.tar_dir.endswith('.calls'):
        tar_list = [args.tar_dir]
    else:
        tar_list = glob.glob(args.tar_dir + '/*.vcf')  # try vcf and
        tar_list += glob.glob(args.tar_dir + '/*.calls')  # g1k format
    if len(tar_list) < 1:
        print('target files were not found')
    tar_samples = {}
    tar_samples = {i.replace('.vcf', '').replace('.calls', '').rsplit('/')[-1].rsplit('_')[0].split('.')[0]: i for i in
                   tar_list}  # file name from path
    print('target files located are:')
    for k in tar_samples: print(k)
else:
    tar_list = []
    tar_samples = samples
    print('target file or directory not found')

# target file or directory following the same usage as bam dir

# collect read depth and length estimates for RD analysis
RD, RL = None, None
if args.read_depth is not None:
    RD = args.read_depth  # CSL returns a list of 1+
else:
    print('read_depth not given')
if args.read_length is not None:
    RL = args.read_length  # CSL returns a list of 1+
else:
    print('read_length not given')
# collect read depth and length estimates for RD analysis
if not (RD is None) and not (RL is None):
    auto_RD_RL = False
else:
    auto_RD_RL = True

if args.output_dir is not None:
    out_dir = args.output_dir
    if not os.path.exists(out_dir): os.makedirs(out_dir)
else:
    print('out put directory not specified')
    raise IOError

if args.wall_time is not None:
    walltime = args.wall_time
else:
    walltime = '24:00:00'

if args.memory is not None:
    ram = args.memory
else:
    ram = '32gb'

if args.processors is not None:
    cpus = args.processors
else:
    cpus = 1

if args.email_address is not None:
    email = args.email_address
else:
    print('email not valid')
    raise AttributeError

# make a temp directory for pbs scripts
pbs_dir = out_dir + 'PBS'
if not os.path.exists(pbs_dir): os.makedirs(pbs_dir)

# write out the .pbs scripts
software_path = os.path.dirname(os.path.abspath(__file__)) + '/../../'
python = software_path + '/anaconda/bin/python'
vp = software_path + '/SVE/scripts/variant_processor.py'
ml = 'module load '  # everything that you need below...
modules = [ml + 'perl/5.16.3', ml + 'gcc/4.9.2', ml + 'R/3.2.1',
           ml + 'Root/v5.34.18', ml + 'samtools/1.2',
           ml + 'pbs-drmaa/1.0.17']  # everything that you need here...
PBS = []
if args.all_at_once:
    print('doing all samples at once...')
    bams = ','.join([samples[k] for k in samples])
    group_pbs = pbs_dir + '/' + 'all_' + '-'.join(stages) + '.pbs'  # name of pbs script for the group
    sub_dir = out_dir + '/' + 'all' + '/'
    PBS += [group_pbs]
    if not os.path.exists(sub_dir): os.makedirs(sub_dir)  # make the sub_directory
    with open(group_pbs, 'w') as pbs:
        if auto_RD_RL:
            print('using the automatic RD and RL estimation step: bam_stats stage added')
            call = [python, vp, '-r', ref_path, '-b', bams, '-s', ','.join(stages),
                    '-o', sub_dir]
            if args.debug: call += ['--debug']
        else:
            call = [python, vp, '-r', ref_path, '-b', bams, '-s', ','.join(stages),
                    '-D', str(RD), '-L', str(RL), '-o', sub_dir]
            if args.debug: call += ['--debug']
        pbs.write('#!/bin/bash\n' + '\n'.join(modules) + '\n' + ' '.join(call))

else:
    print('doing each sample seperately...')
    if len(tar_list) < 1 and set(tar_samples).issubset(set(samples)):
        print('not using targeting files..')
        for k in tar_samples:  # need these keys to match
            sample_pbs = pbs_dir + '/' + k + '-'.join(stages) + '.pbs'  # name of pbs script for sample k
            sub_dir = out_dir + '/' + k + '/'
            PBS += [sample_pbs]
            if not os.path.exists(sub_dir): os.makedirs(sub_dir)  # make the sub_directory
            with open(sample_pbs, 'w') as pbs:
                if auto_RD_RL:
                    call = [python, vp, '-r', ref_path, '-b', samples[k],
                            '-s', ','.join(stages), '-o', sub_dir]
                    if args.debug: call += ['--debug']
                else:
                    call = [python, vp, '-r', ref_path, '-b', samples[k],
                            '-s', ','.join(stages), '-o', sub_dir, '-D', str(RD), '-L', str(RL)]
                    if args.debug: call += ['--debug']
                pbs.write('#!/bin/bash\n' + '\n'.join(modules) + '\n' + ' '.join(call))
    elif set(tar_samples).issubset(set(samples)):
        print('using targeting files..')
        for k in tar_samples:  # need these keys to match
            sample_pbs = pbs_dir + '/' + k + '-'.join(stages) + '.pbs'  # name of pbs script for sample k
            sub_dir = out_dir + '/' + k + '/'
            PBS += [sample_pbs]
            if not os.path.exists(sub_dir): os.makedirs(sub_dir)  # make the sub_directory
            with open(sample_pbs, 'w') as pbs:
                if auto_RD_RL:
                    call = [python, vp, '-r', ref_path, '-b', samples[k], '-t', tar_samples[k],
                            '-s', ','.join(stages), '-o', sub_dir]
                    if args.debug: call += ['--debug']
                else:
                    call = [python, vp, '-r', ref_path, '-b', samples[k], '-t', tar_samples[k],
                            '-s', ','.join(stages), '-o', sub_dir, '-D', str(RD), '-L', str(RL)]
                    if args.debug: call += ['--debug']
                pbs.write('#!/bin/bash\n' + '\n'.join(modules) + '\n' + ' '.join(call))
    elif set(samples).issubset(set(tar_samples)):
        print('using targeting files..')
        for k in samples:  # need these keys to match
            sample_pbs = pbs_dir + '/' + k + '-'.join(stages) + '.pbs'  # name of pbs script for sample k
            sub_dir = out_dir + '/' + k + '/'
            PBS += [sample_pbs]
            if not os.path.exists(sub_dir): os.makedirs(sub_dir)  # make the sub_directory
            with open(sample_pbs, 'w') as pbs:
                if auto_RD_RL:
                    call = [python, vp, '-r', ref_path, '-b', samples[k], '-t', tar_samples[k],
                            '-s', ','.join(stages), '-o', sub_dir]
                    if args.debug: call += ['--debug']
                else:
                    call = [python, vp, '-r', ref_path, '-b', samples[k], '-t', tar_samples[k],
                            '-s', ','.join(stages), '-D', str(RD), '-L', str(RL), '-o', sub_dir]
                    if args.debug: call += ['--debug']
                pbs.write('#!/bin/bash\n' + '\n'.join(modules) + '\n' + ' '.join(call))
# execute qsub with the scripts, getting the jids back (can display these or attach to further monitor progress)
output, err = '', {}
for pbs in PBS:  # test with one of these and a fast caller on a small file...
    print('processing %s' % pbs)
    try:
        job_stub = pbs.rsplit('/')[-1].rsplit('.')[0]
        job_stub = job_stub[:max(len(job_stub),16)]
        command = ['qsub', '-l', 'walltime=%s' % walltime + ',mem=%s' % ram + ',procs=%s' % cpus,
                   '-m', 'e', '-M', email,                   # email prefs
                   '-N',job_stub,                            # name of job
                   '-o',pbs[0:-4] + '.log',                  # output log,
                   '-j oe',pbs]                              # output and pbs script to run
        print(' '.join(command))  # don't run it yet!
        output += subprocess.check_output(' '.join(command), stderr=subprocess.STDOUT, shell=True)
    # catch all errors that arise under normal call behavior
    except subprocess.CalledProcessError as E:
        print('call error: ' + E.output)  # what you would see in the term
        err['output'] = E.output
        # the python exception issues (shouldn't have any...
        print('message: ' + E.message)  # ?? empty
        err['message'] = E.message
        # return codes used for failure....
        print('code: ' + str(E.returncode))  # return 1 for a fail in art?
        err['code'] = E.returncode
    except OSError as E:
        print('os error: ' + E.strerror)  # what you would see in the term
        err['output'] = E.strerror
        # the python exception issues (shouldn't have any...
        print('message: ' + E.message)  # ?? empty
        err['message'] = E.message
        # the error num
        print('code: ' + str(E.errno))
        err['code'] = E.errno
print('output:\n' + output)

# remove/delete the intermediate pbs scripts TO DO...