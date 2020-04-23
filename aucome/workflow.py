#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome workflow --run=ID [-S=STR] [--orthogroups] [--keep-tmp] [--cpu=INT] [-v] [--vv] [--filtering] [--threshold=FLOAT]

options:
    --run=ID    Pathname to the comparison workspace.
    --orthogroups    Use Orthogroups instead of Orthologues after Orthofinder.
    -S=STR    Sequence search program for Orthofinder [Default: diamond].
        Options: blast, mmseqs, blast_gz, diamond
    --keep-tmp    Keep temporary file (especially sequence of predicted gene linked to reaction).
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
    --vv    Very verbose.
    --filtering     Use a filter to limit propagation, by default it is 0.05, if you want to modify the value use --threshold.
    --threshold=FLOAT     Threshold of the filter to limit propagation to use with the --filtering argument.
"""
import aucome
import docopt
import sys
import time


def command_help():
    print(docopt.docopt(__doc__))


def workflow_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    orthogroups = args['--orthogroups']
    sequence_search_prg = args['-S']
    keep_tmp = args['--keep-tmp']
    verbose = args['-v']
    veryverbose = args['--vv']
    filtering = args['--filtering']
    filtering_threshold = args['--threshold']

    if filtering:
        if not filtering_threshold:
            filtering_threshold = 0.05
    else:
        if filtering_threshold:
            sys.exit('--threshold must be used with --filtering.')

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    if veryverbose and not verbose:
        verbose = veryverbose

    run_workflow(run_id, nb_cpu_to_use, orthogroups, sequence_search_prg, filtering_threshold, keep_tmp, verbose, veryverbose)


def run_workflow(run_id, nb_cpu_to_use, orthogroups, sequence_search_prg, filtering_threshold, keep_tmp, verbose, veryverbose=None):
    if verbose:
        print('--- Running workflow ---')
    workflow_start_time = time.time()

    aucome.check.run_check(run_id, nb_cpu_to_use, verbose, veryverbose)

    aucome.reconstruction.run_reconstruction(run_id, nb_cpu_to_use, verbose, veryverbose)

    aucome.orthology.run_orthology(run_id, orthogroups, sequence_search_prg, nb_cpu_to_use, filtering_threshold, verbose, veryverbose)

    aucome.structural.run_structural(run_id, keep_tmp, nb_cpu_to_use, verbose)

    aucome.merge.run_merge(run_id, nb_cpu_to_use, verbose, veryverbose)

    workflow_end_time = (time.time() - workflow_start_time)
    integer_part, decimal_part = str(workflow_end_time).split('.')
    workflow_time = ".".join([integer_part, decimal_part[:3]])

    if verbose:
        print("--- workflow step done in: %ss ---" %workflow_time)
