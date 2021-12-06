#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome workflow --run=ID [--sequence_search_prg=STR] [--keep-tmp] [--cpu=INT] [-v] [--vv] [--filtering] [--threshold=FLOAT] [--union] [--intersection]

options:
    --run=ID    Pathname to the comparison workspace.
    --sequence_search_prg=STR    Sequence search program for Orthofinder [Default: diamond].
        Options: blast, mmseqs, blast_gz, diamond
    --keep-tmp    Keep temporary file (especially sequence of predicted gene linked to reaction).
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
    --vv    Very verbose.
    --filtering     Use a filter to limit propagation, by default it is 0.05, if you want to modify the value use --threshold.
    --threshold=FLOAT     Threshold of the filter to limit propagation to use with the --filtering argument.
    --union          Use the union filter between five threshold values [0.01, 0.05, 0.1, 0.15, 0.2] to limit propagation, to use with the --filtering argument.
    --intersection   Use the intersection filter between five threshold values [0.01, 0.05, 0.1, 0.15, 0.2] to limit propagation, to use with the --filtering argument.
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
    sequence_search_prg = args['--sequence_search_prg']
    cpu = args['--cpu']
    keep_tmp = args['--keep-tmp']
    verbose = args['-v']
    veryverbose = args['--vv']
    filtering = args['--filtering']
    threshold = args['--threshold']
    union = args['--union']
    intersection = args['--intersection']
    filtering_threshold_list = []
    
    if filtering:
        if threshold:
            try:
                filtering_threshold_list.append(float(threshold))
            except:
                sys.exit('filtering_threshold value must be a float between 0 and 1.')       
        else:
            filtering_threshold_list.append(0.05)      
        if union or intersection:
            filtering_threshold_list = [0.01, 0.05, 0.1, 0.15, 0.2]
    else:
        if threshold:
            sys.exit('--threshold must be used with --filtering.')
        if union:
            sys.exit('--union must be used with --filtering.')
        if intersection:
            sys.exit('--intersection must be used with --filtering.')
            
    if cpu:
        nb_cpu_to_use = int(cpu)
    else:
        nb_cpu_to_use = 1

    if veryverbose and not verbose:
        verbose = veryverbose

    run_workflow(run_id, nb_cpu_to_use, sequence_search_prg, filtering_threshold_list, union, intersection, keep_tmp, verbose, veryverbose)


def run_workflow(run_id, nb_cpu_to_use, sequence_search_prg, filtering_threshold_list, union, intersection, keep_tmp, verbose, veryverbose=None):
    if verbose:
        print('--- Running workflow ---')
    workflow_start_time = time.time()

    aucome.check.run_check(run_id, nb_cpu_to_use, verbose, veryverbose)

    aucome.reconstruction.run_reconstruction(run_id, nb_cpu_to_use, verbose, veryverbose)

    aucome.orthology.run_orthology(run_id, sequence_search_prg, nb_cpu_to_use, filtering_threshold_list, union, intersection, verbose, veryverbose)

    aucome.structural.run_structural(run_id, keep_tmp, nb_cpu_to_use, verbose)

    aucome.merge.run_merge(run_id, nb_cpu_to_use, verbose, veryverbose)

    workflow_end_time = (time.time() - workflow_start_time)
    integer_part, decimal_part = str(workflow_end_time).split('.')
    workflow_time = ".".join([integer_part, decimal_part[:3]])

    if verbose:
        print("--- workflow step done in: %ss ---" %workflow_time)
