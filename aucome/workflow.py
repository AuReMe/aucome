#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome workflow --run=ID [-S=STR] [--orthogroups] [--cpu=INT] [-v] [--vv] [--filtering[=FLOAT]]

options:
    --run=ID    Pathname to the comparison workspace.
    --orthogroups    Use Orthogroups instead of Orthologues after Orthofinder.
    -S=STR    Sequence search program for Orthofinder [Default: diamond].
        Options: blast, mmseqs, blast_gz, diamond
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
    --vv    Very verbose.
    --filtering[=FLOAT]     Use a filter to limit propagation [Default: 0.05].
"""

import docopt

import aucome


def command_help():
    print(docopt.docopt(__doc__))


def workflow_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    orthogroups = args['--orthogroups']
    sequence_search_prg = args['-S']
    verbose = args['-v']
    veryverbose = args['--vv']
    filtering = args['--filtering']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    if veryverbose and not verbose:
        verbose = veryverbose

    run_workflow(run_id, nb_cpu_to_use, orthogroups, sequence_search_prg, filtering, verbose, veryverbose)


def run_workflow(run_id, nb_cpu_to_use, orthogroups, sequence_search_prg, filtering, verbose, veryverbose=None):
    aucome.check.run_check(run_id, nb_cpu_to_use, verbose)

    aucome.reconstruction.run_reconstruction(run_id, nb_cpu_to_use, verbose, veryverbose)

    aucome.orthology.run_orthology(run_id, orthogroups, sequence_search_prg, nb_cpu_to_use, filtering, verbose, veryverbose)

    aucome.structural.run_structural(run_id, nb_cpu_to_use, verbose)

    aucome.merge.run_merge(run_id, nb_cpu_to_use, verbose, veryverbose)