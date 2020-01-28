#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome workflow --run=ID [-S=STR] [--orthogroups] [--cpu=INT] [-v] [--filtering]

options:
    --run=ID    Pathname to the comparison workspace.
    --orthogroups    Use Orthogroups instead of Orthologues after Orthofinder.
    -S=STR    Sequence search program for Orthofinder [Default: diamond].
        Options: blast, mmseqs, blast_gz, diamond
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
    --filtering     Use a filter to limit propagation.
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
    filtering = args['--filtering']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_workflow(run_id, nb_cpu_to_use, orthogroups, sequence_search_prg, filtering, verbose)


def run_workflow(run_id, nb_cpu_to_use, orthogroups, sequence_search_prg, filtering, verbose):
    aucome.check.run_check(run_id, nb_cpu_to_use, verbose)

    aucome.reconstruction.run_reconstruction(run_id, nb_cpu_to_use, verbose)

    aucome.orthology.run_orthology(run_id, orthogroups, sequence_search_prg, nb_cpu_to_use, filtering, verbose)

    aucome.structural.run_structural(run_id, nb_cpu_to_use, verbose)

    aucome.merge.run_merge(run_id, nb_cpu_to_use, verbose)