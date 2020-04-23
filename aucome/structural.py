#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome structural --run=ID [--keep-tmp] [--cpu=INT] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --keep-tmp    Keep temporary file (especially sequence of predicted gene linked to reaction).
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
"""

import docopt
import time

from padmet.utils.exploration import prot2genome

from aucome.utils import parse_config_file


def command_help():
    print(docopt.docopt(__doc__))


def structural_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    keep_tmp = args['--keep-tmp']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_structural(run_id, keep_tmp, nb_cpu_to_use, verbose)


def run_structural(run_id, keep_tmp, nb_cpu_to_use, verbose):
    if verbose:
        print('--- Running structural check step ---')
    structural_start_time = time.time()

    config_data = parse_config_file(run_id)
    database_path = config_data['database_path']

    prot2genome.fromAucome(run_id, nb_cpu_to_use, database_path, blastp=True, tblastn=True, exonerate=True, keep_tmp=keep_tmp, debug=False)

    structural_end_time = (time.time() - structural_start_time)
    integer_part, decimal_part = str(structural_end_time).split('.')
    structural_time = ".".join([integer_part, decimal_part[:3]])

    if verbose:
        print("--- structural step done in: %ss ---" %structural_time)