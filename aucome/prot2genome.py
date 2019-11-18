#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome prot2genome --run=ID [--cpu=INT] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.
"""

import docopt
from padmet.utils.exploration import prot2genome

from aucome.utils import parse_config_file


def command_help():
    print(docopt.docopt(__doc__))


def prot2genome_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_prto2genome(run_id, nb_cpu_to_use, verbose)


def run_prto2genome(run_id, nb_cpu_to_use, verbose):

    config_data = parse_config_file(run_id)
    database_path = config_data['database_path']

    prot2genome.fromAucome(run_id, nb_cpu_to_use, database_path, blastp=True, tblastn=True, exonerate=True, debug=False)

