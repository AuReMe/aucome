#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
usage:
    aucome analysis --run=ID [--cpu=INT] [-v]

options:
    --run=ID    Pathname to the comparison workspace.
    --cpu=INT     Number of cpu to use for the multiprocessing (if none use 1 cpu).
    -v     Verbose.

"""

import configparser
import csv
import docopt
import eventlet
import mpwt
import os
import re
import subprocess
import time

from aucome.utils import parse_config_file
from Bio import SeqIO
from multiprocessing import Pool


def command_help():
    print(docopt.docopt(__doc__))


def analysis_parse_args(command_args):
    args = docopt.docopt(__doc__, argv=command_args)
    run_id = args['--run']
    verbose = args['-v']

    if args["--cpu"]:
        nb_cpu_to_use = int(args["--cpu"])
    else:
        nb_cpu_to_use = 1

    run_analysis(run_id, nb_cpu_to_use, verbose)


def run_analysis(run_id, nb_cpu_to_use, verbose):

    config_data = parse_config_file(run_id)

    padmet_utils_path = config_data['padmet_utils_path']
    database_path = config_data['database_path']
    padmet_from_networks_path = config_data['padmet_from_networks_path']
    analysis_path = config_data['analysis_path']

    cmds = ["python3",  padmet_utils_path + "/padmet_utils/exploration/compare_padmet.py", "--padmet", padmet_from_networks_path,
            "--output", analysis_path, "--padmetRef", database_path]

    if verbose:
        cmds.append('-v')

    subprocess.call(cmds)

    cmds = ["python3",  padmet_utils_path + "/padmet_utils/exploration/dendrogram_reactions_distance.py", "--reactions", analysis_path + '/reactions.csv',
            "--output", analysis_path + '/dendrogram_output', "--padmetRef", database_path]

    if verbose:
        cmds.append('-v')

    subprocess.call(cmds)
