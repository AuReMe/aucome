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
import pandas as pa
import re
import subprocess
import sys
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

    aucome_pool = Pool(nb_cpu_to_use)

    config_data = parse_config_file(run_id)

    analysis_path = config_data['analysis_path']
    analysis_group_file_path = config_data['analysis_group_file_path']

    group_data = []
    with open(analysis_group_file_path, 'r') as group_file:
        group_reader = csv.reader(group_file, delimiter='\t')
        for row in group_reader:
            group_name = row[0]
            groups = [org_name for org_name in row[1:] if org_name]
            tmp_data = {'padmet_utils_path': config_data['padmet_utils_path'], 'database_path': config_data['database_path'],
                        'padmet_from_networks_path': config_data['padmet_from_networks_path'], 'analysis_path': analysis_path,
                        'verbose': verbose, 'group_name': group_name,
                        'groups': groups}

            group_data.append(tmp_data)

    aucome_pool.map(analysis_on_group, group_data)

    aucome_pool.close()
    aucome_pool.join()

def analysis_on_group(group_data):
    group_name = group_data['group_name']
    groups = group_data['groups']
    padmet_utils_path = group_data['padmet_utils_path']
    database_path = group_data['database_path']
    padmet_from_networks_path = group_data['padmet_from_networks_path']
    analysis_path = group_data['analysis_path']
    verbose = group_data['verbose']
    group_name = group_data['group_name']
    all_padmet_path = [os.path.join(padmet_from_networks_path,name+".padmet") for name in groups ]
    group_analysis_path = analysis_path + '/' + group_name

    if not os.path.isdir(group_analysis_path):
        if len(groups) == 1:
            sys.exit('A group must contain more than one member.')

        for padmet_path in all_padmet_path:
            if not os.path.exists(padmet_path):
                org_name = os.path.splitext(os.path.basename(padmet_path))[0]
                sys.exit("Padmet file of organism %s from group %s not found in %s" %(org_name, group_name, padmet_from_networks_path))
        
        cmds = ["python3",  padmet_utils_path + "/padmet_utils/exploration/compare_padmet.py", "--padmet", ",".join(all_padmet_path),
                "--output", group_analysis_path, "--padmetRef", database_path]

        if verbose:
            cmds.append('-v')

        subprocess.call(cmds)

        cmds = ["python3",  padmet_utils_path + "/padmet_utils/exploration/dendrogram_reactions_distance.py", "--reactions", group_analysis_path + '/reactions.csv',
                "--output", group_analysis_path + '/dendrogram_output', "--padmetRef", database_path]

        if verbose:
            cmds.append('-v')

        subprocess.call(cmds)

    else:
        print(group_analysis_path + ' already exists. Delete it if you want to relaunch the analysis.')